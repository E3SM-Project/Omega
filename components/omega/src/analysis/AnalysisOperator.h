#ifndef OMEGA_ANALYSISOP_H
#define OMEGA_ANALYSISOP_H

//===-- analysis/AnalysisOperator.h - AnalysisOperator ---*- C++ -*-===//
//
/// \file
/// \brief Defines the AnalysisOperator base class and configuration helpers
///
/// AnalysisOperator is the abstract base class from which all concrete
/// analysis operators are derived. Each operator performs a single,
/// well-defined transformation on input fields (e.g., spatial reduction,
/// temporal averaging). Derived classes are templated on the Kokkos array
/// type of their primary input field, enabling type-safe dispatch based on
/// scalar type, rank, and memory location.
///
/// Operators declare their input field dependencies at construction and
/// create their output Fields in the Field registry. The initialize()
/// method is called after all Fields exist to store mesh/environment
/// pointers needed during compute(). The compute() method retrieves input
/// data from the Field registry and writes results to operator-owned output
/// arrays.
///
/// An output Field inherits its CF metadata from the input it is derived
/// from: the base class provides helpers that collect an input's units,
/// standard name and cell methods, append the cell method describing the new
/// reduction, and create the output Field with them. Operators whose result
/// has different units (an area-weighted sum, a product) derive them with
/// CFUnits before creating the output Field.
///
/// This file also provides helper functions (opParam, makeOpConfig) for
/// constructing Config objects inline when instantiating operators,
/// providing a uniform parameter-passing mechanism whether operators are
/// created from user config or programmatically by bundled AnalysisGroups.
///
//===----------------------------------------------------------------------===//

#include "CFUnits.h"
#include "Config.h"
#include "DataTypes.h"
#include "Dimension.h"
#include "Error.h"
#include "Field.h"
#include "HorzMesh.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"
#include "VertCoord.h"

#include <string>

namespace OMEGA {

/// Base case for variadic template recursion - returns an empty Config.
/// This terminates the recursive construction in makeOpConfig.
inline Config makeOpConfig() { return Config(); }

/// Constructs a Config object from variadic key-value pairs using
/// recursive template expansion. Enables inline operator parameter
/// specification:
/// \code
///   makeOpConfig(opParam("Period", "1day"), opParam("Layer", 10))
/// \endcode
/// This provides a uniform parameter interface whether operators are
/// instantiated from user config or programmatically by AnalysisGroups.
template <typename T, typename... Args>
Config makeOpConfig(const std::pair<std::string, T> &Param, Args... OtherArgs) {
   Config Cfg = makeOpConfig(OtherArgs...); // Recurse from end
   Cfg.add(Param.first, Param.second);
   return Cfg;
}

/// Type alias for operator parameter pairs. Used with opParam() helper
/// to construct key-value pairs for makeOpConfig().
template <typename T>
using OpParamPair = std::pair<std::string, std::decay_t<T>>;

/// Helper function to create operator parameter pairs with perfect
/// forwarding. Usage: opParam("Key", Value) creates a pair suitable for
/// makeOpConfig(). Type decay ensures no reference issues when pairs are
/// passed to Config.
template <typename T>
OpParamPair<T> opParam(std::string Key, ///< [in] parameter name
                       T &&Value        ///< [in] parameter value
) {
   return {std::move(Key), std::forward<T>(Value)};
}

/// AnalysisOperator is the abstract base class for all concrete analysis
/// operators. Derived classes are templated on the Kokkos array type of
/// their primary input field, enabling type-safe factory dispatch. Each
/// operator performs a single, well-defined transformation (spatial
/// reduction, temporal averaging, binary operations, etc.).
///
/// Operators allocate their output data arrays as members and register
/// output Fields in the Field registry during construction. The
/// initialize() method stores mesh/environment pointers after all Fields
/// exist. The compute() method retrieves input data from the Field registry
/// by name and writes results to the operator-owned output arrays.
///
/// Cache validation via isCacheValid() prevents redundant computation when
/// multiple downstream operators share an intermediate result.
class AnalysisOperator {

 public:
   /// Default constructor
   AnalysisOperator();

   /// Constructor with operator type name
   AnalysisOperator(const std::string &OperatorType ///< [in] operator type name
   );

   /// Virtual destructor allows polymorphic deletion of derived classes
   virtual ~AnalysisOperator();

   /// Returns the operator type name (e.g., "SpatialMax", "TimeMean")
   const std::string getOperatorType();

   /// Returns the unique instance name for this operator, derived from the
   /// concatenated upstream field names and operator type. For example,
   /// "Temperature_SpatialMean" for a SpatialMean operator with
   /// Temperature input.
   const std::string getName();

   /// Returns the names of input fields required by this operator. These
   /// names are used during dependency resolution to link operators in the
   /// dependency graph.
   const std::vector<std::string> getInputFieldNames();

   /// Returns the names of output fields produced by this operator. These
   /// are the Field registry names where the operator writes its results.
   const std::vector<std::string> getOutputFieldNames();

   /// Returns true if the operator's output is already valid for the given
   /// timestamp (cache hit). Used by computeRecursive to avoid redundant
   /// computation when multiple downstream operators share this result.
   bool isCacheValid(const TimeInstant &TimeStamp ///< [in] timestamp to check
   );

   /// Initializes the operator after all Fields exist in the registry.
   /// Stores pointers to mesh, environment, and other resources needed
   /// during compute() calls. Called once during Analysis construction
   /// after the dependency graph is built.
   virtual void
   initialize(const MachEnv *Env,      ///< [in] machine environment
              const HorzMesh *Mesh,    ///< [in] horizontal mesh
              const VertCoord *VCoord, ///< [in] vertical coordinate
              Config Options           ///< [in] operator-specific options
   );

   /// Sets the period alarm for temporal reduction operators. The alarm
   /// pointer is used to detect when the accumulation period ends and the
   /// operator should finalize its output (e.g. divide sum by sample count
   /// for TimeMean). Default implementation does nothing; only temporal
   /// operators override.
   virtual void setPeriodAlarm(Alarm *Alarm ///< [in] period alarm ptr
   ) {}

   /// Pure virtual compute method - must be implemented by all derived
   /// classes. Retrieves input field data from the Field registry, performs
   /// the operator's transformation, and writes results to operator-owned
   /// output arrays. Updates LastComputed timestamp and FieldComputed flag
   /// for caching.
   virtual void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                        ) = 0;

 protected:
   /// CF metadata an output Field inherits from the input it is derived
   /// from. Units and standard name are copied and are empty when the input
   /// has none. CellMethods records every reduction applied so far, in
   /// order, as the CF cell_methods attribute (e.g. "area: mean time: mean").
   struct InheritedMetadata {
      std::string Units;       ///< units of the input, empty if unknown
      std::string StdName;     ///< CF standard name, empty if none
      std::string CellMethods; ///< CF cell_methods, empty if none
   };

   /// Collects the CF metadata an output Field inherits from an input Field,
   /// with the cell method of the new reduction (if any) appended to those
   /// the input already records. The units and standard name are copied
   /// unchanged; an operator that changes the units or invalidates the
   /// standard name edits the result before creating its output Field.
   static InheritedMetadata
   inheritMetadata(const std::string &InputName,      ///< [in] input field name
                   const std::string &CellMethod = "" ///< [in] new reduction
   );

   /// Returns the CF cell method for a spatial reduction of an input Field
   /// over all of its owned mesh entities and layers: "area: <Method>" for a
   /// horizontal field and "area: depth: <Method>" for a field with a
   /// vertical dimension. Method is a CF cell method name (mean, minimum,
   /// maximum, standard_deviation, sum).
   static std::string
   spatialCellMethod(const std::string &InputName, ///< [in] input field name
                     const std::string &Method     ///< [in] CF cell method
   );

   /// Returns the value of a string attribute (units, standard_name,
   /// cell_methods, long_name) of a Field, or an empty string if the Field
   /// has no such attribute
   static std::string
   getFieldAttribute(const std::string &FieldName, ///< [in] field name
                     const std::string &AttName    ///< [in] attribute name
   );

   /// Creates an output Field with the given description, valid range and
   /// dimensions, carrying the inherited units, standard name and cell
   /// methods. Empty units or standard name are passed through to
   /// Field::create; an empty cell_methods is not added.
   static std::shared_ptr<Field> createOutputField(
       const std::string &OutputName,  ///< [in] output field name
       const std::string &Description, ///< [in] long name of the output
       const InheritedMetadata &Meta,  ///< [in] inherited CF metadata
       const std::any ValidMin,        ///< [in] min valid value
       const std::any ValidMax,        ///< [in] max valid value
       const int NumDims,              ///< [in] number of dimensions
       const std::vector<std::string> &DimNames ///< [in] dimension names
   );

   const HorzMesh *Mesh;    ///< Horizontal mesh for spatial operations
   const VertCoord *VCoord; ///< Vertical coordinate for vertical ops
   MPI_Comm Comm;           ///< MPI communicator for parallel reductions

   std::string OperatorTypeName;        ///< Operator type (e.g., "SpatialMean")
   std::string InstanceName;            ///< Unique instance name
   std::vector<std::string> InputNames; ///< Required input field names
   std::vector<std::string> OutputNames; ///< Produced output field names

   TimeInstant LastComputed; ///< Timestamp of last compute for caching
   bool FieldComputed;       ///< Flag indicating whether output is valid

}; // end class AnalysisOperator

} // end namespace OMEGA

#endif
