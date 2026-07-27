#ifndef OMEGA_MOC_H
#define OMEGA_MOC_H

//===-- analysis/analysisGroups/MOC.h - MOC AnalysisGroup -------*- C++ -*-===//
//
/// \file
/// \brief Defines the MOC analysis group for computing Meridional Overturning
/// Circulation
///
/// MOC is a bundled AnalysisGroup that computes the Meridional Overturning
/// Circulation (MOC) streamfunction using two methods:
///
/// Latitude-binned MOC (for regions):
/// 1. Converting vertical pseudo-velocity to geometric coordinates
/// 2. Computing vertical flux (VerticalVelocity × AreaCell)
/// 3. Optionally applying a regional mask (ExtractRegion)
/// 4. Accumulating flux into latitude bins
/// 5. Performing horizontal integration (cumulative sum from south to north)
/// 6. Converting to Sverdrups (Sv = 1e6 m³/s)
///
/// Transect-based MOC:
/// 1. Converting PseudoThickness to geometric layer thickness
/// 2. Computing edge transport (LayerThickness × NormalVelocity × DvEdge)
/// 3. Accumulating transport across transect edges
/// 4. Performing vertical integration (cumulative sum from bottom to top)
/// 5. Converting to Sverdrups (Sv = 1e6 m³/s)
///
/// Configuration options:
/// - NumBins: Number of latitude bins (default: 180)
/// - MinLat: Minimum latitude in degrees (default: -90.0)
/// - MaxLat: Maximum latitude in degrees (default: 90.0)
/// - Regions: List of region names for regional MOC (optional)
/// - Transects: List of transect names for transect MOC (optional)
/// - ReductionPeriod: Time-averaging period (e.g., "1month")
/// - SnapshotPeriod: Instantaneous output period (e.g., "1month")
///
/// Example usage in config:
/// \code
/// MOC:
///   NumBins: 180
///   MinLat: -90.0
///   MaxLat: 90.0
///   Regions: [Global, Atlantic, Pacific]
///   Transects: [Drake, Atlantic26N]
///   ReductionPeriod: [1Month]
///   SnapshotPeriod: [1Day]
/// \endcode
///
/// The MOC streamfunction represents the zonally integrated meridional mass
/// transport as a function of latitude and depth. Units are Sverdrups (Sv),
/// where 1 Sv = 10^6 m³/s.
///
//===----------------------------------------------------------------------===//

#include "Analysis.h"
#include "AnalysisGroup.h"
#include "Config.h"
#include "operators/Ops.h"
#include <string>

namespace OMEGA {

/// MOC is a bundled AnalysisGroup that computes the Meridional Overturning
/// Circulation (MOC) streamfunction using two methods. The latitude-binned
/// method converts vertical pseudo-velocity to geometric coordinates, computes
/// vertical flux (velocity × area), accumulates flux into latitude bins, and
/// integrates horizontally (south to north) to produce the streamfunction.
/// The transect-based method converts pseudo layer thickness to geometric
/// coordinates, computes edge transport (LayerThickness × NormalVelocity ×
/// DvEdge), accumulates transport across transect edges, and integrates
/// vertically (bottom to top). Both methods convert results to Sverdrups and
/// support optional temporal averaging.
class MOC : public AnalysisGroup {
 public:
   /// Constructs a MOC analysis group. Reads configuration options (number
   /// of bins, latitude range, regions, output frequency), builds operator
   /// chains for MOC computation, creates IOStreams for output, and
   /// associates operators with streams.
   MOC(const std::string &GroupName, ///< [in] name of this group
       Config &AnalysisGroupOptions, ///< [in] group configuration
       Analysis *AnalysisManager     ///< [in] analysis manager
   );

   /// Default destructor
   ~MOC() = default;

 private:
   /// Builds the complete operator chain string for a single MOC region.
   /// Constructs the chain string encoding all transformation steps:
   /// 1. PseudoToGeometric: Convert vertical pseudo-velocity to geometric
   /// 2. BinaryMultiply: Compute vertical flux (velocity × area)
   /// 3. ExtractRegion: Apply regional mask (optional)
   /// 4. BinnedAccumulator: Accumulate flux into latitude bins
   /// 5. PrefixSum: Horizontal integration (south to north across bins)
   /// 6. ScalarMultiply: Convert to Sverdrups
   ///
   /// The chain string is later parsed by buildTemporalChains() which
   /// appends temporal operators and calls parseChainAndBuildOps().
   ///
   /// Returns the complete chain string (stem) for temporal operator appending.
   std::string buildMOCChain(
       const std::string &RegionName, ///< [in] region name (e.g., "Global")
       const std::string
           &RegionMaskName,      ///< [in] mask field name (empty for global)
       I4 NumBins,               ///< [in] number of latitude bins
       Real MinLat,              ///< [in] minimum latitude
       Real MaxLat,              ///< [in] maximum latitude
       Analysis *AnalysisManager ///< [in] analysis manager
   );

   /// Builds the complete operator chain string for a single transect MOC.
   /// Constructs the chain string encoding all transformation steps:
   /// 1. PseudoToGeometric: Convert PseudoThickness to geometric layer
   /// thickness
   /// 2. BinaryMultiply: Multiply by NormalVelocity (thickness × velocity)
   /// 3. BinaryMultiply: Multiply by DvEdge (transport × edge width)
   /// 4. TransectAccumulator: Accumulate across transect edges
   /// 5. PrefixSum: Vertical integration (bottom to top)
   /// 6. ScalarMultiply: Convert to Sverdrups
   ///
   /// The chain string is later parsed by buildTemporalChains() which
   /// appends temporal operators and calls parseChainAndBuildOps().
   ///
   /// Returns the complete chain string (stem) for temporal operator appending.
   std::string buildTransectMOCChain(
       const std::string &TransectName, ///< [in] transect name
       Analysis *AnalysisManager        ///< [in] analysis manager
   );

   /// Configuration storage for each MOC chain. Vector of Config objects
   /// containing operator-specific parameters (NumBins, Dimension, Reverse,
   /// TransectName). Indexed in the same order as the combined region and
   /// transect chain stems, these configs are passed to buildTemporalChains().
   std::vector<Config> ChainConfigs;

}; // end class MOC

} // end namespace OMEGA

#endif
