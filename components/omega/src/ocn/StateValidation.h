#ifndef OMEGA_STATEVALIDATION_H
#define OMEGA_STATEVALIDATION_H
//===-- ocn/StateValidation.h - ocean state validation ----------*- C++ -*-===//
//
/// \file
/// \brief Declares state validation functions for ocean state validation
///
/// Provides two functions:
///   - checkOceanState: checks for NaN and out-of-bounds conditions and
///     returns the NaN count and out-of-bounds count as a pair without
///     aborting. Suitable for testing.
///   - validateOceanState: calls checkOceanState and conditionally aborts
///     via MPI_Abort based on runtime YAML options AbortOnNan and AbortOnOOB
///     (under Omega::State::Validation), first writing the stream named by
///     WriteStreamOnAbort. When both flags are false a CRITICAL log message
///     is emitted and the run continues.
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "OceanState.h"
#include "TimeMgr.h"
#include "Tracers.h"
#include "VertCoord.h"

#include <string>
#include <vector>

namespace OMEGA {

/// One kind of problem (NaN or out of range) found in one field on this task,
/// with the local cell, layer and value of the first entry that has it
struct ValidationIssue {
   std::string Field; ///< Name of the field, e.g. Temperature
   bool IsNaN;        ///< NaN (true) or out of the valid range (false)
   I4 Count;          ///< Number of entries on this task with the problem
   I4 Cell;           ///< Local index of the cell of the first entry
   I4 Layer;          ///< Vertical layer of the first entry
   Real Value;        ///< Value of the first entry
};

/// Check ocean state fields for NaN values and out-of-bounds conditions
/// without aborting, returning counts of each error type found.
///
/// Only active ocean cells (where CellMask > 0) are checked. A critical log
/// message is emitted for each type of error found, giving the number of
/// entries and the global cell ID, latitude, longitude and layer of the
/// first. If Issues is given, the same problems are appended to it.
///
/// The following fields are validated:
///   - PseudoThickness     : [1e-10, 1000]  (from OceanState)
///   - KineticEnergyCell   : [0, 10]
///                           (from AuxiliaryState::KineticAux)
///   - Temperature tracer  : [-10, 50]      (from Tracers)
///   - Salinity tracer     : [-2, 60]       (from Tracers)
///
/// \param[in] State       Ocean state to validate
/// \param[in] AuxState    Auxiliary state containing KineticEnergyCell
/// \param[in] VCoord      Vertical coordinate containing the CellMask
/// \param[in] TimeLevel   Time level index to validate (typically 0 = current)
/// \param[out] Issues     Optional list to which the problems found are added
/// \return std::pair<I4,I4> {NaN count, out-of-bounds count}; both 0 if valid
std::pair<I4, I4>
checkOceanState(const OceanState *State, const AuxiliaryState *AuxState,
                const VertCoord *VCoord, I4 TimeLevel,
                std::vector<ValidationIssue> *Issues = nullptr);

/// Check ocean state fields for NaN values and out-of-bounds conditions.
///
/// Calls checkOceanState and conditionally aborts via MPI_Abort based on the
/// runtime YAML options AbortOnNan and AbortOnOOB (under
/// Omega::State::Validation) and the counts summed over all tasks. It must
/// therefore be called on every task. Only active ocean cells (where
/// CellMask > 0) are checked.
///
/// The following fields are validated:
///   - PseudoThickness     : [1e-10, 1000]  (from OceanState)
///   - KineticEnergyCell   : [0, 10]
///                           (from AuxiliaryState::KineticAux)
///   - Temperature tracer  : [-10, 50]      (from Tracers)
///   - Salinity tracer     : [-2, 60]       (from Tracers)
///
/// If AbortOnNan is true and NaN values are detected, or AbortOnOOB is true
/// and out-of-bounds values are detected, the stream named by the YAML option
/// WriteStreamOnAbort (default History) is written on every task if it is
/// defined and ModelClock is given, a critical error is logged with a stack
/// backtrace and the run is aborted via MPI_Abort. When both flags are false
/// but errors are found, a CRITICAL log message is emitted and execution
/// continues.
///
/// \param[in] State       Ocean state to validate
/// \param[in] AuxState    Auxiliary state containing KineticEnergyCell
/// \param[in] VCoord      Vertical coordinate containing the CellMask
/// \param[in] TimeLevel   Time level index to validate (typically 0 = current)
/// \param[in] ModelClock  Optional clock for the time stamp of the stream
///                        written before an abort; without it, no stream is
///                        written
void validateOceanState(const OceanState *State, const AuxiliaryState *AuxState,
                        const VertCoord *VCoord, I4 TimeLevel,
                        const Clock *ModelClock = nullptr);

} // namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_STATEVALIDATION_H
