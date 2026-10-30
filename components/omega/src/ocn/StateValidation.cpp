//===-- ocn/StateValidation.cpp - ocean state validation --------*- C++ -*-===//
//
// Validates ocean state fields by checking for NaN values and
// out-of-bounds conditions. Abort behaviour is controlled at runtime via the
// YAML options Omega::State::Validation::AbortOnNan (default true) and
// Omega::State::Validation::AbortOnOOB (default false). The decision is made
// from the counts on all tasks. When an abort flag is set, every task writes
// the stream named by Omega::State::Validation::WriteStreamOnAbort (default
// History), a critical error log with backtrace is emitted and MPI_Abort is
// called; otherwise a CRITICAL log message is written and execution continues.
//
//===----------------------------------------------------------------------===//

#include "StateValidation.h"

#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Eos.h"
#include "Error.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "Reductions.h"
#include "Tracers.h"
#include "VertCoord.h"
#include "mpi.h"

#include <climits>
#include <cmath>
#include <cpptrace/cpptrace.hpp>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace OMEGA {

struct ValidationAbortConfig {
   bool AbortOnOOB                = false;
   bool AbortOnNan                = true;
   std::string WriteStreamOnAbort = "History";
};

static ValidationAbortConfig getValidationAbortConfig() {
   ValidationAbortConfig AbortConfig;

   Config *OmegaConfig = Config::getOmegaConfig();
   if (!OmegaConfig) {
      return AbortConfig;
   }

   Config StateConfig("State");
   if (!OmegaConfig->get(StateConfig).isSuccess()) {
      return AbortConfig;
   }

   Config ValidationConfig("Validation");
   if (!StateConfig.get(ValidationConfig).isSuccess()) {
      return AbortConfig;
   }

   ValidationConfig.get("AbortOnOOB", AbortConfig.AbortOnOOB);
   ValidationConfig.get("AbortOnNan", AbortConfig.AbortOnNan);
   ValidationConfig.get("WriteStreamOnAbort", AbortConfig.WriteStreamOnAbort);

   return AbortConfig;
}

//------------------------------------------------------------------------------
// Helper: abort the run, first writing the given stream if there is one and a
// clock to stamp it with. It must be called on every task, since writing a
// stream is collective. The tasks that found problems have already logged
// them; the master task adds the totals and a backtrace.
static void abortRun(I8 NaNs, I8 OOBs, const std::string &StreamName,
                     const Clock *ModelClock) {

   MachEnv *DefEnv = MachEnv::getDefault();
   MPI_Comm Comm   = DefEnv->getComm();
   bool IsMaster   = DefEnv->isMasterTask();

   if (IsMaster) {
      LOG_CRITICAL("StateValidation: Ocean state validation aborting with {} "
                   "NaNs and {} OOBs on all tasks. The tasks that found them "
                   "give the field and location in critical messages, on "
                   "stderr for tasks that do not log.",
                   NaNs, OOBs);
   }

   if (!StreamName.empty() && ModelClock) {
      if (IOStream::isDefined(StreamName)) {
         if (IsMaster) {
            LOG_CRITICAL("StateValidation: Writing stream {} before aborting, "
                         "with the time at the start of the failing step",
                         StreamName);
         }
         IOStream::write(StreamName, ModelClock, /*ForceWrite=*/true);
      } else if (IsMaster) {
         LOG_CRITICAL("StateValidation: Not writing stream {} before "
                      "aborting because it is not defined",
                      StreamName);
      }
   }

   if (IsMaster) {
      cpptrace::generate_trace().print();
   }
   MPI_Barrier(Comm);
   MPI_Abort(Comm, static_cast<int>(ErrorCode::Critical));
}

//------------------------------------------------------------------------------
// Helper: count the entries over the first NRows rows and NCols columns for
// which IsBad(Row, Col) is true
template <class F>
static I4 countEntries(const std::string &Label, I4 NRows, I4 NCols,
                       const F &IsBad) {
   I4 Count = 0;
   parallelReduce(
       Label, {NRows, NCols},
       KOKKOS_LAMBDA(int Row, int Col, int &Accum) {
          if (IsBad(Row, Col)) {
             ++Accum;
          }
       },
       Count);
   return Count;
}

//------------------------------------------------------------------------------
// Helper: the smallest flat index Row * NCols + Col for which IsBad(Row, Col)
// is true, or -1 if there is none
template <class F>
static I4 firstEntry(const std::string &Label, I4 NRows, I4 NCols,
                     const F &IsBad) {
   I4 First = INT_MAX;
   parallelReduce(
       Label, {NRows, NCols},
       KOKKOS_LAMBDA(int Row, int Col, I4 &Accum) {
          if (IsBad(Row, Col)) {
             Accum = Kokkos::min(Accum, Row * NCols + Col);
          }
       },
       Kokkos::Min<I4>(First));
   return First == INT_MAX ? -1 : First;
}

//------------------------------------------------------------------------------
// NaN and out-of-range counts for one field, with the flat index of the first
// of each (-1 if there is none)
struct FieldCheck {
   I4 NaNs     = 0;
   I4 OOBs     = 0;
   I4 FirstNaN = -1;
   I4 FirstOOB = -1;
};

//------------------------------------------------------------------------------
// Helper: count the NaN and out-of-range entries of a field and, only when
// there are some, find the first of each, so a valid field costs only the
// two counts
template <class FNaN, class FOOB>
static FieldCheck checkEntries(I4 NRows, I4 NCols, const FNaN &IsNaN,
                               const FOOB &IsOOB) {
   FieldCheck Check;
   Check.NaNs = countEntries("CheckNaN", NRows, NCols, IsNaN);
   Check.OOBs = countEntries("CheckBounds", NRows, NCols, IsOOB);
   if (Check.NaNs > 0) {
      Check.FirstNaN = firstEntry("FirstNaN", NRows, NCols, IsNaN);
   }
   if (Check.OOBs > 0) {
      Check.FirstOOB = firstEntry("FirstOOB", NRows, NCols, IsOOB);
   }
   return Check;
}

//------------------------------------------------------------------------------
// Helper: check a 2-D Real device array over the first NCells/NEdges rows and
// NVert columns for NaN and out-of-range entries, restricted to active cells
// (CellMask > 0).
static FieldCheck checkArray2D(const Array2DReal &Arr, I4 NRows, I4 NCols,
                               Real MinVal, Real MaxVal, bool CheckMin,
                               const Array2DReal &CellMask) {
   auto IsNaN = KOKKOS_LAMBDA(int Row, int Col) {
      return CellMask(Row, Col) != 0 && Kokkos::isnan(Arr(Row, Col));
   };
   // NaNs are excluded explicitly, since comparisons with NaN are false
   auto IsOOB = KOKKOS_LAMBDA(int Row, int Col) {
      if (CellMask(Row, Col) == 0) {
         return false;
      }
      Real Val = Arr(Row, Col);
      return !Kokkos::isnan(Val) &&
             (Val > MaxVal || (CheckMin && Val < MinVal));
   };
   return checkEntries(NRows, NCols, IsNaN, IsOOB);
}

//------------------------------------------------------------------------------
// Helper: check a single tracer (row = cell, col = vert) of the 3-D tracer
// array at the given tracer index for NaN and out-of-range entries,
// restricted to active cells (CellMask > 0).
static FieldCheck checkTracerArray(const Array3DReal &Tracers3D, I4 TracerIdx,
                                   I4 NCells, I4 NVert, Real MinVal,
                                   Real MaxVal, const Array2DReal &CellMask) {
   auto IsNaN = KOKKOS_LAMBDA(int Cell, int K) {
      return CellMask(Cell, K) != 0 &&
             Kokkos::isnan(Tracers3D(TracerIdx, Cell, K));
   };
   auto IsOOB = KOKKOS_LAMBDA(int Cell, int K) {
      if (CellMask(Cell, K) == 0) {
         return false;
      }
      Real Val = Tracers3D(TracerIdx, Cell, K);
      return !Kokkos::isnan(Val) && (Val < MinVal || Val > MaxVal);
   };
   return checkEntries(NCells, NVert, IsNaN, IsOOB);
}

//------------------------------------------------------------------------------
// Helpers: copy a single entry of a device array to the host
static Real getEntry(const Array2DReal &Arr, I4 Row, I4 Col) {
   Real Val;
   Kokkos::deep_copy(Val, Kokkos::subview(Arr, Row, Col));
   return Val;
}

static Real getEntry(const Array3DReal &Arr, I4 Idx, I4 Row, I4 Col) {
   Real Val;
   Kokkos::deep_copy(Val, Kokkos::subview(Arr, Idx, Row, Col));
   return Val;
}

//------------------------------------------------------------------------------
// Helper: describe where a cell-layer entry is, for critical messages
static std::string describeLocation(I4 Cell, I4 K) {
   std::string Where = "local cell " + std::to_string(Cell);

   Decomp *DefDecomp = Decomp::getDefault();
   if (DefDecomp) {
      Where = "cell ID " + std::to_string(DefDecomp->CellIDH(Cell));
   }

   HorzMesh *DefMesh = HorzMesh::getDefault();
   if (DefMesh) {
      const Real RadToDeg = 180.0 / Pi;
      Where += fmt::format(" (lat {:.4f}, lon {:.4f})",
                           DefMesh->LatCellH(Cell) * RadToDeg,
                           DefMesh->LonCellH(Cell) * RadToDeg);
   }

   return Where + ", layer " + std::to_string(K);
}

//------------------------------------------------------------------------------
// Helper: log a critical message for each problem found in one field, giving
// where the first one is, record them in Issues if it is given, and add them
// to the totals
static void reportField(const std::string &Field, const std::string &Range,
                        const FieldCheck &Check, I4 NVert,
                        const std::function<Real(I4, I4)> &GetValue,
                        I4 &TotalNaNs, I4 &TotalOOBs,
                        std::vector<ValidationIssue> *Issues) {

   if (Check.NaNs > 0) {
      I4 Cell = Check.FirstNaN / NVert;
      I4 K    = Check.FirstNaN % NVert;
      LOG_CRITICAL("StateValidation: {} contains {} NaN value(s); the first "
                   "is at {}",
                   Field, Check.NaNs, describeLocation(Cell, K));
      if (Issues) {
         Issues->push_back(
             {Field, true, Check.NaNs, Cell, K, GetValue(Cell, K)});
      }
      TotalNaNs += Check.NaNs;
   }
   if (Check.OOBs > 0) {
      I4 Cell  = Check.FirstOOB / NVert;
      I4 K     = Check.FirstOOB % NVert;
      Real Val = GetValue(Cell, K);
      LOG_CRITICAL("StateValidation: {} has {} value(s) outside valid range "
                   "{}; the first is {} at {}",
                   Field, Check.OOBs, Range, Val, describeLocation(Cell, K));
      if (Issues) {
         Issues->push_back({Field, false, Check.OOBs, Cell, K, Val});
      }
      TotalOOBs += Check.OOBs;
   }
}

//------------------------------------------------------------------------------
/// Check ocean state fields for NaN and out-of-bounds conditions.
/// Only active cells (where CellMask > 0) are checked.
/// Returns the total count of errors found; does not abort.
std::pair<I4, I4> checkOceanState(const OceanState *State,
                                  const AuxiliaryState *AuxState,
                                  const VertCoord *VCoord, I4 TimeLevel,
                                  std::vector<ValidationIssue> *Issues) {

   I4 TotalNaNs = 0;
   I4 TotalOOBs = 0;

   const Array2DReal &CellMask = VCoord->CellMask;
   const I4 NCells             = State->NCellsOwned;
   const I4 NVert              = State->NVertLayers;

   // -------------------------------------------------------------------------
   // PseudoThickness: valid range [1e-10, 1000]
   // -------------------------------------------------------------------------
   {
      Array2DReal PseudoThick = State->getPseudoThickness(TimeLevel);
      FieldCheck Check =
          checkArray2D(PseudoThick, NCells, NVert, static_cast<Real>(1e-10),
                       static_cast<Real>(1000.0), /*CheckMin=*/true, CellMask);
      reportField(
          "PseudoThickness", "[1e-10, 1000]", Check, NVert,
          [&](I4 Cell, I4 K) { return getEntry(PseudoThick, Cell, K); },
          TotalNaNs, TotalOOBs, Issues);
   }

   // -------------------------------------------------------------------------
   // KineticEnergyCell: valid range [0, 10]
   // -------------------------------------------------------------------------
   {
      const Array2DReal &KE = AuxState->KineticAux.KineticEnergyCell;
      FieldCheck Check =
          checkArray2D(KE, NCells, NVert, static_cast<Real>(0.0),
                       static_cast<Real>(10.0), /*CheckMin=*/true, CellMask);
      reportField(
          "KineticEnergyCell", "[0, 10]", Check, NVert,
          [&](I4 Cell, I4 K) { return getEntry(KE, Cell, K); }, TotalNaNs,
          TotalOOBs, Issues);
   }

   // -------------------------------------------------------------------------
   // Temperature tracer: valid range [-10, 50]
   // -------------------------------------------------------------------------
   if (Tracers::IndxTemp != -1) {
      Array3DReal AllTracers = Tracers::getAll(TimeLevel);
      const I4 Indx          = Tracers::IndxTemp;
      FieldCheck Check       = checkTracerArray(AllTracers, Indx, NCells, NVert,
                                                static_cast<Real>(-10.0),
                                                static_cast<Real>(50.0), CellMask);
      reportField(
          "Temperature", "[-10, 50]", Check, NVert,
          [&](I4 Cell, I4 K) { return getEntry(AllTracers, Indx, Cell, K); },
          TotalNaNs, TotalOOBs, Issues);
   }

   // -------------------------------------------------------------------------
   // Salinity tracer: valid range [-2, 60]
   // -------------------------------------------------------------------------
   if (Tracers::IndxSalt != -1) {
      Array3DReal AllTracers = Tracers::getAll(TimeLevel);
      const I4 Indx          = Tracers::IndxSalt;
      FieldCheck Check       = checkTracerArray(AllTracers, Indx, NCells, NVert,
                                                static_cast<Real>(-2.0),
                                                static_cast<Real>(60.0), CellMask);
      reportField(
          "Salinity", "[-2, 60]", Check, NVert,
          [&](I4 Cell, I4 K) { return getEntry(AllTracers, Indx, Cell, K); },
          TotalNaNs, TotalOOBs, Issues);
   }

   return {TotalNaNs, TotalOOBs};
}

//------------------------------------------------------------------------------
/// Validate ocean state fields for NaN and out-of-bounds conditions.
/// Only active cells (where CellMask > 0) are checked.
/// Aborts via MPI_Abort on failure.
void validateOceanState(const OceanState *State, const AuxiliaryState *AuxState,
                        const VertCoord *VCoord, I4 TimeLevel,
                        const Clock *ModelClock) {

   auto [LocalNaNs, LocalOOBs] =
       checkOceanState(State, AuxState, VCoord, TimeLevel);
   static const auto AbortConfig = getValidationAbortConfig();

   MachEnv *DefEnv  = MachEnv::getDefault();
   Eos *EosInstance = Eos::getInstance();

   // One reduction gathers everything validation needs from all tasks, so
   // that every task makes the same decision about aborting
   I8 LocalNNegSalinity = EosInstance ? EosInstance->NNegSalinity : 0;
   std::vector<I8> Totals =
       globalSum(std::vector<I8>{LocalNaNs, LocalOOBs, LocalNNegSalinity},
                 DefEnv->getComm());
   const I8 NaNs         = Totals[0];
   const I8 OOBs         = Totals[1];
   const I8 NNegSalinity = Totals[2];

   bool abort = false;

   if (NaNs > 0 && AbortConfig.AbortOnNan) {
      abort = true;
   }
   if (OOBs > 0 && AbortConfig.AbortOnOOB) {
      abort = true;
   }
   if (abort) {
      abortRun(NaNs, OOBs, AbortConfig.WriteStreamOnAbort, ModelClock);
   } else if ((NaNs + OOBs) > 0 && DefEnv->isMasterTask()) {
      // Only report the flags that are actually suppressing an abort, i.e.
      // those corresponding to the issues that were detected
      std::string Flags;
      if (NaNs > 0) {
         Flags = "AbortOnNan is false";
      }
      if (OOBs > 0) {
         if (!Flags.empty())
            Flags += " and ";
         Flags += "AbortOnOOB is false";
      }
      LOG_CRITICAL("StateValidation: Ocean state validation failed with {} "
                   "NaNs and {} OOBs on all tasks. {}; continuing after "
                   "logging.",
                   NaNs, OOBs, Flags);
   }

   if (EosInstance) {
      if (NNegSalinity > 0) {
         LOG_WARN("StateValidation: {} negative salinity value(s) seen by "
                  "the EOS during this time step (summed over all stages "
                  "and calls)",
                  NNegSalinity);
      }
      EosInstance->NNegSalinity = 0;
   }
}

} // namespace OMEGA

//===----------------------------------------------------------------------===//
