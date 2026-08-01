//===-- Test driver for OMEGA Pressure Gradient (PGrad) --------------*-
// C++-*-===/
//
/// \file
/// \brief Test driver for PressureGrad module
//
//===----------------------------------------------------------------------===/

#include "PGrad.h"

#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Eos.h"
#include "Error.h"
#include "Field.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "PGrad.h"
#include "PGradRecon.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertCoord.h"
#include "mpi.h"
#include <limits>
#include <vector>

using namespace OMEGA;

void initPGradTest() {

   Error Err;
   int Err1;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize the Logging system
   initLogging(DefEnv);

   // Read default config if present
   Config("Omega");
   Config::readAll("omega.yml");

   // First step of time stepper initialization needed for IOstream
   TimeStepper::init1();
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   // Initialize the IO system
   IO::init(DefComm);

   // Create the default decomposition (initializes the decomposition)
   Decomp::init();

   // Initialize Field infrastructure
   Field::init(ModelClock);

   // Initialize IOStreams - this does not yet validate the contents
   // of each file, only creates streams from Config
   IOStream::init(ModelClock);

   // Initialize the default halo
   Err1 = Halo::init();
   if (Err1 != 0) {
      LOG_ERROR("PGrad: error initializing default halo");
      Err += Error(ErrorCode::Fail, "PGrad: error initializing default halo");
   }

   // Initialize the default mesh
   HorzMesh::init(ModelClock);

   // Initialize the default vertical coordinate
   VertCoord::init();

   // Initialize the equation of state
   Eos::init();

   // Initialize ocean state
   OceanState::init();

   // Initialize tracers
   Tracers::init();

   CHECK_ERROR_ABORT(Err, "PGrad: error during initialization");
}

// Check the mean-preserving linear reconstruction estimator on a single
// column, without a mesh.
//
// Layer means are sampled from a profile that is exactly linear in pressure,
// so the recovered slope must match the profile's slope to round-off in every
// layer, including the two that use a one-sided difference. The layer
// thicknesses are deliberately non-uniform, which is the case this test exists
// for: a formula that assumed uniform thickness passes the uniform control
// below and fails here. Both are run, so a failure says which of the two it
// is.
//
// Also asserted are the two properties the scheme leans on: that the
// reconstruction reproduces the layer mean at mid-layer pressure, and that its
// deviation integrates to zero over the layer.
int testReconstruction() {

   int Err = 0;

   const Real Eps = std::numeric_limits<Real>::epsilon();
   const Real Tol = 16.0_Real;

   // a profile linear in pressure: Theta(p) = Value0 + Slope0 * p
   const Real Value0 = 12.5_Real;
   const Real Slope0 = -1.5e-6_Real;

   // deliberately non-uniform thicknesses in Pa, then a uniform control
   const std::vector<std::vector<Real>> ThicknessSets = {
       {1.0e5, 3.0e5, 2.0e5, 7.0e5, 1.5e5, 9.0e5, 4.0e5, 2.5e5},
       {3.0e5, 3.0e5, 3.0e5, 3.0e5, 3.0e5, 3.0e5, 3.0e5, 3.0e5}};
   const std::vector<std::string> SetNames = {"non-uniform", "uniform"};

   for (int ISet = 0; ISet < ThicknessSets.size(); ++ISet) {

      const std::vector<Real> &DeltaPress = ThicknessSets[ISet];
      const I4 NLayers                    = DeltaPress.size();
      const I4 KMin                       = 0;
      const I4 KMax                       = NLayers - 1;

      // interface pressures accumulated from the surface, mid-layer pressures
      // as their exact arithmetic midpoint, and layer means sampled from the
      // profile at mid-layer pressure -- which is the exact layer average of a
      // linear profile
      std::vector<Real> PressInterface(NLayers + 1);
      std::vector<Real> PressMid(NLayers);
      std::vector<Real> Value(NLayers);
      PressInterface[0] = 0.0_Real;
      for (int K = 0; K < NLayers; ++K) {
         PressInterface[K + 1] = PressInterface[K] + DeltaPress[K];
         PressMid[K] = 0.5_Real * (PressInterface[K] + PressInterface[K + 1]);
         Value[K]    = Value0 + Slope0 * PressMid[K];
      }

      Real MaxSlopeErr = 0.0_Real;
      Real MaxMeanErr  = 0.0_Real;
      Real MaxIntegral = 0.0_Real;

      for (int K = 0; K < NLayers; ++K) {

         I4 KLo, KHi;
         linearReconStencil(K, KMin, KMax, KLo, KHi);
         const Real Slope = linearReconSlope(Value[KLo], Value[KHi],
                                             PressMid[KLo], PressMid[KHi]);

         // round-off in the difference is set by the size of the layer means,
         // divided by the pressure interval differenced over
         const Real ValueScale =
             std::max(std::abs(Value[KLo]), std::abs(Value[KHi]));
         const Real SlopeBound =
             Eps * ValueScale / std::abs(PressMid[KHi] - PressMid[KLo]);
         MaxSlopeErr =
             std::max(MaxSlopeErr, std::abs(Slope - Slope0) / SlopeBound);

         // the reconstruction reproduces the layer mean at mid-layer pressure
         const Real AtMid =
             linearReconEval(Value[K], Slope, PressMid[K], PressMid[K]);
         MaxMeanErr = std::max(MaxMeanErr,
                               std::abs(AtMid - Value[K]) / (Eps * ValueScale));

         // the deviation integrates to zero over the layer, by two-point
         // Gauss quadrature, which is exact for a linear integrand
         const Real HalfWidth = 0.5_Real * DeltaPress[K];
         const Real Offset    = HalfWidth / std::sqrt(3.0_Real);
         const Real DevLo =
             linearReconDeviation(Slope, PressMid[K], PressMid[K] - Offset);
         const Real DevHi =
             linearReconDeviation(Slope, PressMid[K], PressMid[K] + Offset);
         const Real Integral = HalfWidth * (DevLo + DevHi);
         const Real DevScale = HalfWidth * std::abs(Slope) * HalfWidth;
         if (DevScale > 0.0_Real)
            MaxIntegral =
                std::max(MaxIntegral, std::abs(Integral) / (Eps * DevScale));
      }

      LOG_INFO("PGradTest: reconstruction on {} layers: slope error {} eps, "
               "layer-mean error {} eps, layer integral {} eps",
               SetNames[ISet], MaxSlopeErr, MaxMeanErr, MaxIntegral);

      if (MaxSlopeErr > Tol) {
         LOG_ERROR("PGradTest: reconstruction slope FAIL on {} layers: {} eps",
                   SetNames[ISet], MaxSlopeErr);
         ++Err;
      }
      if (MaxMeanErr > Tol) {
         LOG_ERROR("PGradTest: reconstruction layer mean FAIL on {} layers: "
                   "{} eps",
                   SetNames[ISet], MaxMeanErr);
         ++Err;
      }
      if (MaxIntegral > Tol) {
         LOG_ERROR("PGradTest: reconstruction layer integral FAIL on {} "
                   "layers: {} eps",
                   SetNames[ISet], MaxIntegral);
         ++Err;
      }
   }

   // a column with a single valid layer has nothing to difference against, so
   // a constant is the only mean-preserving reconstruction available
   I4 KLo, KHi;
   linearReconStencil(3, 3, 3, KLo, KHi);
   const Real SingleSlope =
       linearReconSlope(Value0, Value0, 1.0e5_Real, 1.0e5_Real);
   if (KLo == 3 && KHi == 3 && SingleSlope == 0.0_Real) {
      LOG_INFO("PGradTest: reconstruction single-layer column PASS");
   } else {
      LOG_ERROR("PGradTest: reconstruction single-layer column FAIL");
      ++Err;
   }

   if (Err == 0)
      LOG_INFO("PGradTest: reconstruction estimator PASS");

   return Err;

} // end testReconstruction

// Pin the per-column pressure lookup with direct property tests, on a
// fabricated pair of columns and without a mesh.
//
// These tests are not optional, and they are the only thing standing between a
// correct implementation and a silently wrong one. No answer-level check
// anywhere in the test plan can distinguish looking a column's state up by
// pressure from looking it up by layer index: on a profile linear in pressure
// every layer's mean-preserving reconstruction is that same line, so looking
// up the wrong layer costs nothing, and on a curved profile the two rules
// differ by less than a factor of two. An implementation that indexes by layer
// therefore passes every exactness and accuracy gate specified.
//
// Two columns are built sharing a surface and a bottom pressure, with the
// second column's interfaces bulging away from the first in the middle of the
// column so that the two are strongly offset there -- the situation coordinate
// tilt produces. Four properties are asserted:
//
//   - the returned layer's interfaces actually bracket the pressure;
//   - the returned layer differs from the edge layer index by at least two
//     somewhere, so the lookup is demonstrably not an index lookup;
//   - the answer does not depend on the starting hint, so the incremented
//     cursors used in the column scan cannot disagree with a search;
//   - a pressure outside the column clamps to the outermost valid layer, which
//     is the rule where the edge control volume extends past a column's own
//     floor.
int testPressureLookup() {

   int Err = 0;

   const I4 NLayers = 32;
   const I4 KMin    = 0;
   const I4 KMax    = NLayers - 1;

   // Column 0 is uniform in pressure; column 1 has the same surface and
   // bottom pressure but bulges away from it in between, by up to three layer
   // thicknesses. The bulge is small enough that column 1's interfaces stay
   // monotonic.
   const Real DeltaP = 1.25e6_Real;
   const Real Bulge  = 6.0_Real * DeltaP;

   const Real Pi = 3.14159265358979323846_Real;

   HostArray2DReal PressInterface("PressInterface", 2, NLayers + 1);
   for (int K = 0; K <= NLayers; ++K) {
      const Real Uniform   = K * DeltaP;
      PressInterface(0, K) = Uniform;
      PressInterface(1, K) = Uniform + Bulge * std::sin(Pi * K / NLayers);
   }

   // check that column 1 stayed monotonic, or the test is not testing what it
   // means to
   for (int K = 0; K < NLayers; ++K) {
      if (PressInterface(1, K + 1) <= PressInterface(1, K)) {
         LOG_ERROR("PGradTest: pressure lookup setup FAIL: column 1 is not "
                   "monotonic at layer {}",
                   K);
         ++Err;
      }
   }

   I4 MaxOffset     = 0;
   int NBracketFail = 0;
   int NHintFail    = 0;

   for (int K = 0; K < NLayers; ++K) {

      // the edge layer's mid pressure, as the edge average of the two
      // columns' interface pressures
      const Real PressTop =
          0.5_Real * (PressInterface(0, K) + PressInterface(1, K));
      const Real PressBot =
          0.5_Real * (PressInterface(0, K + 1) + PressInterface(1, K + 1));
      const Real Press = 0.5_Real * (PressTop + PressBot);

      for (int ICell = 0; ICell < 2; ++ICell) {

         const I4 KFound =
             findLayerForPress(PressInterface, ICell, KMin, KMax, Press, K);

         // the interfaces of the returned layer must bracket the pressure,
         // unless the search clamped at an end of the column
         const bool Above =
             (KFound == KMin) && (Press < PressInterface(ICell, KMin));
         const bool Below =
             (KFound == KMax) && (Press > PressInterface(ICell, KMax + 1));
         const bool Brackets = Press >= PressInterface(ICell, KFound) &&
                               Press <= PressInterface(ICell, KFound + 1);
         if (!Brackets && !Above && !Below)
            ++NBracketFail;

         const I4 Offset = std::abs(KFound - K);
         if (Offset > MaxOffset)
            MaxOffset = Offset;

         // the answer must not depend on where the search started
         for (I4 KHint = KMin; KHint <= KMax; ++KHint) {
            if (findLayerForPress(PressInterface, ICell, KMin, KMax, Press,
                                  KHint) != KFound)
               ++NHintFail;
         }
      }
   }

   if (NBracketFail == 0) {
      LOG_INFO("PGradTest: pressure lookup bracketing PASS");
   } else {
      LOG_ERROR("PGradTest: pressure lookup bracketing FAIL: {} pressures not "
                "bracketed by the returned layer",
                NBracketFail);
      ++Err;
   }

   // Under tilt the layer containing a pressure is generally not the edge
   // layer index. If this never differed by much, the test configuration
   // would not be exercising the distinction the lookup exists for.
   if (MaxOffset >= 2) {
      LOG_INFO("PGradTest: pressure lookup is not an index lookup PASS: "
               "returned layer differs from the edge layer by up to {}",
               MaxOffset);
   } else {
      LOG_ERROR("PGradTest: pressure lookup is not an index lookup FAIL: "
                "returned layer never differs from the edge layer by more "
                "than {}",
                MaxOffset);
      ++Err;
   }

   if (NHintFail == 0) {
      LOG_INFO("PGradTest: pressure lookup hint independence PASS");
   } else {
      LOG_ERROR("PGradTest: pressure lookup hint independence FAIL: {} "
                "disagreements between starting hints",
                NHintFail);
      ++Err;
   }

   // Pressures outside the column clamp to the outermost valid layer, whose
   // reconstruction is then extrapolated
   const Real AboveTop = PressInterface(0, KMin) - DeltaP;
   const Real BelowBot = PressInterface(0, KMax + 1) + DeltaP;
   const I4 KAbove =
       findLayerForPress(PressInterface, 0, KMin, KMax, AboveTop, KMax);
   const I4 KBelow =
       findLayerForPress(PressInterface, 0, KMin, KMax, BelowBot, KMin);
   if (KAbove == KMin && KBelow == KMax) {
      LOG_INFO("PGradTest: pressure lookup out-of-column clamping PASS");
   } else {
      LOG_ERROR("PGradTest: pressure lookup out-of-column clamping FAIL: got "
                "{} above the column and {} below it",
                KAbove, KBelow);
      ++Err;
   }

   return Err;

} // end testPressureLookup

// Build the two-column test state. Every edge joins cells 0 and 1 and is
// DC apart. Layer interfaces are staggered between the two columns by
// TiltFactor, which tilts them relative to surfaces of constant pressure;
// TiltFactor = 0.5 gives two identical columns and no tilt. Temperature and
// salinity vary linearly with geometric height. Pseudo-thickness, specific
// volume and pressure are iterated to consistency, after which the VertCoord
// pressure and geometric height fields are filled.
void setupTwoColumnState(HorzMesh *Mesh,    ///< [in] Horizontal mesh
                         VertCoord *VCoord, ///< [inout] Vertical coordinate
                         OceanState *State, ///< [inout] Ocean state
                         Eos *EqState,      ///< [inout] Equation of state
                         I4 NVertLayers,    ///< [in] Number of layers
                         Real DC,           ///< [in] Distance between cells
                         Real TiltFactor,   ///< [in] Interface tilt, 0.5 = none
                         Real ZBottom       ///< [in] Depth of the sea floor
) {

   const I4 NCellsAll = Mesh->NCellsAll;
   const I4 NEdgesAll = Mesh->NEdgesAll;

   VCoord->NVertLayers   = NVertLayers;
   VCoord->NVertLayersP1 = NVertLayers + 1;

   Array2DReal SpecVolOld("SpecVolOld", Mesh->NCellsSize, NVertLayers);
   Array1DReal SurfacePressure("SurfacePressure", Mesh->NCellsSize);

   auto &MinLayerCell = VCoord->MinLayerCell;
   auto &MaxLayerCell = VCoord->MaxLayerCell;
   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(int i) {
          MinLayerCell(i) = 0;
          MaxLayerCell(i) = NVertLayers - 1;
       });

   auto &MinLayerEdgeBot = VCoord->MinLayerEdgeBot;
   auto &MaxLayerEdgeTop = VCoord->MaxLayerEdgeTop;
   parallelFor(
       {NEdgesAll}, KOKKOS_LAMBDA(int i) {
          MinLayerEdgeBot(i) = 0;
          MaxLayerEdgeTop(i) = NVertLayers - 1;
       });

   auto &CellsOnEdge = Mesh->CellsOnEdge;
   auto &DcEdge      = Mesh->DcEdge;
   parallelFor(
       {NEdgesAll}, KOKKOS_LAMBDA(int i) {
          CellsOnEdge(i, 0) = 0;
          CellsOnEdge(i, 1) = 1;
          DcEdge(i)         = DC;
       });

   // Fetch reference density from Config
   const Real Density0 = RhoSw;

   const I4 TimeLevel = 0;

   // get state and tracer arrays
   Array2DReal PseudoThick = State->getPseudoThickness(TimeLevel);
   Array2DReal Temp        = Tracers::getByName(TimeLevel, "Temperature");
   Array2DReal Salinity    = Tracers::getByName(TimeLevel, "Salinity");

   // set Z interface and mid-point locations
   const Real DZ         = 2.0_Real * (-ZBottom / NVertLayers);
   auto &BottomGeomDepth = VCoord->BottomGeomDepth;
   auto &GeomZInterface  = VCoord->GeomZInterface;
   auto &GeomZMid        = VCoord->GeomZMid;
   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(int i) {
          GeomZInterface(i, NVertLayers) = ZBottom;
          SurfacePressure(i)             = 0.0_Real;
          BottomGeomDepth(i)             = 0.0_Real;
          for (int k = NVertLayers - 1; k >= 0; --k) {
             Real X = (k + i) % 2;
             Real Dz =
                 (2.0_Real * TiltFactor - 1.0_Real) * X * DZ +
                 (1.0_Real - TiltFactor) * DZ; // staggered pseudo-thickness
             GeomZInterface(i, k) = GeomZInterface(i, k + 1) + Dz;
             PseudoThick(i, k) =
                 GeomZInterface(i, k) - GeomZInterface(i, k + 1);
             GeomZMid(i, k) =
                 0.5_Real * (GeomZInterface(i, k) + GeomZInterface(i, k + 1));
             BottomGeomDepth(i) += Dz;
          }
       });

   // set simple temperature and salinity profiles
   auto &SpecVol = EqState->SpecVol;
   parallelFor(
       {NCellsAll, NVertLayers}, KOKKOS_LAMBDA(int i, int k) {
          Real T0 = 30.0;
          Real TB = 5.0;
          Real S0 = 30.0;
          Real SB = 40.0;

          Real Phi0 = (GeomZMid(i, k) - ZBottom) / (-ZBottom);
          Real PhiB = 1.0_Real - Phi0;

          Temp(i, k)       = T0 * Phi0 + TB * PhiB;
          Salinity(i, k)   = S0 * Phi0 + SB * PhiB;
          SpecVol(i, k)    = 1.0_Real / Density0;
          SpecVolOld(i, k) = SpecVol(i, k);
       });

   // Iterate to converge PseudoThick, SpecVol, PressureMid
   auto &PressureMid = VCoord->PressureMid;
   VCoord->computePressure(PseudoThick, SurfacePressure);
   for (int Iteration = 0; Iteration < 15; ++Iteration) {

      // compute specific volume from EOS
      VCoord->computePressure(PseudoThick, SurfacePressure);
      EqState->computeSpecVol(Temp, Salinity, PressureMid);

      // compute psuedo thickness from specific volume
      parallelFor(
          {NCellsAll, NVertLayers}, KOKKOS_LAMBDA(int i, int k) {
             PseudoThick(i, k) =
                 1.0_Real / (SpecVol(i, k) * Density0) *
                 (GeomZInterface(i, k) - GeomZInterface(i, k + 1));
          });

      // compute difference from previous iteration
      Real MaxValue = 0.0_Real;
      parallelReduce(
          {NCellsAll, NVertLayers},
          KOKKOS_LAMBDA(int i, int k, Real &max) {
             Real Diff = Kokkos::abs(SpecVol(i, k) - SpecVolOld(i, k));
             if (Diff > max)
                max = Diff;
          },
          Kokkos::Max<Real>(MaxValue));

      // check convergence
      if (MaxValue < 1e-12_Real) {
         LOG_INFO("converged: max diff = {}", MaxValue);
         break;
      } else {
         parallelFor(
             {NCellsAll, NVertLayers},
             KOKKOS_LAMBDA(int i, int k) { SpecVolOld(i, k) = SpecVol(i, k); });
      }
   }

   // compute pressure once more with converged PseudoThick
   VCoord->computePressure(PseudoThick, SurfacePressure);

   // compute z levels
   VCoord->computeGeomZHeight(PseudoThick, SpecVol);

} // end setupTwoColumnState

// Assert the identity of design section 3.9 (test section 5.5): with the tidal
// and self-attraction-and-loading potentials zero, PressureGradCentered is
// exactly
//
//    -g / d_e * S_{e,k}
//    S_{e,k} = 1/2 (dZ_k + dZ_{k+1}) + alphaBar / (2 g) (dq_k + dq_{k+1})
//
// where dZ and dq are the cross-edge differences of GeomZInterface and
// PressureInterface and alphaBar is the edge average of SpecVol. S is the
// first-order conversion of a height difference taken at fixed layer index
// into one taken at fixed pressure, so this pins what error the centered
// scheme makes as well as what it computes.
//
// This is a permanent regression test rather than a transitional one. Because
// the two expressions read the mesh, VertCoord and Eos state through
// independently written code, their agreement tests the shared upstream state
// -- edge masks, interface indexing, VertCoord conventions -- and not just the
// pressure gradient arithmetic. Expect agreement to round-off rather than
// bit-for-bit: PressureMid may be formed as zBot + h/2 rather than as the mean
// of the two interface pressures, which is algebraically but not bitwise the
// same.
int testCenteredIdentity(HorzMesh *Mesh,    ///< [in] Horizontal mesh
                         VertCoord *VCoord, ///< [inout] Vertical coordinate
                         OceanState *State, ///< [inout] Ocean state
                         Eos *EqState       ///< [inout] Equation of state
) {

   int Err = 0;

   PressureGrad *DefPGrad = PressureGrad::getDefault();
   if (!DefPGrad || DefPGrad->getType() != PressureGradType::Centered) {
      LOG_ERROR("PGradTest: centered identity needs a Centered default");
      return 1;
   }

   const I4 NVertLayers = 60;
   const Real DC        = 30000.0_Real;
   const Real ZBottom   = -1000.0_Real;
   const Real Eps       = std::numeric_limits<Real>::epsilon();

   // Round-off in the centered functor is set by the size of the hydrostatic
   // terms it forms and cancels, so the tolerance is scaled by those rather
   // than by their differences. Measured at or below half an epsilon; the
   // factor here is headroom for compiler and precision differences.
   const Real Tol = 16.0_Real;

   // A sweep of tilts, with 0.5 the untilted case
   const std::vector<Real> Tilts = {0.5, 0.49, 0.45, 0.4, 0.3, 0.2, 0.1, 0.05};

   Real MaxResid = 0.0_Real;
   Real MinGuard = std::numeric_limits<Real>::max();

   for (Real TiltFactor : Tilts) {

      setupTwoColumnState(Mesh, VCoord, State, EqState, NVertLayers, DC,
                          TiltFactor, ZBottom);

      Array2DReal Tend("TendCentered", Mesh->NEdgesSize, NVertLayers);
      deepCopy(Tend, 0.0_Real);

      Array2DReal PseudoThick = State->getPseudoThickness(0);
      Array2DReal Temp        = Tracers::getByName(0, "Temperature");
      Array2DReal Salinity    = Tracers::getByName(0, "Salinity");

      const auto &PressureMid       = VCoord->PressureMid;
      const auto &PressureInterface = VCoord->PressureInterface;
      const auto &GeomZInterface    = VCoord->GeomZInterface;
      const auto &EdgeMask          = VCoord->EdgeMask;
      const auto &SpecVol           = EqState->SpecVol;
      const auto &CellsOnEdge       = Mesh->CellsOnEdge;
      const auto &DcEdge            = Mesh->DcEdge;

      DefPGrad->computePressureGrad(Tend, PressureMid, PressureInterface,
                                    SpecVol, GeomZInterface, PseudoThick, Temp,
                                    Salinity, EqState);

      // Residual of the identity, in units of eps times the hydrostatic
      // scale, and the same measure for a deliberately broken version of the
      // identity that drops the pressure half of S. The guard is what makes
      // the check non-vacuous: a mis-derivation must miss by far more than
      // the tolerance.
      Real LocResid = 0.0_Real;
      Real LocGuard = std::numeric_limits<Real>::max();
      parallelReduce(
          {Mesh->NEdgesAll, NVertLayers},
          KOKKOS_LAMBDA(int IEdge, int K, Real &MaxR, Real &MinG) {
             const I4 ICell0  = CellsOnEdge(IEdge, 0);
             const I4 ICell1  = CellsOnEdge(IEdge, 1);
             const Real InvDc = 1.0_Real / DcEdge(IEdge);

             const Real Z0K  = GeomZInterface(ICell0, K);
             const Real Z1K  = GeomZInterface(ICell1, K);
             const Real Z0K1 = GeomZInterface(ICell0, K + 1);
             const Real Z1K1 = GeomZInterface(ICell1, K + 1);

             const Real Q0K  = PressureInterface(ICell0, K);
             const Real Q1K  = PressureInterface(ICell1, K);
             const Real Q0K1 = PressureInterface(ICell0, K + 1);
             const Real Q1K1 = PressureInterface(ICell1, K + 1);

             const Real AlphaBar =
                 0.5_Real * (SpecVol(ICell0, K) + SpecVol(ICell1, K));

             const Real ShiftZ = 0.5_Real * ((Z1K - Z0K) + (Z1K1 - Z0K1));
             const Real ShiftQ = AlphaBar / (2.0_Real * Gravity) *
                                 ((Q1K - Q0K) + (Q1K1 - Q0K1));

             const Real Expected =
                 EdgeMask(IEdge, K) * (-Gravity * InvDc) * (ShiftZ + ShiftQ);
             const Real Broken =
                 EdgeMask(IEdge, K) * (-Gravity * InvDc) * ShiftZ;

             // the hydrostatic terms the centered functor forms and cancels
             const Real HydroScale =
                 InvDc * 0.25_Real *
                 (Gravity * (Kokkos::abs(Z0K) + Kokkos::abs(Z1K) +
                             Kokkos::abs(Z0K1) + Kokkos::abs(Z1K1)) +
                  AlphaBar * (Kokkos::abs(Q0K) + Kokkos::abs(Q1K) +
                              Kokkos::abs(Q0K1) + Kokkos::abs(Q1K1)));

             const Real Resid = Kokkos::abs(Tend(IEdge, K) - Expected);
             const Real Guard = Kokkos::abs(Tend(IEdge, K) - Broken);

             // masked edges carry no tendency, so neither expression has
             // anything to say there
             if (EdgeMask(IEdge, K) > 0.0_Real && HydroScale > 0.0_Real) {
                const Real Scale = Eps * HydroScale;
                if (Resid / Scale > MaxR)
                   MaxR = Resid / Scale;
                if (Guard / Scale < MinG)
                   MinG = Guard / Scale;
             }
          },
          Kokkos::Max<Real>(LocResid), Kokkos::Min<Real>(LocGuard));

      LOG_INFO("PGradTest: centered identity at tilt {}: residual {} eps, "
               "dropped-pressure-term guard {} eps",
               TiltFactor, LocResid, LocGuard);

      if (LocResid > MaxResid)
         MaxResid = LocResid;
      // the untilted case has nothing for the guard to detect
      if (TiltFactor != 0.5_Real && LocGuard < MinGuard)
         MinGuard = LocGuard;
   }

   if (MaxResid <= Tol) {
      LOG_INFO("PGradTest: centered identity PASS: max residual {} eps",
               MaxResid);
   } else {
      LOG_ERROR("PGradTest: centered identity FAIL: max residual {} eps "
                "exceeds {} eps",
                MaxResid, Tol);
      ++Err;
   }

   // A mis-derivation that drops the pressure half of S misses by many orders
   // of magnitude; if it does not, the identity is being satisfied trivially.
   //
   // The comparison is against the residual rather than against an absolute
   // count of epsilons. Both quantities are measured in epsilons of the same
   // scale, so their ratio is what the check is really about, and it means
   // the same thing in either precision: an absolute threshold of a thousand
   // epsilons happens to leave a comfortable margin in double precision and
   // only a factor of 1.4 in single, which is too near the edge for a check
   // whose whole purpose is to prove the identity is not vacuous.
   const Real GuardRatio =
       (MaxResid > 0.0_Real) ? MinGuard / MaxResid : 0.0_Real;
   if (GuardRatio > 100.0_Real) {
      LOG_INFO("PGradTest: centered identity guard PASS: dropping the "
               "pressure term misses by {} eps, {} times the residual",
               MinGuard, GuardRatio);
   } else {
      LOG_ERROR("PGradTest: centered identity guard FAIL: dropping the "
                "pressure term changes the answer by only {} eps, {} times "
                "the residual",
                MinGuard, GuardRatio);
      ++Err;
   }

   return Err;

} // end testCenteredIdentity

// Check that the PressureGrad configuration group is parsed as expected and
// that a FiniteVolume instance can be created and dispatched to. The
// FiniteVolume sub-options are read from the same group as PressureGradType
// and default to the Phase 1 values when a key is absent.
int testPGradConfig(const HorzMesh *Mesh,   ///< [in] Horizontal mesh
                    const VertCoord *VCoord ///< [in] Vertical coordinate
) {

   int Err = 0;

   // The default instance is built from the PressureGrad group in the
   // config, which selects the centered scheme
   PressureGrad *DefPGrad = PressureGrad::getDefault();

   if (DefPGrad->getType() == PressureGradType::Centered) {
      LOG_INFO("PGradTest: default PressureGradType parse PASS");
   } else {
      LOG_ERROR("PGradTest: default PressureGradType parse FAIL");
      ++Err;
   }

   // The Phase 1 sub-option values, whether read from the config or taken
   // from the defaults
   if (DefPGrad->getHorzOrder() == 2 &&
       DefPGrad->getVertRecon() == PressureGradVertRecon::Linear &&
       DefPGrad->getQuadraturePoints() >= 1) {
      LOG_INFO("PGradTest: FiniteVolume sub-option parse PASS");
   } else {
      LOG_ERROR("PGradTest: FiniteVolume sub-option parse FAIL: HorzOrder={} "
                "QuadraturePoints={}",
                DefPGrad->getHorzOrder(), DefPGrad->getQuadraturePoints());
      ++Err;
   }

   // Create a second instance that selects the FiniteVolume scheme. The
   // sub-config shares its node with the parent, so resetting the type here
   // is what the new instance sees. Note the explicit std::string: a string
   // literal would select the boolean overload of Config::set.
   Config *Options = Config::getOmegaConfig();
   Config PGradConfig("PressureGrad");
   Err += (Options->get(PGradConfig).isFail() ? 1 : 0);
   PGradConfig.set("PressureGradType", std::string("FiniteVolume"));

   PressureGrad *FVPGrad =
       PressureGrad::create("TestFiniteVolume", Mesh, VCoord, Options);

   if (FVPGrad && FVPGrad->getType() == PressureGradType::FiniteVolume) {
      LOG_INFO("PGradTest: FiniteVolume PressureGradType parse PASS");
   } else {
      LOG_ERROR("PGradTest: FiniteVolume PressureGradType parse FAIL");
      ++Err;
   }

   // Dispatch check: computePressureGrad must take the FiniteVolume branch
   // and leave a finite tendency everywhere it touches
   OceanState *State       = OceanState::getDefault();
   Eos *EqState            = Eos::getInstance();
   Array2DReal PseudoThick = State->getPseudoThickness(0);
   Array2DReal TendFV("TendFiniteVolume", Mesh->NEdgesSize,
                      VCoord->NVertLayers);
   deepCopy(TendFV, 0.0_Real);

   Array2DReal ConservTemp = Tracers::getByName(0, "Temperature");
   Array2DReal AbsSalinity = Tracers::getByName(0, "Salinity");

   FVPGrad->computePressureGrad(
       TendFV, VCoord->PressureMid, VCoord->PressureInterface, EqState->SpecVol,
       VCoord->GeomZInterface, PseudoThick, ConservTemp, AbsSalinity, EqState);

   I4 NBad = 0;
   parallelReduce(
       {Mesh->NEdgesAll, VCoord->NVertLayers},
       KOKKOS_LAMBDA(int IEdge, int K, I4 &LSum) {
          if (!Kokkos::isfinite(TendFV(IEdge, K)))
             ++LSum;
       },
       Kokkos::Sum<I4>(NBad));

   if (NBad == 0) {
      LOG_INFO("PGradTest: FiniteVolume dispatch PASS");
   } else {
      LOG_ERROR("PGradTest: FiniteVolume dispatch FAIL: {} non-finite values",
                NBad);
      ++Err;
   }

   PGradConfig.set("PressureGradType", std::string("Centered"));
   PressureGrad::erase("TestFiniteVolume");

   return Err;

} // end testPGradConfig

int main(int argc, char *argv[]) {
   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   {
      initPGradTest();
      // Initialize default PressureGrad
      PressureGrad::init();

      HorzMesh *DefMesh    = HorzMesh::getDefault();
      VertCoord *VCoord    = VertCoord::getDefault();
      OceanState *DefState = OceanState::getDefault();
      Eos *DefEos          = Eos::getInstance();

      // create arrays: Tend on edges, Pressure/Geopotential/SpecVol on cells
      Array2DReal Tend("Tend", DefMesh->NEdgesSize, VCoord->NVertLayers);

      I4 NEdgesAll = DefMesh->NEdgesAll;

      I4 NVertLayers  = 60;
      Real DC         = 30000.0_Real;
      Real TiltFactor = 0.495_Real;
      Real ZBottom    = -1000.0_Real;
      I4 NRefinements = 4;
      HostArray1DReal Rmse("Rmse", NRefinements);
      for (int Refinement = 0; Refinement < NRefinements; ++Refinement) {

         LOG_INFO("PGradTest: Starting refinement level {}", Refinement);

         setupTwoColumnState(DefMesh, VCoord, DefState, DefEos, NVertLayers, DC,
                             TiltFactor, ZBottom);

         Array2DReal PseudoThick = DefState->getPseudoThickness(0);
         Array2DReal Temp        = Tracers::getByName(0, "Temperature");
         Array2DReal Salinity    = Tracers::getByName(0, "Salinity");
         auto &SpecVol           = DefEos->SpecVol;
         auto &PressureMid       = VCoord->PressureMid;
         auto &GeomZInterface    = VCoord->GeomZInterface;

         LOG_INFO("NVertLayers = {}", NVertLayers);
         LOG_INFO("dC = {}", DC);
         DefState->copyToHost(0);
         HostArray2DReal PseudoThickH = DefState->getPseudoThicknessH(0);
         for (int i = 0; i < 2; ++i) {
            for (int k = 0; k < 2; ++k) {
               LOG_INFO("PseudoThick({}, {}) = {}", i, k, PseudoThickH(i, k));
            }
         }

         // get PressureGrad instance
         PressureGrad *DefPGrad = PressureGrad::getDefault();
         if (!DefPGrad) {
            LOG_INFO("PGrad: default instance not created by init");
         }

         // compute pressure gradient
         deepCopy(Tend, 0.0_Real);

         const auto &PressureInterface = VCoord->PressureInterface;
         DefPGrad->computePressureGrad(Tend, PressureMid, PressureInterface,
                                       SpecVol, GeomZInterface, PseudoThick,
                                       Temp, Salinity, DefEos);

         // compute errors
         Real MaxValue = 0.0_Real;
         parallelReduce(
             {NEdgesAll, NVertLayers - 2},
             KOKKOS_LAMBDA(int i, int k, Real &max) {
                Real Val = Kokkos::abs(Tend(i, k + 1));
                if (Val > max)
                   max = Val;
             },
             Kokkos::Max<Real>(MaxValue));
         Real SumValue = 0.0_Real;
         parallelReduce(
             {NEdgesAll, NVertLayers - 2},
             KOKKOS_LAMBDA(int i, int k, Real &LSum) {
                LSum += Tend(i, k + 1) * Tend(i, k + 1);
             },
             Kokkos::Sum<Real>(SumValue));
         Real RmseVal = std::sqrt(SumValue / (NEdgesAll * (NVertLayers - 2)));
         Rmse(Refinement) = RmseVal;

         LOG_INFO("refinement level {}: max |Tend| = {}, average Tend = {}",
                  Refinement, MaxValue, RmseVal);

         // coarsen for next iteration
         DC          = DC * 2.0_Real;
         NVertLayers = NVertLayers / 2;

      } // refinement loop

      // Test the centered identity of design section 3.9 across a sweep of
      // tilts. This resets the two-column state, so it runs after the
      // refinement loop.
      // Test the reconstruction estimator on its own. It needs no mesh, and
      // it is the most likely place for the exactness gate to fail, so a
      // failure here localizes immediately.
      int ReconErr = testReconstruction();

      // Pin the per-column pressure lookup. No answer-level check anywhere
      // can distinguish this from a layer-index lookup, so these property
      // tests are the only protection against that failure mode.
      int LookupErr = testPressureLookup();

      int IdentityErr = testCenteredIdentity(DefMesh, VCoord, DefState, DefEos);

      // Test parsing and dispatch of the PressureGrad configuration options
      int ConfigErr = testPGradConfig(DefMesh, VCoord);

      // Test for second order convergence
      // resolution (dC) increases in refimenent loop
      if (Rmse(0) < Rmse(NRefinements - 1) / pow(4.0_Real, NRefinements - 1)) {
         RetVal = 0;
      } else {
         RetVal = 1;
      }

      RetVal += ReconErr + LookupErr + IdentityErr + ConfigErr;

      // cleanup
      PressureGrad::clear();
      IOStream::finalize();
      TimeStepper::clear();
      Tracers::clear();
      Eos::destroyInstance();
      OceanState::clear();
      VertCoord::clear();
      Field::clear();
      Dimension::clear();
      HorzMesh::clear();
      Halo::clear();
      Decomp::clear();
      MachEnv::removeAll();
   }
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;
   return RetVal;
}
