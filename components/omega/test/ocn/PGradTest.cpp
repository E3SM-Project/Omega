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
#include "PGradFiniteVolume.h"
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

// The exact layer average over [PressLo, PressHi] of a profile that is
// quadratic in pressure, used to initialize layer means the way the exactness
// gate requires: as the exact layer averages of a prescribed continuous
// profile, which under tilt come out different in the two columns.
KOKKOS_INLINE_FUNCTION Real exactLayerAvg(
    const Real Coeff0,  ///< [in] constant term
    const Real Coeff1,  ///< [in] coefficient of p
    const Real Coeff2,  ///< [in] coefficient of p squared
    const Real PressLo, ///< [in] top interface pressure
    const Real PressHi  ///< [in] bottom interface pressure
) {
   const Real PressMid   = 0.5_Real * (PressLo + PressHi);
   const Real DeltaPress = PressHi - PressLo;
   return Coeff0 + Coeff1 * PressMid +
          Coeff2 * (PressMid * PressMid + DeltaPress * DeltaPress / 12.0_Real);
}

// Check the matched-pressure integrand on a fabricated column pair, without a
// mesh. This exercises the reconstruction, the pressure lookup and the
// edge-shared expansion together, which is the combination the whole scheme
// rests on.
//
// On a profile linear in pressure the integrand must be zero at *every*
// quadrature point, not merely in the integral. That is the sharpest available
// form of the check: it localizes a failure to the layer that caused it, and
// it is what a per-layer-index sharing rule fails. Both columns' layer means
// are the exact layer averages of one prescribed continuous profile, so under
// the offset between the two columns those means genuinely differ -- a
// configuration built by copying identical means into both columns would not
// exercise the property at all.
//
// A quadratic profile runs as a control. The reconstruction does not resolve
// it, so the integrand is nonzero there, and the gap between the two cases is
// what makes the linear result meaningful rather than a statement that
// everything in the test is zero.
int testMatchedPressIntegrand() {

   int Err = 0;

   const Real Eps = std::numeric_limits<Real>::epsilon();
   const Real Tol = 16.0_Real;

   const I4 NLayers = 32;
   const I4 KMin    = 0;
   const I4 KMax    = NLayers - 1;

   const Real DeltaP = 1.25e6_Real;
   const Real Bulge  = 6.0_Real * DeltaP;
   const Real Pi     = 3.14159265358979323846_Real;

   // representative shared edge coefficients; their values cannot affect the
   // exact-set result, since they multiply a quantity that is identically zero
   PGradEdgeEos EdgeEos;
   EdgeEos.SpecVol0    = 9.7e-4_Real;
   EdgeEos.SpecVolDCt  = -2.5e-7_Real;
   EdgeEos.SpecVolDSa  = -7.5e-7_Real;
   EdgeEos.SpecVolDP   = -4.0e-13_Real;
   EdgeEos.ConservTemp = 8.0_Real;
   EdgeEos.AbsSalinity = 34.9_Real;
   EdgeEos.Press       = 2.0e7_Real;

   // the prescribed continuous profile, linear then quadratic in pressure
   const Real Temp0 = 12.0_Real, Temp1 = -2.0e-7_Real;
   const Real Salt0 = 34.5_Real, Salt1 = 2.5e-8_Real;

   Real MaxLinear    = 0.0_Real;
   Real MaxQuadratic = 0.0_Real;

   for (int ICase = 0; ICase < 2; ++ICase) {

      const bool Curved = (ICase == 1);
      const Real Temp2  = Curved ? 3.0e-15_Real : 0.0_Real;
      const Real Salt2  = Curved ? -4.0e-16_Real : 0.0_Real;

      // Column 0 is uniform in pressure; column 1 shares its surface and
      // bottom pressure but is offset from it by up to three layer
      // thicknesses in between, as coordinate tilt would make it.
      HostArray2DReal PressInterface("PressInterface", 2, NLayers + 1);
      HostArray2DReal PressMid("PressMid", 2, NLayers);
      HostArray2DReal ConservTemp("ConservTemp", 2, NLayers);
      HostArray2DReal AbsSalinity("AbsSalinity", 2, NLayers);
      HostArray2DReal SlopeCt("SlopeCt", 2, NLayers);
      HostArray2DReal SlopeSa("SlopeSa", 2, NLayers);

      for (int ICell = 0; ICell < 2; ++ICell) {
         for (int K = 0; K <= NLayers; ++K) {
            const Real Uniform = K * DeltaP;
            PressInterface(ICell, K) =
                (ICell == 1) ? Uniform + Bulge * std::sin(Pi * K / NLayers)
                             : Uniform;
         }
         for (int K = 0; K < NLayers; ++K) {
            PressMid(ICell, K) = 0.5_Real * (PressInterface(ICell, K) +
                                             PressInterface(ICell, K + 1));
            ConservTemp(ICell, K) =
                exactLayerAvg(Temp0, Temp1, Temp2, PressInterface(ICell, K),
                              PressInterface(ICell, K + 1));
            AbsSalinity(ICell, K) =
                exactLayerAvg(Salt0, Salt1, Salt2, PressInterface(ICell, K),
                              PressInterface(ICell, K + 1));
         }
         for (int K = 0; K < NLayers; ++K) {
            I4 KLo, KHi;
            linearReconStencil(K, KMin, KMax, KLo, KHi);
            SlopeCt(ICell, K) = linearReconSlope(
                ConservTemp(ICell, KLo), ConservTemp(ICell, KHi),
                PressMid(ICell, KLo), PressMid(ICell, KHi));
            SlopeSa(ICell, K) = linearReconSlope(
                AbsSalinity(ICell, KLo), AbsSalinity(ICell, KHi),
                PressMid(ICell, KLo), PressMid(ICell, KHi));
         }
      }

      Real Nodes[MaxPGradQuadPoints];
      Real Weights[MaxPGradQuadPoints];
      const I4 NQuad = 2;
      gaussLegendreRule(NQuad, Nodes, Weights);

      Real MaxScaled = 0.0_Real;
      Real MaxAbs    = 0.0_Real;

      for (int K = 0; K < NLayers; ++K) {

         // the edge layer spans the average of the two columns' interfaces
         const Real EdgeTop =
             0.5_Real * (PressInterface(0, K) + PressInterface(1, K));
         const Real EdgeBot =
             0.5_Real * (PressInterface(0, K + 1) + PressInterface(1, K + 1));
         const Real EdgeMid  = 0.5_Real * (EdgeTop + EdgeBot);
         const Real EdgeHalf = 0.5_Real * (EdgeBot - EdgeTop);

         for (int IQuad = 0; IQuad < NQuad; ++IQuad) {

            const Real Press = EdgeMid + EdgeHalf * Nodes[IQuad];

            Real TempAt[2];
            Real SaltAt[2];
            for (int ICell = 0; ICell < 2; ++ICell) {
               // each column's own layer containing this pressure, which
               // under the offset is generally not layer K
               const I4 KFound = findLayerForPress(PressInterface, ICell, KMin,
                                                   KMax, Press, K);
               TempAt[ICell]   = linearReconEval(ConservTemp(ICell, KFound),
                                                 SlopeCt(ICell, KFound),
                                                 PressMid(ICell, KFound), Press);
               SaltAt[ICell]   = linearReconEval(AbsSalinity(ICell, KFound),
                                                 SlopeSa(ICell, KFound),
                                                 PressMid(ICell, KFound), Press);
            }

            const Real SpecVolDiff = matchedPressSpecVolDiff(
                EdgeEos, TempAt[0], SaltAt[0], TempAt[1], SaltAt[1]);

            // scale by the size of the terms that had to cancel, not by an
            // absolute constant
            const Real Scale = std::abs(EdgeEos.SpecVolDCt * TempAt[0]) +
                               std::abs(EdgeEos.SpecVolDSa * SaltAt[0]);

            MaxAbs = std::max(MaxAbs, std::abs(SpecVolDiff));
            MaxScaled =
                std::max(MaxScaled, std::abs(SpecVolDiff) / (Eps * Scale));
         }
      }

      LOG_INFO("PGradTest: matched-pressure integrand, {} profile: max "
               "|dSpecVol| = {} m3/kg = {} eps of the cancelling terms",
               Curved ? "quadratic" : "linear", MaxAbs, MaxScaled);

      if (Curved) {
         MaxQuadratic = MaxScaled;
      } else {
         MaxLinear = MaxScaled;
      }
   }

   if (MaxLinear <= Tol) {
      LOG_INFO("PGradTest: matched-pressure integrand PASS: {} eps on a "
               "profile linear in pressure",
               MaxLinear);
   } else {
      LOG_ERROR("PGradTest: matched-pressure integrand FAIL: {} eps on a "
                "profile linear in pressure exceeds {} eps",
                MaxLinear, Tol);
      ++Err;
   }

   // the linear result means nothing unless the same machinery gives a large
   // answer on a profile the reconstruction does not resolve
   if (MaxQuadratic > 1.0e6_Real) {
      LOG_INFO("PGradTest: matched-pressure integrand control PASS: the "
               "quadratic profile gives {} eps",
               MaxQuadratic);
   } else {
      LOG_ERROR("PGradTest: matched-pressure integrand control FAIL: the "
                "quadratic profile gives only {} eps, so the linear result "
                "may be trivial",
                MaxQuadratic);
      ++Err;
   }

   return Err;

} // end testMatchedPressIntegrand

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

// Build a two-column state whose vertical profile lies in the scheme's exact
// set: temperature and salinity vary linearly with pressure, and each column's
// layer means are the *exact layer averages* of that one continuous profile
// over that column's own layers.
//
// Under the offset between the two columns those averages come out different
// in the two, and that is the point. A configuration built by copying
// identical layer means into both columns would not exercise the property at
// all, because the condition that matters is a property of the reconstructed
// profiles rather than of the layer means.
//
// The two columns share a surface pressure and a bottom pressure and sit on a
// flat floor; only the distribution of pressure between them differs, which is
// what an ALE coordinate does. A horizontally uniform surface pressure is
// required for an exactness gate: where surface pressure varies the state
// carries a real fixed-pressure height difference at the surface and zero is
// the wrong expectation.
//
// Returns the bottom pressure so tests can scale their tolerances by the
// hydrostatic terms.
Real setupExactSetState(HorzMesh *Mesh,    ///< [in] Horizontal mesh
                        VertCoord *VCoord, ///< [inout] Vertical coordinate
                        OceanState *State, ///< [inout] Ocean state
                        Eos *EqState,      ///< [inout] Equation of state
                        I4 NVertLayers,    ///< [in] Number of layers
                        Real DC,           ///< [in] Distance between cells
                        Real Offset        ///< [in] interface offset, in layers
) {

   const I4 NCellsAll = Mesh->NCellsAll;
   const I4 NEdgesAll = Mesh->NEdgesAll;

   VCoord->NVertLayers   = NVertLayers;
   VCoord->NVertLayersP1 = NVertLayers + 1;

   // the prescribed continuous profile, linear in pressure
   const Real PressBot = 4.0e7_Real;
   const Real Temp0 = 12.0_Real, Temp1 = -2.0e-7_Real;
   const Real Salt0 = 34.0_Real, Salt1 = 2.5e-8_Real;
   const Real Pi = 3.14159265358979323846_Real;

   const Real DeltaP = PressBot / NVertLayers;
   const Real Bulge  = Offset * DeltaP;

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

   Array1DReal SurfacePressure("SurfacePressure", Mesh->NCellsSize);
   Array2DReal PseudoThick = State->getPseudoThickness(0);
   Array2DReal Temp        = Tracers::getByName(0, "Temperature");
   Array2DReal Salinity    = Tracers::getByName(0, "Salinity");

   // Column 0 is uniform in pressure; column 1 bulges away from it by up to
   // Offset layer thicknesses in the interior while sharing both end
   // pressures. Pseudo-thickness and pressure thickness are proportional, so
   // setting one sets the other.
   auto &BottomGeomDepth = VCoord->BottomGeomDepth;
   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(int i) {
          SurfacePressure(i) = 0.0_Real;
          // a flat floor: the anchor's height difference is then exact input
          // and vanishes identically
          BottomGeomDepth(i) = 4000.0_Real;
          for (int k = 0; k < NVertLayers; ++k) {
             const Real Shape = (i == 1) ? 1.0_Real : 0.0_Real;
             const Real PTop =
                 k * DeltaP + Shape * Bulge * Kokkos::sin(Pi * k / NVertLayers);
             const Real PBot =
                 (k + 1) * DeltaP +
                 Shape * Bulge * Kokkos::sin(Pi * (k + 1) / NVertLayers);
             PseudoThick(i, k) = (PBot - PTop) / (Gravity * RhoSw);
          }
       });

   // pressure from the pseudo-thickness the model will actually use
   VCoord->computePressure(PseudoThick, SurfacePressure);

   // layer means as the exact layer averages of the prescribed profile, taken
   // over the model's own interface pressures so the two are consistent
   const auto &PressureInterface = VCoord->PressureInterface;
   parallelFor(
       {NCellsAll, NVertLayers}, KOKKOS_LAMBDA(int i, int k) {
          Temp(i, k) =
              exactLayerAvg(Temp0, Temp1, 0.0_Real, PressureInterface(i, k),
                            PressureInterface(i, k + 1));
          Salinity(i, k) =
              exactLayerAvg(Salt0, Salt1, 0.0_Real, PressureInterface(i, k),
                            PressureInterface(i, k + 1));
       });

   // specific volume and its derivatives, from one equation-of-state pass
   EqState->computeSpecVolAndDerivs(Temp, Salinity, VCoord->PressureMid);

   VCoord->computeGeomZHeight(PseudoThick, EqState->SpecVol);

   return PressBot;

} // end setupExactSetState

// Assert that the column scan returns a zero fixed-pressure height difference
// at *every* interface, not merely in the layer mean, for a profile the
// reconstruction resolves exactly.
//
// The shape of a failure is diagnostic and is reported: a residual growing
// with depth points at the recurrence, while one flat with depth points at the
// anchor. Here the two columns share a bottom pressure and a flat floor, so
// the anchor is identically zero by construction and any residual at all is
// the recurrence's.
int testColumnScan(HorzMesh *Mesh,    ///< [in] Horizontal mesh
                   VertCoord *VCoord, ///< [inout] Vertical coordinate
                   OceanState *State, ///< [inout] Ocean state
                   Eos *EqState       ///< [inout] Equation of state
) {

   int Err = 0;

   const I4 NVertLayers = 32;
   const Real DC        = 4000.0_Real;
   const Real Offset    = 3.0_Real;
   const Real Eps       = std::numeric_limits<Real>::epsilon();

   const Real PressBot = setupExactSetState(Mesh, VCoord, State, EqState,
                                            NVertLayers, DC, Offset);

   // A FiniteVolume instance, created after the state so that its working
   // arrays are sized for this number of layers
   Config *Options = Config::getOmegaConfig();
   Config PGradConfig("PressureGrad");
   Options->get(PGradConfig);
   PGradConfig.set("PressureGradType", std::string("FiniteVolume"));
   PressureGrad *FVPGrad =
       PressureGrad::create("TestColumnScan", Mesh, VCoord, Options);
   PGradConfig.set("PressureGradType", std::string("Centered"));

   Array2DReal Tend("TendColumnScan", Mesh->NEdgesSize, NVertLayers);
   deepCopy(Tend, 0.0_Real);

   Array2DReal PseudoThick = State->getPseudoThickness(0);
   Array2DReal Temp        = Tracers::getByName(0, "Temperature");
   Array2DReal Salinity    = Tracers::getByName(0, "Salinity");

   FVPGrad->computePressureGrad(
       Tend, VCoord->PressureMid, VCoord->PressureInterface, EqState->SpecVol,
       VCoord->GeomZInterface, PseudoThick, Temp, Salinity, EqState);

   // Scale by the size of the terms that had to cancel: the hydrostatic
   // contribution of the temperature and salinity structure over the column.
   const auto &SpecVolDCt = EqState->SpecVolDCt;
   const auto &SpecVolDSa = EqState->SpecVolDSa;
   Real CancelScale       = 0.0_Real;
   parallelReduce(
       {2, NVertLayers},
       KOKKOS_LAMBDA(int i, int k, Real &MaxV) {
          const Real V = (Kokkos::abs(SpecVolDCt(i, k) * Temp(i, k)) +
                          Kokkos::abs(SpecVolDSa(i, k) * Salinity(i, k))) *
                         PressBot / Gravity;
          if (V > MaxV)
             MaxV = V;
       },
       Kokkos::Max<Real>(CancelScale));

   const auto &DeltaZFixedP = FVPGrad->getDeltaZFixedP();

   Real MaxTop = 0.0_Real;
   Real MaxBot = 0.0_Real;
   Real MaxAll = 0.0_Real;
   parallelReduce(
       {Mesh->NEdgesAll, NVertLayers + 1},
       KOKKOS_LAMBDA(int IEdge, int K, Real &MaxA, Real &MaxT, Real &MaxB) {
          const Real V = Kokkos::abs(DeltaZFixedP(IEdge, K));
          if (V > MaxA)
             MaxA = V;
          // the shallower and deeper halves of the column, to show whether a
          // residual grows with depth
          if (K <= NVertLayers / 2 && V > MaxT)
             MaxT = V;
          if (K > NVertLayers / 2 && V > MaxB)
             MaxB = V;
       },
       Kokkos::Max<Real>(MaxAll), Kokkos::Max<Real>(MaxTop),
       Kokkos::Max<Real>(MaxBot));

   const Real Scaled = MaxAll / (Eps * CancelScale);

   LOG_INFO("PGradTest: column scan on the exact set: max |DeltaZFixedP| = {} "
            "m ({} eps of the cancelling terms, which are {} m); upper column "
            "{} m, lower column {} m",
            MaxAll, Scaled, CancelScale, MaxTop, MaxBot);

   // Provisional gate. The residual is expected to be a few eps of the
   // cancelling terms: the integrand is zero pointwise to about one eps and
   // the recurrence accumulates that over the column.
   const Real Tol = 100.0_Real;
   if (Scaled <= Tol) {
      LOG_INFO("PGradTest: column scan PASS: {} eps", Scaled);
   } else {
      LOG_ERROR("PGradTest: column scan FAIL: {} eps exceeds {} eps", Scaled,
                Tol);
      ++Err;
   }

   PressureGrad::erase("TestColumnScan");

   return Err;

} // end testColumnScan

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

      // The matched-pressure integrand must be zero at every quadrature
      // point on a resolved profile, not merely in the integral.
      int IntegrandErr = testMatchedPressIntegrand();

      // The column scan must return a zero fixed-pressure height difference
      // at every interface on the exact set, not merely in the layer mean.
      int ScanErr = testColumnScan(DefMesh, VCoord, DefState, DefEos);

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

      RetVal += ReconErr + LookupErr + IntegrandErr + IdentityErr + ConfigErr;

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
