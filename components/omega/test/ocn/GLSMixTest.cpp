//===-- Test driver for OMEGA GLS mixing -------------------------*- C++
//-*-===//
///
/// \file
/// \brief Unit tests for the Generic Length Scale (GLS) turbulence closure
///
//===----------------------------------------------------------------------===//

#include "GLSMix.h"
#include "Config.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Field.h"
#include "FillValues.h"
#include "GLSConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "KPPMix.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeMgr.h"
#include "VertCoord.h"
#include "mpi.h"

#include <memory>
#include <string>

using namespace OMEGA;

namespace {

constexpr int NVertLayers   = 20;
constexpr int NVertLayersP1 = NVertLayers + 1;
constexpr Real RTol         = 1.0e-10_Real;

std::unique_ptr<Clock> TestClock;
Clock *TestClockPtr = nullptr;

void initGLSMixTest(const std::string &TestGroup) {
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();
   initLogging(DefEnv);
   LOG_INFO("------ GLS Mixing Unit Tests ------");

   Config("Omega");
   Config::readAll("omega.yml");

   Config VertMixConfig("VertMix");
   Config GLSConfig("GLS");
   Error Err;
   Err += Config::getOmegaConfig()->get(VertMixConfig);
   Err += VertMixConfig.get(GLSConfig);
   CHECK_ERROR_ABORT(Err, "GLSMixTest: unable to access GLS configuration");

   if (TestGroup == "config-stability") {
      // GLSMix::init must abort on an unrecognized StabilityFunction.
      GLSConfig.set("StabilityFunction", "NotAGLSStabilityFunction");
   }
   if (TestGroup == "config-exclusive") {
      // GLSMix::init must abort when both GLS and KPP are enabled.
      Config KPPConfig("KPP");
      Err += VertMixConfig.get(KPPConfig);
      CHECK_ERROR_ABORT(Err, "GLSMixTest: unable to access KPP configuration");
      KPPConfig.set("Enable", true);
      GLSConfig.set("Enable", true);
   }

   IO::init(DefComm);
   Decomp::init();
   Halo::init();

   Calendar::init("No Leap");
   TimeInstant StartTime(0, 1, 1, 0, 0, 0.0);
   TimeInterval TimeStep(1, TimeUnits::Hours);
   TestClock    = std::make_unique<Clock>(StartTime, TimeStep);
   TestClockPtr = TestClock.get();
   Field::init(TestClockPtr);
   IOStream::init(TestClockPtr);
   HorzMesh::init(TestClockPtr);
   VertCoord::init(false);

   if (TestGroup == "config-exclusive") {
      KPPMix::init();
   }
   GLSMix::init();
}

void finalizeGLSMixTest() {
   GLSMix::destroyInstance();
   KPPMix::destroyInstance();
   IOStream::finalize();
   VertCoord::clear();
   HorzMesh::clear();
   Halo::clear();
   Decomp::clear();
   Field::clear();
   Dimension::clear();
   TestClockPtr = nullptr;
   TestClock.reset();
   MachEnv::removeAll();
}

void checkResult(const char *TestName, int NumErrors) {
   if (NumErrors != 0) {
      ABORT_ERROR("GLSMixTest: {} FAIL with {} errors", TestName, NumErrors);
   }
   LOG_INFO("GLSMixTest: {} PASS", TestName);
}

/// Pure function tests for the stability-function evaluator: Sm and Sh must
/// be finite and non-negative over a representative range of the
/// non-dimensional buoyancy/shear numbers, for every stability-function
/// choice.
void testStabilityFunctions() {
   constexpr int NTypes = 3;
   constexpr int NGh    = 5;
   int NumErrors         = 0;

   GLSStabilityType Types[NTypes] = {GLSStabilityType::CanutoA,
                                     GLSStabilityType::CanutoB,
                                     GLSStabilityType::Galperin};

   for (int ITy = 0; ITy < NTypes; ++ITy) {
      GLSDerivedConstants Coeffs = computeGLSDerivedConstants(
          Types[ITy], 3.0_Real, 1.5_Real, -1.0_Real, 0.5477_Real, 1.44_Real,
          1.92_Real, -0.4_Real, 1.0_Real, 1.0_Real, 1.3_Real);
      GLSStabilityType Type = Types[ITy];

      int TypeErrors = 0;
      parallelReduce(
          "GLSMixTest-StabilityFunctions", {NGh},
          KOKKOS_LAMBDA(int IGh, int &ErrorCount) {
             const Real Gh = -0.2_Real + 0.05_Real * IGh; // spans Ghmin..Gh0
             const Real GmBound =
                 computeCanutoGmBound(Type, Coeffs, Gh);
             const Real Gm = Kokkos::min(GmBound, 0.01_Real);

             Real Sm = 0.0_Real, Sh = 0.0_Real;
             computeStabilityFunctions(Type, Coeffs, Gh, Gm, Sm, Sh);

             if (Kokkos::isnan(Sm) || Kokkos::isnan(Sh) || Sm < 0.0_Real ||
                 Sh < 0.0_Real) {
                ++ErrorCount;
             }
          },
          TypeErrors);
      NumErrors += TypeErrors;
   }

   checkResult("testStabilityFunctions", NumErrors);
}

/// Regression test for the REMORA/ROMS length-scale-limiter bug: under
/// neutral/unstable stratification (where the Galperin bound does not
/// apply), the unconditional hard cap must still keep the turbulence length
/// scale -- and hence the resulting viscosity -- bounded. Before the fix,
/// an unbounded length scale here could drive viscosity to unphysical
/// values (REMORA observed ~4e6 m^2/s in a 23 m column).
void testLengthScaleCap() {
   const auto Mesh       = HorzMesh::getDefault();
   const auto VCoord     = VertCoord::getDefault();
   VCoord->NVertLayers   = NVertLayers;
   VCoord->NVertLayersP1 = NVertLayersP1;
   I4 NCellsSize         = Mesh->NCellsSize;
   I4 NCellsAll          = Mesh->NCellsAll;
   I4 NEdgesAll          = Mesh->NEdgesAll;

   // Uniform full-depth columns, 1 m layers (23 m total, matching the scale
   // of the documented REMORA failure case).
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   parallelFor(
       "setMinMax", {NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          MinLayerCell(ICell) = 0;
          MaxLayerCell(ICell) = NVertLayers - 1;
       });
   Halo *MeshHalo = Halo::getDefault();
   MeshHalo->exchangeFullArrayHalo(VCoord->MaxLayerCell, OnCell);
   VCoord->minMaxLayerEdge(MeshHalo);

   OMEGA_SCOPE(GeomZInterface, VCoord->GeomZInterface);
   parallelFor(
       "setGeom", {NCellsAll, NVertLayersP1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          GeomZInterface(ICell, K) = -1.0_Real * K; // 1 m layers
       });

   GLSMix *TestGLSMix  = GLSMix::getInstance();
   TestGLSMix->Enabled = true;

   auto NormalVelEdge = Array2DReal("NormalVelEdge", NEdgesAll, NVertLayers);
   auto TangVelEdge   = Array2DReal("TangVelEdge", NEdgesAll, NVertLayers);
   auto BVF = Array2DReal("BVF", NCellsSize, NVertLayersP1);
   auto PseudoThickCell =
       Array2DReal("PseudoThickCell", NCellsSize, NVertLayers);
   auto KineticEnergyCell =
       Array2DReal("KineticEnergyCell", NCellsSize, NVertLayers);
   auto ZonalStressCell = Array1DReal("ZonalStressCell", NCellsSize);
   auto MeridStressCell = Array1DReal("MeridStressCell", NCellsSize);

   // Strong shear, neutral-to-unstable stratification (BVF <= 0): the
   // Galperin bound does not engage here, so only the unconditional hard
   // cap can keep the length scale (and viscosity) bounded.
   deepCopy(NormalVelEdge, 0.3_Real);
   deepCopy(TangVelEdge, 0.0_Real);
   deepCopy(BVF, -1.0e-4_Real);
   deepCopy(PseudoThickCell, 1.0_Real);
   deepCopy(KineticEnergyCell, 0.05_Real);
   deepCopy(ZonalStressCell, 0.1_Real);
   deepCopy(MeridStressCell, 0.0_Real);

   const Real DT = 300.0_Real;
   for (int Step = 0; Step < 5; ++Step) {
      TestGLSMix->advance(NormalVelEdge, TangVelEdge, BVF, PseudoThickCell,
                         KineticEnergyCell, ZonalStressCell, MeridStressCell,
                         DT);
   }

   // The cap bounds L <= LengthScaleCapFraction * Depth; Sm/Sh are O(1), so
   // viscosity/diffusivity should stay well under a generous bound even
   // with a 23 m water column -- in contrast to the ~1e6 m^2/s blow-up
   // documented for the unfixed REMORA code.
   constexpr Real SanityBound = 50.0_Real;

   int NumErrors = 0;
   OMEGA_SCOPE(LocVertVisc, TestGLSMix->VertVisc);
   OMEGA_SCOPE(LocVertDiff, TestGLSMix->VertDiff);
   OMEGA_SCOPE(LocLscale, TestGLSMix->Lscale);
   const Real CapFraction = TestGLSMix->LengthScaleCapFraction;

   parallelReduce(
       "CheckLengthScaleCap", {NCellsAll, NVertLayersP1},
       KOKKOS_LAMBDA(int ICell, int K, int &ErrorCount) {
          const int KMin = MinLayerCell(ICell);
          const int KMax = MaxLayerCell(ICell) + 1;
          if (K < KMin || K > KMax)
             return;

          if (Kokkos::isnan(LocVertVisc(ICell, K)) ||
              Kokkos::isnan(LocVertDiff(ICell, K)) ||
              LocVertVisc(ICell, K) > SanityBound ||
              LocVertDiff(ICell, K) > SanityBound) {
             ++ErrorCount;
          }

          const Real Depth = Kokkos::abs(GeomZInterface(ICell, KMin) -
                                        GeomZInterface(ICell, KMax));
          if (LocLscale(ICell, K) >
              CapFraction * Depth + 1.0e-6_Real * Depth) {
             ++ErrorCount;
          }
       },
       NumErrors);

   checkResult("testLengthScaleCap", NumErrors);
}

} // namespace

int main(int argc, char *argv[]) {
   const std::string TestGroup = argc > 1 ? argv[1] : "all";

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   initGLSMixTest(TestGroup);

   if (TestGroup == "stability" || TestGroup == "all") {
      testStabilityFunctions();
   }
   if (TestGroup == "lengthcap" || TestGroup == "all") {
      testLengthScaleCap();
   }
   if (TestGroup != "stability" && TestGroup != "lengthcap" &&
       TestGroup != "all" && TestGroup != "config-stability" &&
       TestGroup != "config-exclusive") {
      ABORT_ERROR("GLSMixTest: unknown test group '{}'", TestGroup);
   }
   // config-stability and config-exclusive are expected to abort inside
   // initGLSMixTest(); reaching here means they did not, which is itself
   // a failure for those two groups (they are registered WILL_FAIL).

   finalizeGLSMixTest();

   LOG_INFO("------ GLS {} Tests Successful ------", TestGroup);
   Kokkos::finalize();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return 0;
}
