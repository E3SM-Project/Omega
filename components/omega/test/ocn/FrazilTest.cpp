//===-- Test driver for OMEGA FrazilFormation ---------------------------*- C++
//-*-===//
//
/// \file
/// \brief Minimal test driver for OMEGA frazil formation functor
//
//===-----------------------------------------------------------------------===/

#include "Frazil.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeMgr.h"
#include "VertCoord.h"
#include "mpi.h"

using namespace OMEGA;

constexpr int NVertLayers = 60;

void initFrazilTest(const std::string &mesh) {
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);
   LOG_INFO("------ Frazil Unit Tests ------");

   Config("Omega");
   Config::readAll("omega.yml");

   Config *OmegaConfig = Config::getOmegaConfig();
   Config TendConfig("Tendencies");
   OmegaConfig->get(TendConfig);
   TendConfig.set("FrazilTendencyEnable", true);

   Calendar::init("No Leap");
   TimeInstant StartTime(0, 1, 1, 0, 0, 0.0);
   TimeInterval TimeStep(1, TimeUnits::Hours);
   Clock ModelClockTmp(StartTime, TimeStep);
   Clock *ModelClock = &ModelClockTmp;

   IO::init(DefComm);
   Decomp::init(mesh);
   Field::init(ModelClock);
   IOStream::init(ModelClock);
   Halo::init();
   HorzMesh::init(ModelClock);
   VertCoord::init(false);
   Frazil::init();
}

void finalizeFrazilTest() {
   Frazil::clear();
   VertCoord::clear();
   HorzMesh::clear();
   Halo::clear();
   Decomp::clear();
   Field::clear();
   Dimension::clear();
   IOStream::finalize();
   MachEnv::removeAll();
}

// this test only excercises the frazil formation functor (no melt)
// in a cold case: the frazil terms should be
// - strictly positive for ice, liquid, and salt mass
//- strictly negative for ice and liquid energy
// - positive for T tendency and negative for S, H tendencies
void testFrazilFormationCold() {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   VCoord->NVertLayers = NVertLayers;

   const Real SAIn = 35.0_Real;
   const Real CTIn = -2.0_Real;
   const Real PIn  = 100.0_Real;
   const Real h    = 10.0_Real;
   const Real RTol = 1e-10_Real;

   (void)Mesh;

   FrazilFormation ComputeFrazilFormation;
   ComputeFrazilFormation.phi       = 0.75_Real;
   ComputeFrazilFormation.massLimit = 0.1_Real;

   Real AccMIce  = 0.0_Real;
   Real AccMLiq  = 0.0_Real;
   Real AccMSalt = 0.0_Real;
   Real AccELiq  = 0.0_Real;
   Real AccEIce  = 0.0_Real;

   Real HTend = 0.0_Real;
   Real TTend = 0.0_Real;
   Real STend = 0.0_Real;

   ComputeFrazilFormation(SAIn, CTIn, PIn, h, AccMIce, AccMLiq, AccMSalt,
                          AccELiq, AccEIce, HTend, TTend, STend);

   if (AccMIce <= 0.0_Real) {
      ABORT_ERROR(
          "FrazilFormationTestCold: accumulated ice mass is non-positive: {}",
          AccMIce);
   }
   if (AccMLiq <= 0.0_Real) {
      ABORT_ERROR("FrazilFormationTestCold: accumulated liquid mass is "
                  "non-positive: {}",
                  AccMLiq);
   }
   if (AccMSalt <= 0.0_Real) {
      ABORT_ERROR(
          "FrazilFormationTestCold: accumulated salt mass is non-positive: {}",
          AccMSalt);
   }

   if (AccELiq >= 0.0_Real) {
      ABORT_ERROR("FrazilFormationTestCold: accumulated liquid energy is "
                  "positive (exp. negative): {}",
                  AccELiq);
   }
   if (AccEIce >= 0.0_Real) {
      ABORT_ERROR("FrazilFormationTestCold: accumulated ice energy is positive "
                  "(exp. negative): {}",
                  AccEIce);
   }
   if (HTend >= 0.0_Real) {
      ABORT_ERROR(
          "FrazilFormationTestCold: HTend is positive (exp. negative): {}",
          HTend);
   }
   if (TTend <= 0.0_Real) {
      ABORT_ERROR(
          "FrazilFormationTestCold: TTend is negative (exp. positive): {}",
          TTend);
   }
   if (STend >= 0.0_Real) {
      ABORT_ERROR(
          "FrazilFormationTestCold: STend is positive (exp. negative): {}",
          STend);
   }
   LOG_INFO(
       "FrazilFormationTestCold: AccMIce = {}, AccMLiq = {}, AccMSalt = {}, "
       "AccELiq = {}, AccEIce = {}, HTend = {}, TTend = {}, STend = {}",
       AccMIce, AccMLiq, AccMSalt, AccELiq, AccEIce, HTend, TTend, STend);
}

// this test only excercises the frazil formation functor (no melt)
// in a warm case: the frazil FORMATION terms should all be zero
void testFrazilFormationWarm() {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   VCoord->NVertLayers = NVertLayers;

   const Real SAIn = 35.0_Real;
   const Real CTIn = 10.0_Real;
   const Real PIn  = 100.0_Real;
   const Real h    = 10.0_Real;
   const Real RTol = 1e-10_Real;

   (void)Mesh;

   FrazilFormation ComputeFrazilFormation;
   ComputeFrazilFormation.phi       = 0.75_Real;
   ComputeFrazilFormation.massLimit = 0.1_Real;

   Real AccMIce  = 0.0_Real;
   Real AccMLiq  = 0.0_Real;
   Real AccMSalt = 0.0_Real;
   Real AccELiq  = 0.0_Real;
   Real AccEIce  = 0.0_Real;

   Real HTend = 0.0_Real;
   Real TTend = 0.0_Real;
   Real STend = 0.0_Real;

   ComputeFrazilFormation(SAIn, CTIn, PIn, h, AccMIce, AccMLiq, AccMSalt,
                          AccELiq, AccEIce, HTend, TTend, STend);

   if (!isApprox(AccMIce, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFormationTest warm: expected zero AccMIce, got {}",
                  AccMIce);
   }

   if (!isApprox(AccMSalt, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFormationTest warm: expected zero AccMSalt, got {}",
                  AccMSalt);
   }
   if (!isApprox(AccMLiq, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFormationTest warm: expected zero AccMLiq, got {}",
                  AccMLiq);
   }

   if (!isApprox(AccELiq, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFormationTest warm: expected zero AccELiq, got {}",
                  AccELiq);
   }
   if (!isApprox(AccEIce, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFormationTest warm: expected zero AccEIce, got {}",
                  AccEIce);
   }
   if (!isApprox(HTend, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFormationTest warm: expected zero HTend, got {}",
                  HTend);
   }

   if (!isApprox(TTend, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFormationTest warm: expected zero TTend, got {}",
                  TTend);
   }

   if (!isApprox(STend, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFormationTest warm: expected zero STend, got {}",
                  STend);
   }
   LOG_INFO(
       "FrazilFormationTestWarm: AccMIce = {}, AccMLiq = {}, AccMSalt = {}, "
       "AccELiq = {}, AccEIce = {}, HTend = {}, TTend = {}, STend = {}",
       AccMIce, AccMLiq, AccMSalt, AccELiq, AccEIce, HTend, TTend, STend);
}

// this test only excercises the fixed-property frazil formation functor (no
// melt) in a warm case: the frazil FORMATION terms should all be zero
void testFixedPropertyFrazilFormationWarm() {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   VCoord->NVertLayers = NVertLayers;

   const Real SAIn = 35.0_Real;
   const Real CTIn = 10.0_Real;
   const Real PIn  = 100.0_Real;
   const Real h    = 10.0_Real;
   const Real RTol = 1e-10_Real;

   (void)Mesh;

   FixedPropertyFrazilFormation ComputeFrazilFormation;

   Real AccMIce  = 0.0_Real;
   Real AccMLiq  = 0.0_Real;
   Real AccMSalt = 0.0_Real;
   Real AccELiq  = 0.0_Real;
   Real AccEIce  = 0.0_Real;

   Real HTend = 0.0_Real;
   Real TTend = 0.0_Real;
   Real STend = 0.0_Real;

   Real CTfrz = gsw_ct_freezing_poly(SAIn, PIn, 0.0_Real);

   ComputeFrazilFormation(SAIn, CTIn, PIn, h, AccMIce, AccMSalt, AccEIce, HTend,
                          TTend, STend, CTfrz);

   if (!isApprox(AccMIce, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFixedPropertyFormationTest warm: expected zero "
                  "AccMIce, got {}",
                  AccMIce);
   }

   if (!isApprox(AccMSalt, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFixedPropertyFormationTest warm: expected zero "
                  "AccMSalt, got {}",
                  AccMSalt);
   }

   if (!isApprox(AccEIce, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFixedPropertyFormationTest warm: expected zero "
                  "AccEIce, got {}",
                  AccEIce);
   }
   if (!isApprox(HTend, 0.0_Real, RTol)) {
      ABORT_ERROR(
          "FrazilFixedPropertyFormationTest warm: expected zero HTend, got {}",
          HTend);
   }

   if (!isApprox(TTend, 0.0_Real, RTol)) {
      ABORT_ERROR(
          "FrazilFixedPropertyFormationTest warm: expected zero TTend, got {}",
          TTend);
   }

   if (!isApprox(STend, 0.0_Real, RTol)) {
      ABORT_ERROR(
          "FrazilFixedPropertyFormationTest warm: expected zero STend, got {}",
          STend);
   }
   LOG_INFO("FrazilFixedPropertyFormationTestWarm: AccMIce = {}, AccMLiq = {}, "
            "AccMSalt = {}, "
            "AccELiq = {}, AccEIce = {}, HTend = {}, TTend = {}, STend = {}",
            AccMIce, AccMLiq, AccMSalt, AccELiq, AccEIce, HTend, TTend, STend);
}

// this test only excercises the frazil formation functor (no melt)
// in a cold case: the frazil terms should be
// - strictly positive for ice, liquid, and salt mass
//- strictly negative for ice and liquid energy
// - positive for T tendency and negative for S, H tendencies
void testFixedPropertyFrazilFormationCold() {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   VCoord->NVertLayers = NVertLayers;

   const Real SAIn = 35.0_Real;
   const Real CTIn = -2.0_Real;
   const Real PIn  = 100.0_Real;
   const Real h    = 10.0_Real;
   const Real RTol = 1e-10_Real;

   (void)Mesh;

   FixedPropertyFrazilFormation ComputeFrazilFormation;

   Real AccMIce  = 0.0_Real;
   Real AccMLiq  = 0.0_Real;
   Real AccMSalt = 0.0_Real;
   Real AccELiq  = 0.0_Real;
   Real AccEIce  = 0.0_Real;

   Real HTend = 0.0_Real;
   Real TTend = 0.0_Real;
   Real STend = 0.0_Real;

   Real CTfrz = gsw_ct_freezing_poly(SAIn, PIn, 0.0_Real);

   ComputeFrazilFormation(SAIn, CTIn, PIn, h, AccMIce, AccMSalt, AccEIce, HTend,
                          TTend, STend, CTfrz);

   if (AccMIce <= 0.0_Real) {
      ABORT_ERROR(
          "FrazilFixedPropertyFormationTestCold: accumulated ice mass is "
          "non-positive: {}",
          AccMIce);
   }

   if (AccMSalt <= 0.0_Real) {
      ABORT_ERROR(
          "FrazilFixedPropertyFormationTestCold: accumulated salt mass is "
          "non-positive: {}",
          AccMSalt);
   }

   if (AccEIce >= 0.0_Real) {
      ABORT_ERROR("FrazilFixedPropertyFormationTestCold: accumulated ice "
                  "energy is positive "
                  "(exp. negative): {}",
                  AccEIce);
   }
   if (HTend >= 0.0_Real) {
      ABORT_ERROR("FrazilFixedPropertyFormationTestCold: HTend is positive "
                  "(exp. negative): {}",
                  HTend);
   }
   if (TTend <= 0.0_Real) {
      ABORT_ERROR("FrazilFixedPropertyFormationTestCold: TTend is negative "
                  "(exp. positive): {}",
                  TTend);
   }
   if (STend >= 0.0_Real) {
      ABORT_ERROR("FrazilFixedPropertyFormationTestCold: STend is positive "
                  "(exp. negative): {}",
                  STend);
   }
   LOG_INFO("FrazilFixedPropertyFormationTestCold: AccMIce = {}, AccMLiq = {}, "
            "AccMSalt = {}, "
            "AccELiq = {}, AccEIce = {}, HTend = {}, TTend = {}, STend = {}",
            AccMIce, AccMLiq, AccMSalt, AccELiq, AccEIce, HTend, TTend, STend);
}

// this test exercises the frazil formation and melt functors
// in a column of water with both cold and warm layers.
// It turns to frazil column conservation check.
// In dev, there is extra verbose logging in the frazil code (TBRemoved).
void testComputeFrazilColumn() {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();
   auto *TestFrazil  = Frazil::getDefault();

   if (!TestFrazil) {
      ABORT_ERROR("FrazilTestColumn: default frazil object is null");
   }

   const Real RTol    = 1e-10_Real;
   const Real SACold  = 35.0_Real;
   const Real PRef    = 100.0_Real;
   const Real HRef    = 10.0_Real;
   const Real CTCold  = -2.0_Real;
   const Real CTWarm  = 0.0_Real;
   const Real CTWarm2 = -1.9_Real;

   Array2DReal SA("SA", Mesh->NCellsSize, NVertLayers);
   Array2DReal CT("CT", Mesh->NCellsSize, NVertLayers);
   Array2DReal P("P", Mesh->NCellsSize, NVertLayers);
   Array2DReal H("H", Mesh->NCellsSize, NVertLayers);

   deepCopy(SA, SACold);
   deepCopy(CT, CTWarm);
   deepCopy(P, PRef);
   deepCopy(H, HRef);

   deepCopy(TestFrazil->AccMIce, 0.0_Real);
   deepCopy(TestFrazil->AccMLiq, 0.0_Real);
   deepCopy(TestFrazil->AccMSalt, 0.0_Real);
   deepCopy(TestFrazil->AccELiq, 0.0_Real);
   deepCopy(TestFrazil->AccEIce, 0.0_Real);
   deepCopy(TestFrazil->FrazilHTend, 0.0_Real);
   deepCopy(TestFrazil->FrazilTTend, 0.0_Real);
   deepCopy(TestFrazil->FrazilSTend, 0.0_Real);

   auto MinLayerCellH = createHostMirrorCopy(VCoord->MinLayerCell);
   auto MaxLayerCellH = createHostMirrorCopy(VCoord->MaxLayerCell);

   const I4 ICell = 0;
   const I4 KMin  = MinLayerCellH(ICell);
   const I4 KMax  = MaxLayerCellH(ICell);
   if ((KMax - KMin + 1) < 4) {
      ABORT_ERROR("FrazilTestColumn: cell {} has fewer than 4 active layers",
                  ICell);
   }

   const I4 KBottom0  = KMax;
   const I4 KBottom1  = KMax - 1;
   const I4 KWarm     = KMax - 2;
   const I4 KTopCold  = KMax - 3;
   const I4 KCold2    = KMin + 3;
   const I4 KCold3    = KMin + 2;
   const I4 KWarm2    = KMin + 1;
   const I4 KTopCold2 = KMin;

   auto CTH               = createHostMirrorCopy(CT);
   CTH(ICell, KBottom0)   = CTCold;
   CTH(ICell, KBottom1)   = CTCold;
   CTH(ICell, KWarm)      = CTWarm;
   CTH(ICell, KTopCold)   = CTCold;
   CTH(ICell, KCold2 + 2) = CTCold;
   CTH(ICell, KCold2 + 1) = CTCold - .5_Real;
   CTH(ICell, KCold2)     = CTCold;
   CTH(ICell, KCold3)     = CTCold;
   CTH(ICell, KWarm2)     = CTWarm2;
   CTH(ICell, KTopCold2)  = CTCold;
   deepCopy(CT, CTH);

   const bool SavedConservationCheck = TestFrazil->conservationCheck;
   // TestFrazil->conservationCheck     = true;

   TestFrazil->computeFrazil(CT, SA, P, H);
   TestFrazil->conservationCheck = SavedConservationCheck;

   auto HTendH = createHostMirrorCopy(TestFrazil->FrazilHTend);
   auto TTendH = createHostMirrorCopy(TestFrazil->FrazilTTend);
   auto STendH = createHostMirrorCopy(TestFrazil->FrazilSTend);

   if (HTendH(ICell, KBottom0) >= 0.0_Real ||
       TTendH(ICell, KBottom0) <= 0.0_Real ||
       STendH(ICell, KBottom0) >= 0.0_Real) {
      ABORT_ERROR(
          "FrazilTestColumn: bottom cold layer sign check failed (HTend<0, "
          "TTend>0, STend<0 expected)");
   }

   if (HTendH(ICell, KBottom1) >= 0.0_Real ||
       TTendH(ICell, KBottom1) <= 0.0_Real ||
       STendH(ICell, KBottom1) >= 0.0_Real) {
      ABORT_ERROR(
          "FrazilTestColumn: second cold layer sign check failed (HTend<0, "
          "TTend>0, STend<0 expected)");
   }

   if (HTendH(ICell, KWarm) < 0.0_Real || TTendH(ICell, KWarm) > 0.0_Real ||
       STendH(ICell, KWarm) < 0.0_Real) {
      ABORT_ERROR("FrazilTestColumn: warm layer sign check failed (HTend>=0, "
                  "TTend<=0, STend>=0 expected)");
   }

   if (HTendH(ICell, KTopCold) >= 0.0_Real ||
       TTendH(ICell, KTopCold) <= 0.0_Real ||
       STendH(ICell, KTopCold) >= 0.0_Real) {
      ABORT_ERROR(
          "FrazilTestColumn: top cold layer sign check failed (HTend<0, "
          "TTend>0, STend<0 expected)");
   }

   LOG_INFO(
       "FrazilTestColumn: XTend branch-switching checks passed for ICell={}",
       ICell);
}

// this test exercises the frazil formation and melt functors
// with a depth limit set. Layers deeper than the depth limit
// should have zero frazil tendencies.
void testComputeFrazilDepthLimit() {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();
   auto *TestFrazil  = Frazil::getDefault();

   if (!TestFrazil) {
      ABORT_ERROR("FrazilTestColumn: default frazil object is null");
   }

   const Real RTol    = 1e-12_Real;
   const Real SACold  = 35.0_Real;
   const Real PRef    = 100.0_Real;
   const Real HRef    = 10.0_Real;
   const Real CTCold  = -2.0_Real;
   const Real CTWarm  = 0.0_Real;
   const Real CTWarm2 = -1.9_Real;

   Array2DReal SA("SA", Mesh->NCellsSize, NVertLayers);
   Array2DReal CT("CT", Mesh->NCellsSize, NVertLayers);
   Array2DReal P("P", Mesh->NCellsSize, NVertLayers);
   Array2DReal H("H", Mesh->NCellsSize, NVertLayers);

   deepCopy(SA, SACold);
   deepCopy(CT, CTWarm);
   deepCopy(P, PRef);
   deepCopy(H, HRef);

   deepCopy(TestFrazil->AccMIce, 0.0_Real);
   deepCopy(TestFrazil->AccMLiq, 0.0_Real);
   deepCopy(TestFrazil->AccMSalt, 0.0_Real);
   deepCopy(TestFrazil->AccELiq, 0.0_Real);
   deepCopy(TestFrazil->AccEIce, 0.0_Real);
   deepCopy(TestFrazil->FrazilHTend, 0.0_Real);
   deepCopy(TestFrazil->FrazilTTend, 0.0_Real);
   deepCopy(TestFrazil->FrazilSTend, 0.0_Real);

   auto MinLayerCellH = createHostMirrorCopy(VCoord->MinLayerCell);
   auto MaxLayerCellH = createHostMirrorCopy(VCoord->MaxLayerCell);

   const I4 ICell = 0;
   const I4 KMin  = MinLayerCellH(ICell);
   const I4 KMax  = MaxLayerCellH(ICell);
   if ((KMax - KMin + 1) < 10) {
      ABORT_ERROR("FrazilTestColumn: cell {} has fewer than 10 active layers",
                  ICell);
   }

   const I4 KBottom0  = KMax;
   const I4 KBottom1  = KMax - 1;
   const I4 KWarm     = KMax - 2;
   const I4 KTopCold  = KMax - 3;
   const I4 KCold2    = KMin + 3;
   const I4 KCold3    = KMin + 2;
   const I4 KWarm2    = KMin + 1;
   const I4 KTopCold2 = KMin;

   auto CTH               = createHostMirrorCopy(CT);
   CTH(ICell, KBottom0)   = CTCold;
   CTH(ICell, KBottom1)   = CTCold;
   CTH(ICell, KWarm)      = CTWarm;
   CTH(ICell, KTopCold)   = CTCold;
   CTH(ICell, KCold2 + 2) = CTCold;
   CTH(ICell, KCold2 + 1) = CTCold - .5_Real;
   CTH(ICell, KCold2)     = CTCold;
   CTH(ICell, KCold3)     = CTCold;
   CTH(ICell, KWarm2)     = CTWarm2;
   CTH(ICell, KTopCold2)  = CTCold;
   deepCopy(CT, CTH);

   const bool SavedConservationCheck = TestFrazil->conservationCheck;
   const bool SavedDepthLimit        = TestFrazil->depthLimit;
   const Real TestDepthLimit         = 35.0_Real; // this needs to be positive
   // if TestDepthLimit is negative, test will fail:
   // - the code assume depthlimit < 0 mean no limit (i.e. full depth frazil)
   // - the test below will exclude all layers and fail because Tend !=0.

   // Populate GeomZMid explicitly for the test column.
   auto GeomZMidH = createHostMirrorCopy(VCoord->GeomZMid);
   for (I4 K = KMin; K <= KMax; ++K) {
      GeomZMidH(ICell, K) = -10.0_Real * (K - KMin + 1);
   }
   deepCopy(VCoord->GeomZMid, GeomZMidH);

   TestFrazil->conservationCheck = true;
   TestFrazil->depthLimit        = TestDepthLimit;
   TestFrazil->computeFrazil(CT, SA, P, H);
   TestFrazil->conservationCheck = SavedConservationCheck;
   TestFrazil->depthLimit        = SavedDepthLimit;

   auto HTendH = createHostMirrorCopy(TestFrazil->FrazilHTend);
   auto TTendH = createHostMirrorCopy(TestFrazil->FrazilTTend);
   auto STendH = createHostMirrorCopy(TestFrazil->FrazilSTend);

   bool FoundExcludedLayer = false;
   for (I4 K = KMin; K <= KMax; ++K) {
      const Real Depth    = GeomZMidH(ICell, K);
      const Real AbsDepth = Depth < 0.0_Real ? -Depth : Depth;

      if (AbsDepth > TestDepthLimit) {
         FoundExcludedLayer = true;
         if (!isApprox(HTendH(ICell, K), 0.0_Real, RTol) ||
             !isApprox(TTendH(ICell, K), 0.0_Real, RTol) ||
             !isApprox(STendH(ICell, K), 0.0_Real, RTol)) {
            ABORT_ERROR("FrazilDepthLimitTest: excluded layer K={} has "
                        "non-zero tendencies (HTend={}, TTend={}, STend={})",
                        K, HTendH(ICell, K), TTendH(ICell, K),
                        STendH(ICell, K));
         }
      }
   }
   if (!FoundExcludedLayer) {
      ABORT_ERROR("FrazilDepthLimitTest: no layers were excluded for ICell={} "
                  "with depthLimit={}",
                  ICell, TestDepthLimit);
   }

   LOG_INFO("FrazilDepthLimitTest: depthLimit={} exclusion check passed for "
            "ICell={}",
            TestDepthLimit, ICell);
}

void frazilTest(const std::string &MeshFile = "OmegaMesh.nc") {
   initFrazilTest(MeshFile);
   testFrazilFormationCold();
   testFrazilFormationWarm();
   testFixedPropertyFrazilFormationCold();
   testFixedPropertyFrazilFormationWarm();
   testComputeFrazilColumn();
   testComputeFrazilDepthLimit();
   finalizeFrazilTest();
}

int main(int argc, char *argv[]) {
   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   frazilTest();

   LOG_INFO("------ Frazil Unit Tests Successful ------");

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   return 0;
}
