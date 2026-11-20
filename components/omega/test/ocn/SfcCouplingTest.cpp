#include "SfcCoupling.h"
#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Eos.h"
#include "Field.h"
#include "Forcing.h"
#include "Frazil.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "VertAdv.h"
#include "VertCoord.h"
#include "mpi.h"

#include <algorithm>

using namespace OMEGA;

struct TestSetup {

   // Coupler buffer widths (fields per cell)
   int NImports = 17;
   int NExports = 12;

   std::map<std::string, int> ImportIdxMap = {
       {"Foxx_taux", 3},  {"Foxx_tauy", 8},  {"Foxx_swnet", 0},
       {"Foxx_sen", 1},   {"Foxx_lat", 2},   {"Foxx_lwup", 4},
       {"Faxa_lwdn", 5},  {"Fioi_salt", 6},  {"Fioi_melth", 7},
       {"Fioi_meltw", 9}, {"Faxa_snow", 10}, {"Faxa_rain", 11},
       {"Foxx_evap", 12}, {"Foxx_rofl", 13}, {"Foxx_rofi", 14},
       {"Si_bpress", 15}, {"Sa_pslv", 16}};
   std::map<std::string, int> ExportIdxMap = {
       {"So_t", 2},         {"So_s", 4},          {"So_u", 6},
       {"So_v", 9},         {"So_dhdx", 5},       {"So_dhdy", 3},
       {"So_ssh", 1},       {"Faoo_h2otemp", 0},  {"Fioo_q", 8},
       {"Fioo_frazilm", 7}, {"Fioo_frazils", 11}, {"Fioo_frazilh", 10}};
};

CouplingInitParams mockCouplingInitParams(
    CouplingLayout Layout                        = CouplingLayout::MCT,
    std::optional<TimeInterval> CouplingTimeStep = std::nullopt) {

   TestSetup Setup;

   auto *DefTimeStepper = TimeStepper::getDefault();
   TimeInterval CouplingTimeStep_ =
       CouplingTimeStep.value_or(DefTimeStepper->getTimeStep());

   CouplingInitParams CouplingParams{.NImportFields    = Setup.NImports,
                                     .NExportFields    = Setup.NExports,
                                     .ImportIdxMap     = Setup.ImportIdxMap,
                                     .ExportIdxMap     = Setup.ExportIdxMap,
                                     .CouplingTimeStep = CouplingTimeStep_,
                                     .Layout           = Layout};

   return CouplingParams;
}

std::string toString(const CouplingLayout &Layout) {
   switch (Layout) {
   case CouplingLayout::MCT:
      return "MCT";
   case CouplingLayout::MOAB:
      return "MOAB";
   default:
      return "Unknown";
   }
}

HostArray1DReal makeCellVarryingArray(const std::string &Name, const int NCells,
                                      const Real Base) {
   HostArray1DReal Array(Name, NCells);
   for (int Cell = 0; Cell < NCells; Cell++) {
      Array(Cell) = Base + Cell;
   }
   return Array;
}

// Shared index formula for packing/unpacking from raw coupler buffer
int flatIdx(const CouplingLayout Layout, const int Cell, const int Field,
            const int NCells, const int NFields) {
   if (Layout == CouplingLayout::MCT)
      return Cell * NFields + Field;
   else // MOAB
      return Field * NCells + Cell;
}

// Set surface-layer temperature and salinity to Base + Cell in tracer level 0
void setSfcTracers(const int NCells, const Real TempBase,
                   const Real SalinBase) {
   VertCoord *DefVertCoord = VertCoord::getDefault();

   I4 TempIdx, SalinIdx;
   Tracers::getIndex(TempIdx, "Temperature");
   Tracers::getIndex(SalinIdx, "Salinity");

   HostArray2DReal TempH  = Tracers::getHostByIndex(0, TempIdx);
   HostArray2DReal SalinH = Tracers::getHostByIndex(0, SalinIdx);

   for (int Cell = 0; Cell < NCells; Cell++) {
      const int KSfc     = DefVertCoord->MinLayerCellH(Cell);
      TempH(Cell, KSfc)  = TempBase + Cell;
      SalinH(Cell, KSfc) = SalinBase + Cell;
   }
   Tracers::copyToDevice(0);
}

int initSfcCouplingTest(const std::string &MeshFile) {

   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);
   LOG_INFO("------ Surface Coupling unit tests ------");

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   // Force ConstantEos so calcPtFromCt() is a no-op, letting tests below use
   // simple mocked Temperature/Salinity without the TEOS-10 polynomial
   Config *OmegaConfig = Config::getOmegaConfig();
   Config EosConfig("Eos");
   OmegaConfig->get(EosConfig);
   EosConfig.set("EosType", std::string("constant"));

   // Initialize time stepper and get model clock
   TimeStepper::init1();
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   IO::init(DefComm);
   Decomp::init(MeshFile);

   // Initialize streams
   Field::init(ModelClock);
   IOStream::init(ModelClock);

   // Initialize the halo used by the pressure import fields.
   Halo::init();
   HorzMesh::init(ModelClock);

   VertCoord::init();

   int StateErr = OceanState::init();
   if (StateErr != 0) {
      Err++;
      LOG_ERROR("SfcCouplingTest: Error initializing defualt ocean state");
   }

   Eos::init();
   Frazil::init();

   Forcing::init();
   Tracers::init();

   // Needed by SfcCoupling::updateExportFields, which recomputes the
   // momentum/vertical aux variables so SshCell is current
   VertAdv::init();
   AuxiliaryState::init();

   return Err;
}

// Import: fill a raw x2o buffer, import it, and check the owned cells of the
// Forcing and VertCoord arrays it is meant to populate
int testImportFromCoupler(const CouplingLayout Layout) {

   int Err = 0;

   TestSetup Setup;
   auto CouplingParams = mockCouplingInitParams(Layout);
   Err += SfcCoupling::init(CouplingParams);

   SfcCoupling *DefCoupling = SfcCoupling::getDefault();
   Forcing *DefForcing      = Forcing::getDefault();
   VertCoord *DefVertCoord  = VertCoord::getDefault();

   const int NCells   = HorzMesh::getDefault()->NCellsOwned;
   const int NImports = Setup.NImports;
   const int NExports = Setup.NExports;

   std::vector<Real> CplToOcnData(NCells * NImports, 0.0);
   std::vector<Real> OcnToCplData(NCells * NExports, 0.0);

   // Give every field a distinct, cell-varying value: FieldIdx + Cell
   for (const auto &[Name, FieldIdx] : Setup.ImportIdxMap) {
      for (int Cell = 0; Cell < NCells; Cell++) {
         CplToOcnData[flatIdx(Layout, Cell, FieldIdx, NCells, NImports)] =
             static_cast<Real>(FieldIdx + Cell);
      }
   }

   // Push one cell's sea ice pressure above the 5 m limit so the clamp is
   // exercised
   const int BPressIdx = Setup.ImportIdxMap.at("Si_bpress");
   const int PslvIdx   = Setup.ImportIdxMap.at("Sa_pslv");
   CplToOcnData[flatIdx(Layout, 0, BPressIdx, NCells, NImports)] =
       2.0_Real * MaxSeaIcePressure;

   DefCoupling->attachData(CplToOcnData.data(), OcnToCplData.data());
   DefCoupling->importFromCoupler();
   DefCoupling->applyImportFields(DefForcing, DefVertCoord);

   auto checkImportField = [&](const Array1DReal &Field,
                               const std::string &Name) {
      const int FieldIdx = Setup.ImportIdxMap.at(Name);
      HostArray1DReal Expected =
          makeCellVarryingArray("Expected" + Name, NCells, Real(FieldIdx));
      auto Owned = Kokkos::subview(Field, std::pair(0, NCells));
      if (!arraysEqual(Owned, Expected)) {
         LOG_ERROR("SfcCouplingTest: import of {} FAIL", Name);
         return false;
      }
      return true;
   };

   HostArray1DReal ExpectedSfcPress("ExpectedSurfacePressure", NCells);
   for (int Cell = 0; Cell < NCells; Cell++) {
      const Real BPress = std::min(
          CplToOcnData[flatIdx(Layout, Cell, BPressIdx, NCells, NImports)],
          MaxSeaIcePressure);
      ExpectedSfcPress(Cell) = BPress + Real(PslvIdx + Cell) - AtmRefP;
   }
   auto OwnedSfcPress =
       Kokkos::subview(DefVertCoord->SurfacePressure, std::pair(0, NCells));
   const bool SfcPressPass = arraysEqual(OwnedSfcPress, ExpectedSfcPress);
   if (!SfcPressPass)
      LOG_ERROR("SfcCouplingTest: import of SurfacePressure FAIL");

   const auto &Stress = DefForcing->SfcStressForcing;
   const auto &Tracer = DefForcing->TracerForcing;

   // Not short-circuited, so every failing field is logged
   bool ImportPass = SfcPressPass;
   ImportPass &= checkImportField(Stress.ZonalStressCell, "Foxx_taux");
   ImportPass &= checkImportField(Stress.MeridStressCell, "Foxx_tauy");
   ImportPass &= checkImportField(Tracer.ShortWaveHeatFluxCell, "Foxx_swnet");
   ImportPass &= checkImportField(Tracer.SensibleHeatFluxCell, "Foxx_sen");
   ImportPass &= checkImportField(Tracer.LatentHeatFluxEvapCell, "Foxx_lat");
   ImportPass &= checkImportField(Tracer.LongWaveHeatFluxUpCell, "Foxx_lwup");
   ImportPass &= checkImportField(Tracer.LongWaveHeatFluxDownCell, "Faxa_lwdn");
   ImportPass &= checkImportField(Tracer.SeaIceSaltFluxCell, "Fioi_salt");
   ImportPass &= checkImportField(Tracer.SeaIceHeatFluxCell, "Fioi_melth");
   ImportPass &=
       checkImportField(Tracer.SeaIceFreshWaterFluxCell, "Fioi_meltw");
   ImportPass &= checkImportField(Tracer.SnowFluxCell, "Faxa_snow");
   ImportPass &= checkImportField(Tracer.RainFluxCell, "Faxa_rain");
   ImportPass &= checkImportField(Tracer.EvaporationFluxCell, "Foxx_evap");
   ImportPass &= checkImportField(Tracer.RiverRunoffFluxCell, "Foxx_rofl");
   ImportPass &= checkImportField(Tracer.IceRunoffFluxCell, "Foxx_rofi");

   if (ImportPass) {
      LOG_INFO("SfcCouplingTest: import with {} layout PASS", toString(Layout));
   } else {
      Err++;
      LOG_ERROR("SfcCouplingTest: import with {} layout FAIL",
                toString(Layout));
   }

   SfcCoupling::clear();

   return Err;
}

// Averaging: accumulate NSteps of varying surface tracers, export, and check
// the raw o2x buffer holds their mean (converted to coupler units)
int testUpdateExportFields(const I4 NSteps) {

   int Err = 0;

   TestSetup Setup;
   auto *DefStepper            = TimeStepper::getDefault();
   const CouplingLayout Layout = CouplingLayout::MCT;

   // Coupling interval spans NSteps ocean timesteps
   auto CouplingParams =
       mockCouplingInitParams(Layout, DefStepper->getTimeStep() * NSteps);

   Err += SfcCoupling::init(CouplingParams);

   SfcCoupling *DefCoupling = SfcCoupling::getDefault();
   OceanState *DefState     = OceanState::getDefault();
   Frazil *DefFrazil        = Frazil::getDefault();

   const int NCells   = HorzMesh::getDefault()->NCellsOwned;
   const int NImports = Setup.NImports;
   const int NExports = Setup.NExports;

   std::vector<Real> CplToOcnData(NCells * NImports, 0.0);
   std::vector<Real> OcnToCplData(NCells * NExports, 0.0);
   DefCoupling->attachData(CplToOcnData.data(), OcnToCplData.data());

   const int TempIdx  = Setup.ExportIdxMap.at("So_t");
   const int SalinIdx = Setup.ExportIdxMap.at("So_s");

   const Real TempBase  = static_cast<Real>(TempIdx);
   const Real SalinBase = static_cast<Real>(SalinIdx);

   // Clock and alarm behavior is covered by TimeMgrTest; keeping this test
   // focused on field accumulation avoids advancing a shared clock with
   // alarms from cleared test objects.
   for (I4 Step = 0; Step < NSteps; ++Step) {
      setSfcTracers(NCells, TempBase + Step, SalinBase + Step);
      DefCoupling->updateExportFields(DefState, Tracers::getAll(0), DefFrazil);
   }

   DefCoupling->exportToCoupler();

   // Only accumulated round-off error (no discretization/truncation error),
   // so use a tighter tolerance than tests comparing against analytic
   // solutions of numerical operators. ConstantEos makes CT->PT the identity.
   const Real RTol       = sizeof(Real) == 4 ? 1e-5 : 1e-10;
   const Real StepOffset = static_cast<Real>(NSteps - 1) / 2.0;

   I4 TempErr  = 0;
   I4 SalinErr = 0;

   for (int Cell = 0; Cell < NCells; Cell++) {
      const Real ExpectedTemp  = TempBase + Cell + StepOffset + TkFrz;
      const Real ExpectedSalin = (SalinBase + Cell + StepOffset) / Psu2Gpkg;

      const Real Temp =
          OcnToCplData[flatIdx(Layout, Cell, TempIdx, NCells, NExports)];
      const Real Salin =
          OcnToCplData[flatIdx(Layout, Cell, SalinIdx, NCells, NExports)];

      if (!isApprox(Temp, ExpectedTemp, RTol))
         TempErr++;
      if (!isApprox(Salin, ExpectedSalin, RTol))
         SalinErr++;
   }

   if (TempErr == 0) {
      LOG_INFO("SfcCouplingTest: {}-step average of So_t PASS", NSteps);
   } else {
      Err += TempErr;
      LOG_ERROR("SfcCouplingTest: {}-step average of So_t FAIL - outside "
                "relative tolerance of {} for {} cells",
                NSteps, RTol, TempErr);
   }

   if (SalinErr == 0) {
      LOG_INFO("SfcCouplingTest: {}-step average of So_s PASS", NSteps);
   } else {
      Err += SalinErr;
      LOG_ERROR("SfcCouplingTest: {}-step average of So_s FAIL - outside "
                "relative tolerance of {} for {} cells",
                NSteps, RTol, SalinErr);
   }

   SfcCoupling::clear();
   return Err;
}

// Export: check fields land in the right slots of the raw o2x buffer for the
// given layout, then check a second interval doesn't see the first one
int testExportToCoupler(const CouplingLayout Layout) {

   int Err = 0;

   TestSetup Setup;
   auto CouplingParams = mockCouplingInitParams(Layout);
   Err += SfcCoupling::init(CouplingParams);

   SfcCoupling *DefCoupling = SfcCoupling::getDefault();
   OceanState *DefState     = OceanState::getDefault();
   VertCoord *DefVertCoord  = VertCoord::getDefault();
   Frazil *DefFrazil        = Frazil::getDefault();

   const int NCells   = HorzMesh::getDefault()->NCellsOwned;
   const int NImports = Setup.NImports;
   const int NExports = Setup.NExports;

   std::vector<Real> CplToOcnData(NCells * NImports, 0.0);
   std::vector<Real> OcnToCplData(NCells * NExports, 0.0);
   DefCoupling->attachData(CplToOcnData.data(), OcnToCplData.data());

   const int TempIdx  = Setup.ExportIdxMap.at("So_t");
   const int SalinIdx = Setup.ExportIdxMap.at("So_s");
   const int SshIdx   = Setup.ExportIdxMap.at("So_ssh");
   const int UIdx     = Setup.ExportIdxMap.at("So_u");
   const int VIdx     = Setup.ExportIdxMap.at("So_v");

   auto exported = [&](const int Cell, const int FieldIdx) {
      return OcnToCplData[flatIdx(Layout, Cell, FieldIdx, NCells, NExports)];
   };

   // Save the shared normal velocity so it can be restored at the end
   Array2DReal NormVel = DefState->getNormalVelocity(0);
   Array2DReal SavedNormVel("SavedNormVel", NormVel.extent(0),
                            NormVel.extent(1));
   deepCopy(SavedNormVel, NormVel);

   // Interval 1: a single step, so the average is exactly the input.
   // Non-zero velocity so interval 2 can detect a missing reset.
   const Real TempBase  = static_cast<Real>(TempIdx);
   const Real SalinBase = static_cast<Real>(SalinIdx);
   setSfcTracers(NCells, TempBase, SalinBase);
   deepCopy(NormVel, 1.0_Real);

   DefCoupling->updateExportFields(DefState, Tracers::getAll(0), DefFrazil);

   // SSH is instantaneous: whatever SshCell holds at export is sent
   HostArray1DReal ExpectedSsh =
       makeCellVarryingArray("ExpectedSsh", NCells, Real(SshIdx));
   auto SshCellOwned =
       Kokkos::subview(DefVertCoord->SshCell, std::pair(0, NCells));
   deepCopy(SshCellOwned, ExpectedSsh);

   DefCoupling->exportToCoupler();

   // Check 1: values are packed into the right slots and converted to coupler
   // units (ConstantEos makes CT->PT the identity).
   // Velocity recon correctness is tested in HorzOperatorsTest.
   // TODO: Add end-to-end SSH-gradient accumulation/export coverage if we
   // decide to continue passing ssh grad (cf. ssh directly) to mpas-si
   int PackErr = 0;
   for (int Cell = 0; Cell < NCells; Cell++) {
      if (exported(Cell, TempIdx) != TempBase + Cell + TkFrz)
         PackErr++;
      if (exported(Cell, SalinIdx) != (SalinBase + Cell) / Psu2Gpkg)
         PackErr++;
      if (exported(Cell, SshIdx) != ExpectedSsh(Cell))
         PackErr++;
   }

   if (PackErr == 0) {
      LOG_INFO("SfcCouplingTest: exportToCoupler with {} layout PASS",
               toString(Layout));
   } else {
      Err += PackErr;
      LOG_ERROR("SfcCouplingTest: exportToCoupler with {} layout FAIL - "
                "{} packing errors",
                toString(Layout), PackErr);
   }

   // Interval 2: one step with new values. If the accumulators were reset by
   // the first export, the export is exactly the new values; otherwise it is
   // mixed with interval 1.
   const Real Shift = 100.0_Real;
   setSfcTracers(NCells, TempBase + Shift, SalinBase + Shift);
   deepCopy(NormVel, 0.0_Real);

   DefCoupling->updateExportFields(DefState, Tracers::getAll(0), DefFrazil);
   DefCoupling->exportToCoupler();

   // Check 2: nothing from interval 1 leaks into interval 2
   int ResetErr = 0;
   for (int Cell = 0; Cell < NCells; Cell++) {
      if (exported(Cell, TempIdx) != TempBase + Shift + Cell + TkFrz)
         ResetErr++;
      if (exported(Cell, SalinIdx) != (SalinBase + Shift + Cell) / Psu2Gpkg)
         ResetErr++;
      if (exported(Cell, UIdx) != 0.0_Real)
         ResetErr++;
      if (exported(Cell, VIdx) != 0.0_Real)
         ResetErr++;
   }

   if (ResetErr == 0) {
      LOG_INFO("SfcCouplingTest: exportToCoupler reset with {} layout PASS",
               toString(Layout));
   } else {
      Err += ResetErr;
      LOG_ERROR("SfcCouplingTest: exportToCoupler reset with {} layout FAIL "
                "- {} errors",
                toString(Layout), ResetErr);
   }

   deepCopy(NormVel, SavedNormVel);

   SfcCoupling::clear();

   return Err;
}

int testEraseAndGet() {

   int Err = 0;

   TestSetup Setup;
   auto *DefMesh         = HorzMesh::getDefault();
   auto *DefStepper      = TimeStepper::getDefault();
   TimeInterval TimeStep = DefStepper->getTimeStep();

   // test creation of a non-default, named surface coupling object
   SfcCoupling::create("AnotherSfcCoupling", DefMesh, Setup.NImports,
                       Setup.NExports, Setup.ImportIdxMap, Setup.ExportIdxMap,
                       DefStepper, TimeStep, CouplingLayout::MCT);

   if (SfcCoupling::get("AnotherSfcCoupling")) {
      LOG_INFO("SfcCouplingTest: Non-default SfcCoupling retrieval PASS");
   } else {
      Err++;
      LOG_ERROR("SfcCouplingTest: Non-default SfcCoupling retrieval FAIL");
   }

   // test erase
   SfcCoupling::erase("AnotherSfcCoupling");

   // get() logs an expected error here; wrap with spdlog::set_level to
   // silence it if the message becomes confusing in test logs
   if (SfcCoupling::get("AnotherSfcCoupling")) {
      Err++;
      LOG_ERROR("SfcCouplingTest: erased SfcCoupling retrieval FAIL");
   } else {
      LOG_INFO("SfcCouplingTest: erased SfcCoupling retrieval PASS");
   }

   return Err;
}

void finalizeSfcCouplingTest() {

   AuxiliaryState::clear();
   VertAdv::clear();
   Tracers::clear();
   Eos::destroyInstance();
   Frazil::clear();
   Forcing::clear();
   OceanState::clear();
   VertCoord::clear();
   HorzMesh::clear();
   Field::clear();
   Dimension::clear();
   TimeStepper::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

int sfcCouplingTest(const std::string &MeshFile = "OmegaMesh.nc") {

   int Err = initSfcCouplingTest(MeshFile);

   if (Err != 0) {
      LOG_CRITICAL("SfcCouplingTest: Error initializing");
   }

   Err += testImportFromCoupler(CouplingLayout::MCT);
   Err += testImportFromCoupler(CouplingLayout::MOAB);

   Err += testUpdateExportFields(1);
   Err += testUpdateExportFields(5);

   Err += testExportToCoupler(CouplingLayout::MCT);
   Err += testExportToCoupler(CouplingLayout::MOAB);

   Err += testEraseAndGet();

   if (Err == 0) {
      LOG_INFO("SfcCouplingTest: Successful completion");
   }

   finalizeSfcCouplingTest();

   return Err;
}

int main(int argc, char *argv[]) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   RetVal += sfcCouplingTest();

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;
   if (RetVal == 0)
      LOG_INFO("------ SufaceCoupling unit tests successful ------");

   return RetVal;
} // end of main
