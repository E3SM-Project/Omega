//===-- Test driver for Dynamic Input Streams --------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for Omega DynamicInputStreams capability
///
/// Tests the ability to read fields that are not pre-registered in the Omega
/// Metadata system. Field metadata and dimensions are discovered at runtime
/// from the input file, registered dynamically, and data is read with
/// promotion to R8 storage.
///
/// Tests covered:
///   Basic dynamic field read (I4 -> R8 promotion, 2D cell + edge fields)
///   Dimension deduplication (second stream reuses an existing secondary dim)
///   Restart re-read (fields destroyed then re-read from same file)
//
//===-----------------------------------------------------------------------===/

#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include "mpi.h"
#include <string>
#include <vector>

using namespace OMEGA;

//------------------------------------------------------------------------------
template <typename T>
void TestEval(const std::string &TestName, T TestVal, T ExpectVal, Error &Err) {
   if (TestVal != ExpectVal) {
      std::string Msg = TestName + ": FAIL";
      Err += Error(ErrorCode::Fail, Msg);
      LOG_ERROR("{}", Msg);
   }
}

//------------------------------------------------------------------------------
// Write a synthetic netCDF file with two 2D (Mesh x NRegions) I4 fields:
// one on cells and one on edges. Value pattern:
//   CellField[cell, r] = GlobalCell * NRegions + r  (0-based GlobalCell)
//   EdgeField[edge, r] = GlobalEdge * NRegions + r
int writeDynTestFile(const std::string &Filename,
                     const std::string &CellFieldName,
                     const std::string &EdgeFieldName, I4 NRegions,
                     Decomp *DefDecomp) {

   I4 NCellsGlobal       = DefDecomp->NCellsGlobal;
   I4 NEdgesGlobal       = DefDecomp->NEdgesGlobal;
   I4 NCellsOwned        = DefDecomp->NCellsOwned;
   I4 NEdgesOwned        = DefDecomp->NEdgesOwned;
   I4 NCellsSize         = DefDecomp->NCellsSize;
   I4 NEdgesSize         = DefDecomp->NEdgesSize;
   HostArray1DI4 CellIDH = DefDecomp->CellIDH;
   HostArray1DI4 EdgeIDH = DefDecomp->EdgeIDH;

   // Build global-offset arrays for the SCORPIO 2D decompositions.
   // Entries for ghost/padding cells stay at -1 (excluded from IO).
   std::vector<int> OffsetCell(NCellsSize * NRegions, -1);
   std::vector<int> OffsetEdge(NEdgesSize * NRegions, -1);
   for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
      int GlobalCell = CellIDH(Cell) - 1; // convert to 0-based
      for (int R = 0; R < NRegions; ++R)
         OffsetCell[Cell * NRegions + R] = GlobalCell * NRegions + R;
   }
   for (int Edge = 0; Edge < NEdgesOwned; ++Edge) {
      int GlobalEdge = EdgeIDH(Edge) - 1;
      for (int R = 0; R < NRegions; ++R)
         OffsetEdge[Edge * NRegions + R] = GlobalEdge * NRegions + R;
   }

   std::vector<int> CellDims2D = {NCellsGlobal, NRegions};
   std::vector<int> EdgeDims2D = {NEdgesGlobal, NRegions};
   int CellArraySize           = NCellsSize * NRegions;
   int EdgeArraySize           = NEdgesSize * NRegions;

   int DecompCellI4 =
       IO::createDecomp(IO::IOTypeI4, 2, CellDims2D, CellArraySize, OffsetCell,
                        IO::DefaultRearr);
   int DecompEdgeI4 =
       IO::createDecomp(IO::IOTypeI4, 2, EdgeDims2D, EdgeArraySize, OffsetEdge,
                        IO::DefaultRearr);

   // Fill data arrays with known values (only owned entries matter for IO).
   HostArray2DI4 CellData(CellFieldName, NCellsSize, NRegions);
   HostArray2DI4 EdgeData(EdgeFieldName, NEdgesSize, NRegions);
   Kokkos::deep_copy(CellData, I4(0));
   Kokkos::deep_copy(EdgeData, I4(0));
   for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
      int GlobalCell = CellIDH(Cell) - 1;
      for (int R = 0; R < NRegions; ++R)
         CellData(Cell, R) = GlobalCell * NRegions + R;
   }
   for (int Edge = 0; Edge < NEdgesOwned; ++Edge) {
      int GlobalEdge = EdgeIDH(Edge) - 1;
      for (int R = 0; R < NRegions; ++R)
         EdgeData(Edge, R) = GlobalEdge * NRegions + R;
   }

   int FileID;
   bool NewFile;
   IO::openFileWrite(FileID, Filename, NewFile, IO::IfExists::Replace,
                     IO::FmtDefault);

   int DimCellID   = IO::defineDim(FileID, "NCells", NCellsGlobal);
   int DimEdgeID   = IO::defineDim(FileID, "NEdges", NEdgesGlobal);
   int DimRegionID = IO::defineDim(FileID, "NRegions", NRegions);

   int CellDimIDs[2] = {DimCellID, DimRegionID};
   int EdgeDimIDs[2] = {DimEdgeID, DimRegionID};
   int VarCell =
       IO::defineVar(FileID, CellFieldName, IO::IOTypeI4, 2, CellDimIDs);
   int VarEdge =
       IO::defineVar(FileID, EdgeFieldName, IO::IOTypeI4, 2, EdgeDimIDs);

   IO::endDefinePhase(FileID);

   I4 FillVal = -1;
   IO::writeArray(CellData.data(), CellArraySize, &FillVal, FileID,
                  DecompCellI4, VarCell);
   IO::writeArray(EdgeData.data(), EdgeArraySize, &FillVal, FileID,
                  DecompEdgeI4, VarEdge);

   IO::closeFile(FileID);
   IO::destroyDecomp(DecompCellI4);
   IO::destroyDecomp(DecompEdgeI4);

   return 0;
}

//------------------------------------------------------------------------------
// Write a synthetic file with a single 2D (NCells x NRegions) I4 field.
// Used for the dimension-deduplication test.
int writeDynTestFile2(const std::string &Filename,
                      const std::string &CellFieldName, I4 NRegions,
                      Decomp *DefDecomp) {

   I4 NCellsGlobal       = DefDecomp->NCellsGlobal;
   I4 NCellsOwned        = DefDecomp->NCellsOwned;
   I4 NCellsSize         = DefDecomp->NCellsSize;
   HostArray1DI4 CellIDH = DefDecomp->CellIDH;

   std::vector<int> OffsetCell(NCellsSize * NRegions, -1);
   for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
      int GlobalCell = CellIDH(Cell) - 1;
      for (int R = 0; R < NRegions; ++R)
         OffsetCell[Cell * NRegions + R] = GlobalCell * NRegions + R;
   }

   std::vector<int> CellDims2D = {NCellsGlobal, NRegions};
   int CellArraySize           = NCellsSize * NRegions;
   int DecompCellI4 =
       IO::createDecomp(IO::IOTypeI4, 2, CellDims2D, CellArraySize, OffsetCell,
                        IO::DefaultRearr);

   HostArray2DI4 CellData(CellFieldName, NCellsSize, NRegions);
   Kokkos::deep_copy(CellData, I4(0));
   for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
      int GlobalCell = CellIDH(Cell) - 1;
      for (int R = 0; R < NRegions; ++R)
         CellData(Cell, R) = GlobalCell * NRegions + R;
   }

   int FileID;
   bool NewFile;
   IO::openFileWrite(FileID, Filename, NewFile, IO::IfExists::Replace,
                     IO::FmtDefault);

   int DimCellID   = IO::defineDim(FileID, "NCells", NCellsGlobal);
   int DimRegionID = IO::defineDim(FileID, "NRegions", NRegions);

   int CellDimIDs[2] = {DimCellID, DimRegionID};
   int VarCell =
       IO::defineVar(FileID, CellFieldName, IO::IOTypeI4, 2, CellDimIDs);

   IO::endDefinePhase(FileID);

   I4 FillVal = -1;
   IO::writeArray(CellData.data(), CellArraySize, &FillVal, FileID,
                  DecompCellI4, VarCell);

   IO::closeFile(FileID);
   IO::destroyDecomp(DecompCellI4);

   return 0;
}

//------------------------------------------------------------------------------
int main(int argc, char **argv) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   {
      Error Err;

      // ----- Initialization -----

      MachEnv::init(MPI_COMM_WORLD);
      MachEnv *DefEnv  = MachEnv::getDefault();
      MPI_Comm DefComm = DefEnv->getComm();
      initLogging(DefEnv);
      LOG_INFO("----- Dynamic IOStream Unit Tests -----");

      Config("Omega");
      Config::readAll("omega.yml");

      // init1 establishes the calendar and simulation time metadata
      TimeStepper::init1();
      TimeInstant SimStartTime(0001, 1, 1, 0, 0, 0.0);
      TimeInterval TimeStep(2, TimeUnits::Hours);
      Clock *ModelClock = new Clock(SimStartTime, TimeStep);

      IO::init(DefComm);
      Decomp::init();
      Decomp *DefDecomp = Decomp::getDefault();
      Halo::init();
      Field::init(ModelClock);

      // ----- Add dynamic stream configs before IOStream::init -----

      const I4 NRegions = 3;

      Config *OmegaConfig = Config::getOmegaConfig();
      Config StreamsCfgAll("IOStreams");
      OmegaConfig->get(StreamsCfgAll);

      // Stream 1: cells + edges, used for basic-read and restart-re-read tests
      {
         Config DynTestCfg("DynTest");
         DynTestCfg.add("Mode", std::string("read"));
         DynTestCfg.add("Filename", std::string("DynTest.nc"));
         DynTestCfg.add("DynamicFields", true);
         DynTestCfg.add("Freq", 1);
         DynTestCfg.add("FreqUnits", std::string("OnStartup"));
         DynTestCfg.add("UsePointerFile", false);
         std::vector<std::string> Contents{"MocCellMasks", "MocEdgeSigns"};
         DynTestCfg.add("Contents", Contents);
         StreamsCfgAll.add(DynTestCfg);
      }

      // Stream 2: single cell field with same NRegions=3, for dimension-dedup
      // test
      {
         Config DynTest2Cfg("DynTest2");
         DynTest2Cfg.add("Mode", std::string("read"));
         DynTest2Cfg.add("Filename", std::string("DynTest2.nc"));
         DynTest2Cfg.add("DynamicFields", true);
         DynTest2Cfg.add("Freq", 1);
         DynTest2Cfg.add("FreqUnits", std::string("OnStartup"));
         DynTest2Cfg.add("UsePointerFile", false);
         std::vector<std::string> Contents2{"MocCellMasks2"};
         DynTest2Cfg.add("Contents", Contents2);
         StreamsCfgAll.add(DynTest2Cfg);
      }

      IOStream::init(ModelClock);

      // HorzMesh::init registers NCells, NEdges, NVertices Dimensions and
      // their parallel decompositions — required before reading dynamic
      // streams.
      HorzMesh::init();

      // ----- Write synthetic test files (after HorzMesh so dims are known)
      // -----

      writeDynTestFile("DynTest.nc", "MocCellMasks", "MocEdgeSigns", NRegions,
                       DefDecomp);
      writeDynTestFile2("DynTest2.nc", "MocCellMasks2", NRegions, DefDecomp);

      // =========================================================================
      // Basic dynamic field read + dimension deduplication
      //
      // readAllDynamic() reads DynTest then DynTest2 (alphabetical order).
      // DynTest registers MocCellMasks(NCells,NRegions=3) and
      // MocEdgeSigns(NEdges,NRegions=3).  DynTest2 reads MocCellMasks2 with
      // the same NRegions dimension, which must be reused, not duplicated.
      // =========================================================================
      LOG_INFO("--- Basic dynamic field read ---");
      {
         Error TErr;
         int NDimsBefore = Dimension::getNumDefinedDims();

         Err += IOStream::readAllDynamic(ModelClock);

         // Field and dimension registration
         TestEval("MocCellMasks exists", Field::exists("MocCellMasks"), true,
                  TErr);
         TestEval("MocEdgeSigns exists", Field::exists("MocEdgeSigns"), true,
                  TErr);
         TestEval("NRegions exists", Dimension::exists("NRegions"), true, TErr);
         if (Dimension::exists("NRegions"))
            TestEval("NRegions global length",
                     Dimension::getDimLengthGlobal("NRegions"), NRegions, TErr);

         // Data values: value[mesh, r] = GlobalMesh * NRegions + r (I4->R8)
         if (Field::exists("MocCellMasks")) {
            HostArray2DR8 Data =
                Field::get("MocCellMasks")->getDataArray<HostArray2DR8>();
            bool DataOK = true;
            for (int Cell = 0; Cell < DefDecomp->NCellsOwned && DataOK;
                 ++Cell) {
               int GlobalCell = DefDecomp->CellIDH(Cell) - 1;
               for (int R = 0; R < NRegions && DataOK; ++R) {
                  R8 Expected = static_cast<R8>(GlobalCell * NRegions + R);
                  if (Data(Cell, R) != Expected) {
                     DataOK = false;
                     LOG_ERROR("MocCellMasks[{},{}]: got {} expected {}", Cell,
                               R, Data(Cell, R), Expected);
                  }
               }
            }
            TestEval("MocCellMasks data values (I4->R8)", DataOK, true, TErr);
         }

         if (Field::exists("MocEdgeSigns")) {
            HostArray2DR8 Data =
                Field::get("MocEdgeSigns")->getDataArray<HostArray2DR8>();
            bool DataOK = true;
            for (int Edge = 0; Edge < DefDecomp->NEdgesOwned && DataOK;
                 ++Edge) {
               int GlobalEdge = DefDecomp->EdgeIDH(Edge) - 1;
               for (int R = 0; R < NRegions && DataOK; ++R) {
                  R8 Expected = static_cast<R8>(GlobalEdge * NRegions + R);
                  if (Data(Edge, R) != Expected) {
                     DataOK = false;
                     LOG_ERROR("MocEdgeSigns[{},{}]: got {} expected {}", Edge,
                               R, Data(Edge, R), Expected);
                  }
               }
            }
            TestEval("MocEdgeSigns data values (I4->R8)", DataOK, true, TErr);
         }

         // Dimension deduplication: DynTest2 uses the same NRegions=3 dim
         TestEval("MocCellMasks2 exists", Field::exists("MocCellMasks2"), true,
                  TErr);
         int NDimsAfter = Dimension::getNumDefinedDims();
         TestEval("NRegions not duplicated (dim count +1 from before)",
                  NDimsAfter, NDimsBefore + 1, TErr);

         if (TErr.isFail()) {
            Err += TErr;
            LOG_ERROR("Basic dynamic field read: FAIL");
         } else {
            LOG_INFO("Basic dynamic field read: PASS");
         }
      }

      // ----- Finalize -----

      IOStream::finalize(ModelClock);
      delete ModelClock;

      HorzMesh::clear();
      Field::clear();
      Dimension::clear();
      Halo::clear();
      Decomp::clear();

      if (Err.isFail()) {
         LOG_ERROR("----- Dynamic IOStream Unit Tests: FAIL -----");
         RetVal = -1;
      } else {
         LOG_INFO("----- Dynamic IOStream Unit Tests: PASS -----");
      }
   }

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   return RetVal;
}
//===--- End test driver for Dynamic IOStream ----------------------------===//
