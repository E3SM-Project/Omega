//===-- Test: dynamic field name collision -------------------------*- C++
//-*-===/
//
/// \file
/// \brief Verifies that reading a dynamic stream aborts when the requested
/// field name is already registered in the global field registry.
///
/// Two streams both list "MocCellMasks".  After the first stream reads and
/// registers the field, the second stream attempts to register the same name
/// and must abort with a clear error message.  CTest marks this test
/// WILL_FAIL, so the non-zero exit code counts as a pass.
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
// Write a test file with one 2D (NCells x NRegions) I4 field.
static void writeTestFile(const std::string &Filename,
                          const std::string &FieldName, I4 NRegions,
                          Decomp *DefDecomp) {

   I4 NCellsGlobal       = DefDecomp->NCellsGlobal;
   I4 NCellsOwned        = DefDecomp->NCellsOwned;
   I4 NCellsSize         = DefDecomp->NCellsSize;
   HostArray1DI4 CellIDH = DefDecomp->CellIDH;

   std::vector<int> Offset(NCellsSize * NRegions, -1);
   for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
      int GCell = CellIDH(Cell) - 1;
      for (int R = 0; R < NRegions; ++R)
         Offset[Cell * NRegions + R] = GCell * NRegions + R;
   }

   std::vector<int> Dims2D = {NCellsGlobal, NRegions};
   int ArraySize           = NCellsSize * NRegions;
   int DecompID = IO::createDecomp(IO::IOTypeI4, 2, Dims2D, ArraySize, Offset,
                                   IO::DefaultRearr);

   HostArray2DI4 Data(FieldName, NCellsSize, NRegions);
   Kokkos::deep_copy(Data, I4(0));
   for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
      int GCell = CellIDH(Cell) - 1;
      for (int R = 0; R < NRegions; ++R)
         Data(Cell, R) = GCell * NRegions + R;
   }

   int FileID;
   bool NewFile;
   IO::openFileWrite(FileID, Filename, NewFile, IO::IfExists::Replace,
                     IO::FmtDefault);
   int DimCell   = IO::defineDim(FileID, "NCells", NCellsGlobal);
   int DimRegion = IO::defineDim(FileID, "NRegions", NRegions);
   int DimIDs[2] = {DimCell, DimRegion};
   int VarID     = IO::defineVar(FileID, FieldName, IO::IOTypeI4, 2, DimIDs);
   IO::endDefinePhase(FileID);
   I4 FillVal = -1;
   IO::writeArray(Data.data(), ArraySize, &FillVal, FileID, DecompID, VarID);
   IO::closeFile(FileID);
   IO::destroyDecomp(DecompID);
}

//------------------------------------------------------------------------------
int main(int argc, char **argv) {

   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   {
      MachEnv::init(MPI_COMM_WORLD);
      MachEnv *DefEnv  = MachEnv::getDefault();
      MPI_Comm DefComm = DefEnv->getComm();
      initLogging(DefEnv);
      LOG_INFO("----- Dynamic IOStream Collision Test -----");

      Config("Omega");
      Config::readAll("omega.yml");
      TimeStepper::init1();

      TimeInstant SimStartTime(0001, 1, 1, 0, 0, 0.0);
      TimeInterval TimeStep(2, TimeUnits::Hours);
      Clock *ModelClock = new Clock(SimStartTime, TimeStep);

      IO::init(DefComm);
      Decomp::init();
      Decomp *DefDecomp = Decomp::getDefault();
      Halo::init();
      Field::init(ModelClock);

      // Add two streams that both list the same field name.
      Config *OmegaConfig = Config::getOmegaConfig();
      Config StreamsCfg("IOStreams");
      OmegaConfig->get(StreamsCfg);

      auto addStream = [&](const std::string &Name,
                           const std::string &Filename) {
         Config Cfg(Name);
         Cfg.add("Mode", std::string("read"));
         Cfg.add("Filename", Filename);
         Cfg.add("DynamicFields", true);
         Cfg.add("Freq", 1);
         Cfg.add("FreqUnits", std::string("OnStartup"));
         Cfg.add("UsePointerFile", false);
         std::vector<std::string> Contents{"MocCellMasks"};
         Cfg.add("Contents", Contents);
         StreamsCfg.add(Cfg);
      };

      addStream("CollisionStream1", std::string("DynCollision.nc"));
      addStream("CollisionStream2", std::string("DynCollision.nc"));

      IOStream::init(ModelClock);
      HorzMesh::init();

      writeTestFile("DynCollision.nc", "MocCellMasks", 3, DefDecomp);

      // readAllDynamic reads CollisionStream1 then CollisionStream2
      // (alphabetical order). CollisionStream1 succeeds and registers
      // "MocCellMasks". CollisionStream2 must abort: "MocCellMasks" is already
      // registered.
      Error DynErr = IOStream::readAllDynamic(ModelClock);
      CHECK_ERROR_ABORT(DynErr, "readAllDynamic unexpectedly succeeded");

      // Unreachable if the test behaves correctly.
      delete ModelClock;
   }

   Kokkos::finalize();
   MPI_Finalize();
   return 0;
}
//===--- End collision test
//------------------------------------------------===//
