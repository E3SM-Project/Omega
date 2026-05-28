//===-- Test: dynamic field with no mesh dimension ----------------*- C++
//-*-===/
//
/// \file
/// \brief Verifies that reading a dynamic stream aborts when the requested
/// field has no mesh dimension (NCells, NEdges, or NVertices).
///
/// The test file contains a 1D field indexed only on a non-mesh dimension
/// (NRegions).  The dynamic registration must detect the missing mesh
/// dimension and abort.  CTest marks this test WILL_FAIL.
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
// Write a test file with one 1D (NRegions) R8 field — no mesh dimension.
static void writeTestFile(const std::string &Filename,
                          const std::string &FieldName, I4 NRegions) {

   int FileID;
   bool NewFile;
   IO::openFileWrite(FileID, Filename, NewFile, IO::IfExists::Replace,
                     IO::FmtDefault);

   int DimRegion    = IO::defineDim(FileID, "NRegions", NRegions);
   int RegionDim[1] = {DimRegion};
   int VarID = IO::defineVar(FileID, FieldName, IO::IOTypeR8, 1, RegionDim);
   IO::endDefinePhase(FileID);

   // Non-distributed write: all ranks write the same small array.
   std::vector<R8> Data(NRegions);
   for (int R = 0; R < NRegions; ++R)
      Data[R] = static_cast<R8>(R + 1);
   IO::writeNDVar(Data.data(), FileID, VarID);

   IO::closeFile(FileID);
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
      LOG_INFO("----- Dynamic IOStream Bad Field Test -----");

      Config("Omega");
      Config::readAll("omega.yml");
      TimeStepper::init1();

      TimeInstant SimStartTime(0001, 1, 1, 0, 0, 0.0);
      TimeInterval TimeStep(2, TimeUnits::Hours);
      Clock *ModelClock = new Clock(SimStartTime, TimeStep);

      IO::init(DefComm);
      Decomp::init();
      Halo::init();
      Field::init(ModelClock);

      Config *OmegaConfig = Config::getOmegaConfig();
      Config StreamsCfg("IOStreams");
      OmegaConfig->get(StreamsCfg);

      Config BadCfg("BadFieldStream");
      BadCfg.add("Mode", std::string("read"));
      BadCfg.add("Filename", std::string("DynBadField.nc"));
      BadCfg.add("DynamicFields", true);
      BadCfg.add("Freq", 1);
      BadCfg.add("FreqUnits", std::string("OnStartup"));
      BadCfg.add("UsePointerFile", false);
      std::vector<std::string> Contents{"RegionAreas"};
      BadCfg.add("Contents", Contents);
      StreamsCfg.add(BadCfg);

      IOStream::init(ModelClock);

      // HorzMesh registers NCells/NEdges/NVertices as distributed dimensions,
      // making the absence of a mesh dim detectable.
      HorzMesh::init();

      writeTestFile("DynBadField.nc", "RegionAreas", 3);

      // Must abort: "RegionAreas" has no NCells, NEdges, or NVertices
      // dimension.
      Error DynErr = IOStream::readAllDynamic(ModelClock);
      CHECK_ERROR_ABORT(DynErr, "readAllDynamic unexpectedly succeeded");

      // Unreachable if the test behaves correctly.
      delete ModelClock;
   }

   Kokkos::finalize();
   MPI_Finalize();
   return 0;
}
//===--- End bad field test
//------------------------------------------------===//
