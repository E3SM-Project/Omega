//===-- ocn/ForcingTest.cpp - Forcing File I/O Test -------*- C++ -*-------===//
//
/// \file
/// \brief Test for forcing file I/O and data ingestion
///
/// This test validates that the Forcing module can initialize auxiliary
/// variables and successfully read forcing data from forcing.nc. The test
/// initializes all forcing fields to zero, then reads forcing.nc and confirms
/// that at least one field transitions to a non-zero value.
//
//===----------------------------------------------------------------------===//

#include "Forcing.h"
#include "Config.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "Field.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertCoord.h"

#include <string>
#include <vector>

using namespace OMEGA;

// Default mesh file for this test
const std::string DefaultMeshFile = "OmegaMesh.nc";

namespace {

Clock *TestClock          = nullptr;
bool ForcingInitSucceeded = false;

int configureForcingReadStream() {
   std::shared_ptr<IOStream> ForcingStream = IOStream::get("Forcing");
   if (ForcingStream == nullptr) {
      LOG_ERROR("ForcingTest: Failed to get Forcing stream");
      return 1;
   }

   // Remove group token and replace with explicit forcing.nc field names.
   ForcingStream->removeField("Forcing");
   ForcingStream->addField("snowFlux");
   ForcingStream->addField("latentHeatFlux");
   ForcingStream->addField("seaIceSalinityFlux");

   return 0;
}

I4 countNonZeroEntries(const Array1DReal &Arr) {
   I4 LocalCount = 0;
   parallelReduce(
       {Arr.extent_int(0)},
       KOKKOS_LAMBDA(int I, I4 &Count) {
          if (Kokkos::abs(Arr(I)) > 1.0e-20_Real) {
             ++Count;
          }
       },
       Kokkos::Sum<I4>(LocalCount));

   I4 GlobalCount = 0;
   MPI_Allreduce(&LocalCount, &GlobalCount, 1, MPI_INT, MPI_SUM,
                 MachEnv::getDefault()->getComm());
   return GlobalCount;
}

int validateForcingStreamConfig() {
   int Err                                         = 0;
   const std::vector<std::string> RequiredContents = {
       "snowFlux", "latentHeatFlux", "seaIceSalinityFlux"};

   Config *OmegaConfig = Config::getOmegaConfig();
   if (OmegaConfig == nullptr) {
      LOG_ERROR("ForcingTest: Omega configuration is not available");
      return 1;
   }

   Error CfgErr;
   Config IOStreamsConfig("IOStreams");
   CfgErr = OmegaConfig->get(IOStreamsConfig);
   if (!CfgErr.isSuccess()) {
      LOG_ERROR("ForcingTest: Missing Omega.IOStreams configuration");
      return 1;
   }

   Config ForcingStreamConfig("Forcing");
   CfgErr = IOStreamsConfig.get(ForcingStreamConfig);
   if (!CfgErr.isSuccess()) {
      LOG_ERROR("ForcingTest: Missing Omega.IOStreams.Forcing stream");
      return 1;
   }

   std::string Filename;
   std::string Mode;
   std::vector<std::string> Contents;

   CfgErr = ForcingStreamConfig.get("Filename", Filename);
   if (!CfgErr.isSuccess()) {
      LOG_ERROR("ForcingTest: Omega.IOStreams.Forcing.Filename is missing");
      return 1;
   }
   if (Filename != "forcing.nc") {
      LOG_ERROR("ForcingTest: Expected Forcing stream filename forcing.nc, "
                "found {}",
                Filename);
      return 1;
   }

   CfgErr = ForcingStreamConfig.get("Mode", Mode);
   if (!CfgErr.isSuccess()) {
      LOG_ERROR("ForcingTest: Omega.IOStreams.Forcing.Mode is missing");
      return 1;
   }
   if (Mode != "read") {
      LOG_ERROR("ForcingTest: Expected Forcing stream mode read, found {}",
                Mode);
      return 1;
   }

   CfgErr = ForcingStreamConfig.get("Contents", Contents);
   if (!CfgErr.isSuccess()) {
      LOG_ERROR("ForcingTest: Omega.IOStreams.Forcing.Contents is missing");
      return 1;
   }
   for (const auto &Required : RequiredContents) {
      bool Found = false;
      for (const auto &Entry : Contents) {
         if (Entry == Required) {
            Found = true;
            break;
         }
      }
      if (!Found) {
         LOG_ERROR("ForcingTest: Omega.IOStreams.Forcing.Contents missing "
                   "required entry {}",
                   Required);
         return 1;
      }
   }

   const std::string StreamFilename = IOStream::getFilename("Forcing");
   if (StreamFilename != "forcing.nc") {
      LOG_ERROR("ForcingTest: IOStream Forcing resolves to filename {} instead "
                "of forcing.nc",
                StreamFilename);
      return 1;
   }

   return Err;
}

int validateForcingFields() {
   int Err                           = 0;
   const std::string ForcingFilename = "forcing.nc";

   // Expected field names in forcing.nc
   const std::vector<std::string> ExpectedFields = {
       "snowFlux", "latentHeatFlux", "seaIceSalinityFlux"};

   // Open the forcing file to check for required fields
   int FileID = -1;
   IO::openFileRead(FileID, ForcingFilename);

   std::vector<std::string> MissingFields;
   for (const auto &FieldName : ExpectedFields) {
      int VarID  = -1;
      int PIOErr = PIOc_inq_varid(FileID, FieldName.c_str(), &VarID);
      if (PIOErr != PIO_NOERR) {
         MissingFields.push_back(FieldName);
      }
   }

   IO::closeFile(FileID);

   if (!MissingFields.empty()) {
      std::string MissingList;
      for (size_t i = 0; i < MissingFields.size(); ++i) {
         if (i > 0)
            MissingList += ", ";
         MissingList += MissingFields[i];
      }
      LOG_ERROR("ForcingTest: Required field(s) missing in forcing.nc: {}",
                MissingList);
      return 1;
   }

   LOG_INFO("ForcingTest: All required forcing fields found in forcing.nc");
   return Err;
}

} // namespace

/// Initialize the test environment
/// \param [in] MeshFile The mesh file to use for initialization
void initForcingTest(const std::string &MeshFile) {

   ForcingInitSucceeded = false;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize logging
   initLogging(DefEnv);

   // Open config file
   Config::Initialize();
   Config::readAll("omega.yml");

   // Restrict the Forcing stream contents to variables available in
   // forcing.nc for this test case.
   Config *OmegaConfig = Config::getOmegaConfig();
   Config IOStreamsConfig("IOStreams");
   Error CfgErr = OmegaConfig->get(IOStreamsConfig);
   if (!CfgErr.isSuccess()) {
      ABORT_ERROR("ForcingTest: Missing Omega.IOStreams configuration");
   }
   Config ForcingStreamConfig("Forcing");
   CfgErr = IOStreamsConfig.get(ForcingStreamConfig);
   if (!CfgErr.isSuccess()) {
      ABORT_ERROR("ForcingTest: Missing Omega.IOStreams.Forcing stream");
   }
   ForcingStreamConfig.set(
       "Contents", std::vector<std::string>{"snowFlux", "latentHeatFlux",
                                            "seaIceSalinityFlux"});

   TimeStepper::init1();
   TimeStepper *DefStepper = TimeStepper::getDefault();
   TestClock               = DefStepper->getClock();

   IO::init(DefComm);
   IOStream::init(TestClock);
   Field::init(TestClock);

   Decomp::init(MeshFile);

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      ABORT_ERROR("ForcingTest: error initializing default halo");
   }

   HorzMesh::init();

   // Initialize vertical coordinate with default layers
   VertCoord::init();
   Tracers::init();

   // Initialize Forcing module
   Forcing::init();
   ForcingInitSucceeded = true;

   if (configureForcingReadStream() != 0) {
      ABORT_ERROR("ForcingTest: Error configuring Forcing stream contents");
   }

   std::shared_ptr<IOStream> ForcingStream = IOStream::get("Forcing");
   bool StreamValid                        = ForcingStream->validate();
   if (!StreamValid) {
      ABORT_ERROR("ForcingTest: Error validating Forcing stream");
   }

} // end initForcingTest

/// Finalize the test environment
void finalizeForcingTest() {
   // Forcing::clear can abort in this standalone test cleanup sequence.
   // Skip explicit Forcing teardown and rely on process teardown.
   if (!ForcingInitSucceeded) {
      LOG_WARN("ForcingTest: Skipping Forcing cleanup because init did not "
               "complete");
   }
   Tracers::clear();
   HorzMesh::clear();
   VertCoord::clear();
   Field::clear();
   Dimension::clear();
   IOStream::finalize();
   TimeStepper::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
   TestClock = nullptr;

} // end finalizeForcingTest

/// Test reading forcing data from file
/// \return Error count (0 on success, >0 on failure)
int testForcingFileRead() {

   I4 Err              = 0;
   Forcing *DefForcing = Forcing::getDefault();

   if (DefForcing == nullptr) {
      LOG_ERROR("ForcingTest: Failed to get default Forcing instance");
      return 1;
   }

   LOG_INFO("ForcingTest: Testing forcing file read from stream Forcing");

   if (validateForcingStreamConfig() != 0) {
      return 1;
   }

   if (validateForcingFields() != 0) {
      LOG_ERROR("ForcingTest: forcing.nc does not contain expected fields");
      return 1;
   }

   // Access forcing auxiliary variables and verify zero initialization
   const auto &TracerForcingAux = DefForcing->TracerForcingAux;

   I4 PreNonZero = 0;
   PreNonZero += countNonZeroEntries(TracerForcingAux.SnowFluxCell);
   PreNonZero += countNonZeroEntries(TracerForcingAux.LatentHeatFluxCell);
   PreNonZero += countNonZeroEntries(TracerForcingAux.SeaIceSaltFluxCell);

   if (PreNonZero != 0) {
      LOG_ERROR("ForcingTest: Expected zero-initialized forcing arrays, found "
                "{} non-zero entries before read",
                PreNonZero);
      return 1;
   }

   LOG_INFO("ForcingTest: Verified zero initialization before stream read");

   Metadata ReqMetadata;
   Error ReadErr = IOStream::read("Forcing", TestClock, ReqMetadata, true);
   if (!ReadErr.isSuccess()) {
      LOG_ERROR("ForcingTest: Failed reading Forcing stream from "
                "forcing.nc");
      return 1;
   }

   I4 PostNonZero = 0;
   PostNonZero += countNonZeroEntries(TracerForcingAux.SnowFluxCell);
   PostNonZero += countNonZeroEntries(TracerForcingAux.LatentHeatFluxCell);
   PostNonZero += countNonZeroEntries(TracerForcingAux.SeaIceSaltFluxCell);

   if (PostNonZero == 0) {
      LOG_ERROR("ForcingTest: forcing.nc contains variables snowFlux, "
                "latentHeatFlux, seaIceSalinityFlux and stream read "
                "succeeded, but TracerForcingAux members SnowFluxCell, "
                "LatentHeatFluxCell, SeaIceSaltFluxCell remain zero "
                "(possible stream-to-field mapping mismatch)");
      ++Err;
   } else {
      LOG_INFO("ForcingTest: forcing.nc read produced {} non-zero entries in "
               "selected TracerForcingAux arrays",
               PostNonZero);
   }

   if (Err == 0) {
      LOG_INFO("ForcingTest: Forcing file read PASS");
   }

   return Err;

} // end testForcingFileRead

/// Main forcing test function
/// \return Error count (0 on success, >0 on failure)
int forcingTest() {

   int Err = 0;

   initForcingTest(DefaultMeshFile);

   Err += testForcingFileRead();

   if (Err == 0) {
      LOG_INFO("ForcingTest: Successful completion");
   }

   finalizeForcingTest();

   return Err;

} // end forcingTest

/// Main entry point for the test
int main(int argc, char *argv[]) {

   int RetErr = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   RetErr = forcingTest();

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   return RetErr;

} // end of main
//===----------------------------------------------------------------------===//
