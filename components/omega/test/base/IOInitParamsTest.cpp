//===-- Test driver for OMEGA driver-owned IO parameters ---------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA driver-owned IO parameters
///
/// This driver tests the IO::init overload that takes IOInitParams. In a
/// coupled run the IO base task and rearranger belong to the driver, which
/// supplies them from shr_pio_getioroot and shr_pio_getrearranger, while the
/// remaining IO settings still come from the Omega configuration. This test
/// passes values that deliberately differ from the ones in the config file and
/// checks that the supplied values are the ones actually used, so that a
/// regression which silently fell back to the config values would fail here.
///
//
//===-----------------------------------------------------------------------===/

#include "Config.h"
#include "DataTypes.h"
#include "Error.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "Pacer.h"
#include "mpi.h"

#include <string>

using namespace OMEGA;

//------------------------------------------------------------------------------
// The test driver for driver-owned IO parameters.
//
int main(int argc, char *argv[]) {

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   {
      Error Err;

      // Initialize the Machine Environment class and retrieve the default
      // environment and communicator
      MachEnv::init(MPI_COMM_WORLD);
      MachEnv *DefEnv  = MachEnv::getDefault();
      MPI_Comm DefComm = DefEnv->getComm();
      I4 MyTask        = DefEnv->getMyTask();
      I4 NumTasks      = DefEnv->getNumTasks();

      // Initialize the Logging system
      initLogging(DefEnv);
      LOG_INFO("----- Driver-owned IO Parameters Unit Testing -----");

      // This test needs a task the base task can be moved to
      if (NumTasks < 2)
         ABORT_ERROR("IOInitParamsTest: FAIL test requires at least 2 tasks, "
                     "got {}",
                     NumTasks);

      // Open config file
      Config("Omega");
      Config::readAll("omega.yml");
      Config *OmegaConfig = Config::getOmegaConfig();

      // Read the IO settings the configuration file supplies. These are the
      // values IO::init would use if the supplied parameters were ignored.
      Config IOConfig("IO");
      Err = OmegaConfig->get(IOConfig);
      CHECK_ERROR_ABORT(Err, "IOInitParamsTest: FAIL IO group not found in "
                             "config");

      I4 ConfigBaseTask = 0;
      Err               = IOConfig.get("IOBaseTask", ConfigBaseTask);
      CHECK_ERROR_ABORT(Err, "IOInitParamsTest: FAIL IOBaseTask not found in "
                             "config");

      std::string ConfigRearrName = "box";
      Err = IOConfig.get("IORearranger", ConfigRearrName);
      CHECK_ERROR_ABORT(Err, "IOInitParamsTest: FAIL IORearranger not found in "
                             "config");
      IO::Rearranger ConfigRearr = IO::RearrFromString(ConfigRearrName);

      // Choose parameters that differ from the config values, so that the two
      // sources can be told apart. Without this the test would pass whether or
      // not the supplied values were honored.
      I4 DriverBaseTask = (ConfigBaseTask == 0) ? 1 : 0;
      IO::Rearranger DriverRearr =
          (ConfigRearr == IO::RearrSubset) ? IO::RearrBox : IO::RearrSubset;

      LOG_INFO(
          "IOInitParamsTest: config supplies IOBaseTask={} IORearranger={}",
          ConfigBaseTask, static_cast<int>(ConfigRearr));
      LOG_INFO(
          "IOInitParamsTest: driver supplies IOBaseTask={} IORearranger={}",
          DriverBaseTask, static_cast<int>(DriverRearr));

      // Initialize IO with the driver-owned parameters
      IO::IOInitParams DriverParams{DriverBaseTask, DriverRearr};
      IO::init(DefComm, DriverParams);

      // The rearranger actually in use is recorded in DefaultRearr
      if (IO::DefaultRearr != DriverRearr)
         ABORT_ERROR("IOInitParamsTest: FAIL rearranger is {} but the driver "
                     "supplied {}",
                     static_cast<int>(IO::DefaultRearr),
                     static_cast<int>(DriverRearr));

      if (IO::DefaultRearr == ConfigRearr)
         ABORT_ERROR("IOInitParamsTest: FAIL rearranger fell back to the "
                     "config value {}",
                     static_cast<int>(ConfigRearr));

      // Ask SCORPIO which tasks it made IO tasks. The lowest of them is the
      // base task, whatever the task count and stride happen to be, so this
      // checks the supplied base task reached PIOc_Init_Intracomm rather than
      // only being stored somewhere.
      bool IsIOTask = false;
      int PIOErr    = PIOc_iam_iotask(IO::SysID, &IsIOTask);
      if (PIOErr != PIO_NOERR)
         ABORT_ERROR("IOInitParamsTest: FAIL could not query SCORPIO IO tasks");

      I4 MyIORank  = IsIOTask ? MyTask : NumTasks;
      I4 MinIORank = NumTasks;
      MPI_Allreduce(&MyIORank, &MinIORank, 1, MPI_INT, MPI_MIN, DefComm);

      if (MinIORank == NumTasks)
         ABORT_ERROR("IOInitParamsTest: FAIL SCORPIO reports no IO tasks");

      if (MinIORank != DriverBaseTask)
         ABORT_ERROR("IOInitParamsTest: FAIL lowest SCORPIO IO task is {} but "
                     "the driver supplied base task {}",
                     MinIORank, DriverBaseTask);

      LOG_INFO("IOInitParamsTest: lowest SCORPIO IO task is {} as supplied",
               MinIORank);

      // Exit environments
      MachEnv::removeAll();

      LOG_INFO("IOInitParamsTest: Successful completion");
   }

   LOG_INFO("----- Driver-owned IO Parameters Unit Tests Successful -----");
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return 0; // if we made it here, return successfully

} // end of main
//===-----------------------------------------------------------------------===/
