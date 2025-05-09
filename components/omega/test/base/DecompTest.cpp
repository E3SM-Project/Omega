//===-- Test driver for OMEGA Decomp -----------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA domain decomposition (Decomp)
///
/// This driver tests the OMEGA domain decomposition, decomposing the
/// horizontal domain and creating a number of index-space arrays for
/// locating and describing mesh locations within a parallel distributed
/// memory.
///
//
//===-----------------------------------------------------------------------===/

#include "Decomp.h"
#include "Config.h"
#include "DataTypes.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "Pacer.h"
#include "mpi.h"

#include <iostream>

using namespace OMEGA;

//------------------------------------------------------------------------------
// The initialization routine for Decomp testing. It calls various
// init routines, including the creation of the default decomposition.

int initDecompTest() {

   int Err = 0;

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv. Then retrieve the default environment and
   // some needed data members.
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv = MachEnv::getDefault();
   MPI_Comm DefComm       = DefEnv->getComm();

   initLogging(DefEnv);

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   // Initialize the IO system
   Err = IO::init(DefComm);
   if (Err != 0)
      LOG_ERROR("DecompTest: error initializing parallel IO");

   // Create the default decomposition (initializes the decomposition)
   Decomp::init();

   return Err;
}

//------------------------------------------------------------------------------
// The test driver for Decomp. This tests the decomposition of a sample
// horizontal domain and verifies the mesh is decomposed correctly.
//
int main(int argc, char *argv[]) {

   int RetVal = 0;

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
   {
      // Call initialization routine to create the default decomposition
      int Err = initDecompTest();
      if (Err != 0)
         LOG_CRITICAL("DecompTest: Error initializing");

      // Get MPI vars if needed
      MachEnv *DefEnv = MachEnv::getDefault();
      MPI_Comm Comm   = DefEnv->getComm();
      I4 MyTask       = DefEnv->getMyTask();
      I4 NumTasks     = DefEnv->getNumTasks();
      bool IsMaster   = DefEnv->isMasterTask();

      // Test retrieval of the default decomposition
      Decomp *DefDecomp = Decomp::getDefault();
      if (DefDecomp) { // true if non-null ptr
         LOG_INFO("DecompTest: Default decomp retrieval PASS");
      } else {
         LOG_INFO("DecompTest: Default decomp retrieval FAIL");
         return -1;
      }

      // Test that all Cells, Edges, Vertices are accounted for by
      // summing the list of owned values by all tasks. The result should
      // be the sum of the integers from 1 to NCellsGlobal (or edges, vertices).
      I4 RefSumCells    = 0;
      I4 RefSumEdges    = 0;
      I4 RefSumVertices = 0;
      for (int n = 0; n < DefDecomp->NCellsGlobal; ++n)
         RefSumCells += n + 1;
      for (int n = 0; n < DefDecomp->NEdgesGlobal; ++n)
         RefSumEdges += n + 1;
      for (int n = 0; n < DefDecomp->NVerticesGlobal; ++n)
         RefSumVertices += n + 1;
      I4 LocSumCells          = 0;
      I4 LocSumEdges          = 0;
      I4 LocSumVertices       = 0;
      HostArray1DI4 CellIDH   = DefDecomp->CellIDH;
      HostArray1DI4 EdgeIDH   = DefDecomp->EdgeIDH;
      HostArray1DI4 VertexIDH = DefDecomp->VertexIDH;
      for (int n = 0; n < DefDecomp->NCellsOwned; ++n)
         LocSumCells += CellIDH(n);
      for (int n = 0; n < DefDecomp->NEdgesOwned; ++n)
         LocSumEdges += EdgeIDH(n);
      for (int n = 0; n < DefDecomp->NVerticesOwned; ++n)
         LocSumVertices += VertexIDH(n);
      I4 SumCells    = 0;
      I4 SumEdges    = 0;
      I4 SumVertices = 0;
      Err =
          MPI_Allreduce(&LocSumCells, &SumCells, 1, MPI_INT32_T, MPI_SUM, Comm);
      Err =
          MPI_Allreduce(&LocSumEdges, &SumEdges, 1, MPI_INT32_T, MPI_SUM, Comm);
      Err = MPI_Allreduce(&LocSumVertices, &SumVertices, 1, MPI_INT32_T,
                          MPI_SUM, Comm);

      if (SumCells == RefSumCells) {
         LOG_INFO("DecompTest: Sum cell ID test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("DecompTest: Sum cell ID test FAIL {} {}", SumCells,
                  RefSumCells);
      }
      if (SumEdges == RefSumEdges) {
         LOG_INFO("DecompTest: Sum edge ID test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("DecompTest: Sum edge ID test FAIL {} {}", SumEdges,
                  RefSumEdges);
      }
      if (SumVertices == RefSumVertices) {
         LOG_INFO("DecompTest: Sum vertex ID test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("DecompTest: Sum vertex ID test FAIL {} {}", SumVertices,
                  RefSumVertices);
      }

      // Clean up
      Decomp::clear();
      MachEnv::removeAll();

      if (Err == 0)
         LOG_INFO("DecompTest: Successful completion");
   }
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
