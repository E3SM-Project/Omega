//===-- Test driver for OMEGA Halo -------------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA Halo class
///
/// This driver tests the OMEGA model Halo class, which collects and stores
/// everything needed to perform halo exchanges on any supported YAKL array
/// defined on a mesh in OMEGA with a given parallel decomposition. This
/// unit test driver tests functionality by creating YAKL arrays of every
/// type and dimensionality supported in OMEGA, initializing the array based
/// on global IDs of the mesh elememts, performing a halo exchange, and
/// confirming the exchanged array is identical to the initial array.
///
//
//===-----------------------------------------------------------------------===/

#include "DataTypes.h"
#include "Decomp.h"
#include "Halo.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "mpi.h"

#include <iostream>
#include <iomanip>

//------------------------------------------------------------------------------
// This function template performs a single test on a YAKL array type in a
// given index space. Two YAKL arrays of the same type and size are input,
// one with the array initialized based on the global IDs of the mesh elements,
// and the other uninitialized. The Halo class object, a label describing the
// test for output, an integer to accumulate errors, and optionally the index
// space of the input arrays (default is OnCell) are also input. The InitArray
// is first copied into the TestArray, and then a halo exchange is performed on
// the TestArray. Afterwards, InitArray is compared to TestArray and if any
// elements differ the test is a failure and an error is returned.

template <typename T>
void haloExchangeTest(OMEGA::Halo MyHalo,
     T InitArray,    /// Array initialized based on global IDs of mesh elements
     T TestArray,    /// Array to be exchanged and compared
     const char *Label,       /// Unique label for test
     OMEGA::I4 &TotErr,       /// Integer to track errors
     OMEGA::MeshElement ThisElem = OMEGA::OnCell   /// index space
) {

   OMEGA::I4 IErr{0};      // error code

   // Set total array size and ensure arrays are of same size
   OMEGA::I4 NTot = InitArray.totElems();
   if (NTot != TestArray.totElems()) {
      LOG_ERROR("HaloTest: {} arrays must be of same size", Label);
      TotErr += -1;
      return;
   }

   // Copy InitArray to TestArray
   InitArray.deep_copy_to(TestArray);

   // Perform halo exchange
   IErr = MyHalo.exchangeFullArrayHalo(TestArray, ThisElem);
   if (IErr != 0) {
      LOG_ERROR("HaloTest: Error during {} halo exchange", Label);
      TotErr += -1;
      return;
   }

   // Collapse arrays to 1D for easy iteration
   auto CollapsedInit = InitArray.collapse();
   auto CollapsedTest = TestArray.collapse();

   // Confirm all elements are identical, if not set error code
   // and break out of loop
   for (int N = 0; N < NTot; ++N) {
      if (CollapsedInit(N) != CollapsedTest(N)) {
         IErr = -1;
         break;
      }
   }

   if (IErr == 0) {
      LOG_INFO("HaloTest: {} exchange test PASS", Label);
   } else {
      LOG_INFO("HaloTest: {} exchange test FAIL", Label);
      TotErr += -1;
   }

   return;
} // end haloExchangeTest


//------------------------------------------------------------------------------
// The test driver. Performs halo exchange tests of all index spaces and all
// supported YAKL array types.

int main(int argc, char *argv[]) {

//   FILE *ofptr;

   // Error tracking variables
   OMEGA::I4 TotErr{0};
   OMEGA::I4 IErr{0};

   // Initialize global MPI environment and YAKL
   MPI_Init(&argc, &argv);
   yakl::init();

   // Initialize the machine environment and fetch the default environment
   // pointer, the MPI communicator and the task ID of the local task
   OMEGA::MachEnv::init(MPI_COMM_WORLD);
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefaultEnv();
   MPI_Comm DefComm       = DefEnv->getComm();
   OMEGA::I4 MyTask       = DefEnv->getMyTask();

   // Initialize the IO system
   IErr = OMEGA::IOInit(DefComm);
   if (IErr != 0)
      LOG_ERROR("DecompTest: error initializing parallel IO");

   // Create the default decomposition (initializes the decomposition)
   IErr = OMEGA::Decomp::init();
   if (IErr != 0)
      LOG_ERROR("DecompTest: error initializing default decomposition");

   // Retrieve the default decomposition
   OMEGA::Decomp *DefDecomp = OMEGA::Decomp::getDefault();

   // Create the halo exchange object for the given MachEnv and Decomp
   OMEGA::Halo MyHalo(DefEnv, DefDecomp);


   // 1DI4 array tests for each index space (cell, edge, and vertex)
   OMEGA::ArrayHost1DI4 Init1DI4Cell("Init1DI4Cell", DefDecomp->NCellsAll+1);
   OMEGA::ArrayHost1DI4 Test1DI4Cell("Test1DI4Cell", DefDecomp->NCellsAll+1);
   Init1DI4Cell = DefDecomp->CellIDH;
   haloExchangeTest(MyHalo, Init1DI4Cell, Test1DI4Cell,
                                  "1DI4 Cell", TotErr);

   OMEGA::ArrayHost1DI4 Init1DI4Edge("Init1DI4Edge", DefDecomp->NEdgesAll+1);
   OMEGA::ArrayHost1DI4 Test1DI4Edge("Test1DI4Edge", DefDecomp->NEdgesAll+1);
   Init1DI4Edge = DefDecomp->EdgeIDH;
   haloExchangeTest(MyHalo, Init1DI4Edge, Test1DI4Edge,
                                  "1DI4 Edge", TotErr, OMEGA::OnEdge);

   OMEGA::ArrayHost1DI4 Init1DI4Vertex("Init1DI4Vertex",
                                       DefDecomp->NVerticesAll+1);
   OMEGA::ArrayHost1DI4 Test1DI4Vertex("Test1DI4Vertex",
                                        DefDecomp->NVerticesAll+1);
   Init1DI4Vertex = DefDecomp->VertexIDH;
   haloExchangeTest(MyHalo, Init1DI4Vertex, Test1DI4Vertex,
                                  "1DI4 Vertex", TotErr, OMEGA::OnVertex);

   // Declaration of variables for remaining tests

   // Random dimension sizes
   OMEGA::I4 N2{20};
   OMEGA::I4 N3{10};
   OMEGA::I4 N4{3};
   OMEGA::I4 N5{2};

   // Number of owned + halo cells
   OMEGA::I4 NCellsAll = DefDecomp->NCellsAll;

   // Declare init and test arrays for all the remaining array types
   OMEGA::ArrayHost1DI8 Init1DI8("Init1DI8", NCellsAll);
   OMEGA::ArrayHost1DR4 Init1DR4("Init1DR4", NCellsAll);
   OMEGA::ArrayHost1DR8 Init1DR8("Init1DR8", NCellsAll);
   OMEGA::ArrayHost2DI4 Init2DI4("Init2DI4", NCellsAll, N2);
   OMEGA::ArrayHost2DI8 Init2DI8("Init2DI8", NCellsAll, N2);
   OMEGA::ArrayHost2DR4 Init2DR4("Init2DR4", NCellsAll, N2);
   OMEGA::ArrayHost2DR8 Init2DR8("Init2DR8", NCellsAll, N2);
   OMEGA::ArrayHost3DI4 Init3DI4("Init3DI4", N3, NCellsAll, N2);
   OMEGA::ArrayHost3DI8 Init3DI8("Init3DI8", N3, NCellsAll, N2);
   OMEGA::ArrayHost3DR4 Init3DR4("Init3DR4", N3, NCellsAll, N2);
   OMEGA::ArrayHost3DR8 Init3DR8("Init3DR8", N3, NCellsAll, N2);
   OMEGA::ArrayHost4DI4 Init4DI4("Init4DI4", N4, N3, NCellsAll, N2);
   OMEGA::ArrayHost4DI8 Init4DI8("Init4DI8", N4, N3, NCellsAll, N2);
   OMEGA::ArrayHost4DR4 Init4DR4("Init4DR4", N4, N3, NCellsAll, N2);
   OMEGA::ArrayHost4DR8 Init4DR8("Init4DR8", N4, N3, NCellsAll, N2);
   OMEGA::ArrayHost5DI4 Init5DI4("Init5DI4", N5, N4, N3, NCellsAll, N2);
   OMEGA::ArrayHost5DI8 Init5DI8("Init5DI8", N5, N4, N3, NCellsAll, N2);
   OMEGA::ArrayHost5DR4 Init5DR4("Init5DR4", N5, N4, N3, NCellsAll, N2);
   OMEGA::ArrayHost5DR8 Init5DR8("Init5DR8", N5, N4, N3, NCellsAll, N2);
   OMEGA::ArrayHost1DI8 Test1DI8("Test1DI8", NCellsAll);
   OMEGA::ArrayHost1DR4 Test1DR4("Test1DR4", NCellsAll);
   OMEGA::ArrayHost1DR8 Test1DR8("Test1DR8", NCellsAll);
   OMEGA::ArrayHost2DI4 Test2DI4("Test2DI4", NCellsAll, N2);
   OMEGA::ArrayHost2DI8 Test2DI8("Test2DI8", NCellsAll, N2);
   OMEGA::ArrayHost2DR4 Test2DR4("Test2DR4", NCellsAll, N2);
   OMEGA::ArrayHost2DR8 Test2DR8("Test2DR8", NCellsAll, N2);
   OMEGA::ArrayHost3DI4 Test3DI4("Test3DI4", N3, NCellsAll, N2);
   OMEGA::ArrayHost3DI8 Test3DI8("Test3DI8", N3, NCellsAll, N2);
   OMEGA::ArrayHost3DR4 Test3DR4("Test3DR4", N3, NCellsAll, N2);
   OMEGA::ArrayHost3DR8 Test3DR8("Test3DR8", N3, NCellsAll, N2);
   OMEGA::ArrayHost4DI4 Test4DI4("Test4DI4", N4, N3, NCellsAll, N2);
   OMEGA::ArrayHost4DI8 Test4DI8("Test4DI8", N4, N3, NCellsAll, N2);
   OMEGA::ArrayHost4DR4 Test4DR4("Test4DR4", N4, N3, NCellsAll, N2);
   OMEGA::ArrayHost4DR8 Test4DR8("Test4DR8", N4, N3, NCellsAll, N2);
   OMEGA::ArrayHost5DI4 Test5DI4("Test5DI4", N5, N4, N3, NCellsAll, N2);
   OMEGA::ArrayHost5DI8 Test5DI8("Test5DI8", N5, N4, N3, NCellsAll, N2);
   OMEGA::ArrayHost5DR4 Test5DR4("Test5DR4", N5, N4, N3, NCellsAll, N2);
   OMEGA::ArrayHost5DR8 Test5DR8("Test5DR8", N5, N4, N3, NCellsAll, N2);

   // Initialize and run remaining 1D tests
   for (int ICell = 0; ICell < NCellsAll; ++ICell) {
      OMEGA::I4 NewVal = DefDecomp->CellIDH(ICell);
      Init1DI8(ICell) = static_cast<OMEGA::I8>(NewVal);
      Init1DR4(ICell) = static_cast<OMEGA::R4>(NewVal);
      Init1DR8(ICell) = static_cast<OMEGA::R8>(NewVal);
   }

   haloExchangeTest(MyHalo, Init1DI8, Test1DI8, "1DI8", TotErr);
   haloExchangeTest(MyHalo, Init1DR4, Test1DR4, "1DR4", TotErr);
   haloExchangeTest(MyHalo, Init1DR8, Test1DR8, "1DR8", TotErr);

   // Initialize and run 2D tests
   for (int ICell = 0; ICell < NCellsAll; ++ICell) {
      for(int J = 0; J < N2; ++J) {
         OMEGA::I4 NewVal = (J + 1)*DefDecomp->CellIDH(ICell);
         Init2DI4(ICell, J) = NewVal;
         Init2DI8(ICell, J) = static_cast<OMEGA::I8>(NewVal);
         Init2DR4(ICell, J) = static_cast<OMEGA::R4>(NewVal);
         Init2DR8(ICell, J) = static_cast<OMEGA::R8>(NewVal);
      }
   }

   haloExchangeTest(MyHalo, Init2DI4, Test2DI4, "2DI4", TotErr);
   haloExchangeTest(MyHalo, Init2DI8, Test2DI8, "2DI8", TotErr);
   haloExchangeTest(MyHalo, Init2DR4, Test2DR4, "2DR4", TotErr);
   haloExchangeTest(MyHalo, Init2DR8, Test2DR8, "2DR8", TotErr);

   // Initialize and run 3D tests
   for (int K = 0; K < N3; ++K) {
      for (int ICell = 0; ICell < NCellsAll; ++ICell) {
         for (int J = 0; J < N2; ++J) {
            OMEGA::I4 NewVal = (K + 1)*(J + 1)*DefDecomp->CellIDH(ICell);
            Init3DI4(K, ICell, J) = NewVal;
            Init3DI8(K, ICell, J) = static_cast<OMEGA::I8>(NewVal);
            Init3DR4(K, ICell, J) = static_cast<OMEGA::R4>(NewVal);
            Init3DR8(K, ICell, J) = static_cast<OMEGA::R8>(NewVal);
         }
      }
   }

   haloExchangeTest(MyHalo, Init3DI4, Test3DI4, "3DI4", TotErr);
   haloExchangeTest(MyHalo, Init3DI8, Test3DI8, "3DI8", TotErr);
   haloExchangeTest(MyHalo, Init3DR4, Test3DR4, "3DR4", TotErr);
   haloExchangeTest(MyHalo, Init3DR8, Test3DR8, "3DR8", TotErr);

   // Initialize and run 4D tests
   for (int L = 0; L < N4; ++L) {
      for (int K = 0; K < N3; ++K) {
         for (int ICell = 0; ICell < NCellsAll; ++ICell) {
            for (int J = 0; J < N2; ++J) {
               OMEGA::I4 NewVal =
                  (L + 1)*(K + 1)*(J + 1)*DefDecomp->CellIDH(ICell);
               Init4DI4(L, K, ICell, J) = NewVal;
               Init4DI8(L, K, ICell, J) = static_cast<OMEGA::I8>(NewVal);
               Init4DR4(L, K, ICell, J) = static_cast<OMEGA::R4>(NewVal);
               Init4DR8(L, K, ICell, J) = static_cast<OMEGA::R8>(NewVal);
            }
         }
      }
   }

   haloExchangeTest(MyHalo, Init4DI4, Test4DI4, "4DI4", TotErr);
   haloExchangeTest(MyHalo, Init4DI8, Test4DI8, "4DI8", TotErr);
   haloExchangeTest(MyHalo, Init4DR4, Test4DR4, "4DR4", TotErr);
   haloExchangeTest(MyHalo, Init4DR8, Test4DR8, "4DR8", TotErr);

   // Initialize and run 5D tests
   for (int M = 0; M < N5; ++M) {
      for (int L = 0; L < N4; ++L) {
         for (int K = 0; K < N3; ++K) {
            for (int ICell = 0; ICell < NCellsAll; ++ICell) {
               for (int J = 0; J < N2; ++J) {
                  OMEGA::I4 NewVal =
                     (M + 1)*(L + 1)*(K + 1)*(J + 1)*DefDecomp->CellIDH(ICell);
                  Init5DI4(M ,L, K, ICell, J) = NewVal;
                  Init5DI8(M ,L, K, ICell, J) = static_cast<OMEGA::I8>(NewVal);
                  Init5DR4(M ,L, K, ICell, J) = static_cast<OMEGA::R4>(NewVal);
                  Init5DR8(M ,L, K, ICell, J) = static_cast<OMEGA::R8>(NewVal);
               }
            }
         }
      }
   }

   haloExchangeTest(MyHalo, Init5DI4, Test5DI4, "5DI4", TotErr);
   haloExchangeTest(MyHalo, Init5DI8, Test5DI8, "5DI8", TotErr);
   haloExchangeTest(MyHalo, Init5DR4, Test5DR4, "5DR4", TotErr);
   haloExchangeTest(MyHalo, Init5DR8, Test5DR8, "5DR8", TotErr);

   if (TotErr == 0) {
      LOG_INFO("HaloTest: Successful completion");
   } else {
      LOG_INFO("HaloTest: Failed");
      return -1;
   }
   yakl::finalize();
   MPI_Finalize();


} // end of main
//===-----------------------------------------------------------------------===/
