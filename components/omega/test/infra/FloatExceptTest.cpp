//===-- Test driver for Omega floating-point exception control ---*- C++ -*-===/
//
/// \file
/// \brief Test driver for Omega floating-point exception control
///
/// This driver tests the capabilities of Omega to unmask floating-point
/// exceptions. It intentionally causes an exception and is expected to fail.
//
//===-----------------------------------------------------------------------===/

#include "FloatExcept.h"
#include "MachEnv.h"

using namespace OMEGA;

void divByZero() {
   volatile double D0   = 0.0;
   volatile double D1   = 1.0;
   volatile double DInf = D1 / D0;
}

int main(int argc, char **argv) {
   MPI_Init(&argc, &argv);
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv = OMEGA::MachEnv::getDefault();
   initLogging(DefEnv);

   FloatExceptStatus Status = enableFloatExceptions();

   divByZero();

   CHECK_ERROR_ABORT(Status.Err, "Floating-point Exceptions Tests FAIL");

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return 0;
}
