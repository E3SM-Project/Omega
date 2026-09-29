//===-- FloatExcept.cpp - Omega floating-point exception handling implementation
//-*- C++ -*-===//
//
/// \file
/// \brief Utilities for handling floating-point exceptions
///
/// This file contains the implementation of the floating-point exception
/// handling functions for Omega.
//
//===---------------------------------------------------------------------------------------===//

#include "FloatExcept.h"
#include "OmegaKokkos.h"
#include <Kokkos_Core.hpp>
#include <cfenv>

namespace OMEGA {

// Floating-point exception handling using the feenableexcept GNU extension
#ifdef OMEGA_HAVE_FEENABLEEXCEPT
static FloatExceptStatus enableFloatExceptionsGnuExtension(int Exceptions) {
   FloatExceptStatus Status;

   // Clear previous exceptions
   int Clear = std::feclearexcept(FE_ALL_EXCEPT);
   if (Clear != 0) {
      Status.Err +=
          Error(ErrorCode::Fail, "Failed to clear previous exceptions");
      return Status;
   }

   int EnableVal = feenableexcept(Exceptions);
   if (EnableVal == -1) {
      Status.Err +=
          Error(ErrorCode::Fail, "Failed to enable floating-point exceptions");
   } else {
      Status.OldExceptions = EnableVal;
   }
   return Status;
}

static FloatExceptStatus disableFloatExceptionsGnuExtension(int Exceptions) {
   FloatExceptStatus Status;

   int DisableVal = fedisableexcept(Exceptions);

   if (DisableVal == -1) {
      Status.Err +=
          Error(ErrorCode::Fail, "Failed to disable floating-point exceptions");
   } else {
      Status.OldExceptions = DisableVal;
   }

   return Status;
}
#endif

// Floating-point exception handling for ARM64
#ifdef __aarch64__

// Exception bits in fpcr are shifted by this amount
static constexpr int FPCRShift = 8;

static FloatExceptStatus enableFloatExceptionsArm64(int Exceptions) {
   FloatExceptStatus Status;

   std::fenv_t env;
   if (fegetenv(&env) != 0) {
      Status.Err += Error(ErrorCode::Fail, "Failed to get floating-point env");
      return Status;
   }

   // Get currently enabled exceptions
   Status.OldExceptions = (env.__fpcr >> FPCRShift) & FE_ALL_EXCEPT;

   // Clear all exceptions
   env.__fpsr &= ~FE_ALL_EXCEPT;

   // Enable selected exceptions
   env.__fpcr |= (Exceptions << FPCRShift);

   if (fesetenv(&env) != 0) {
      Status.Err += Error(ErrorCode::Fail, "Failed to set floating-point env");
   }

   return Status;
}

static FloatExceptStatus disableFloatExceptionsArm64(int Exceptions) {
   FloatExceptStatus Status;

   std::fenv_t env;
   if (fegetenv(&env) != 0) {
      Status.Err += Error(ErrorCode::Fail, "Failed to get floating-point env");
      return Status;
   }

   // Get currently enabled exceptions
   Status.OldExceptions = (env.__fpcr >> FPCRShift) & FE_ALL_EXCEPT;

   // Disable selected exceptions
   env.__fpcr &= ~(Exceptions << FPCRShift);

   if (fesetenv(&env) != 0) {
      Status.Err += Error(ErrorCode::Fail, "Failed to set floating-point env");
   }

   return Status;
}
#endif

// Enable selected floating-point exceptions
// If feenableexcept exists then use it, else
// dispatch based on architecture
FloatExceptStatus enableFloatExceptions(int Exceptions) {
   auto EnableLambda = [=]() {
#ifdef OMEGA_HAVE_FEENABLEEXCEPT
      return enableFloatExceptionsGnuExtension(Exceptions);
#elif defined __aarch64__
      return enableFloatExceptionsArm64(Exceptions);
#else // unsupported arch
      return FloatExceptStatus{Error(ErrorCode::Fail,
                                     "Omega doesn't support enabling floating "
                                     "point exceptions on this architecture")};
#endif
   };

   // Enable exceptions on the main thread
   FloatExceptStatus Status = EnableLambda();

   // Enable exceptions on threads spawned by Kokkos
#if defined(KOKKOS_ENABLE_OPENMP) || defined(KOKKOS_ENABLE_THREADS)
   if (Kokkos::is_initialized()) {
      const auto NumThreads = ExecSpace().concurrency();
      // Not using Omege wrappers because we want to always use a static
      // schedule
      const auto Policy =
          Kokkos::RangePolicy<ExecSpace, Kokkos::Schedule<Kokkos::Static>>(
              0, NumThreads);
      Kokkos::parallel_reduce(
          Policy, KOKKOS_LAMBDA(int, Error &Err) { Err += EnableLambda().Err; },
          Status.Err);
   }
#endif

   return Status;
}

// Disable selected floating-point exceptions
// If feenableexcept exists then use it, else
// dispatch based on architecture
FloatExceptStatus disableFloatExceptions(int Exceptions) {
   auto DisableLambda = [=]() {
#ifdef OMEGA_HAVE_FEENABLEEXCEPT
      return disableFloatExceptionsGnuExtension(Exceptions);
#elif defined __aarch64__
      return disableFloatExceptionsArm64(Exceptions);
#else // unsupported arch
      return FloatExceptStatus{Error(ErrorCode::Fail,
                                     "Omega doesn't support disabling floating "
                                     "point exceptions on this architecture")};
#endif
   };

   // Disable exceptions on the main thread
   FloatExceptStatus Status = DisableLambda();

   // Disable exceptions on threads spawned by Kokkos
#if defined(KOKKOS_ENABLE_OPENMP) || defined(KOKKOS_ENABLE_THREADS)

   if (Kokkos::is_initialized()) {
      const auto NumThreads = ExecSpace().concurrency();
      // Not using Omega wrappers because we want to always use a static
      // schedule
      const auto Policy =
          Kokkos::RangePolicy<ExecSpace, Kokkos::Schedule<Kokkos::Static>>(
              0, NumThreads);
      Kokkos::parallel_reduce(
          Policy,
          KOKKOS_LAMBDA(int, Error &Err) { Err += DisableLambda().Err; },
          Status.Err);
   }
#endif

   return Status;
}

// A helper function to conditionally enable floating-point exceptions in Omega
// unit tests Floating-point exceptions are enabled when explicitly requested
// using the OMEGA_FPEXCEPTS_IN_TESTS CMake option or when using the GNU
// compiler
FloatExceptStatus enableFloatExceptionsInTests(int Exceptions) {
   FloatExceptStatus Status;
// Enable if explicitly requested
#ifdef OMEGA_FPEXCEPTS_IN_TESTS
   Status = enableFloatExceptions(Exceptions);
#else
   // Enable if using the GNU compiler
#if defined(KOKKOS_COMPILER_GNU)
   Status = enableFloatExceptions(Exceptions);
#endif
#endif
   return Status;
}

} // namespace OMEGA
//===----------------------------------------------------------------------===//
