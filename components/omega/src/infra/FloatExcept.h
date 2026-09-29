#ifndef OMEGA_FLOATEXCEPT_H
#define OMEGA_FLOATEXCEPT_H
//===-- FloatExcept.h - Omega floating-point exception control --*- C++ -*-===//
//
/// \file
/// \brief Utilities for handling floating-point exceptions
///
/// This header defines functions for enabling and disabling floating-point
/// exceptions in Omega
//
//===----------------------------------------------------------------------===//

#include "Error.h"
#include <cfenv>

namespace OMEGA {

// Return value for the functions in this module combining error code with
// old exceptions mask
struct FloatExceptStatus {
   Error Err{};
   int OldExceptions{};
};

// Enable selected floating-point exceptions
// Returns FloatExceptStatus, which on success contains exceptions mask before
// the call to this function
FloatExceptStatus enableFloatExceptions(int Exceptions = (FE_DIVBYZERO |
                                                          FE_INVALID |
                                                          FE_OVERFLOW));

// Disable selected floating-point exceptions
// Returns FloatExceptStatus, which on success contains exceptions mask before
// the call to this function
FloatExceptStatus disableFloatExceptions(int Exceptions = FE_ALL_EXCEPT);

// A helper function to conditionally enable floating-point exceptions in Omega
// unit tests Floating-point exceptions are enabled when explicitly requested
// using the OMEGA_FPEXCEPTS_IN_TESTS CMake option or when using the GNU
// compiler
FloatExceptStatus enableFloatExceptionsInTests(int Exceptions = (FE_DIVBYZERO |
                                                                 FE_INVALID |
                                                                 FE_OVERFLOW));
} // namespace OMEGA
//===----------------------------------------------------------------------===//
#endif // OMEGA_FLOATEXCEPT_H
