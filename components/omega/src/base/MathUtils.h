#ifndef OMEGA_MATHUTILS_H
#define OMEGA_MATHUTILS_H

//===-- base/MathUtils.h - mathematical utilities  --------------*- C++ -*-===//
//
/// \file
/// \brief Defines mathematical functions and utilities
///
/// This header defines mathematical functions (pow, exp, sin, etc) and other
/// general mathematical utilities needed in Omega.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include <concepts>
#include <type_traits>

namespace OMEGA::Math {

// Specialized power functions for compile-time constant integer exponent

// Non-negative power case
template <int N, class T>
KOKKOS_INLINE_FUNCTION constexpr T pow(T X)
   requires(N >= 0)
{
   if constexpr (N == 0) {
      return static_cast<T>(1);
   } else if constexpr (N % 2 == 0) {
      const T Y = pow<N / 2>(X);
      return Y * Y;
   } else if constexpr (N % 2 == 1) {
      const T Y = pow<N / 2>(X);
      return Y * Y * X;
   }
}

// Negative power case
// Only for floating-point types
template <int N, std::floating_point T>
KOKKOS_INLINE_FUNCTION constexpr T pow(T X)
   requires(N < 0)
{
   return static_cast<T>(1) / pow<-N>(X);
}

// C++ and Kokkos don't provide min/max mathematical functions
// They only provide fmin/fmax for floating-point arguments
// There are min/max *algorithms*, but they take their arguments by reference,
// which can lead to compilation issues on GPU, and to sub-optimal code
// generation Hence, we roll our own min/max functions that call fmin/fmax when
// appropriate We also provide versions with more than two arguments

// Min of two integral types
template <std::integral T> KOKKOS_INLINE_FUNCTION constexpr T min(T A, T B) {
   return B < A ? B : A;
}

// Min of two floating-point types
template <std::floating_point T>
KOKKOS_INLINE_FUNCTION constexpr T min(T A, T B) {
   return Kokkos::fmin(A, B);
}

// Min of two arbitrary types
template <class T1, class T2>
KOKKOS_INLINE_FUNCTION constexpr auto min(T1 A, T2 B) {
   using CT = std::common_type_t<T1, T2>;
   return min(static_cast<CT>(A), static_cast<CT>(B));
}

// Min of three or more arguments
template <class T, class... Ts>
   requires(sizeof...(Ts) >= 2)
KOKKOS_INLINE_FUNCTION constexpr auto min(T X0, Ts... Xs) {
   return min(X0, min(Xs...));
}

// Max of two integral types
template <std::integral T> KOKKOS_INLINE_FUNCTION constexpr T max(T A, T B) {
   return A < B ? B : A;
}

// Max of two floating-point types
template <std::floating_point T>
KOKKOS_INLINE_FUNCTION constexpr T max(T A, T B) {
   return Kokkos::fmax(A, B);
}

// Max of two arbitrary types
template <class T1, class T2>
KOKKOS_INLINE_FUNCTION constexpr auto max(T1 A, T2 B) {
   using CT = std::common_type_t<T1, T2>;
   return max(static_cast<CT>(A), static_cast<CT>(B));
}

// Max of three or more arguments
template <class T, class... Ts>
   requires(sizeof...(Ts) >= 2)
KOKKOS_INLINE_FUNCTION constexpr auto max(T X0, Ts... Xs) {
   return max(X0, max(Xs...));
}

KOKKOS_INLINE_FUNCTION
bool isApprox(Real X, Real Y, Real RTol, Real ATol = 0) {
   if (Kokkos::isnan(X) || Kokkos::isnan(Y) || Kokkos::isinf(X) ||
       Kokkos::isinf(Y)) {
      return false; // Treat NaN or Inf as failure
   }

   return Kokkos::abs(X - Y) <=
          max(ATol, RTol * max(Kokkos::abs(X), Kokkos::abs(Y)));
}

// Below is a list of mathematical functions that we take from Kokkos without
// modification See
// https://kokkos.org/kokkos-core-wiki/API/core/numerics/mathematical-functions.html

// Basic operations
using Kokkos::abs;
using Kokkos::fdim;
using Kokkos::fma;
using Kokkos::nan;
using Kokkos::remainder;
// using Kokkos::remquo;

// Exponential functions
using Kokkos::exp;
using Kokkos::exp2;
using Kokkos::expm1;
using Kokkos::log;
using Kokkos::log10;
using Kokkos::log1p;
using Kokkos::log2;

// Power functions
using Kokkos::cbrt;
using Kokkos::hypot;
using Kokkos::pow;
using Kokkos::sqrt;

// Trigonometric functions
using Kokkos::acos;
using Kokkos::asin;
using Kokkos::atan;
using Kokkos::atan2;
using Kokkos::cos;
using Kokkos::sin;
using Kokkos::tan;

// Hyperbolic functions
using Kokkos::acosh;
using Kokkos::asinh;
using Kokkos::atanh;
using Kokkos::cosh;
using Kokkos::sinh;
using Kokkos::tanh;

// Error and gamma functions
using Kokkos::erf;
using Kokkos::erfc;
using Kokkos::lgamma;
using Kokkos::tgamma;

// Nearest integer floating point operations
using Kokkos::ceil;
using Kokkos::floor;
// using Kokkos::llrint;
// using Kokkos::llround;
// using Kokkos::lrint;
// using Kokkos::lround;
using Kokkos::nearbyint;
// using Kokkos::rint;
using Kokkos::round;
using Kokkos::trunc;

// Floating-point manipulation functions
using Kokkos::copysign;
// using Kokkos::frexp;
// using Kokkos::ilogb;
// using Kokkos::ldexp;
using Kokkos::logb;
// using Kokkos::modf;
using Kokkos::nextafter;
// using Kokkos::nexttoward;
// using Kokkos::scalbln;
// using Kokkos::scalbn;

// Classification and comparison
// using Kokkos::fpclassify;
using Kokkos::isfinite;
// using Kokkos::isgreater;
// using Kokkos::isgreaterequal;
using Kokkos::isinf;
// using Kokkos::isless;
// using Kokkos::islessequal;
// using Kokkos::islessgreater;
using Kokkos::isnan;
// using Kokkos::isnormal;
// using Kokkos::isunordered;
using Kokkos::signbit;

// Non-standard
// using Kokkos::rcp;
using Kokkos::rsqrt;
} // namespace OMEGA::Math

#endif
