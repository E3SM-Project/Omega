(omega-dev-math-utils)=

# Mathematical Utilities

`MathUtils.h` provides Omega-specific mathematical utilities and makes a
selected set of Kokkos mathematical functions available in the
`OMEGA::Math` namespace. Include `MathUtils.h` and use the `Math::` prefix
for these functions.

## Omega utilities

### Compile-time integer powers

`Math::pow<N>(X)` computes an integer power using the exponent `N` as a
compile-time parameter. Non-negative exponents are supported for any type
that supports multiplication. Negative exponents are supported for
floating-point types and compute the reciprocal of the corresponding
positive power.

```c++
const Real Cube = Math::pow<3>(X);
const Real Inverse = Math::pow<-2>(X);
```

The header also imports `Kokkos::pow`, so runtime-exponent overloads
available in Kokkos can also be called as `Math::pow`.

### Minimum and maximum

`Math::min(A, B)` and `Math::max(A, B)` provide minimum and maximum
operations for integral and floating-point values. Mixed input types are
converted to their common type. Each function also accepts three or more
arguments:

```c++
const Real Smallest = Math::min(A, B, C);
const int Largest = Math::max(I, J, K, L);
```

These functions are Omega utilities; they are distinct from the algorithms
`std::min` and `std::max`.

### Clamping

`Math::clamp(X, Lo, Hi)` returns `X` limited to the inclusive interval
`[Lo, Hi]`. The lower bound must not exceed the upper bound. Its arguments can
have different types and are compared using the common type of those arguments:

```c++
const Real BoundedDepth = Math::clamp(Depth, 0.0_Real, MaxDepth);
```

### Linear interpolation

`Math::linearInterp(X, Y1, X1, Y2, X2)` evaluates the line through
`(X1, Y1)` and `(X2, Y2)` at `X`, including extrapolation outside the two
points. When `X1` and `X2` are equal, it returns `Y1`:

```c++
const Real Value = Math::linearInterp(X, Y1, X1, Y2, X2);
```

### Approximate equality

`Math::isApprox(X, Y, RTol, ATol = 0)` returns whether two `Real` values
are within combined relative and absolute tolerances:

```text
abs(X - Y) <= max(ATol, RTol * max(abs(X), abs(Y)))
```

It returns `false` if either value is infinite or NaN.

## Mathematical functions from Kokkos

The following functions are imported from Kokkos without modification and
are available in `OMEGA::Math`:

| Category | Functions |
| --- | --- |
| Basic operations | `abs`, `fdim`, `fma`, `nan`, `remainder` |
| Exponential | `exp`, `exp2`, `expm1`, `log`, `log10`, `log1p`, `log2` |
| Power | `cbrt`, `hypot`, `pow`, `sqrt` |
| Trigonometric | `acos`, `asin`, `atan`, `atan2`, `cos`, `sin`, `tan` |
| Hyperbolic | `acosh`, `asinh`, `atanh`, `cosh`, `sinh`, `tanh` |
| Error and gamma | `erf`, `erfc`, `lgamma`, `tgamma` |
| Nearest-integer operations | `ceil`, `floor`, `round`, `trunc` |
| Floating-point manipulation | `copysign`, `logb`, `nextafter` |
| Classification and comparison | `isfinite`, `isinf`, `isnan`, `signbit` |
| Non-standard | `rsqrt` |

For example, call `Math::sqrt(X)` or `Math::isfinite(X)`. These functions
retain the behavior and supported types provided by the Kokkos version used
to build Omega.
