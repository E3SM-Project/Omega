#ifndef OMEGA_GLS_CONSTANTS_H
#define OMEGA_GLS_CONSTANTS_H
//===-- ocn/GLSConstants.h - GLS Constants and Stability Functions --*- C++
//-*-===//
//
/// \file
/// \brief Generic Length Scale (GLS) closure constants and stability
/// functions
///
/// This header defines the closure-family exponents, stability-function
/// coefficient tables, and the device-callable stability-function evaluator
/// for the GLS turbulence closure (Umlauf and Burchard, 2003; Warner et al.,
/// 2005). Follows the ROMS/REMORA notation and default parameter tables.
//
//===----------------------------------------------------------------------===//

#include "GlobalConstants.h"
#include "OmegaKokkos.h"

namespace OMEGA {

/// @brief Selects the stability-function formulation used to relate the
/// non-dimensional shear/buoyancy numbers (Gm, Gh) to the momentum and
/// tracer stability functions (Sm, Sh).
enum class GLSStabilityType : I4 {
   CanutoA  = 0, ///< Canuto et al. (2001) stability functions, variant A
   CanutoB  = 1, ///< Canuto et al. (2001) stability functions, variant B
   Galperin = 2  ///< Galperin et al. (1988) / Mellor-Yamada quasi-equilibrium
};

/// @brief Precomputed, run-constant quantities derived from the GLS
/// configuration. Computed once in GLSMix::init() and captured by value in
/// device kernels so the hot loop never calls pow()/exp() on configuration
/// values.
struct GLSDerivedConstants {
   // Exponent groupings (Umlauf and Burchard 2003 notation)
   Real CmuFac1  = 0.0_Real; ///< cmu0^(-p/n)
   Real CmuFac2  = 0.0_Real; ///< cmu0^(3+p/n)
   Real CmuFac3  = 0.0_Real; ///< 1 / cmu0^2
   Real GlsFac2  = 0.0_Real; ///< cmu0^p * n * kappa^n (bottom psi flux)
   Real GlsFac3  = 0.0_Real; ///< cmu0^p * n (surface psi flux)
   Real GlsFac4  = 0.0_Real; ///< cmu0^p (bottom psi Dirichlet BC)
   Real GlsFac5  = 0.0_Real; ///< 0.56^(n/2) * cmu0^p (psi comparison bound)
   Real GlsFac6  = 0.0_Real; ///< 8 / cmu0^6
   Real GlsExp1  = 0.0_Real; ///< 1/n
   Real TkeExp1  = 0.0_Real; ///< m/n
   Real TkeExp2  = 0.0_Real; ///< 0.5 + m/n
   Real TkeExp4  = 0.0_Real; ///< m + 0.5*n
   Real Cmu0ExpP = 0.0_Real; ///< cmu0^p
   Real Cmu0Cube = 0.0_Real; ///< cmu0^3
   Real Sqrt2    = 0.0_Real; ///< sqrt(2)

   // Stability-function limits
   Real Gh0   = 0.0_Real; ///< Steady-state Gh upper limit
   Real Ghcri = 0.0_Real; ///< Critical Gh for the smooth cutoff
   Real Ghmin = -0.28_Real; ///< Lower bound on Gh

   // Canuto et al. (2001) quartic rational-function coefficients
   Real S0 = 0.0_Real, S1 = 0.0_Real, S2 = 0.0_Real;
   Real S4 = 0.0_Real, S5 = 0.0_Real, S6 = 0.0_Real;
   Real B0 = 0.0_Real, B1 = 0.0_Real, B2 = 0.0_Real;
   Real B3 = 0.0_Real, B4 = 0.0_Real, B5 = 0.0_Real;

   // Galperin / Mellor-Yamada quasi-equilibrium coefficients
   Real MySm2 = 0.0_Real, MySm3 = 0.0_Real, MySm4 = 0.0_Real;
   Real MySh1 = 0.0_Real, MySh2 = 0.0_Real;
};

/// @brief Compute the run-constant GLS derived quantities from the
/// configured closure exponents and stability-function choice. Intended to
/// be called once, on the host, from GLSMix::init().
///
/// The Canuto L1-L8 coefficients and Gh0/Ghcri limits follow the
/// stability-function-specific defaults tabulated in Umlauf and Burchard
/// (2003) Table 2 and Warner et al. (2005), reproduced here from REMORA's
/// GLS implementation (itself ported from ROMS gls_corstep.F).
inline GLSDerivedConstants
computeGLSDerivedConstants(GLSStabilityType Type, Real P, Real M, Real N,
                           Real Cmu0, Real C1, Real C2, Real C3Minus,
                           Real C3Plus, Real SigmaK, Real SigmaPsi) {
   GLSDerivedConstants C;

   C.CmuFac1 = std::pow(Cmu0, -P / N);
   C.CmuFac2 = std::pow(Cmu0, 3.0_Real + P / N);
   C.CmuFac3 = 1.0_Real / (Cmu0 * Cmu0);

   C.GlsFac2 = std::pow(Cmu0, P) * N * std::pow(VonKar, N);
   C.GlsFac3 = std::pow(Cmu0, P) * N;
   C.GlsFac4 = std::pow(Cmu0, P);
   C.GlsFac5 = std::pow(0.56_Real, 0.5_Real * N) * std::pow(Cmu0, P);
   C.GlsFac6 = 8.0_Real / std::pow(Cmu0, 6.0_Real);

   C.GlsExp1 = 1.0_Real / N;
   C.TkeExp1 = M / N;
   C.TkeExp2 = 0.5_Real + M / N;
   C.TkeExp4 = M + 0.5_Real * N;

   C.Cmu0ExpP = std::pow(Cmu0, P);
   C.Cmu0Cube = Cmu0 * Cmu0 * Cmu0;
   C.Sqrt2    = std::sqrt(2.0_Real);

   Real L1 = 0.0_Real, L2 = 0.0_Real, L3 = 0.0_Real, L4 = 0.0_Real;
   Real L5 = 0.0_Real, L6 = 0.0_Real, L7 = 0.0_Real, L8 = 0.0_Real;

   switch (Type) {
   case GLSStabilityType::CanutoA:
      C.Gh0   = 0.0329_Real;
      C.Ghcri = 0.03_Real;
      L1 = 0.107_Real;
      L2 = 0.0032_Real;
      L3 = 0.0864_Real;
      L4 = 0.12_Real;
      L5 = 11.9_Real;
      L6 = 0.4_Real;
      L7 = 0.0_Real;
      L8 = 0.48_Real;
      break;
   case GLSStabilityType::CanutoB:
      C.Gh0   = 0.0444_Real;
      C.Ghcri = 0.0414_Real;
      L1 = 0.127_Real;
      L2 = 0.00336_Real;
      L3 = 0.0906_Real;
      L4 = 0.101_Real;
      L5 = 11.2_Real;
      L6 = 0.4_Real;
      L7 = 0.0_Real;
      L8 = 0.318_Real;
      break;
   case GLSStabilityType::Galperin:
   default:
      C.Gh0   = 0.028_Real;
      C.Ghcri = 0.02_Real;
      break;
   }

   if (Type == GLSStabilityType::CanutoA || Type == GLSStabilityType::CanutoB) {
      C.S0 = 1.5_Real * L1 * L5 * L5;
      C.S1 = -L4 * (L6 + L7) +
             2.0_Real * L4 * L5 * (L1 - L2 / 3.0_Real - L3) +
             1.5_Real * L1 * L5 * L8;
      C.S2 = -0.375_Real * L1 * (L6 * L6 - L7 * L7);
      C.S4 = 2.0_Real * L5;
      C.S5 = 2.0_Real * L4;
      C.S6 = (2.0_Real / 3.0_Real) * L5 * (3.0_Real * L3 * L3 - L2 * L2) -
             0.5_Real * L5 * L1 * (3.0_Real * L3 - L2) +
             0.75_Real * L1 * (L6 - L7);
      C.B0 = 3.0_Real * L5 * L5;
      C.B1 = L5 * (7.0_Real * L4 + 3.0_Real * L8);
      C.B2 = L5 * L5 * (3.0_Real * L3 * L3 - L2 * L2) -
             0.75_Real * (L6 * L6 - L7 * L7);
      C.B3 = L4 * (4.0_Real * L4 + 3.0_Real * L8);
      C.B5 = 0.25_Real * (L2 * L2 - 3.0_Real * L3 * L3) *
             (L6 * L6 - L7 * L7);
      C.B4 = L4 * (L2 * L6 - 3.0_Real * L3 * L7 - L5 * (L2 * L2 - L3 * L3)) +
             L5 * L8 * (3.0_Real * L3 * L3 - L2 * L2);
   } else {
      // Mellor-Yamada quasi-equilibrium coefficients (Galperin et al. 1988),
      // using the standard Mellor-Yamada level-2 constants.
      constexpr Real MyA1 = 0.92_Real;
      constexpr Real MyA2 = 0.74_Real;
      constexpr Real MyB1 = 16.6_Real;
      constexpr Real MyB2 = 10.1_Real;
      constexpr Real MyC1 = 0.08_Real;

      C.MySm2 = 9.0_Real * MyA1 * MyA2;
      C.MySm3 = MyA1 * (1.0_Real - 3.0_Real * MyC1 - 6.0_Real * MyA1 / MyB1);
      C.MySm4 = 18.0_Real * MyA1 * MyA1 + 9.0_Real * MyA1 * MyA2;
      C.MySh1 = MyA2 * (1.0_Real - 6.0_Real * MyA1 / MyB1);
      C.MySh2 = 3.0_Real * MyA2 * (6.0_Real * MyA1 + MyB2);
   }

   return C;
}

/// @brief Canuto stability functions require an upper bound on the shear
/// number Gm (set by the requirement that the denominator of the stability
/// functions stays positive), evaluated from Gh before Gm is clipped to its
/// physical shear-based value. No-op (returns a large number) for
/// non-Canuto stability types.
KOKKOS_INLINE_FUNCTION
Real computeCanutoGmBound(GLSStabilityType Type, const GLSDerivedConstants &C,
                          Real Gh) {
   if (Type != GLSStabilityType::CanutoA && Type != GLSStabilityType::CanutoB)
      return 1.0e30_Real;

   return (C.B0 / C.GlsFac6 - C.B1 * Gh + C.B3 * C.GlsFac6 * Gh * Gh) /
          (C.B2 - C.B4 * C.GlsFac6 * Gh);
}

/// @brief Evaluate the non-dimensional stability functions Sm (momentum) and
/// Sh (tracers) from the shear and buoyancy numbers Gm, Gh.
///
/// @param Type Stability-function formulation
/// @param C Run-constant derived quantities from computeGLSDerivedConstants
/// @param Gh Non-dimensional buoyancy number (-N^2 L^2 / 2k, bounded)
/// @param Gm Non-dimensional shear number (shear^2 L^2 / 2k, Canuto only)
/// @param[out] Sm Momentum stability function
/// @param[out] Sh Tracer stability function
KOKKOS_INLINE_FUNCTION
void computeStabilityFunctions(GLSStabilityType Type,
                               const GLSDerivedConstants &C, Real Gh, Real Gm,
                               Real &Sm, Real &Sh) {
   if (Type == GLSStabilityType::CanutoA || Type == GLSStabilityType::CanutoB) {
      const Real StabCff = C.B0 - C.B1 * C.GlsFac6 * Gh +
                           C.B2 * C.GlsFac6 * Gm +
                           C.B3 * C.GlsFac6 * C.GlsFac6 * Gh * Gh -
                           C.B4 * C.GlsFac6 * C.GlsFac6 * Gh * Gm +
                           C.B5 * C.GlsFac6 * C.GlsFac6 * Gm * Gm;
      Sm = (C.S0 - C.S1 * C.GlsFac6 * Gh + C.S2 * C.GlsFac6 * Gm) / StabCff;
      Sh = (C.S4 - C.S5 * C.GlsFac6 * Gh + C.S6 * C.GlsFac6 * Gm) / StabCff;
      Sm = Kokkos::max(Sm, 0.0_Real);
      Sh = Kokkos::max(Sh, 0.0_Real);

      // Convert from Canuto to ROMS/GLS normalization
      Sm = Sm * C.Sqrt2 / C.Cmu0Cube;
      Sh = Sh * C.Sqrt2 / C.Cmu0Cube;
   } else if (Type == GLSStabilityType::Galperin) {
      const Real CffGalperin = 1.0_Real - C.MySh2 * Gh;
      Sh                     = C.MySh1 / CffGalperin;
      Sm = (C.MySm3 + Sh * Gh * C.MySm4) / (1.0_Real - C.MySm2 * Gh);
   } else {
      Sm = 0.0_Real;
      Sh = 0.0_Real;
   }
}

} // namespace OMEGA

#endif // OMEGA_GLS_CONSTANTS_H
