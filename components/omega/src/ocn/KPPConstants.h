#ifndef OMEGA_KPP_CONSTANTS_H
#define OMEGA_KPP_CONSTANTS_H
//===-- ocn/KPPConstants.h - KPP Constants and Profiles --------*- C++ -*-===//
//
/// \file
/// \brief KPP-specific constants, parameters, and profile functions
///
/// This header defines constants, parameters, and inline functions for the
/// K-Profile Parameterization boundary layer mixing scheme. Includes
/// non-dimensional profile functions used in KPP coefficient calculations.
//
//===----------------------------------------------------------------------===//

#include "GlobalConstants.h"
#include "OmegaKokkos.h"

namespace OMEGA::KPP {

// ==========================================================================
// Physical Parameters for KPP
// ==========================================================================

/// Critical bulk Richardson number for defining OBL (Large et al. 1994)
constexpr Real RICRIT = 0.3;

/// Parameter for smoothing velocity shear profiles
constexpr Real ZETA_M_SCALE = VonKar; // Momentum scale (normalized)
constexpr Real ZETA_S_SCALE = 0.16;   // Tracer/salt scale
constexpr Real ZETA_T_SCALE = 0.16;   // Temperature scale

// ==========================================================================
// Monin-Obukhov stability function parameters (Large et al. 1994, App. B)
// Transition thresholds between weakly and strongly unstable regimes
// ==========================================================================

/// Transition zeta for momentum: below this value, strongly-unstable formula
/// is used. Default: -0.2 (CVMix default).
constexpr Real ZETA_M = -0.2_Real;

/// Transition zeta for scalars: below this value, strongly-unstable formula
/// is used. Default: -1.0 (CVMix default).
constexpr Real ZETA_S = -1.0_Real;

/// Derived constants for phi_m^{-1} strongly-unstable branch (momentum).
/// a_m = (1-16*ZETA_M)^{-0.25} * (1 - 4*ZETA_M)
constexpr Real A_MO_M = 1.2573615702_Real;
/// c_m = (1-16*ZETA_M)^{-0.25} * 12
constexpr Real C_MO_M = 8.3824104679_Real;

/// Derived constants for phi_s^{-1} strongly-unstable branch (scalar).
/// a_s = sqrt(1-16*ZETA_S) * (1 + 8*ZETA_S)  (can be negative)
constexpr Real A_MO_S = -28.8617393793_Real;
/// c_s = 24 * sqrt(1-16*ZETA_S)
constexpr Real C_MO_S = 98.9545350148_Real;

/// Surface mixing coefficients
constexpr Real HUON = 0.03;  // Surface momentum mixing parameter
constexpr Real BD   = 1.0;   // Buoyancy parameter (dimensionless)
constexpr Real C1   = 0.112; // Langmuir circulation parameter

/// Langmuir enhancement factor parameters
constexpr Real PEC_LANGMUIR = 0.5; // Peclet number for Langmuir

/// Minimum/maximum bounds on friction velocity values
constexpr Real MIN_USTAR = 1.0e-4; // Minimum friction velocity (m/s)
constexpr Real MAX_USTAR = 1.0;    // Upper limit to clamp rare extremes

/// Maximum vertical levels (for static allocations if needed)
constexpr I4 NLEV_MAX = 500;

// ==========================================================================
// OBL Depth Computation Parameters
// ==========================================================================

/// Safety multiplier for OBL search (prevents searching too deep)
/// Default: 1.0 (search to 1.0 * Ri_crit threshold)
constexpr Real STOP_OBL_SEARCH_MULT = 1.0;

/// Minimum OBL depth (m)
constexpr Real MIN_OBL_DEPTH = 2.0;

/// Minimum OBL under sea ice (m) when ice fraction > 0.15
constexpr Real MIN_OBL_UNDER_ICE = 5.0;

/// Ice fraction threshold below which OBL is fully computed
constexpr Real ICE_FRACTION_THRESHOLD = 0.05;

/// Ice fraction for triggering minimum OBL enforcement
constexpr Real ICE_SUPPRESSION_THRESHOLD = 0.15;

// ==========================================================================
// Surface Layer Parameters
// ==========================================================================

/// Surface layer extent (fraction of OBL depth)
/// Used for averaging turbulent scales near surface
constexpr Real SURFACE_LAYER_EXTENT = 0.1;

/// Number of smoothing passes for Richardson number (reduce noise)
constexpr I4 RI_SMOOTH_LOOPS = 2;

/// Prandtl number for converting momentum viscosity to tracer diffusivity
constexpr Real PRANDTL_NUMBER = 1.0;

// ==========================================================================
// KPP Profile Functions
// ==========================================================================

/// @brief G(sigma) - Non-local flux profile function
/// Non-zero only within the OBL. Applied to surface tracer fluxes.
/// REFERENCES: Large et al. (1994) Eq. (12)-(13), Large et al. (1997)
///
/// @param sigma Normalized vertical position (-z/h), 0 at surface, -1 at base
/// @return G(sigma) dimensionless profile value
KOKKOS_INLINE_FUNCTION
Real KPPProfileG(Real sigma) {
   // Omega uses sigma in [-1,0]. Convert to CVMix sigma_mu in [0,1]
   // where sigma_mu=0 at surface and sigma_mu=1 at OBL base.
   sigma = Kokkos::fmax(-1.0, Kokkos::fmin(0.0, sigma));

   const Real sigma_mu = -sigma;
   return sigma_mu * (1.0 - sigma_mu) * (1.0 - sigma_mu);
}

/// @brief G_pnl(sigma) - Parabolic non-local flux profile
/// Used only for non-local tracer flux when matching option is
/// "ParabolicNonLocal". This must NOT be used for KPP viscosity/diffusivity.
///
/// @param sigma Normalized vertical position (-z/h), 0 at surface, -1 at base
/// @return Dimensionless non-local profile value
KOKKOS_INLINE_FUNCTION
Real KPPProfileGParabolicNonLocal(Real sigma) {
   sigma = Kokkos::fmax(-1.0, Kokkos::fmin(0.0, sigma));

   // Parabolic non-local option in CVMix is (1-sigma_mu)^2.
   const Real sigma_mu  = -sigma;
   const Real one_minus = 1.0 - sigma_mu;
   return one_minus * one_minus;
}

/// @brief G_matchboth(sigma) - Cubic LMD-style non-local profile
/// Used for MatchBoth to distinguish it from SimpleShapes without relying on
/// external CVMix calls. In sigma_mu coordinates this is:
///   G = (1 - sigma_mu)^2 * (1 + 2*sigma_mu)
/// where sigma_mu in [0,1] is 0 at surface and 1 at OBL base.
///
/// @param sigma Normalized vertical position (-z/h), 0 at surface, -1 at base
/// @return Dimensionless non-local profile value
KOKKOS_INLINE_FUNCTION
Real KPPProfileGMatchBoth(Real sigma) {
   sigma = Kokkos::fmax(-1.0, Kokkos::fmin(0.0, sigma));

   const Real sigma_mu  = -sigma;
   const Real one_minus = 1.0 - sigma_mu;
   return one_minus * one_minus * (1.0 + 2.0 * sigma_mu);
}

/// @brief M1(sigma) - Momentum mixing profile function
/// Multiplies friction velocity and turbulent velocity scale
/// REFERENCES: Large et al. (1994) Eq. (11)
///
/// @param sigma Normalized vertical position (-z/h)
/// @return K*w_s profile multiplier (dimensionless)
KOKKOS_INLINE_FUNCTION
Real KPPProfileM1(Real sigma) {
   // CVMix simple gradient shape: sigma_mu*(1-sigma_mu)^2.
   sigma = Kokkos::fmax(-1.0, Kokkos::fmin(0.0, sigma));

   const Real sigma_mu = -sigma;
   return sigma_mu * (1.0 - sigma_mu) * (1.0 - sigma_mu);
}

/// @brief Matched KPP gradient profile shape.
///
/// Uses the SimpleShapes gradient profile plus a smooth correction that is
/// zero at the surface and equals ShapeAtBase at the OBL base. This lets
/// MatchBoth profiles meet pre-existing interior mixing at the BLD base while
/// preserving SimpleShapes behavior when ShapeAtBase is zero.
KOKKOS_INLINE_FUNCTION
Real KPPProfileMatched(Real sigma, Real ShapeAtBase) {
   sigma = Kokkos::fmax(-1.0, Kokkos::fmin(0.0, sigma));

   const Real sigma_mu = -sigma;
   const Real simple   = sigma_mu * (1.0 - sigma_mu) * (1.0 - sigma_mu);
   const Real smooth   = sigma_mu * sigma_mu * (3.0 - 2.0 * sigma_mu);
   return simple + ShapeAtBase * smooth;
}

/// @brief phi_m^{-1}(zeta) - Inverse momentum Monin-Obukhov stability function
/// Multiplied by von Karman constant and friction velocity to give turbulent
/// momentum velocity scale: w_m = kappa * u* * phi_m^{-1}(zeta)
/// Three-regime formulation per Large et al. (1994) Appendix B and CVMix.
///
/// @param zeta Monin-Obukhov stability coordinate (dimensionless)
/// @return phi_m^{-1} (dimensionless, > 0)
KOKKOS_INLINE_FUNCTION
Real KPPProfileM2(Real zeta) {
   if (zeta >= 0.0_Real) {
      // Stable regime
      return 1.0_Real / (1.0_Real + 5.0_Real * zeta);
   } else if (zeta >= ZETA_M) {
      // Weakly unstable: (1 - 16*zeta)^{1/4}
      return Kokkos::pow(1.0_Real - 16.0_Real * zeta, 0.25_Real);
   } else {
      // Strongly unstable: (a_m - c_m*zeta)^{1/3}
      return Kokkos::pow(A_MO_M - C_MO_M * zeta, 1.0_Real / 3.0_Real);
   }
}

/// @brief S1(sigma) - Tracer/scalar mixing profile function
/// Similar to momentum but computed separately for tracers
/// REFERENCES: Large et al. (1994) Eq. (11)
///
/// @param sigma Normalized vertical position
/// @return K*w_s profile multiplier for tracers (dimensionless)
KOKKOS_INLINE_FUNCTION
Real KPPProfileS1(Real sigma) {
   // CVMix simple gradient shape for tracers: sigma_mu*(1-sigma_mu)^2.
   sigma = Kokkos::fmax(-1.0, Kokkos::fmin(0.0, sigma));

   const Real sigma_mu = -sigma;
   return sigma_mu * (1.0 - sigma_mu) * (1.0 - sigma_mu);
}

/// @brief phi_s^{-1}(zeta) - Inverse scalar Monin-Obukhov stability function
/// Multiplied by von Karman constant and friction velocity to give turbulent
/// scalar velocity scale: w_s = kappa * u* * phi_s^{-1}(zeta)
/// Three-regime formulation per Large et al. (1994) Appendix B and CVMix.
/// Note: scalar and momentum exponents differ in the weakly-unstable regime.
///
/// @param zeta Monin-Obukhov stability coordinate (dimensionless)
/// @return phi_s^{-1} (dimensionless, > 0)
KOKKOS_INLINE_FUNCTION
Real KPPProfileS2(Real zeta) {
   if (zeta >= 0.0_Real) {
      // Stable regime
      return 1.0_Real / (1.0_Real + 5.0_Real * zeta);
   } else if (zeta >= ZETA_S) {
      // Weakly unstable: (1 - 16*zeta)^{1/2}  (scalar uses 1/2, not 1/4)
      return Kokkos::sqrt(1.0_Real - 16.0_Real * zeta);
   } else {
      // Strongly unstable: (a_s - c_s*zeta)^{1/3}
      return Kokkos::pow(A_MO_S - C_MO_S * zeta, 1.0_Real / 3.0_Real);
   }
}

/// @brief Hu(sigma) - Momentum surface value scaling
/// Sets surface boundary condition for momentum mixing
/// REFERENCES: Large et al. (1994)
///
/// @param sigma Normalized vertical position
/// @return Normalized profile value
KOKKOS_INLINE_FUNCTION
Real KPPHu(Real sigma) {
   // At sigma=0 (surface), should return HUON=0.03
   // Simple linear decay: Hu(sigma) = HUON * (1 + sigma)
   return HUON * (1.0 + sigma);
}

// ==========================================================================
// Langmuir Enhancement Factor (Theory-based Wave Model)
// ==========================================================================

/// @brief Estimate Stokes drift velocity scale from wind speed
/// Theory-based approach (no active wave data needed)
/// REFERENCES: Li et al. 2016, cvmix_kpp_ustokes_SL_model
///
/// @param wind10m Wind speed at 10 m height (m/s)
/// @param h_bl Boundary layer depth (m)
/// @return Stokes drift velocity scale (m/s)
KOKKOS_INLINE_FUNCTION
Real EstokesSLModel(Real wind10m, Real h_bl) {
   // u_s,BL = (u_10/362.0) * sqrt(2*alpha*cd) * lambda/h_bl
   // Simplified: alpha=0.84, cd ~ 1.2e-3, lambda ~ 2pi*g/w^2
   wind10m = Kokkos::fmax(0.0, wind10m);
   h_bl    = Kokkos::fmax(1.0, h_bl);

   // Approximate relation from Li et al.
   const Real C_drag     = 1.2e-3;
   const Real alpha_wave = 0.84;

   // Typical Stokes drift scale
   Real u_s = 0.016 * wind10m; // Simplified; sqrt(2*alpha*C_d)*(wind/g)

   return Kokkos::fmax(0.0, u_s);
}

/// @brief Langmuir number from friction velocity and Stokes drift
/// REFERENCES: Large et al. 2015 Eq. 6
///
/// @param u_star Friction velocity (m/s)
/// @param u_stokes Stokes drift at surface (m/s)
/// @return Langmuir number (dimensionless)
KOKKOS_INLINE_FUNCTION
Real ComputeLangmuirNumber(Real u_star, Real u_stokes) {
   u_star   = Kokkos::fmax(MIN_USTAR, u_star);
   u_stokes = Kokkos::fmax(1.0e-8, u_stokes);

   // La = sqrt(u_star / u_stokes)
   return Kokkos::sqrt(u_star / u_stokes);
}

/// @brief Langmuir enhancement factor for KPP from wind
/// Theory-based approach: depends on Langmuir number
/// REFERENCES: Li et al. (2016) Eq. 1-3
///
/// @param wind10m Wind speed at 10 m (m/s)
/// @param u_star Friction velocity (m/s)
/// @param h_bl Boundary layer depth (m)
/// @return Enhancement factor R_L (dimensionless, > 1.0 enhances mixing)
KOKKOS_INLINE_FUNCTION
Real ComputeEnhancementFactor(Real wind10m, Real u_star, Real h_bl) {
   u_star  = Kokkos::fmax(MIN_USTAR, u_star);
   wind10m = Kokkos::fmax(0.0, wind10m);
   h_bl    = Kokkos::fmax(1.0, h_bl);

   // Estimate Stokes drift from wind
   Real u_stokes = EstokesSLModel(wind10m, h_bl);

   // Compute Langmuir number
   Real la = ComputeLangmuirNumber(u_star, u_stokes);

   // Enhancement factor: R_L = sqrt(1 + 0.5 * (u_stokes/u_star)^2)
   // Alternative form based on Langmuir number:
   // R_L = sqrt(1 + 0.5 / La^2) for La > 0.5
   Real la_inv = 1.0 / Kokkos::fmax(0.5, la);
   Real r_l    = Kokkos::sqrt(1.0 + 0.5 * la_inv * la_inv);

   // Clamp to reasonable range [1.0, 2.0]
   return Kokkos::fmin(2.0, Kokkos::fmax(1.0, r_l));
}

// ==========================================================================
// Utility Functions for OBL Depth Computation
// ==========================================================================

/// @brief Check if a point should be suppressed (e.g., under ice)
/// Sets OBL to minimum if ice coverage or land ice present
///
/// @param ice_fraction Sea ice coverage (0-1)
/// @param land_ice_mask Land ice mask (0=ocean, non-zero=ice)
/// @return True if suppression applies
KOKKOS_INLINE_FUNCTION
bool ShouldSuppressOBL(Real ice_fraction, I4 land_ice_mask) {
   return (land_ice_mask != 0) || (ice_fraction > ICE_SUPPRESSION_THRESHOLD);
}

/// @brief Apply OBL depth constraints based on column properties
///
/// @param h_obl Current OBL depth (m)
/// @param layer_thickness Surface layer thickness (m)
/// @param water_depth Total water depth (m)
/// @param ice_fraction Sea ice coverage (0-1)
/// @return Constrained OBL depth (m)
KOKKOS_INLINE_FUNCTION
Real ConstrainOBLDepth(Real h_obl, Real layer_thickness, Real water_depth,
                       Real ice_fraction) {
   // Lower bound: at least half the surface layer thickness
   h_obl = Kokkos::fmax(h_obl, layer_thickness * 0.5);

   // Enforce minimum under ice
   if (ice_fraction > ICE_SUPPRESSION_THRESHOLD) {
      h_obl = Kokkos::fmax(h_obl, MIN_OBL_UNDER_ICE);
   }

   // Upper bound: cannot exceed water depth
   h_obl = Kokkos::fmin(h_obl, water_depth * 0.95);

   return h_obl;
}

// ==========================================================================
// Turbulent Velocity Scale Computation
// ==========================================================================

/// @brief Compute turbulent velocity scale (w_s)
/// Combined velocity scale for momentum and buoyancy
/// REFERENCES: Large et al. (1994) Eq. (9)-(10)
///
/// @param u_star Friction velocity (m/s)
/// @param b0 Surface buoyancy flux (m²/s³)
/// @param h_obl Boundary layer depth (m)
/// @return Turbulent velocity scale w_s (m/s)
KOKKOS_INLINE_FUNCTION
Real ComputeTurbulentVelocityScale(Real u_star, Real b0, Real h_obl) {
   u_star = Kokkos::fmax(0.0_Real, u_star);
   h_obl  = Kokkos::fmax(0.0_Real, h_obl);

   // w_s = (u_star^3 + 0.35 * b0 * h_obl)^(1/3)
   // Note: v_t = 0.35 in KPP (empirical constant)
   const Real v_t = 0.35;

   // Momentum contribution
   const Real w_m = u_star * u_star * u_star;

   // Buoyancy contribution for unstable (cooling/densifying) forcing.
   // In this sign convention, free convection corresponds to b0 < 0.
   const Real w_b = v_t * Kokkos::fmax(0.0_Real, -b0) * h_obl;

   // Combined scale
   const Real w_s =
       Kokkos::pow(Kokkos::fmax(0.0_Real, w_m + w_b), 1.0_Real / 3.0_Real);

   return w_s;
}

} // namespace OMEGA::KPP

#endif // OMEGA_KPP_CONSTANTS_H
