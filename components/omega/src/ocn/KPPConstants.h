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
// Monin-Obukhov stability function parameters (Large et al. 1994, App. B)
//
// The stability coordinate is zeta = d/L, where d is depth below the surface
// and L is the Monin-Obukhov length. zeta > 0 is stable (surface warming or
// salinification), zeta < 0 is unstable (convective). The values below set
// where the weakly-unstable branch hands off to the strongly-unstable branch.
// ==========================================================================

/// Transition zeta for momentum: below this value, strongly-unstable formula
/// is used. Default: -0.2 (CVMix default).
constexpr Real ZetaM = -0.2_Real;

/// Transition zeta for scalars: below this value, strongly-unstable formula
/// is used. Default: -1.0 (CVMix default).
constexpr Real ZetaS = -1.0_Real;

// The four constants below are fixed by requiring the strongly-unstable
// branch to match the weakly-unstable branch in value at the transition zeta.

/// a_m = (1-16*ZetaM)^{-0.25} * (1 - 4*ZetaM)
constexpr Real AMoM = 1.2573615702_Real;
/// c_m = (1-16*ZetaM)^{-0.25} * 12
constexpr Real CMoM = 8.3824104679_Real;

/// a_s = sqrt(1-16*ZetaS) * (1 + 8*ZetaS)  (can be negative)
constexpr Real AMoS = -28.8617393793_Real;
/// c_s = 24 * sqrt(1-16*ZetaS)
constexpr Real CMoS = 98.9545350148_Real;

/// Surface value of the momentum shape function (dimensionless)
constexpr Real HuOn = 0.03;

/// Floor on friction velocity (m/s). Keeps the turbulent velocity scales and
/// the Monin-Obukhov length finite in near-calm conditions.
constexpr Real MinUStar = 1.0e-4;

// ==========================================================================
// OBL Depth Computation Parameters
// ==========================================================================

/// Safety multiplier for OBL search (prevents searching too deep)
/// Default: 1.0 (search to 1.0 * Ri_crit threshold)
constexpr Real StopOBLSearchMult = 1.0;

/// Minimum OBL depth under sea ice (m), applied above IceSuppressThresh
constexpr Real MinOBLUnderIce = 5.0;

/// Ice fraction above which Langmuir enhancement is disabled
constexpr Real IceFracThresh = 0.05;

/// Ice fraction above which the minimum OBL depth is enforced
constexpr Real IceSuppressThresh = 0.15;

// ==========================================================================
// Surface Layer Parameters
// ==========================================================================

/// Surface layer extent as a fraction of the trial OBL depth (epsilon in
/// Large et al. 1994). Reference values entering the bulk Richardson number
/// are averaged over the top SurfaceLayerExtent * d of the column.
constexpr Real SurfaceLayerExtent = 0.1;

/// Empirical convective velocity coefficient in the turbulent velocity scale
/// w_s = (u*^3 + ConvectiveVelFac * max(-B_0,0) * h)^(1/3)
constexpr Real ConvectiveVelFac = 0.35;

// ==========================================================================
// KPP Shape and Stability Functions
//
// Vertical position inside the OBL is expressed two ways:
//   Sigma    in [-1,0], the Omega convention: 0 at the surface, -1 at the
//            OBL base (it follows the sign of the z coordinate).
//   SigmaMu  in [ 0,1], the CVMix/Large et al. convention: SigmaMu = -Sigma,
//            so 0 at the surface and 1 at the OBL base.
// The shape functions below take Sigma and convert internally, so callers
// never need to flip signs.
// ==========================================================================

/// @brief Momentum gradient shape function, multiplied by h and the turbulent
/// velocity scale to give the KPP viscosity: Kx = h * w_m * G(sigma).
/// REFERENCES: Large et al. (1994) Eq. (11)
///
/// @param Sigma Normalized vertical position (-z/h)
/// @return Dimensionless shape value
KOKKOS_INLINE_FUNCTION
Real kppShapeMomentum(Real Sigma) {
   Sigma = Kokkos::fmax(-1.0, Kokkos::fmin(0.0, Sigma));

   const Real SigmaMu = -Sigma;
   return SigmaMu * (1.0 - SigmaMu) * (1.0 - SigmaMu);
}

/// @brief Matched KPP gradient shape function.
///
/// Uses the SimpleShapes gradient shape plus a smooth correction that is
/// zero at the surface and equals ShapeAtBase at the OBL base. This lets
/// MatchBoth profiles meet pre-existing interior mixing at the BLD base while
/// preserving SimpleShapes behavior when ShapeAtBase is zero.
KOKKOS_INLINE_FUNCTION
Real kppShapeMatched(Real Sigma, Real ShapeAtBase) {
   Sigma = Kokkos::fmax(-1.0, Kokkos::fmin(0.0, Sigma));

   const Real SigmaMu = -Sigma;
   const Real Simple  = SigmaMu * (1.0 - SigmaMu) * (1.0 - SigmaMu);
   const Real Smooth  = SigmaMu * SigmaMu * (3.0 - 2.0 * SigmaMu);
   return Simple + ShapeAtBase * Smooth;
}

/// @brief phi_m^{-1}(zeta) - Inverse momentum Monin-Obukhov stability function
/// Multiplied by von Karman constant and friction velocity to give turbulent
/// momentum velocity scale: w_m = kappa * u* * phi_m^{-1}(zeta)
/// Three-regime formulation per Large et al. (1994) Appendix B and CVMix.
///
/// @param Zeta Monin-Obukhov stability coordinate d/L (dimensionless)
/// @return phi_m^{-1} (dimensionless, > 0)
KOKKOS_INLINE_FUNCTION
Real kppPhiInvMomentum(Real Zeta) {
   if (Zeta >= 0.0_Real) {
      // Stable regime
      return 1.0_Real / (1.0_Real + 5.0_Real * Zeta);
   } else if (Zeta >= ZetaM) {
      // Weakly unstable: (1 - 16*zeta)^{1/4}
      return Kokkos::pow(1.0_Real - 16.0_Real * Zeta, 0.25_Real);
   } else {
      // Strongly unstable (convective): (a_m - c_m*zeta)^{1/3}
      return Kokkos::pow(AMoM - CMoM * Zeta, 1.0_Real / 3.0_Real);
   }
}

/// @brief Scalar gradient shape function, multiplied by h and the turbulent
/// velocity scale to give the KPP diffusivity: Kx = h * w_s * G(sigma).
/// The non-local flux reuses this same shape, scaled by C_s instead of h*w_s.
/// REFERENCES: Large et al. (1994) Eq. (11), Eq. (12)-(13), Large et al. (1997)
///
/// @param Sigma Normalized vertical position (-z/h)
/// @return Dimensionless shape value
KOKKOS_INLINE_FUNCTION
Real kppShapeScalar(Real Sigma) {
   Sigma = Kokkos::fmax(-1.0, Kokkos::fmin(0.0, Sigma));

   const Real SigmaMu = -Sigma;
   return SigmaMu * (1.0 - SigmaMu) * (1.0 - SigmaMu);
}

/// @brief phi_s^{-1}(zeta) - Inverse scalar Monin-Obukhov stability function
/// Multiplied by von Karman constant and friction velocity to give turbulent
/// scalar velocity scale: w_s = kappa * u* * phi_s^{-1}(zeta)
/// Three-regime formulation per Large et al. (1994) Appendix B and CVMix.
/// Scalars mix more efficiently than momentum in unstable conditions, which
/// is why the weakly-unstable exponent is 1/2 here and 1/4 for momentum.
///
/// @param Zeta Monin-Obukhov stability coordinate d/L (dimensionless)
/// @return phi_s^{-1} (dimensionless, > 0)
KOKKOS_INLINE_FUNCTION
Real kppPhiInvScalar(Real Zeta) {
   if (Zeta >= 0.0_Real) {
      // Stable regime
      return 1.0_Real / (1.0_Real + 5.0_Real * Zeta);
   } else if (Zeta >= ZetaS) {
      // Weakly unstable: (1 - 16*zeta)^{1/2}
      return Kokkos::sqrt(1.0_Real - 16.0_Real * Zeta);
   } else {
      // Strongly unstable (convective): (a_s - c_s*zeta)^{1/3}
      return Kokkos::pow(AMoS - CMoS * Zeta, 1.0_Real / 3.0_Real);
   }
}

/// @brief Hu(sigma) - Momentum surface value scaling
/// Sets surface boundary condition for momentum mixing
/// REFERENCES: Large et al. (1994)
///
/// @param Sigma Normalized vertical position
/// @return Normalized shape value
KOKKOS_INLINE_FUNCTION
Real kppSurfaceMomentumScale(Real Sigma) {
   // Linear decay from HuOn at the surface to zero at the OBL base
   return HuOn * (1.0 + Sigma);
}

// ==========================================================================
// Langmuir Enhancement Factor (Theory-based Wave Model)
//
// Langmuir circulations are wind-and-wave driven counter-rotating vortices
// that deepen and strengthen boundary layer mixing beyond what wind stress
// alone produces. Their strength is measured by the turbulent Langmuir
// number La = sqrt(u* / u_stokes): small La means wave forcing dominates.
// With no wave model coupled, the Stokes drift is estimated from the 10 m
// wind following Li et al. (2016).
// ==========================================================================

/// @brief Estimate the surface-layer Stokes drift velocity scale from wind
/// REFERENCES: Li et al. 2016, cvmix_kpp_ustokes_SL_model
///
/// TODO: this is a placeholder empirical fit, not the full surface-layer
/// averaged Stokes drift of the reference model, which needs the wave
/// spectrum (or a wave component) to evaluate
///   u_s,BL = (U10/362) * sqrt(2*alpha*Cd) * lambda/HBL
/// with alpha = 0.84, Cd ~ 1.2e-3 and lambda ~ 2*pi*g/omega^2. Until those
/// inputs are available, HBL is accepted but unused.
///
/// @param Wind10m Wind speed at 10 m height (m/s)
/// @param HBL Boundary layer depth (m)
/// @return Stokes drift velocity scale (m/s)
KOKKOS_INLINE_FUNCTION
Real estimateStokesDriftSL(Real Wind10m, Real HBL) {
   Wind10m = Kokkos::fmax(0.0, Wind10m);
   HBL     = Kokkos::fmax(1.0, HBL);

   const Real UStokes = 0.016 * Wind10m;

   return Kokkos::fmax(0.0, UStokes);
}

/// @brief Turbulent Langmuir number from friction velocity and Stokes drift
/// REFERENCES: Large et al. 2015 Eq. 6
///
/// @param UStar Friction velocity (m/s)
/// @param UStokes Stokes drift at surface (m/s)
/// @return Langmuir number (dimensionless)
KOKKOS_INLINE_FUNCTION
Real computeLangmuirNumber(Real UStar, Real UStokes) {
   UStar   = Kokkos::fmax(MinUStar, UStar);
   UStokes = Kokkos::fmax(1.0e-8, UStokes);

   return Kokkos::sqrt(UStar / UStokes);
}

/// @brief Langmuir enhancement factor applied to the KPP velocity scales
/// REFERENCES: Li et al. (2016) Eq. 1-3
///
/// @param Wind10m Wind speed at 10 m (m/s)
/// @param UStar Friction velocity (m/s)
/// @param HBL Boundary layer depth (m)
/// @return Enhancement factor R_L (dimensionless, > 1.0 enhances mixing)
KOKKOS_INLINE_FUNCTION
Real computeLangmuirEnhancement(Real Wind10m, Real UStar, Real HBL) {
   UStar   = Kokkos::fmax(MinUStar, UStar);
   Wind10m = Kokkos::fmax(0.0, Wind10m);
   HBL     = Kokkos::fmax(1.0, HBL);

   const Real UStokes = estimateStokesDriftSL(Wind10m, HBL);
   const Real La      = computeLangmuirNumber(UStar, UStokes);

   // R_L = sqrt(1 + 0.5/La^2); La is floored at 0.5 so that the weak-wave
   // limit returns the unenhanced scales rather than diverging.
   const Real LaInv = 1.0 / Kokkos::fmax(0.5, La);
   const Real RL    = Kokkos::sqrt(1.0 + 0.5 * LaInv * LaInv);

   return Kokkos::fmin(2.0, Kokkos::fmax(1.0, RL));
}

// ==========================================================================
// Utility Functions for OBL Depth Computation
// ==========================================================================

/// @brief Check if a point should be suppressed (e.g., under ice)
/// Sets OBL to minimum if ice coverage or land ice present
///
/// @param IceFrac Sea ice coverage (0-1)
/// @param LandIceMask Land ice mask (0=ocean, non-zero=ice)
/// @return True if suppression applies
KOKKOS_INLINE_FUNCTION
bool shouldSuppressOBL(Real IceFrac, I4 LandIceMask) {
   return (LandIceMask != 0) || (IceFrac > IceSuppressThresh);
}

/// @brief Apply OBL depth constraints based on column properties
///
/// @param HOBL Current OBL depth (m)
/// @param LayerThickness Surface layer thickness (m)
/// @param WaterDepth Total water depth (m)
/// @param IceFrac Sea ice coverage (0-1)
/// @return Constrained OBL depth (m)
KOKKOS_INLINE_FUNCTION
Real constrainOBLDepth(Real HOBL, Real LayerThickness, Real WaterDepth,
                       Real IceFrac) {
   // Lower bound: at least half the surface layer thickness
   HOBL = Kokkos::fmax(HOBL, LayerThickness * 0.5);

   // Enforce minimum under ice
   if (IceFrac > IceSuppressThresh) {
      HOBL = Kokkos::fmax(HOBL, MinOBLUnderIce);
   }

   // Upper bound: cannot exceed water depth
   HOBL = Kokkos::fmin(HOBL, WaterDepth * 0.95);

   return HOBL;
}

// ==========================================================================
// Turbulent Velocity Scale Computation
// ==========================================================================

/// @brief Compute the depth-independent turbulent velocity scale
/// Blends the shear-driven scale u* with the convective scale so that the
/// result stays finite in both the wind-driven and free-convection limits.
/// REFERENCES: Large et al. (1994) Eq. (9)-(10)
///
/// @param UStar Friction velocity (m/s)
/// @param BuoyFlux Surface buoyancy flux (m^2/s^3), negative when convective
/// @param HOBL Boundary layer depth (m)
/// @return Turbulent velocity scale w_s (m/s)
KOKKOS_INLINE_FUNCTION
Real computeTurbVelocityScale(Real UStar, Real BuoyFlux, Real HOBL) {
   UStar = Kokkos::fmax(0.0_Real, UStar);
   HOBL  = Kokkos::fmax(0.0_Real, HOBL);

   // Momentum contribution
   const Real WMom = UStar * UStar * UStar;

   // Buoyancy contribution for unstable (cooling/densifying) forcing.
   // In this sign convention, free convection corresponds to BuoyFlux < 0.
   const Real WBuoy =
       ConvectiveVelFac * Kokkos::fmax(0.0_Real, -BuoyFlux) * HOBL;

   return Kokkos::pow(Kokkos::fmax(0.0_Real, WMom + WBuoy),
                      1.0_Real / 3.0_Real);
}

} // namespace OMEGA::KPP

#endif // OMEGA_KPP_CONSTANTS_H
