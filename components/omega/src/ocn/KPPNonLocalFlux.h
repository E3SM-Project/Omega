#ifndef OMEGA_KPP_NONLOCAL_FLUX_H
#define OMEGA_KPP_NONLOCAL_FLUX_H
//===-- ocn/KPPNonLocalFlux.h - Non-local Flux Computation -----*- C++ -*-===//
//
/// \file
/// \brief Compute KPP non-local flux profiles for tracers
///
/// This header defines functors for computing the non-local flux coefficient
/// G(σ) which is applied to surface tracer fluxes to produce vertical mixing
/// of tracers. The non-local flux represents transport by coherent plumes
/// within the OBL.
//
//===----------------------------------------------------------------------===//

#include "KPPConstants.h"
#include "OmegaKokkos.h"

namespace OMEGA::KPP {

/// @brief Non-local flux profile functor
/// Computes G(σ) applied to surface tracer fluxes
///
/// The non-local flux produces vertical transport:
///   flux(z) = G(σ) × Q_surf
/// where σ = -z/h_OBL (normalized depth)
///
/// REFERENCES: Large et al. (1994) Eq. (12)-(13), Large et al. (1997)
class KPPComputeNonLocalFlux {

 public:
   Array1DReal ZInterface; ///< Depth at interfaces (m, negative down)
   Array1DReal ZCenter;    ///< Depth at cell centers (m)
   Array1DI4 MinLayerCell; ///< Min layer index per cell
   Array1DI4 MaxLayerCell; ///< Max layer index per cell

   // OBL depth information
   Real OBLDepth; ///< Current OBL depth (m)
   I4 OBLIndex;   ///< Layer index of OBL base

   // Reference profiles for shear stability correction
   Array1DReal GradientRichardsonNum; ///< Ri_g for stability correction

   // Output
   Array1DReal NonLocalFluxProfile; ///< G(σ) values at interfaces

   /// @brief Constructor
   KPPComputeNonLocalFlux(const Array1DReal &ZInterface_in,
                          const Array1DReal &ZCenter_in,
                          const Array1DI4 &MinLayerCell_in,
                          const Array1DI4 &MaxLayerCell_in, Real obl_depth,
                          I4 obl_index, const Array1DReal &RiGrad_in,
                          const Array1DReal &G_profile_out)
       : ZInterface(ZInterface_in), ZCenter(ZCenter_in),
         MinLayerCell(MinLayerCell_in), MaxLayerCell(MaxLayerCell_in),
         OBLDepth(obl_depth), OBLIndex(obl_index),
         GradientRichardsonNum(RiGrad_in), NonLocalFluxProfile(G_profile_out) {}

   /// @brief Compute non-local flux profile G(σ)
   ///
   /// Algorithm:
   /// 1. For each layer k from surface to OBL base:
   ///    a. Compute normalized depth σ = -z/h_OBL
   ///    b. Evaluate G(σ) profile function
   ///    c. Apply stability correction if needed
   /// 2. Set G(σ) = 0 below OBL base
   ///
   KOKKOS_FUNCTION
   void computeNonLocalFlux(I4 ICell) const {

      const I4 KMin = MinLayerCell(ICell);
      const I4 KMax = MaxLayerCell(ICell);

      // Clamp OBL depth to reasonable bounds
      Real h_obl = Kokkos::fmax(1.0, OBLDepth);

      // =======================================================================
      // Compute G(σ) at each interface
      // =======================================================================
      for (I4 k = KMin; k <= KMax + 1; ++k) {

         Real z_interface = Kokkos::abs(ZInterface(k));

         // Check if point is within OBL
         if (z_interface <= h_obl) {

            // Normalized depth: σ = -z/h (negative in ocean convention)
            Real sigma = -(z_interface / h_obl); // -1 <= sigma <= 0

            // Evaluate G(σ) profile
            Real g_sigma = KPPProfileG(sigma);

            NonLocalFluxProfile(k) = g_sigma;

         } else {
            // Below OBL base: no non-local flux
            NonLocalFluxProfile(k) = 0.0;
         }
      }
   }

}; // class KPPComputeNonLocalFlux

} // namespace OMEGA::KPP

#endif // OMEGA_KPP_NONLOCAL_FLUX_H
