#ifndef OMEGA_PGRAD_FINITE_VOLUME_H
#define OMEGA_PGRAD_FINITE_VOLUME_H
//===-- ocn/PGradFiniteVolume.h - Finite-Volume PGrad ---------*- C++ -*-===//
///
/// The equation-of-state expansion shared across an edge, used by the
/// FiniteVolume pressure gradient.
///
/// Specific volume is never integrated directly. For each cell and layer it is
/// expanded to first order about a reference state,
///
///    alpha(Ct, Sa, p) = alpha0 + alphaCt (Ct - Ct0) + alphaSa (Sa - Sa0)
///                              + alphaP  (p  - p0)
///
/// with the four coefficients coming from a single equation-of-state
/// evaluation per cell per layer -- the SpecVol, SpecVolDCt, SpecVolDSa and
/// SpecVolDP fields Eos already computes. Because temperature and salinity are
/// reconstructed as low-order polynomials in pressure and the expansion is
/// linear in them, the expanded specific volume is a low-order polynomial in
/// pressure and every integral the scheme needs is available in closed form.
/// No equation-of-state evaluation occurs inside any integral, which is what
/// bounds the cost.
///
/// The expansion point and the coefficients are shared across each edge:
/// averaged from the edge's two cells and used for *both* columns. This is
/// what the robustness property rests on. Give each column its own expansion
/// point and the two columns describe the same water with two slightly
/// different approximations to the same equation of state, their alpha0 and
/// alphaP terms no longer cancel in the matched-pressure difference, and the
/// scheme generates spurious flow out of nothing but its own equation-of-state
/// approximation.
///
/// Note what is and is not load-bearing here. That *one* set multiplies both
/// columns is essential. *Which* set it is -- selected by edge layer, or by
/// some other rule -- is an ordinary accuracy question and cannot break
/// exactness, because the coefficients end up multiplying a quantity that is
/// identically zero on the exact set.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "OmegaKokkos.h"

namespace OMEGA {

/// The equation-of-state expansion shared by both columns of an edge layer.
/// One instance multiplies both columns' contributions.
struct PGradEdgeEos {
   Real SpecVol0;    ///< specific volume at the shared expansion state
   Real SpecVolDCt;  ///< d(alpha)/d(ConservTemp) at the shared state
   Real SpecVolDSa;  ///< d(alpha)/d(AbsSalinity) at the shared state
   Real SpecVolDP;   ///< d(alpha)/d(Pressure) at the shared state
   Real ConservTemp; ///< shared expansion state, conservative temperature
   Real AbsSalinity; ///< shared expansion state, absolute salinity
   Real Press;       ///< shared expansion state, pressure
};

/// Build the shared expansion for one edge layer by averaging the two
/// adjacent cells' coefficients and reference states. The pressure of the
/// shared state is the edge average of the two columns' mid-layer pressures,
/// which is exactly the mid-layer pressure of the edge control volume when
/// that volume is taken as the average of the two columns' interface
/// pressures.
KOKKOS_INLINE_FUNCTION PGradEdgeEos buildEdgeEos(
    const Array2DReal &SpecVol,     ///< [in] specific volume
    const Array2DReal &SpecVolDCt,  ///< [in] d(alpha)/d(ConservTemp)
    const Array2DReal &SpecVolDSa,  ///< [in] d(alpha)/d(AbsSalinity)
    const Array2DReal &SpecVolDP,   ///< [in] d(alpha)/d(Pressure)
    const Array2DReal &ConservTemp, ///< [in] layer-mean temperature
    const Array2DReal &AbsSalinity, ///< [in] layer-mean salinity
    const Array2DReal &PressureMid, ///< [in] mid-layer pressure
    const I4 ICell0,                ///< [in] first cell on the edge
    const I4 ICell1,                ///< [in] second cell on the edge
    const I4 K                      ///< [in] edge layer
) {
   PGradEdgeEos Eos;
   Eos.SpecVol0   = 0.5_Real * (SpecVol(ICell0, K) + SpecVol(ICell1, K));
   Eos.SpecVolDCt = 0.5_Real * (SpecVolDCt(ICell0, K) + SpecVolDCt(ICell1, K));
   Eos.SpecVolDSa = 0.5_Real * (SpecVolDSa(ICell0, K) + SpecVolDSa(ICell1, K));
   Eos.SpecVolDP  = 0.5_Real * (SpecVolDP(ICell0, K) + SpecVolDP(ICell1, K));
   Eos.ConservTemp =
       0.5_Real * (ConservTemp(ICell0, K) + ConservTemp(ICell1, K));
   Eos.AbsSalinity =
       0.5_Real * (AbsSalinity(ICell0, K) + AbsSalinity(ICell1, K));
   Eos.Press = 0.5_Real * (PressureMid(ICell0, K) + PressureMid(ICell1, K));
   return Eos;
}

/// One column's expanded specific volume at a pressure, using the edge-shared
/// expansion. Needed only where a single column's specific volume is
/// integrated on its own, which in Phase 1 is only the anchor of the column
/// scan; the interior of the scheme differences two columns at matched
/// pressure, where the SpecVol0 and SpecVolDP terms cancel and never have to
/// be formed.
KOKKOS_INLINE_FUNCTION Real edgeSpecVol(
    const PGradEdgeEos &Eos, ///< [in] edge-shared expansion
    const Real ConservTemp,  ///< [in] this column's temperature at p
    const Real AbsSalinity,  ///< [in] this column's salinity at p
    const Real Press         ///< [in] pressure
) {
   return Eos.SpecVol0 + Eos.SpecVolDCt * (ConservTemp - Eos.ConservTemp) +
          Eos.SpecVolDSa * (AbsSalinity - Eos.AbsSalinity) +
          Eos.SpecVolDP * (Press - Eos.Press);
}

} // namespace OMEGA

//===----------------------------------------------------------------------===//
#endif
