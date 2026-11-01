#ifndef OMEGA_PGRAD_RECON_H
#define OMEGA_PGRAD_RECON_H
//===-- ocn/PGradRecon.h - Pressure Gradient Reconstruction ---*- C++ -*-===//
///
/// Mean-preserving vertical reconstruction of conservative temperature and
/// absolute salinity in pressure, used by the FiniteVolume pressure gradient.
///
/// Within layer K of a column the prognostic layer mean is supplemented by a
/// deviation that integrates to zero over the layer,
///
///    Theta(p) = Theta_k + Slope_k * (p - PressureMid_k)
///
/// The linear deviations make the pressure gradient exact for profiles that
/// vary linearly with pressure; the parabolic ones further down this file add a
/// curvature term and make it exact for profiles that vary quadratically.
/// Neither has a limiter: the reconstruction feeds an integral rather than an
/// advective flux, and a limiter active on smooth data would break the
/// cancellation the scheme rests on precisely where the profile is well
/// resolved.
///
/// Two properties carry the weight, described here for the linear
/// reconstruction.
///
/// Mean-preserving. PressureMid is the exact arithmetic midpoint of the two
/// interface pressures, so the integral of (p - PressureMid) over the layer is
/// zero and the deviation cannot move the layer mean. This is also what makes
/// VertCoord's midpoint rule the exact layer integral of a linear
/// reconstruction, and hence what spares VertCoord any change.
///
/// Exact on a non-uniform grid. For a profile linear in pressure the layer
/// means lie exactly on the line as a function of mid-layer pressure, because
/// the mean of a linear function over a layer is its value at the layer's
/// midpoint. A centered difference of the layer means with respect to
/// PressureMid therefore recovers the slope exactly, for any distribution of
/// layer thicknesses. This is why the estimator differences against the actual
/// mid-layer pressures and not against layer index: a formula that assumed
/// uniform thickness would pass a uniform-thickness test and fail on Omega's
/// grid.
///
/// The functions here take scalars, or are templated on the array type, so
/// that they can be unit tested without a mesh.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "OmegaKokkos.h"

namespace OMEGA {

/// The two layers the linear slope estimator differences, for layer K of a
/// column whose valid layers run from KMin to KMax. The difference is centered
/// in the interior and one-sided in the shallowest and deepest valid layer.
/// The one-sided branches are not an edge case: the deepest valid layer is a
/// partial cell, which carries the whole signal where bathymetry steps. Where
/// a column has a single valid layer both indices are K, which the slope
/// estimator turns into a zero slope -- a constant being the only
/// mean-preserving reconstruction available there.
KOKKOS_INLINE_FUNCTION void
linearReconStencil(const I4 K,    ///< [in] layer to reconstruct
                   const I4 KMin, ///< [in] shallowest valid layer
                   const I4 KMax, ///< [in] deepest valid layer
                   I4 &KLo,       ///< [out] shallower layer to difference
                   I4 &KHi        ///< [out] deeper layer to difference
) {
   KLo = (K > KMin) ? K - 1 : K;
   KHi = (K < KMax) ? K + 1 : K;
}

/// The slope, per Pa, of the mean-preserving linear reconstruction, from the
/// layer means and mid-layer pressures of the two layers linearReconStencil
/// selects. Exact for a profile linear in pressure on an arbitrary grid and
/// second-order accurate otherwise. Returns zero where the two mid-layer
/// pressures coincide, which happens only when the column has a single valid
/// layer.
KOKKOS_INLINE_FUNCTION Real linearReconSlope(
    const Real ValueLo,    ///< [in] layer mean at KLo
    const Real ValueHi,    ///< [in] layer mean at KHi
    const Real PressMidLo, ///< [in] mid-layer pressure at KLo
    const Real PressMidHi  ///< [in] mid-layer pressure at KHi
) {
   const Real DeltaPress = PressMidHi - PressMidLo;
   return (DeltaPress != 0.0_Real) ? (ValueHi - ValueLo) / DeltaPress
                                   : 0.0_Real;
}

/// The deviation of the reconstruction from the layer mean at a pressure.
/// Integrates to zero over the layer because PressMid is the exact arithmetic
/// midpoint of the two interface pressures. The pressure need not lie within
/// the layer: near the sea floor a column's deepest reconstruction is
/// evaluated slightly below its own floor, where this extrapolates.
KOKKOS_INLINE_FUNCTION Real linearReconDeviation(
    const Real Slope,    ///< [in] slope from
                         ///< linearReconSlope
    const Real PressMid, ///< [in] mid-layer pressure
    const Real Press     ///< [in] pressure to evaluate at
) {
   return Slope * (Press - PressMid);
}

/// The reconstruction evaluated at a pressure. Equal to the layer mean at
/// PressMid, by construction.
KOKKOS_INLINE_FUNCTION Real linearReconEval(
    const Real Value,    ///< [in] prognostic layer mean
    const Real Slope,    ///< [in] slope from linearReconSlope
    const Real PressMid, ///< [in] mid-layer pressure
    const Real Press     ///< [in] pressure to evaluate at
) {
   return Value + linearReconDeviation(Slope, PressMid, Press);
}

/// The layer of column ICell whose interfaces bracket the pressure Press,
/// searched from the hint KHint and clamped to the column's valid range
/// [KMin, KMax].
///
/// This is the lookup that makes the scheme a fixed-pressure comparison rather
/// than a fixed-layer-index one, and it is the single place where that
/// distinction lives. At a pressure in edge layer K, each column supplies its
/// own temperature and salinity from whichever of *its own* layers contains
/// that pressure, which under tilt is generally not layer K -- at a coordinate
/// tilt of 50 m/km with 64 m layers the two columns' layer K are offset by
/// nearly three layer thicknesses and do not overlap in pressure at all.
/// Looking a column's state up by layer index instead would silently
/// reintroduce the very error the scheme exists to remove, and would still
/// pass every exactness and accuracy check available, which is why the
/// property tests on this function are not optional.
///
/// Interface pressures increase downward, so the search walks up while Press
/// lies above the layer's top interface and down while it lies below the
/// bottom one. Within the column scan the answer advances monotonically with
/// K, so passing the previous answer as KHint makes this a pair of incremented
/// cursors rather than a search; the result does not depend on the hint.
///
/// Where Press lies outside the column altogether the outermost valid layer is
/// returned and its reconstruction is extrapolated. This is the rule at the
/// top and bottom of the column: the edge control volume's pressure range is
/// in general neither column's own, so near the sea floor a column's deepest
/// reconstruction must be evaluated slightly below its own floor. Exactness
/// survives extrapolation, because on the exact set an extrapolated
/// reconstruction still reproduces the true profile.
///
/// Templated on the array type so that the same code is exercised by the
/// device kernels and by the host-side property tests that pin it.
template <class ArrayType>
KOKKOS_INLINE_FUNCTION I4
findLayerForPress(const ArrayType &PressInterface, ///< [in] interface pressures
                  const I4 ICell,                  ///< [in] cell to search
                  const I4 KMin,    ///< [in] shallowest valid layer
                  const I4 KMax,    ///< [in] deepest valid layer
                  const Real Press, ///< [in] pressure to locate
                  const I4 KHint    ///< [in] starting guess
) {
   I4 K = KHint;
   if (K < KMin)
      K = KMin;
   if (K > KMax)
      K = KMax;

   while (K > KMin && Press < PressInterface(ICell, K))
      --K;
   while (K < KMax && Press > PressInterface(ICell, K + 1))
      ++K;

   return K;
}

//------------------------------------------------------------------------------
// The parabolic reconstruction
//
// The piecewise parabolic method of Colella and Woodward (1984) without its
// limiter, as design section 3.4 specifies it. Within layer K the deviation is
//
//    Slope_k * (p - PressureMid_k)
//       + Curv_k * ((p - PressureMid_k)^2 - Thick_k^2 / 12)
//
// which integrates to zero over the layer for any Slope and Curv, again because
// PressureMid is the exact midpoint. The parabola keeps the layer mean and
// takes, at each of the layer's two interfaces, the derivative there of the
// polynomial P through the cumulative integral of the profile at five
// consecutive interfaces centred on it. The layer means are exact layer
// averages, so P reproduces the cumulative integral of any profile cubic in
// pressure whatever the layer thicknesses, the interface values are then exact,
// and the reconstruction is exact for profiles quadratic in pressure.

/// The most interfaces the parabolic estimator interpolates through: two
/// layers on either side of the interface.
constexpr int MaxPpmWindow = 5;

/// The first and last interface, JLo and JHi, of the window the parabolic
/// estimator interpolates through for the value at interface J of a column
/// whose valid layers run from KMin to KMax: five consecutive interfaces
/// centred on J, shifted inward near the ends of the column so that they stay
/// among its valid layers, or all of the column's interfaces where it has
/// fewer than four valid layers. P then has as high a degree as the column's
/// layers allow -- in a column of two layers it is the polynomial whose
/// derivative is the linear reconstruction, and in one of a single layer the
/// constant.
KOKKOS_INLINE_FUNCTION void ppmReconWindow(const I4 J,    ///< [in] interface
                                           const I4 KMin, ///< [in] shallowest
                                                          ///< valid layer
                                           const I4 KMax, ///< [in] deepest
                                                          ///< valid layer
                                           I4 &JLo, ///< [out] first interface
                                           I4 &JHi  ///< [out] last interface
) {
   if (KMax - KMin + 1 < MaxPpmWindow - 1) {
      JLo = KMin;
      JHi = KMax + 1;
      return;
   }
   JLo = J - 2;
   if (JLo < KMin)
      JLo = KMin;
   if (JLo > KMax + 2 - MaxPpmWindow)
      JLo = KMax + 2 - MaxPpmWindow;
   JHi = JLo + MaxPpmWindow - 1;
}

/// The first four derivatives at the pressure Press of the polynomial through
/// the cumulative integral of column ICell's layer means, less RefValue, at
/// the interfaces JLo to JHi.
///
/// No cumulative integral is ever summed. Its first divided differences
/// between neighbouring interfaces are the layer means themselves, so the
/// Newton form of the polynomial is built from them directly. Subtracting
/// RefValue shifts the first derivative by that constant and leaves the others
/// unchanged; taking it to be the mean of the layer being reconstructed makes
/// the round-off scale with the variation of the profile across the window
/// rather than with the values themselves or with the depth of the column.
template <class ArrayType>
KOKKOS_INLINE_FUNCTION void
ppmReconDerivs(const ArrayType &Values,         ///< [in] layer means
               const ArrayType &PressInterface, ///< [in] interface pressures
               const I4 ICell,                  ///< [in] cell
               const I4 JLo,                    ///< [in] first interface
               const I4 JHi,                    ///< [in] last interface
               const Real Press,    ///< [in] pressure to evaluate at
               const Real RefValue, ///< [in] value subtracted from the means
               Real Derivs[4]       ///< [out] first to fourth derivatives
) {
   const int NLayers = JHi - JLo;

   // the interfaces relative to the evaluation point, and the
   // divided-difference table, whose first column is the layer means
   Real Nodes[MaxPpmWindow];
   Real Table[MaxPpmWindow - 1];
   for (int I = 0; I <= NLayers; ++I)
      Nodes[I] = PressInterface(ICell, JLo + I) - Press;
   for (int I = 0; I < NLayers; ++I)
      Table[I] = Values(ICell, JLo + I) - RefValue;

   // the coefficients, in powers of s = p - Press, of the Newton basis
   // polynomial of the current degree, the product of (s - Nodes[I]) over the
   // first Degree nodes
   Real Basis[MaxPpmWindow] = {1.0_Real, 0.0_Real, 0.0_Real, 0.0_Real,
                               0.0_Real};
   for (int Q = 0; Q < 4; ++Q)
      Derivs[Q] = 0.0_Real;

   for (int Degree = 1; Degree <= NLayers; ++Degree) {
      if (Degree > 1) {
         for (int I = 0; I <= NLayers - Degree; ++I)
            Table[I] =
                (Table[I + 1] - Table[I]) / (Nodes[I + Degree] - Nodes[I]);
      }
      for (int Q = Degree; Q > 0; --Q)
         Basis[Q] = Basis[Q - 1] - Nodes[Degree - 1] * Basis[Q];
      Basis[0] = -Nodes[Degree - 1] * Basis[0];

      // the Q-th derivative at s = 0 is Q! times the coefficient of s^Q
      Real Factorial = 1.0_Real;
      for (int Q = 1; Q <= Degree && Q <= 4; ++Q) {
         Factorial *= Q;
         Derivs[Q - 1] += Table[0] * Factorial * Basis[Q];
      }
   }
}

/// The slope and curvature of the parabolic reconstruction of layer K of
/// column ICell, whose valid layers run from KMin to KMax.
///
/// Where the layer's two interfaces take their values from different windows
/// -- the interior of a column -- the parabola is the one through the two
/// interface values that keeps the layer mean,
///
///    Slope = (ABot - ATop) / Thick,
///    Curv  = 3 * (ATop + ABot - 2 * Value) / Thick^2.
///
/// Where both come from the same window -- the two outermost valid layers at
/// each end of a column, and every layer of a column of four or fewer -- the
/// coefficients are the polynomial's own derivatives at the layer midpoint,
///
///    Slope = P''(PressMid) + Thick^2 * P''''(PressMid) / 24,
///    Curv  = P'''(PressMid) / 2,
///
/// which is the same thing rearranged. The two agree in exact arithmetic, but
/// in a thin layer the first forms the curvature from a small difference of
/// nearly equal values and divides it by the square of the thickness. That does
/// no harm inside the layer and a great deal beyond it, and the outermost
/// layers are the ones extrapolated past the ends of the column: built the
/// first way, a single-precision bottom cell a ten-thousandth of a layer thick
/// extrapolates with an error of order 10 degrees C.
template <class ArrayType>
KOKKOS_INLINE_FUNCTION void
ppmReconCoeffs(const ArrayType &Values,         ///< [in] layer means
               const ArrayType &PressInterface, ///< [in] interface pressures
               const I4 ICell,                  ///< [in] cell
               const I4 K,                      ///< [in] layer to reconstruct
               const I4 KMin,                   ///< [in] shallowest valid layer
               const I4 KMax,                   ///< [in] deepest valid layer
               Real &Slope,                     ///< [out] slope, per Pa
               Real &Curv ///< [out] curvature, per Pa squared
) {
   const Real Thick = PressInterface(ICell, K + 1) - PressInterface(ICell, K);
   if (!(Thick > 0.0_Real)) {
      Slope = 0.0_Real;
      Curv  = 0.0_Real;
      return;
   }

   I4 JLoTop, JHiTop, JLoBot, JHiBot;
   ppmReconWindow(K, KMin, KMax, JLoTop, JHiTop);
   ppmReconWindow(K + 1, KMin, KMax, JLoBot, JHiBot);
   const Real RefValue = Values(ICell, K);

   Real Derivs[4];
   if (JLoTop == JLoBot) {
      const Real PressMid =
          0.5_Real * (PressInterface(ICell, K) + PressInterface(ICell, K + 1));
      ppmReconDerivs(Values, PressInterface, ICell, JLoTop, JHiTop, PressMid,
                     RefValue, Derivs);
      Slope = Derivs[1] + Thick * Thick * Derivs[3] / 24.0_Real;
      Curv  = 0.5_Real * Derivs[2];
   } else {
      // the interface values as departures from this layer's mean
      ppmReconDerivs(Values, PressInterface, ICell, JLoTop, JHiTop,
                     PressInterface(ICell, K), RefValue, Derivs);
      const Real Top = Derivs[0];
      ppmReconDerivs(Values, PressInterface, ICell, JLoBot, JHiBot,
                     PressInterface(ICell, K + 1), RefValue, Derivs);
      const Real Bot = Derivs[0];
      Slope          = (Bot - Top) / Thick;
      Curv           = 3.0_Real * (Top + Bot) / (Thick * Thick);
   }
}

/// The parabolic reconstruction evaluated at a pressure. Its mean over the
/// layer is Value, by construction, and like the linear one it extrapolates
/// where the pressure lies outside the layer.
KOKKOS_INLINE_FUNCTION Real ppmReconEval(
    const Real Value,    ///< [in] prognostic layer mean
    const Real Slope,    ///< [in] slope from ppmReconCoeffs
    const Real Curv,     ///< [in] curvature from ppmReconCoeffs
    const Real PressMid, ///< [in] mid-layer pressure
    const Real Thick,    ///< [in] layer thickness in pressure
    const Real Press     ///< [in] pressure to evaluate at
) {
   const Real Offset = Press - PressMid;
   return Value + Slope * Offset +
          Curv * (Offset * Offset - Thick * Thick / 12.0_Real);
}

/// Layer K of column ICell's reconstruction at a pressure, linear or
/// parabolic.
///
/// This is the one function through which the column scan's three evaluation
/// sites -- the integrand of each edge layer, the anchor's two short integrals,
/// and the anchor's carry-down -- and the unit test's host-side assembly
/// evaluate a column's temperature or salinity, so that the anchor cannot
/// describe different water from the scan. Curv and the interface pressures
/// are read only for the parabolic reconstruction, and the linear one is
/// evaluated exactly as linearReconEval does.
template <class ArrayType>
KOKKOS_INLINE_FUNCTION Real
reconEval(const bool Parabolic,            ///< [in] parabolic, not linear
          const ArrayType &Values,         ///< [in] layer means
          const ArrayType &Slope,          ///< [in] slopes
          const ArrayType &Curv,           ///< [in] curvatures, if parabolic
          const ArrayType &PressMid,       ///< [in] mid-layer pressures
          const ArrayType &PressInterface, ///< [in] interface pressures
          const I4 ICell,                  ///< [in] cell
          const I4 K,                      ///< [in] the cell's layer
          const Real Press                 ///< [in] pressure to evaluate at
) {
   if (Parabolic)
      return ppmReconEval(
          Values(ICell, K), Slope(ICell, K), Curv(ICell, K), PressMid(ICell, K),
          PressInterface(ICell, K + 1) - PressInterface(ICell, K), Press);
   return linearReconEval(Values(ICell, K), Slope(ICell, K), PressMid(ICell, K),
                          Press);
}

} // namespace OMEGA

//===----------------------------------------------------------------------===//
#endif
