#ifndef OMEGA_KPP_MIX_H
#define OMEGA_KPP_MIX_H
//===-- ocn/KPPMix.h - K-Profile Parameterization --------*- C++ -*-===//
//
/// \file
/// \brief K-Profile Parameterization (KPP) boundary layer mixing scheme
///
/// This header defines the KPPMix class for computing ocean boundary layer
/// mixing coefficients using the K-Profile Parameterization scheme.
/// Follows Large et al. (1994) formulation with optional Langmuir circulation
/// enhancement.
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "HorzOperators.h"
#include "KPPConstants.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"
#include "VertCoord.h"
#include <string>

namespace OMEGA {

/// @brief enum sets which matching criterion is used at the BLD base.
/// SimpleShapes does no matching, MatchBoth matches KPP diagnosed
/// diffusivity to the interior mixing coefficient.
/// Also selects the shape used for the non-local flux, which follows the
/// scalar diffusivity profile.
enum class KPPMatchType : I4 {
   SimpleShapes = 0, ///< Unmatched Large et al. (1994) cubic shapes
   MatchBoth    = 1  ///< Match the interior coefficient at the OBL base
};

/// @brief Langmuir circulation enhancement factor applied to buoyancy forcing
class KPPLangmuirFactor {
 public:
   bool UseLangmuirCirculation = true; ///< Apply wave enhancement
   /// Disable Langmuir above this ice fraction for theory wave model
   /// This should be disabled for active wave configurations
   Real IceFracThresholdForLangmuir = KPP::IceFracThresh;

   KOKKOS_FUNCTION void operator()(const Array1DReal &LangmuirFactor, I4 ICell,
                                   const Array1DReal &IceFraction,
                                   const Array1DReal &SurfaceFrictionVelocity,
                                   const Array1DReal &WindSpeed10m) const {

      const Real IceFrac = IceFraction(ICell);
      if (UseLangmuirCirculation && IceFrac < IceFracThresholdForLangmuir) {
         const Real UStar = SurfaceFrictionVelocity(ICell);
         const Real Wind10m =
             (WindSpeed10m.extent(0) > 0) ? WindSpeed10m(ICell) : 0.0_Real;
         LangmuirFactor(ICell) =
             KPP::computeLangmuirEnhancement(Wind10m, UStar, 50.0);
      } else {
         LangmuirFactor(ICell) = 1.0_Real;
      }
   }
};

/// @brief Stage 1 kernel: Compute OBL depth from the bulk Richardson number
/// criterion
///
/// The OBL depth is determined by searching for the depth at which the bulk
/// Richardson number exceeds the critical value.  In KPP the bulk Richardson
/// Number is augmented by a turbulent velocity which is designed to
/// guarantee the entrainment heat flux is -0.2 times the surface buoyancy flux.
///
/// Velocities live on edges (C-grid), so the shear entering Ri is a
/// kite-area weighted average over the edges of each cell.
class KPPOBLDepthSearch {
 public:
   ///< Ri_crit for determining OBL base
   Real CriticalRichardson = KPP::CriticalRi;
   Real SurfaceLayerExtent = KPP::SurfaceLayerExtent; ///< Frac of OBL depth
   /// Apply minimum OBL depth above this ice fraction
   Real IceFracThresholdForMinimumOBL = KPP::IceSuppressThresh;
   /// Min OBL depth under sea ice (m)
   Real MinimumOBLUnderSeaIce = KPP::MinOBLUnderIce;
   /// Continue the search below the crossing to fill the Ri diagnostics
   bool FullRiProfile = false;

   /// Constructor for KPPOBLDepthSearch
   KPPOBLDepthSearch(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void operator()(const Array1DReal &BoundaryLayerDepth,
                                   const Array1DI4 &IndexBoundaryLayerDepth,
                                   const Array2DReal &BulkRichardsonNumber,
                                   const Array2DReal &BulkRichardsonShear,
                                   const Array2DReal &UnresolvedShear,
                                   const Array2DReal &BuoyancyJump, I4 ICell,
                                   const Array2DReal &PotentialDensity,
                                   const Array2DReal &NormalVelocity,
                                   const Array2DReal &TangentialVelocity,
                                   const Array1DReal &SurfaceFrictionVelocity,
                                   const Array1DReal &SurfaceBuoyancyFlux,
                                   const Array2DReal &BruntVaisalaFreqSq,
                                   const Array1DReal &IceFraction,
                                   const Array1DReal &LangmuirFactor) const {

      const Real UStar = SurfaceFrictionVelocity(ICell);
      // The buoyancy flux defined here is the direct flux only, it does not
      // include any enthalpy from mass fluxes.
      const Real BuoyFlux = SurfaceBuoyancyFlux(ICell);

      const I4 KMin     = MinLayerCell(ICell);
      const I4 KMax     = MaxLayerCell(ICell);
      const I4 KIntTop  = Kokkos::min(KMin + 1, NVertLayers);
      const I4 KIntDeep = Kokkos::min(KMax + 1, NVertLayers);

      const Real IceFrac = IceFraction(ICell);

      // KPP depths are measured below the free surface, so geometric
      // heights must be offset by the sea surface height.
      const Real Ssh = SshCell(ICell);

      // Default to the full water column; overwritten if Ri crosses.
      Real OBLDepth         = Ssh - ZInterface(ICell, KIntDeep);
      I4 KCross             = -1;
      const Real RiCritical = CriticalRichardson;
      // Bulk Ri is evaluated at cell centers while the reference average spans
      // the top epsilon*d; this factor corrects for that offset.
      const Real RiScaling   = 1.0_Real - 0.5_Real * SurfaceLayerExtent;
      const Real BuoyFluxEff = BuoyFlux * LangmuirFactor(ICell);

      // Unresolved shear coefficient, Large et al. (1994) Eq. (23):
      // Vt^2 = Cv * sqrt(-beta_T/(c_s*eps)) / (kappa^2 * Ri_crit) * d*N*w_s
      // CSUnres is c_s for the strongly-unstable scalar branch and VtCoef
      // collects the constant prefactor.
      const Real CSUnres = 24.0_Real * Kokkos::sqrt(17.0_Real);
      const Real VtCoef =
          Kokkos::sqrt(0.2_Real / Kokkos::max(KPP::NumericalTolerance,
                                              CSUnres * SurfaceLayerExtent)) /
          (VonKar * VonKar);

      // ----------------------------------------------------------------
      // Edge weights are the MPAS kite areas (0.25*dc*dv) normalized by
      // the cell area.
      // ----------------------------------------------------------------
      const I4 NEdges                 = NEdgesOnCell(ICell);
      const I4 NEdgesEff              = Kokkos::min(NEdges, MaxEdgesBound);
      bool EdgeValid[MaxEdgesBound]   = {};
      Real EdgeWeights[MaxEdgesBound] = {};
      const Real InvAreaCell =
          1.0_Real / Kokkos::max(AreaCell(ICell), KPP::Tiny);
      for (I4 J = 0; J < NEdgesEff; ++J) {
         const I4 IEdge = EdgesOnCell(ICell, J);
         const I4 KEMin = MinLayerEdgeBot(IEdge);
         const I4 KEMax = MaxLayerEdgeTop(IEdge);
         EdgeValid[J]   = (KEMax >= KEMin && KEMin >= 0 && KEMin < NVertLayers);
         if (EdgeValid[J]) {
            EdgeWeights[J] =
                0.25_Real * DcEdge(IEdge) * DvEdge(IEdge) * InvAreaCell;
         }
      }
      if (NEdgesEff > 0) {
         Real SumW = 0.0_Real;
         for (I4 J = 0; J < NEdgesEff; ++J) {
            if (EdgeValid[J]) {
               SumW += EdgeWeights[J];
            }
         }
         if (SumW < KPP::Tiny) {
            // Degenerate kite areas: fall back to equal weighting.
            I4 NEdgesValid = 0;
            for (I4 J = 0; J < NEdgesEff; ++J) {
               if (EdgeValid[J]) {
                  ++NEdgesValid;
               }
            }
            if (NEdgesValid > 0) {
               const Real EqualW = 1.0_Real / static_cast<Real>(NEdgesValid);
               for (I4 J = 0; J < NEdgesEff; ++J) {
                  EdgeWeights[J] = EdgeValid[J] ? EqualW : 0.0_Real;
               }
            }
         } else {
            const Real InvSumW = 1.0_Real / SumW;
            for (I4 J = 0; J < NEdgesEff; ++J) {
               if (EdgeValid[J]) {
                  EdgeWeights[J] *= InvSumW;
               }
            }
         }
      }

      // ----------------------------------------------------------------
      // Bulk Richardson search, Large et al. (1994) Eq. (21):
      //   Ri_b(d) = (B_r - B(d)) d / (|V_r - V(d)|^2 + Vt^2(d))
      // where the reference values B_r, V_r are averaged over the top
      // epsilon*d (d is the candidate OBL depth) of the column. Since
      // epsilon*d grows monotonically with
      // the trial depth, the averaging window only ever extends downward,
      // so the running sums are carried across trial depths rather than
      // rebuilt from the surface at each one.
      // The OBL base is the first d at which Ri_b reaches RiCritical, and
      // the search stops there unless the full profile is wanted for the
      // Ri diagnostics.
      // ----------------------------------------------------------------

      // Cell surface-layer density average
      I4 KSurfaceAvg = KMin;
      const Real ThickTop =
          Kokkos::abs(ZInterface(ICell, KMin + 1) - ZInterface(ICell, KMin));
      Real SumThickness = Kokkos::max(ThickTop, KPP::NumericalTolerance);
      Real SumRho       = PotentialDensity(ICell, KMin) * SumThickness;

      // Per-edge surface-layer velocity averages
      // These are updated as the surface layer depth grows with increasing
      // trial depth
      I4 KSurfE[MaxEdgesBound]      = {};
      Real SumThickE[MaxEdgesBound] = {};
      Real SumUnE[MaxEdgesBound]    = {};
      Real SumVtE[MaxEdgesBound]    = {};

      for (I4 J = 0; J < NEdgesEff; ++J) {
         if (!EdgeValid[J]) {
            continue;
         }
         const I4 IEdge  = EdgesOnCell(ICell, J);
         const I4 JCell  = CellsOnCell(ICell, J);
         const I4 KEMin  = MinLayerEdgeBot(IEdge);
         KSurfE[J]       = KEMin;
         const I4 KIntE0 = Kokkos::min(KEMin + 1, NVertLayers);
         const Real ZEdgeTop =
             0.5_Real * (ZInterface(ICell, KEMin) + ZInterface(JCell, KEMin));
         const Real ZEdgeBot =
             0.5_Real * (ZInterface(ICell, KIntE0) + ZInterface(JCell, KIntE0));
         const Real Thick0 = Kokkos::abs(ZEdgeBot - ZEdgeTop);
         SumThickE[J]      = Kokkos::max(Thick0, KPP::NumericalTolerance);
         const I4 KE0      = Kokkos::min(KEMin, NVertLayers - 1);
         SumUnE[J]         = NormalVelocity(IEdge, KE0) * SumThickE[J];
         SumVtE[J]         = TangentialVelocity(IEdge, KE0) * SumThickE[J];
      }

      for (I4 K = KMin; K <= KMax; ++K) {
         const I4 KCell     = Kokkos::min(K, NVertLayers - 1);
         const I4 KInt      = Kokkos::min(K + 1, NVertLayers);
         const Real ZDepth  = Ssh - ZInterface(ICell, KInt);
         const Real ZCenter = Ssh - ZMid(ICell, KCell);
         if (ZDepth < KPP::NumericalTolerance)
            continue;

         const Real SurfLayerDepth = SurfaceLayerExtent * ZDepth;

         // Advance cell surface average for density
         while (KSurfaceAvg < K &&
                (Ssh - ZInterface(ICell, KSurfaceAvg + 1)) < SurfLayerDepth) {
            ++KSurfaceAvg;
            const I4 KSA      = Kokkos::min(KSurfaceAvg, NVertLayers - 1);
            const Real DZ     = Kokkos::abs(ZInterface(ICell, KSurfaceAvg + 1) -
                                            ZInterface(ICell, KSurfaceAvg));
            const Real ThickK = Kokkos::max(DZ, KPP::NumericalTolerance);
            SumThickness += ThickK;
            SumRho += PotentialDensity(ICell, KSA) * ThickK;
         }

         // Advance per-edge surface averages for velocity
         for (I4 J = 0; J < NEdgesEff; ++J) {
            if (!EdgeValid[J]) {
               continue;
            }
            const I4 IEdge     = EdgesOnCell(ICell, J);
            const I4 JCell     = CellsOnCell(ICell, J);
            const I4 KEMax     = MaxLayerEdgeTop(IEdge);
            const Real SshEdge = 0.5_Real * (SshCell(ICell) + SshCell(JCell));
            while (KSurfE[J] < K &&
                   (SshEdge - 0.5_Real * (ZInterface(ICell, KSurfE[J] + 1) +
                                          ZInterface(JCell, KSurfE[J] + 1))) <
                       SurfLayerDepth) {
               ++KSurfE[J];
               const I4 KE = Kokkos::min(
                   Kokkos::max(KSurfE[J], MinLayerEdgeBot(IEdge)), KEMax);
               const Real ZEdgeTop = 0.5_Real * (ZInterface(ICell, KSurfE[J]) +
                                                 ZInterface(JCell, KSurfE[J]));
               const Real ZEdgeBot =
                   0.5_Real * (ZInterface(ICell, KSurfE[J] + 1) +
                               ZInterface(JCell, KSurfE[J] + 1));
               const Real DZ     = Kokkos::abs(ZEdgeBot - ZEdgeTop);
               const Real ThickK = Kokkos::max(DZ, KPP::NumericalTolerance);
               SumThickE[J] += ThickK;
               SumUnE[J] += NormalVelocity(IEdge, KE) * ThickK;
               SumVtE[J] += TangentialVelocity(IEdge, KE) * ThickK;
            }
         }

         // Compute the average density in the surface layer
         const Real InvSumThickness =
             1.0_Real / Kokkos::max(SumThickness, KPP::NumericalTolerance);
         const Real RhoAvgSurf = SumRho * InvSumThickness;

         // Buoyancy jump B_r - B(d), positive for stable stratification
         const Real RhoK     = PotentialDensity(ICell, KCell);
         const Real DeltaRho = RhoK - RhoAvgSurf;
         const Real DeltaB   = Gravity * DeltaRho / RhoSw;
         // The following field is for diagnostic purposes only
         // and does not affect the computation
         BuoyancyJump(ICell, KInt) = DeltaB;

         // Resolved shear |V_r - V(d)|^2, averaged over the cell edges
         Real DeltaVSq = 0.0_Real;
         if (NEdges > 0) {
            for (I4 J = 0; J < NEdgesEff; ++J) {
               if (!EdgeValid[J]) {
                  continue;
               }
               const I4 IEdge = EdgesOnCell(ICell, J);
               const I4 KEMin = MinLayerEdgeBot(IEdge);
               const I4 KEMax = MaxLayerEdgeTop(IEdge);
               const I4 KE    = Kokkos::min(Kokkos::max(K, KEMin), KEMax);
               const Real InvThickE =
                   1.0_Real /
                   Kokkos::max(SumThickE[J], KPP::NumericalTolerance);
               const Real UnAvg = SumUnE[J] * InvThickE;
               const Real VtAvg = SumVtE[J] * InvThickE;
               const Real UnK   = NormalVelocity(IEdge, KE);
               const Real VtK   = TangentialVelocity(IEdge, KE);
               const Real DUn   = UnK - UnAvg;
               const Real DVt   = VtK - VtAvg;
               DeltaVSq += EdgeWeights[J] * (DUn * DUn + DVt * DVt);
            }
         }
         BulkRichardsonShear(ICell, KInt) = DeltaVSq;

         const Real SigmaLoc =
             Kokkos::fmin(1.0_Real, Kokkos::fmax(0.0_Real, SurfaceLayerExtent));

         // Turbulent scalar velocity scale w_s at the surface-layer depth
         Real WTurb = 0.0_Real;
         if (UStar > KPP::NumericalTolerance) {
            const Real U3   = UStar * UStar * UStar;
            const Real Zeta = SigmaLoc * ZDepth * VonKar * BuoyFluxEff /
                              Kokkos::max(U3, KPP::Tiny);
            const Real PhiInvS = KPP::kppPhiInvScalar(Zeta);
            WTurb = VonKar * UStar * Kokkos::max(PhiInvS, 0.0_Real);
         } else if (BuoyFluxEff < 0.0_Real) {
            // Free convection limit: u* drops out and w_s ~ (c_s d B_0)^1/3
            const Real CS  = KPP::CMoS;
            const Real WS3 = -CS * SigmaLoc * ZDepth * VonKar * BuoyFluxEff;
            WTurb          = VonKar * Kokkos::pow(Kokkos::max(WS3, 0.0_Real),
                                                  1.0_Real / 3.0_Real);
         }

         // Unresolved turbulent shear Vt^2 (m^2/s^2), Large et al. Eq.
         // (23). Cv ramps from 2.1 to 1.7 as stratification strengthens.
         const Real NCntr = Kokkos::sqrt(
             Kokkos::max(0.0_Real, BruntVaisalaFreqSq(ICell, KInt)));
         const Real Cv =
             (NCntr < 0.002_Real) ? (2.1_Real - 200.0_Real * NCntr) : 1.7_Real;
         const Real Vt2 =
             Kokkos::max(KPP::MinUnresolvedShearSq,
                         Cv * VtCoef * ZCenter * NCntr * WTurb /
                             Kokkos::max(RiCritical, KPP::NumericalTolerance));
         UnresolvedShear(ICell, KInt) = Vt2;

         const Real VelScaleSq = DeltaVSq + Vt2;

         const Real RiBulk = RiScaling * DeltaB * ZCenter /
                             Kokkos::max(VelScaleSq, KPP::NumericalTolerance);
         BulkRichardsonNumber(ICell, KInt) = RiBulk;

         if (KCross < 0 && RiBulk > RiCritical) {
            KCross = K;
            // Levels below the crossing only feed the Ri diagnostics
            if (!FullRiProfile)
               break;
         }
      }

      if (KCross >= KMin) {
         if (KCross > KMin) {
            // Ri values are defined at cell centers, so interpolate on
            // center depths to keep the abscissa consistent.
            // Define quantities to conduct the interpolation of Ri
            // between the two model levels to determine the final
            // boundary layer depth.
            const I4 KAbove    = Kokkos::max(KMin, KCross - 1);
            const I4 KBelow    = Kokkos::min(KCross, NVertLayers - 1);
            const I4 KAboveRi  = Kokkos::min(KAbove + 1, NVertLayers);
            const I4 KBelowRi  = Kokkos::min(KBelow + 1, NVertLayers);
            const Real ZAbove  = Ssh - ZMid(ICell, KAbove);
            const Real ZBelow  = Ssh - ZMid(ICell, KBelow);
            const Real RiAbove = BulkRichardsonNumber(ICell, KAboveRi);
            const Real RiBelow = BulkRichardsonNumber(ICell, KBelowRi);

            const Real H = ZBelow - ZAbove;
            if (H > KPP::NumericalTolerance) {
               // CVMix-style QUAD interpolation for OBL crossing:
               // - first interior crossing uses zero slope at top point
               // - deeper crossings use upstream slope.
               // Quadratic interpolation is recommended by Danabosglu et al
               // (2006)
               // https://journals.ametsoc.org/view/journals/clim/19/11/jcli3739.1.pdf
               Real SlopeAbove = 0.0_Real;
               if (KCross > KMin + 1) {
                  const I4 KPrev    = Kokkos::max(KMin, KAbove - 1);
                  const I4 KPrevRi  = Kokkos::min(KPrev + 1, NVertLayers);
                  const Real ZPrev  = Ssh - ZMid(ICell, KPrev);
                  const Real RiPrev = BulkRichardsonNumber(ICell, KPrevRi);
                  const Real DZPrev = ZAbove - ZPrev;
                  if (Kokkos::abs(DZPrev) > KPP::NumericalTolerance) {
                     SlopeAbove = (RiAbove - RiPrev) / DZPrev;
                  }
               }

               // In local coordinate T = z - ZAbove:
               // Ri(T) = QuadA T^2 + SlopeAbove T + RiAbove, with QuadA
               // fixed by requiring Ri(H) = RiBelow. The OBL base is the
               // root of Ri(T) = RiCritical.
               const Real QuadA =
                   (RiBelow - RiAbove - SlopeAbove * H) / (H * H);
               const Real QuadC = RiAbove - RiCritical;

               Real TCross = H;
               if (Kokkos::abs(QuadA) < 1.0e-14_Real) {
                  // Degenerate quadratic -> linear fallback.
                  const Real DRi = RiBelow - RiAbove;
                  if (Kokkos::abs(DRi) > KPP::NumericalTolerance) {
                     const Real Frac = Kokkos::fmax(
                         0.0_Real,
                         Kokkos::fmin(1.0_Real, (RiCritical - RiAbove) / DRi));
                     TCross = Frac * H;
                  }
               } else {
                  const Real Disc =
                      SlopeAbove * SlopeAbove - 4.0_Real * QuadA * QuadC;
                  if (Disc >= 0.0_Real) {
                     const Real SqrtDisc = Kokkos::sqrt(Disc);
                     const Real T1 =
                         (-SlopeAbove + SqrtDisc) / (2.0_Real * QuadA);
                     const Real T2 =
                         (-SlopeAbove - SqrtDisc) / (2.0_Real * QuadA);

                     const bool T1Ok = (T1 >= 0.0_Real && T1 <= H);
                     const bool T2Ok = (T2 >= 0.0_Real && T2 <= H);
                     if (T1Ok && T2Ok) {
                        // Both roots lie in the interval; prefer the one
                        // nearest mid-interval, as CVMix does.
                        const Real Mid = 0.5_Real * H;
                        TCross =
                            (Kokkos::abs(T1 - Mid) <= Kokkos::abs(T2 - Mid))
                                ? T1
                                : T2;
                     } else if (T1Ok) {
                        TCross = T1;
                     } else if (T2Ok) {
                        TCross = T2;
                     } else {
                        TCross = H;
                     }
                  }
               }

               TCross   = Kokkos::fmax(0.0_Real, Kokkos::fmin(H, TCross));
               OBLDepth = ZAbove + TCross;
            } else {
               OBLDepth = ZBelow;
            }
         } else {
            // Match center-based OBL convention when crossing occurs in
            // the top interval.
            OBLDepth = Ssh - ZMid(ICell, KMin);
         }
      } else {
         OBLDepth = Ssh - ZInterface(ICell, KIntDeep);
      }

      const Real TopLayerThickness =
          Kokkos::abs(ZInterface(ICell, KIntTop) - ZInterface(ICell, KMin));
      const Real MinOBLDepth = 0.5_Real * TopLayerThickness;
      const Real MaxOBLDepth = Ssh - ZMid(ICell, KMax);
      // Impose chosen limits on the depth of the OBL
      OBLDepth = KPP::kppClampOBLDepth(OBLDepth, MinOBLDepth, MaxOBLDepth,
                                       IceFrac > IceFracThresholdForMinimumOBL,
                                       MinimumOBLUnderSeaIce);

      const I4 KFinal =
          KPP::kppOBLIndex(ZInterface, ICell, KMin, KMax, Ssh, OBLDepth);

      BoundaryLayerDepth(ICell)      = OBLDepth;
      IndexBoundaryLayerDepth(ICell) = KFinal;
   }

 private:
   /// Compile-time bound for the per-thread edge scratch arrays
   static constexpr I4 MaxEdgesBound = HorzMesh::MaxEdgesBound;

   I4 NVertLayers;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
   Array2DReal ZInterface;
   Array2DReal ZMid;
   Array1DReal SshCell;
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnCell;
   Array1DReal AreaCell;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
};

/// @brief Area-weighted smoothing of the boundary layer depth over each cell
/// and its neighbors. This is a laplacian smoothing operation. This suppresses
/// the grid-scale noise that the discrete Ri crossing search and no time
/// history in the OBL depth calculation can introduce.
class KPPBLDSmooth {
 public:
   /// Constructor for KPPBLDSmooth
   KPPBLDSmooth(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   operator()(const Array1DReal &BoundaryLayerDepthSmooth, I4 ICell,
              const Array1DReal &BoundaryLayerDepth) const {

      const I4 KMin = MinLayerCell(ICell);
      if (KMin < 0 || KMin >= NVertLayers) {
         BoundaryLayerDepthSmooth(ICell) = BoundaryLayerDepth(ICell);
         return;
      }

      const I4 NEdges = NEdgesOnCell(ICell);
      Real AreaSum    = 0.0_Real;
      Real BLDSum     = 0.0_Real;
      I4 EdgeCount    = 0;

      for (I4 J = 0; J < NEdges; ++J) {
         const I4 INeighbor = CellsOnCell(ICell, J);
         if (INeighbor == NCellsAll) {
            continue;
         }

         const I4 KMinNbr = MinLayerCell(INeighbor);
         if (KMinNbr < 0 || KMinNbr >= NVertLayers) {
            continue;
         }

         const Real NbrArea = AreaCell(INeighbor);
         BLDSum += 2.0_Real * NbrArea * BoundaryLayerDepth(INeighbor);
         AreaSum += 2.0_Real * NbrArea;
         ++EdgeCount;
      }

      if (EdgeCount > 0) {
         const Real SelfArea = AreaCell(ICell);
         BLDSum += BoundaryLayerDepth(ICell) * static_cast<Real>(EdgeCount) *
                   SelfArea;
         AreaSum += static_cast<Real>(EdgeCount) * SelfArea;
      }

      if (AreaSum > 0.0_Real) {
         BoundaryLayerDepthSmooth(ICell) = BLDSum / AreaSum;
      } else {
         BoundaryLayerDepthSmooth(ICell) = BoundaryLayerDepth(ICell);
      }
   }

 private:
   I4 NVertLayers;
   I4 NCellsAll;
   Array1DI4 MinLayerCell;
   Array1DI4 NEdgesOnCell;
   Array2DI4 CellsOnCell;
   Array1DReal AreaCell;
};

/// @brief Re-clamp and re-index the boundary layer depth after smoothing
class KPPBLDCommit {
 public:
   /// Constructor for KPPBLDCommit
   KPPBLDCommit(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   operator()(const Array1DReal &BoundaryLayerDepth,
              const Array1DI4 &IndexBoundaryLayerDepth, I4 ICell,
              const Array1DReal &BoundaryLayerDepthSmooth) const {

      const I4 KMin = MinLayerCell(ICell);
      const I4 KMax = MaxLayerCell(ICell);
      if (KMin < 0 || KMax < KMin || KMin >= NVertLayers) {
         return;
      }

      const Real Ssh = SshCell(ICell);

      const I4 KIntTop = Kokkos::min(KMin + 1, NVertLayers);
      const Real TopLayerThickness =
          Kokkos::abs(ZInterface(ICell, KIntTop) - ZInterface(ICell, KMin));
      const Real MinOBLDepth = 0.5_Real * TopLayerThickness;
      const Real MaxOBLDepth = Ssh - ZMid(ICell, KMax);

      // The sea-ice minimum is deliberately not reapplied here; it was
      // already enforced before smoothing.
      const Real OBLDepth =
          KPP::kppClampOBLDepth(BoundaryLayerDepthSmooth(ICell), MinOBLDepth,
                                MaxOBLDepth, false, 0.0_Real);

      const I4 KFinal =
          KPP::kppOBLIndex(ZInterface, ICell, KMin, KMax, Ssh, OBLDepth);

      BoundaryLayerDepth(ICell)      = OBLDepth;
      IndexBoundaryLayerDepth(ICell) = KFinal;
   }

 private:
   I4 NVertLayers;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array2DReal ZInterface;
   Array2DReal ZMid;
   Array1DReal SshCell;
};

/// @brief Initialize the coefficient arrays with zero KPP contribution, or
/// with precomputed interior mixing for matched-coefficient construction.
class KPPCoeffsInit {
 public:
   bool UseInteriorMix = false; ///< Interior coefficients were supplied

   KOKKOS_FUNCTION void operator()(const Array2DReal &VertDiff,
                                   const Array2DReal &VertVisc,
                                   const Array2DReal &VertNonLocalFlux,
                                   const Array2DReal &TurbulentVelocityScale,
                                   I4 ICell, I4 K,
                                   const Array2DReal &InteriorVertDiff,
                                   const Array2DReal &InteriorVertVisc) const {

      VertDiff(ICell, K) =
          UseInteriorMix ? InteriorVertDiff(ICell, K) : 0.0_Real;
      VertVisc(ICell, K) =
          UseInteriorMix ? InteriorVertVisc(ICell, K) : 0.0_Real;
      VertNonLocalFlux(ICell, K)       = 0.0_Real;
      TurbulentVelocityScale(ICell, K) = 0.0_Real;
   }
};

/// @brief Stage 2 kernel: KPP profile-based diffusivity, viscosity and
/// non-local flux, with optional enhanced mixing at the OBL base.
class KPPMixingCoeffs {
 public:
   Real SurfaceLayerExtent   = KPP::SurfaceLayerExtent; ///< Frac of OBL depth
   bool UseEnhancedDiffusion = true;  ///< Apply enhanced mixing at OBL base
   bool UseInteriorMix       = false; ///< Interior coefficients were supplied
   bool UseMatchedShapes     = false; ///< Match interior value at the OBL base

   /// Constructor for KPPMixingCoeffs
   KPPMixingCoeffs(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void operator()(const Array2DReal &VertDiff,
                                   const Array2DReal &VertVisc,
                                   const Array2DReal &VertNonLocalFlux,
                                   const Array2DReal &TurbulentVelocityScale,
                                   I4 ICell,
                                   const Array1DReal &BoundaryLayerDepth,
                                   const Array1DI4 &IndexBoundaryLayerDepth,
                                   const Array1DReal &SurfaceFrictionVelocity,
                                   const Array1DReal &SurfaceBuoyancyFlux,
                                   const Array2DReal &InteriorVertDiff,
                                   const Array2DReal &InteriorVertVisc) const {

      const Real HOBL = BoundaryLayerDepth(ICell);

      const I4 KMin = MinLayerCell(ICell);
      const I4 KMax = MaxLayerCell(ICell);
      if (KMin < 0 || KMin >= NVertLayers || KMax < KMin) {
         return;
      }
      const I4 KMatch =
          Kokkos::min(KMax + 1, IndexBoundaryLayerDepth(ICell) + 1);

      // KPP depths are measured below the free surface, so geometric
      // heights must be offset by the sea surface height.
      const Real Ssh = SshCell(ICell);

      const Real UStar      = SurfaceFrictionVelocity(ICell);
      const Real BuoyFlux   = SurfaceBuoyancyFlux(ICell);
      const Real NonLocalCs = KPP::kppNonLocalCs(VonKar, SurfaceLayerExtent);

      for (I4 K = KMin; K <= KMax + 1; ++K) {
         const I4 KIface   = Kokkos::min(Kokkos::max(K, 0), NVertLayers);
         const Real ZDepth = Ssh - ZInterface(ICell, KIface);

         // Check if within OBL using depth below the free surface.
         if (ZDepth <= HOBL && HOBL > 0.0_Real) {
            // Normalized depth in Omega sign convention: sigma in [-1,0].
            Real Sigma = -ZDepth / HOBL;
            Sigma      = Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, Sigma));

            // CVMix-style turbulent scales: w = kappa*u*/phi in general,
            // with explicit free-convection limits when u*=0. The scales
            // are frozen at the surface-layer depth below the surface
            // layer, so SigmaLoc is capped at SurfaceLayerExtent.
            const Real SigmaCoord = -Sigma; // [0,1]
            const Real SigmaLoc   = Kokkos::fmin(
                SurfaceLayerExtent, Kokkos::fmax(0.0_Real, SigmaCoord));

            Real WMTurb = 0.0_Real;
            Real WSTurb = 0.0_Real;
            KPP::kppTurbScales(UStar, BuoyFlux, HOBL, SigmaLoc, VonKar, WMTurb,
                               WSTurb);

            // For MatchBoth, the shape value the KPP profile must reach at
            // the OBL base so that it joins the interior coefficient there.
            const Real MatchViscShape =
                UseMatchedShapes
                    ? KPP::kppMatchShape(InteriorVertVisc(ICell, KMatch), HOBL,
                                         WMTurb)
                    : 0.0_Real;
            const Real MatchDiffShape =
                UseMatchedShapes
                    ? KPP::kppMatchShape(InteriorVertDiff(ICell, KMatch), HOBL,
                                         WSTurb)
                    : 0.0_Real;

            // Momentum mixing contribution.
            const Real ShapeM =
                UseMatchedShapes ? KPP::kppShapeMatched(Sigma, MatchViscShape)
                                 : KPP::kppShapeMomentum(Sigma);
            VertVisc(ICell, K) = HOBL * WMTurb * ShapeM;

            // Tracer mixing contribution.
            const Real ShapeS =
                UseMatchedShapes ? KPP::kppShapeMatched(Sigma, MatchDiffShape)
                                 : KPP::kppShapeScalar(Sigma);
            VertDiff(ICell, K)               = HOBL * WSTurb * ShapeS;
            TurbulentVelocityScale(ICell, K) = WSTurb;

            // Non-local flux: C_s * G(sigma).
            // C_s = C* * kappa * (c_s * kappa * epsilon)^(1/3)
            // per Large et al. (1994) eq. 20 (~6.33 with default constants).
            // The non-local shape is always the unmatched scalar shape,
            // independent of MatchTechnique, so that gamma vanishes at the
            // OBL base. The matched shape is non-zero there by construction,
            // which would leave a non-local flux at the base that jumps to
            // zero just below it. CVMix likewise keeps matching (a
            // diffusivity choice) separate from the non-local shape.
            const Real NonLocalShape = KPP::kppShapeScalar(Sigma);

            // Match CVMix behavior: apply non-local term only when
            // surface buoyancy forcing is unstable/neutral.
            if (BuoyFlux <= 0.0_Real) {
               VertNonLocalFlux(ICell, K) = NonLocalCs * NonLocalShape;
            } else {
               VertNonLocalFlux(ICell, K) = 0.0_Real;
            }

         } else {
            // Below OBL: preserve interior values for MatchBoth, otherwise
            // no KPP contribution.
            VertDiff(ICell, K) =
                UseInteriorMix ? InteriorVertDiff(ICell, K) : 0.0_Real;
            VertVisc(ICell, K) =
                UseInteriorMix ? InteriorVertVisc(ICell, K) : 0.0_Real;
            VertNonLocalFlux(ICell, K)       = 0.0_Real;
            TurbulentVelocityScale(ICell, K) = 0.0_Real;
         }
      }

      // Optional enhanced diffusion/viscosity treatment at OBL base.
      // Match CVMix Appendix D weighting at the interface nearest HOBL:
      // the OBL base rarely lands on an interface, so the coefficient at
      // the neighboring interface KTarget is replaced by a quadratic blend
      // of the KPP value extrapolated from KKtup and the value already
      // there, weighted by where HOBL falls between the two cell centers.
      if (UseEnhancedDiffusion && HOBL > 0.0_Real) {
         const I4 KOBL = Kokkos::max(
             KMin, Kokkos::min(IndexBoundaryLayerDepth(ICell), KMax));
         const Real ZMidOBL = Ssh - ZMid(ICell, KOBL);

         const bool TargetOutsideOBL = HOBL >= ZMidOBL;
         const I4 KKtup = TargetOutsideOBL ? KOBL : Kokkos::max(KMin, KOBL - 1);
         const I4 KTarget = TargetOutsideOBL ? Kokkos::min(KOBL + 1, KMax + 1)
                                             : Kokkos::max(KMin + 1, KOBL);

         const Real ZKtup = Ssh - ZMid(ICell, KKtup);
         const Real ZNext = (KKtup < KMax)
                                ? (Ssh - ZMid(ICell, KKtup + 1))
                                : (Ssh - ZInterface(ICell, KKtup + 1));
         const Real Delta = Kokkos::fmax(
             0.0_Real,
             Kokkos::fmin(1.0_Real, (HOBL - ZKtup) /
                                        Kokkos::max(ZNext - ZKtup,
                                                    KPP::NumericalTolerance)));
         const Real OneMinusDelta = 1.0_Real - Delta;

         Real SigmaKtup = -ZKtup / HOBL;
         SigmaKtup = Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, SigmaKtup));
         const Real SigmaCoord = -SigmaKtup;
         const Real SigmaLoc   = Kokkos::fmin(SurfaceLayerExtent,
                                              Kokkos::fmax(0.0_Real, SigmaCoord));

         Real WMKtup = 0.0_Real;
         Real WSKtup = 0.0_Real;
         KPP::kppTurbScales(UStar, BuoyFlux, HOBL, SigmaLoc, VonKar, WMKtup,
                            WSKtup);

         const Real MatchViscShape =
             UseMatchedShapes
                 ? KPP::kppMatchShape(InteriorVertVisc(ICell, KMatch), HOBL,
                                      WMKtup)
                 : 0.0_Real;
         const Real MatchDiffShape =
             UseMatchedShapes
                 ? KPP::kppMatchShape(InteriorVertDiff(ICell, KMatch), HOBL,
                                      WSKtup)
                 : 0.0_Real;

         const Real ViscKtup =
             HOBL * WMKtup *
             (UseMatchedShapes ? KPP::kppShapeMatched(SigmaKtup, MatchViscShape)
                               : KPP::kppShapeMomentum(SigmaKtup));
         const Real DiffKtup =
             HOBL * WSKtup *
             (UseMatchedShapes ? KPP::kppShapeMatched(SigmaKtup, MatchDiffShape)
                               : KPP::kppShapeScalar(SigmaKtup));

         const Real ViscProfile = VertVisc(ICell, KTarget);
         const Real DiffProfile = VertDiff(ICell, KTarget);

         const Real EnhVisc = OneMinusDelta * OneMinusDelta * ViscKtup +
                              Delta * Delta * ViscProfile;
         const Real EnhDiff = OneMinusDelta * OneMinusDelta * DiffKtup +
                              Delta * Delta * DiffProfile;

         const Real OldVisc =
             UseInteriorMix ? InteriorVertVisc(ICell, KTarget) : 0.0_Real;
         const Real OldDiff =
             UseInteriorMix ? InteriorVertDiff(ICell, KTarget) : 0.0_Real;
         const Real NewVisc = OneMinusDelta * OldVisc + Delta * EnhVisc;
         const Real NewDiff = OneMinusDelta * OldDiff + Delta * EnhDiff;

         VertVisc(ICell, KTarget) = NewVisc;
         VertDiff(ICell, KTarget) = NewDiff;

         // Keep the non-local term consistent with the rescaled diffusivity
         if (!TargetOutsideOBL && DiffProfile != 0.0_Real) {
            VertNonLocalFlux(ICell, KTarget) =
                VertNonLocalFlux(ICell, KTarget) * NewDiff / DiffProfile;
         } else if (!TargetOutsideOBL) {
            VertNonLocalFlux(ICell, KTarget) = 0.0_Real;
         }
      }
   }

 private:
   I4 NVertLayers;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array2DReal ZInterface;
   Array2DReal ZMid;
   Array1DReal SshCell;
};

/// @brief KPP Boundary Layer Mixing Scheme
///
/// Implements the K-Profile Parameterization following Large et al. (1994)
/// with optional Langmuir circulation enhancement. Computes vertical
/// diffusivity, viscosity, and non-local flux coefficients for the OBL.
///
/// Two-stage computation:
/// 1. Stage 1: Compute OBL depth from bulk Richardson criterion
/// 2. Stage 2: Compute mixing coefficients within and below OBL
///
class KPPMix {

 public:
   /// @brief Singleton instance management
   static KPPMix *getInstance();
   static void init();
   static void destroyInstance();

   /// @brief Main computation routine
   /// Calls Stage 1 and Stage 2 computation in sequence
   ///
   /// Input arrays should be pre-populated with current state.
   /// NormalVelocity and TangentialVelocity are edge-based quantities
   /// (C-grid convention): dimensions [NEdges × NVertLayers].
   /// Output arrays are computed in-place.
   void computeKPPMix(
       const Array2DReal
           &PotentialDensity, ///< Density (kg/m³) [NCells×NLevels]
       const Array2DReal &NormalVelocity,     ///< Normal vel on edges (m/s)
       const Array2DReal &TangentialVelocity, ///< Tangential vel on edges (m/s)
       const Array1DReal &SurfaceFrictionVelocity, ///< u* (m/s)
       const Array1DReal &SurfaceBuoyancyFlux,     ///< B_0 (m²/s³)
       const Array2DReal &BruntVaisalaFreqSq,      ///< N² (s⁻²)
       const Array1DReal &IceFraction,             ///< Sea ice cover (0-1)
       const Array1DReal &WindSpeed10m =
           Array1DReal() ///< Wind for Langmuir (m/s)
   );

   // =======================================================================
   // Output Fields
   // =======================================================================

   /// @brief Vertical diffusivity at layer interfaces (m²/s)
   /// Size: [nCells][nLevels+1]
   Array2DReal VertDiff;

   /// @brief Vertical viscosity at layer interfaces (m²/s)
   /// Size: [nCells][nLevels+1]
   Array2DReal VertVisc;

   /// @brief Boundary layer depth (m)
   /// Size: [nCells]
   Array1DReal BoundaryLayerDepth;

   /// @brief OBL depth as layer index
   /// Size: [nCells]
   Array1DI4 IndexBoundaryLayerDepth;

   /// @brief Non-local flux coefficient profile G(σ) (dimensionless)
   /// Size: [nCells][nLevels+1]
   /// Applied to surface tracer fluxes to compute vertical transport
   Array2DReal VertNonLocalFlux;

   /// @brief Bulk Richardson number profile used in OBL search (dimensionless)
   /// Size: [nCells][nLevels+1]
   Array2DReal BulkRichardsonNumber;

   /// @brief Shear contribution to bulk Richardson denominator (m^2/s^2)
   /// Size: [nCells][nLevels+1]
   Array2DReal BulkRichardsonShear;

   /// @brief Unresolved shear contribution Vt^2 (m^2/s^2)
   /// Size: [nCells][nLevels+1]
   Array2DReal UnresolvedShear;

   /// @brief Buoyancy jump (density anomaly converted to buoyancy) (m/s²)
   /// Size: [nCells][nLevels+1]
   /// Captures delta_b = g * delta_rho / rho_sw at each layer during OBL search
   Array2DReal BuoyancyJump;

   /// @brief Turbulent velocity scale profile (m/s), tracer branch
   /// Size: [nCells][nLevels+1]
   Array2DReal TurbulentVelocityScale;

   /// @brief Potential density used by KPP OBL search (kg/m^3)
   /// Size: [nCells][nLevels]
   Array2DReal PotentialDensity;

   /// @brief Surface friction velocity u* (m/s)
   /// Size: [nCells]
   Array1DReal SurfaceFrictionVelocity;

   /// @brief Surface buoyancy flux B_0 (m²/s³)
   /// Size: [nCells]
   Array1DReal SurfaceBuoyancyFlux;

   // =======================================================================
   // Configuration Parameters
   // =======================================================================

   bool Enabled = true; ///< Enable/disable KPP mixing

   // Defaults below may be overridden from the Config file; where a value also
   // appears in KPPConstants.h, that is the authoritative default.
   Real CriticalRichardson = KPP::CriticalRi;         ///< Ri_crit for OBL base
   Real SurfaceLayerExtent = KPP::SurfaceLayerExtent; ///< Frac of OBL depth

   bool UseLangmuirCirculation = true;  ///< Apply wave enhancement
   bool DebugDiagnostics       = false; ///< Print per-step KPP diagnostics

   // Ice/Langmuir controls (kept configurable to match reference semantics)
   /// Disable Langmuir above this ice fraction
   Real IceFractionThresholdForLangmuir = KPP::IceFracThresh;
   /// Apply minimum OBL depth above this ice fraction
   Real IceFractionThresholdForMinimumOBL = KPP::IceSuppressThresh;
   /// Min OBL depth under sea ice (m)
   Real MinimumOBLUnderSeaIce = KPP::MinOBLUnderIce;

   Real BackgroundVisc = 1.0e-4; ///< Background viscosity below OBL (m²/s)
   Real BackgroundDiff = 1.0e-5; ///< Background diffusivity below OBL (m²/s)

   // KPP matching/profile controls (CVMix-style semantics)
   KPPMatchType MatchTechnique = KPPMatchType::SimpleShapes;
   std::string InterpType2Str  = "LMD94"; ///< Linear, Quadratic, Cubic, LMD94
   bool UseEnhancedDiffusion   = true;    ///< Apply enhanced mixing at OBL base
   bool UseBLDSmoothing = true; ///< Apply MPAS-style BLD horizontal smoothing

   // Field names for I/O
   std::string BuoyancyJumpFldName;
   std::string VertDiffFldName;
   std::string VertViscFldName;
   std::string OBLDepthFldName;
   std::string OBLDepthIndexFldName;
   std::string NonLocalFluxFldName;
   std::string BulkRichardsonFldName;
   std::string BulkRichardsonShearFldName;
   std::string UnresolvedShearFldName;
   std::string TurbulentVelScaleFldName;
   std::string PotentialDensityFldName;
   std::string SurfFricVelFldName;
   std::string SurfBuoyFluxFldName;
   std::string Name;

 private:
   /// @brief Private constructor for singleton pattern
   KPPMix(const std::string &InName, const HorzMesh *InMesh,
          const VertCoord *InVCoord);

   /// @brief Private destructor
   ~KPPMix();

   /// @brief Static singleton instance
   static KPPMix *Instance;

   /// @brief Mesh and coordinate references
   const HorzMesh *Mesh;
   const VertCoord *VCoord;

 public:
   /// @brief Stage 1: Compute OBL depth using edge-based velocity shear
   void computeOBLDepth(const Array2DReal &PotentialDensity,
                        const Array2DReal &NormalVelocity,
                        const Array2DReal &TangentialVelocity,
                        const Array1DReal &SurfaceFrictionVelocity,
                        const Array1DReal &SurfaceBuoyancyFlux,
                        const Array2DReal &BruntVaisalaFreqSq,
                        const Array1DReal &IceFraction,
                        const Array1DReal &WindSpeed10m);

   /// @brief Stage 2: Compute KPP mixing contribution or matched coefficients
   void computeMixingCoefficients(
       const Array2DReal &PotentialDensity,
       const Array1DReal &SurfaceFrictionVelocity,
       const Array1DReal &SurfaceBuoyancyFlux,
       const Array2DReal &InteriorVertDiff = Array2DReal(),
       const Array2DReal &InteriorVertVisc = Array2DReal());

 private:
   /// @brief Print targeted diagnostics for KPP troubleshooting
   void logDiagnostics(const Array2DReal &PotentialDensity,
                       const Array2DReal &NormalVelocity,
                       const Array2DReal &TangentialVelocity,
                       const Array1DReal &SurfaceFrictionVelocity,
                       const Array1DReal &SurfaceBuoyancyFlux,
                       const Array1DReal &WindSpeed10m);

   /// @brief Register fields with I/O system
   void defineFields();

   // Delete copy and move constructors/assignment
   KPPMix(const KPPMix &)            = delete;
   KPPMix &operator=(const KPPMix &) = delete;
   KPPMix(KPPMix &&)                 = delete;
   KPPMix &operator=(KPPMix &&)      = delete;

}; // class KPPMix

} // namespace OMEGA

#endif // OMEGA_KPP_MIX_H
