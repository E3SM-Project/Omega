//===-- ocn/KPPMix.cpp - KPP Boundary Layer Mixing Implementation --*- C++
//-*-===//
//
/// \file
/// \brief Implementation of KPP boundary layer mixing scheme
///
/// This file implements the KPPMix class for computing ocean boundary layer
/// mixing following Large et al. (1994) with optional Langmuir enhancement.
//
//===----------------------------------------------------------------------===//

#include "KPPMix.h"
#include "DataTypes.h"
#include "Error.h"
#include "FillValues.h"
#include "GlobalConstants.h"
#include "KPPConstants.h"
#include "Logging.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"
#include <limits>

namespace OMEGA {

// Singleton instance
KPPMix *KPPMix::Instance = nullptr;

namespace {

bool matchTypeFromString(const std::string &Name, KPPMatchType &Type) {
   if (Name == "SimpleShapes") {
      Type = KPPMatchType::SimpleShapes;
      return true;
   }
   if (Name == "MatchBoth") {
      Type = KPPMatchType::MatchBoth;
      return true;
   }
   return false;
}

const char *matchTypeName(KPPMatchType Type) {
   switch (Type) {
   case KPPMatchType::MatchBoth:
      return "MatchBoth";
   default:
      return "SimpleShapes";
   }
}

} // anonymous namespace

/// Constructor for KPPMix
KPPMix::KPPMix(const std::string &InName, const HorzMesh *InMesh,
               const VertCoord *InVCoord)
    : Name(InName), Mesh(InMesh), VCoord(InVCoord) {

   // Allocate output arrays
   VertDiff = Array2DReal("VertDiff", Mesh->NCellsSize, VCoord->NVertLayersP1);
   VertVisc = Array2DReal("VertVisc", Mesh->NCellsSize, VCoord->NVertLayersP1);
   BoundaryLayerDepth = Array1DReal("BoundaryLayerDepth", Mesh->NCellsSize);
   IndexBoundaryLayerDepth =
       Array1DI4("IndexBoundaryLayerDepth", Mesh->NCellsSize);
   VertNonLocalFlux =
       Array2DReal("VertNonLocalFlux", Mesh->NCellsSize, VCoord->NVertLayersP1);
   BulkRichardsonNumber = Array2DReal("BulkRichardsonNumber", Mesh->NCellsSize,
                                      VCoord->NVertLayersP1);
   BulkRichardsonShear  = Array2DReal("BulkRichardsonShear", Mesh->NCellsSize,
                                      VCoord->NVertLayersP1);
   UnresolvedShear =
       Array2DReal("UnresolvedShear", Mesh->NCellsSize, VCoord->NVertLayersP1);
   BuoyancyJump =
       Array2DReal("BuoyancyJump", Mesh->NCellsSize, VCoord->NVertLayersP1);
   TurbulentVelocityScale = Array2DReal(
       "TurbulentVelocityScale", Mesh->NCellsSize, VCoord->NVertLayersP1);
   PotentialDensity =
       Array2DReal("PotentialDensity", Mesh->NCellsSize, VCoord->NVertLayers);
   SurfaceFrictionVelocity =
       Array1DReal("SurfaceFrictionVelocity", Mesh->NCellsSize);
   SurfaceBuoyancyFlux = Array1DReal("SurfaceBuoyancyFlux", Mesh->NCellsSize);

   // Set field names
   VertDiffFldName            = "VertDiff";
   VertViscFldName            = "VertVisc";
   OBLDepthFldName            = "BoundaryLayerDepth";
   OBLDepthIndexFldName       = "BoundaryLayerDepthIndex";
   NonLocalFluxFldName        = "VertNonLocalFlux";
   BulkRichardsonFldName      = "BulkRichardsonNumber";
   BulkRichardsonShearFldName = "BulkRichardsonShear";
   UnresolvedShearFldName     = "UnresolvedShear";
   BuoyancyJumpFldName        = "BuoyancyJump";
   TurbulentVelScaleFldName   = "TurbulentVelocityScale";
   PotentialDensityFldName    = "PotentialDensity";
   SurfFricVelFldName         = "SurfaceFrictionVelocity";
   SurfBuoyFluxFldName        = "SurfaceBuoyancyFlux";

   if (Name != "Default") {
      VertDiffFldName.append(Name);
      VertViscFldName.append(Name);
      OBLDepthFldName.append(Name);
      OBLDepthIndexFldName.append(Name);
      NonLocalFluxFldName.append(Name);
      BulkRichardsonFldName.append(Name);
      BulkRichardsonShearFldName.append(Name);
      UnresolvedShearFldName.append(Name);
      BuoyancyJumpFldName.append(Name);
      TurbulentVelScaleFldName.append(Name);
      PotentialDensityFldName.append(Name);
      SurfFricVelFldName.append(Name);
      SurfBuoyFluxFldName.append(Name);
   }

   defineFields();
}

/// Destructor for KPPMix
KPPMix::~KPPMix() {}

/// Get singleton instance
KPPMix *KPPMix::getInstance() { return Instance; }

/// Destroy singleton instance
void KPPMix::destroyInstance() {
   delete Instance;
   Instance = nullptr;
}

/// Initialize KPPMix from configuration
void KPPMix::init() {
   if (!Instance) {
      Instance = new KPPMix("Default", HorzMesh::getDefault(),
                            VertCoord::getDefault());
   }

   Error Err;
   KPPMix *DefKPPMix   = KPPMix::getInstance();
   Config *OmegaConfig = Config::getOmegaConfig();

   // Get VertMix config group
   Config VertMixConfig("VertMix");
   Err += OmegaConfig->get(VertMixConfig);
   CHECK_ERROR_ABORT(Err, "KPPMix::init: VertMix group not found in Config");

   // Get KPP config subgroup
   Config KPPConfig("KPP");
   Err += VertMixConfig.get(KPPConfig);
   if (Err.isFail()) {
      LOG_WARN("KPPMix::init: KPP subgroup not found, using defaults");
      return; // Continue with defaults
   }

   // Read KPP parameters
   bool Enable = true;
   Err += KPPConfig.get("Enable", Enable);
   DefKPPMix->Enabled = Enable;

   Err += KPPConfig.get("CriticalBulkRichardsonNumber",
                        DefKPPMix->CriticalRichardson);
   Err += KPPConfig.get("SurfaceLayerExtent", DefKPPMix->SurfaceLayerExtent);

   // KPP matching/profile semantics.
   std::string MatchStr = "SimpleShapes";
   Error MatchErr       = KPPConfig.get("MatchTechnique", MatchStr);
   if (!MatchErr.isSuccess()) {
      MatchErr.reset();
   }
   Error InterpErr = KPPConfig.get("InterpType2", DefKPPMix->InterpType2Str);
   if (!InterpErr.isSuccess()) {
      InterpErr.reset();
   }
   Error EnhancedErr =
       KPPConfig.get("UseEnhancedDiffusion", DefKPPMix->UseEnhancedDiffusion);
   if (!EnhancedErr.isSuccess()) {
      EnhancedErr.reset();
   }
   Error BLDSmoothErr =
       KPPConfig.get("UseBLDSmoothing", DefKPPMix->UseBLDSmoothing);
   if (!BLDSmoothErr.isSuccess()) {
      BLDSmoothErr.reset();
   }

   if (!matchTypeFromString(MatchStr, DefKPPMix->MatchTechnique)) {
      ABORT_ERROR("KPPMix::init: Invalid MatchTechnique='{}', must be "
                  "SimpleShapes or MatchBoth",
                  MatchStr);
   }

   // Wave and flux options
   Err += KPPConfig.get("UseLangmuirCirculation",
                        DefKPPMix->UseLangmuirCirculation);
   Err += KPPConfig.get("IceFractionThresholdForLangmuir",
                        DefKPPMix->IceFractionThresholdForLangmuir);
   Err += KPPConfig.get("IceFractionThresholdForMinimumOBL",
                        DefKPPMix->IceFractionThresholdForMinimumOBL);
   Err +=
       KPPConfig.get("MinimumOBLUnderSeaIce", DefKPPMix->MinimumOBLUnderSeaIce);
   Error DebugErr =
       KPPConfig.get("DebugDiagnostics", DefKPPMix->DebugDiagnostics);
   if (!DebugErr.isSuccess()) {
      DebugErr.reset();
      DefKPPMix->DebugDiagnostics = false;
   }
   if (DefKPPMix->DebugDiagnostics) {
      LOG_WARN("KPP debug diagnostics enabled");
   }

   // Background mixing
   Err += KPPConfig.get("BackgroundViscosity", DefKPPMix->BackgroundVisc);
   Err += KPPConfig.get("BackgroundDiffusivity", DefKPPMix->BackgroundDiff);

   LOG_WARN("KPPMix::init: KPP initialized enabled={} debugDiagnostics={} "
            "match={}",
            DefKPPMix->Enabled, DefKPPMix->DebugDiagnostics,
            matchTypeName(DefKPPMix->MatchTechnique));
}

/// Main computation routine
void KPPMix::computeKPPMix(const Array2DReal &PotentialDensity,
                           const Array2DReal &NormalVelocity,
                           const Array2DReal &TangentialVelocity,
                           const Array1DReal &SurfaceFrictionVelocity,
                           const Array1DReal &SurfaceBuoyancyFlux,
                           const Array2DReal &BruntVaisalaFreqSq,
                           const Array1DReal &IceFraction,
                           const Array1DReal &WindSpeed10m) {

   if (!Enabled) {
      return; // Skip if disabled
   }

   // Retain PotentialDensity for diagnostics/stream output.
   deepCopy(this->PotentialDensity, PotentialDensity);

   // =======================================================================
   // Stage 1: Compute OBL Depth
   // =======================================================================
   computeOBLDepth(PotentialDensity, NormalVelocity, TangentialVelocity,
                   SurfaceFrictionVelocity, SurfaceBuoyancyFlux,
                   BruntVaisalaFreqSq, IceFraction, WindSpeed10m);

   // =======================================================================
   // Stage 2: Compute Mixing Coefficients
   // =======================================================================
   computeMixingCoefficients(PotentialDensity, SurfaceFrictionVelocity,
                             SurfaceBuoyancyFlux);

   if (DebugDiagnostics) {
      logDiagnostics(PotentialDensity, NormalVelocity, TangentialVelocity,
                     SurfaceFrictionVelocity, SurfaceBuoyancyFlux,
                     WindSpeed10m);
   }
}

void KPPMix::logDiagnostics(const Array2DReal &PotentialDensity,
                            const Array2DReal &NormalVelocity,
                            const Array2DReal &TangentialVelocity,
                            const Array1DReal &SurfaceFrictionVelocity,
                            const Array1DReal &SurfaceBuoyancyFlux,
                            const Array1DReal &WindSpeed10m) {

   using namespace KPP;

   const auto MinLayerCellH = createHostMirrorCopy(VCoord->MinLayerCell);
   const auto MaxLayerCellH = createHostMirrorCopy(VCoord->MaxLayerCell);
   const auto ZInterfaceH   = createHostMirrorCopy(VCoord->GeomZInterface);
   const auto SshCellH      = createHostMirrorCopy(VCoord->SshCell);
   const auto DensityH      = createHostMirrorCopy(PotentialDensity);
   const auto UStarH        = createHostMirrorCopy(SurfaceFrictionVelocity);
   const auto B0H           = createHostMirrorCopy(SurfaceBuoyancyFlux);
   const auto OBLDepthH     = createHostMirrorCopy(BoundaryLayerDepth);
   const auto OBLIndexH     = createHostMirrorCopy(IndexBoundaryLayerDepth);
   const auto VertDiffH     = createHostMirrorCopy(VertDiff);
   const auto VertViscH     = createHostMirrorCopy(VertVisc);

   // NormalVelocity and TangentialVelocity are edge-based; not accessed here
   (void)NormalVelocity;
   (void)TangentialVelocity;

   const int NCellsAll = Mesh->NCellsAll;
   if (NCellsAll <= 0) {
      return;
   }

   // Domain-wide diagnostic to avoid misleading single-cell checks.
   Real MaxAbsB0       = 0.0_Real;
   Real MaxAbsUStar    = 0.0_Real;
   Real MaxVertDiff    = 0.0_Real;
   Real MaxVertVisc    = 0.0_Real;
   int MaxAbsB0Cell    = -1;
   int MaxAbsUStarCell = -1;
   int MaxVertDiffCell = -1;
   int MaxVertDiffK    = -1;
   int MaxVertViscCell = -1;
   int MaxVertViscK    = -1;
   for (int C = 0; C < NCellsAll; ++C) {
      const Real AbsB0    = Kokkos::abs(B0H(C));
      const Real AbsUStar = Kokkos::abs(UStarH(C));
      if (AbsB0 > MaxAbsB0) {
         MaxAbsB0     = AbsB0;
         MaxAbsB0Cell = C;
      }
      if (AbsUStar > MaxAbsUStar) {
         MaxAbsUStar     = AbsUStar;
         MaxAbsUStarCell = C;
      }
      const int KCMin = MinLayerCellH(C);
      const int KCMax = MaxLayerCellH(C) + 1;
      for (int K = KCMin; K <= KCMax; ++K) {
         const Real Diff = VertDiffH(C, K);
         const Real Visc = VertViscH(C, K);
         if (Diff > MaxVertDiff) {
            MaxVertDiff     = Diff;
            MaxVertDiffCell = C;
            MaxVertDiffK    = K;
         }
         if (Visc > MaxVertVisc) {
            MaxVertVisc     = Visc;
            MaxVertViscCell = C;
            MaxVertViscK    = K;
         }
      }
   }
   LOG_WARN("KPP debug domain post-coeff: max|b0|={} at cell={} max|u*|={} "
            "at cell={} maxKPPDiff={} at cell={},k={} maxKPPVisc={} at "
            "cell={},k={}",
            MaxAbsB0, MaxAbsB0Cell, MaxAbsUStar, MaxAbsUStarCell, MaxVertDiff,
            MaxVertDiffCell, MaxVertDiffK, MaxVertVisc, MaxVertViscCell,
            MaxVertViscK);

   const int ICell       = 0;
   const int KMin        = MinLayerCellH(ICell);
   const int KMax        = MaxLayerCellH(ICell);
   const int NVertLayers = VCoord->NVertLayers;

   if (KMin > KMax) {
      return;
   }

   const int KSurf     = Kokkos::min(KMin, NVertLayers - 1);
   const Real RhoSurf  = DensityH(ICell, KSurf);
   const Real UStar    = UStarH(ICell);
   const Real UStarEff = Kokkos::fmax(KPP::MinUStar, UStar);
   const Real BuoyFlux = B0H(ICell);
   Real Wind10m        = 0.0_Real;
   if (WindSpeed10m.extent(0) > 0) {
      const auto Wind10mH = createHostMirrorCopy(WindSpeed10m);
      Wind10m             = Wind10mH(ICell);
   }
   const Real LangmuirFactor =
       UseLangmuirCirculation
           ? computeLangmuirEnhancement(Wind10m, UStarEff, 50.0)
           : 1.0_Real;
   const Real BuoyFluxEff = BuoyFlux * LangmuirFactor;

   LOG_WARN("KPP debug: cell={} h_obl={} m k_obl={} u*={} b0={} b0_eff={} "
            "langmuir={}",
            ICell, OBLDepthH(ICell), OBLIndexH(ICell), UStar, BuoyFlux,
            BuoyFluxEff, LangmuirFactor);

   const int KOblIface = Kokkos::min(
       NVertLayers, Kokkos::max(KMin, static_cast<int>(OBLIndexH(ICell)) + 1));
   LOG_WARN("KPP debug coeff target: cell={} k_obl={} iface={} diff={} "
            "visc={}",
            ICell, OBLIndexH(ICell), KOblIface, VertDiffH(ICell, KOblIface),
            VertViscH(ICell, KOblIface));

   const int KTop  = Kokkos::min(KMax, KMin + 3);
   const int KOBL  = OBLIndexH(ICell);
   const Real HOBL = OBLDepthH(ICell);

   for (int K = KMin; K <= KTop; ++K) {
      const int KCell   = Kokkos::min(K, NVertLayers - 1);
      const int KIface  = Kokkos::min(K + 1, NVertLayers);
      const Real ZDepth = SshCellH(ICell) - ZInterfaceH(ICell, KIface);

      const Real RhoK     = DensityH(ICell, KCell);
      const Real DeltaRho = RhoK - RhoSurf;
      const Real DeltaB   = Gravity * DeltaRho / RhoSw;
      const Real WTurb =
          computeTurbVelocityScale(UStarEff, BuoyFluxEff, ZDepth);
      const Real RiBulk = DeltaB * ZDepth / (WTurb * WTurb + 1.0e-12_Real);

      Real Sigma = 0.0_Real;
      if (K <= KOBL) {
         Sigma = -1.0_Real * static_cast<Real>(K - KMin) /
                 static_cast<Real>(KOBL - KMin + 1);
         Sigma = Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, Sigma));
      }

      // Monin-Obukhov coordinate zeta = d/L at this depth
      const Real ZLocal = -Sigma * HOBL;
      Real Zeta         = 0.0_Real;
      const Real Denom  = VonKar * BuoyFlux;
      if (Kokkos::abs(Denom) > 1.0e-16_Real) {
         const Real LMoninObukhov = (UStarEff * UStarEff * UStarEff) / Denom;
         if (Kokkos::abs(LMoninObukhov) > 1.0e-16_Real) {
            Zeta = ZLocal / LMoninObukhov;
         }
      }

      const Real PhiInvM = KPP::kppPhiInvMomentum(Zeta);
      const Real PhiInvS = KPP::kppPhiInvScalar(Zeta);

      LOG_WARN(
          "KPP debug top: cell={} k={} z={} ri_b={} zeta={} phi_m={} phi_s={}",
          ICell, K, ZDepth, RiBulk, Zeta, PhiInvM, PhiInvS);
   }
}

/// Stage 1: Compute OBL depth using bulk Richardson search with edge-based
/// velocity shear following the MPAS CVMix reference implementation
/// (mpas_ocn_vmix_cvmix.F)
void KPPMix::computeOBLDepth(const Array2DReal &PotentialDensity,
                             const Array2DReal &NormalVelocity,
                             const Array2DReal &TangentialVelocity,
                             const Array1DReal &SurfaceFrictionVelocity,
                             const Array1DReal &SurfaceBuoyancyFlux,
                             const Array2DReal &BruntVaisalaFreqSq,
                             const Array1DReal &IceFraction,
                             const Array1DReal &WindSpeed10m) {

   using namespace KPP;

   I4 NVertLayers = VCoord->NVertLayers;

   // =======================================================================
   // Compute Langmuir enhancement factors if wind speed is available
   // =======================================================================
   Array1DReal LangmuirFactor("LangmuirFactor", Mesh->NCellsSize);
   const bool LocUseLangmuirCirculation      = UseLangmuirCirculation;
   const Real LocSurfaceLayerExtent          = SurfaceLayerExtent;
   const Real LocCriticalRichardson          = CriticalRichardson;
   const Real LocIceFracThresholdForLangmuir = IceFractionThresholdForLangmuir;
   parallelFor(
       "KPP-Langmuir", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          const Real IceFrac = IceFraction(ICell);
          if (LocUseLangmuirCirculation &&
              IceFrac < LocIceFracThresholdForLangmuir) {
             const Real UStar = SurfaceFrictionVelocity(ICell);
             const Real Wind10m =
                 (WindSpeed10m.extent(0) > 0) ? WindSpeed10m(ICell) : 0.0_Real;
             LangmuirFactor(ICell) =
                 computeLangmuirEnhancement(Wind10m, UStar, 50.0);
          } else {
             LangmuirFactor(ICell) = 1.0;
          }
       });

   // =======================================================================
   // Stage 1: Compute OBL depth using edge-based velocity shear
   // =======================================================================

   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   OMEGA_SCOPE(ZInterface, VCoord->GeomZInterface);
   OMEGA_SCOPE(ZMid, VCoord->GeomZMid);
   OMEGA_SCOPE(LocSshCell, VCoord->SshCell);
   OMEGA_SCOPE(NEdgesOnCell, Mesh->NEdgesOnCell);
   OMEGA_SCOPE(EdgesOnCell, Mesh->EdgesOnCell);
   OMEGA_SCOPE(CellsOnCell, Mesh->CellsOnCell);
   OMEGA_SCOPE(AreaCell, Mesh->AreaCell);
   OMEGA_SCOPE(DcEdge, Mesh->DcEdge);
   OMEGA_SCOPE(DvEdge, Mesh->DvEdge);
   OMEGA_SCOPE(LocPotentialDensity, PotentialDensity);
   OMEGA_SCOPE(LocNormalVelocity, NormalVelocity);
   OMEGA_SCOPE(LocTangentialVelocity, TangentialVelocity);
   OMEGA_SCOPE(LocBruntVaisalaFreqSq, BruntVaisalaFreqSq);
   OMEGA_SCOPE(LocIceFraction, IceFraction);
   OMEGA_SCOPE(LocLangmuirFactor, LangmuirFactor);
   OMEGA_SCOPE(LocBoundaryLayerDepth, BoundaryLayerDepth);
   OMEGA_SCOPE(LocIndexBoundaryLayerDepth, IndexBoundaryLayerDepth);
   OMEGA_SCOPE(LocBulkRichardson, BulkRichardsonNumber);
   OMEGA_SCOPE(LocBulkRichardsonShear, BulkRichardsonShear);
   OMEGA_SCOPE(LocUnresolvedShear, UnresolvedShear);
   OMEGA_SCOPE(LocBuoyancyJump, BuoyancyJump);
   const bool LocUseBLDSmoothing           = UseBLDSmoothing;
   const bool LocFullRiProfile             = DebugDiagnostics;
   const Real LocIceFracThresholdForMinOBL = IceFractionThresholdForMinimumOBL;
   const Real LocMinimumOBLUnderSeaIce     = MinimumOBLUnderSeaIce;

   deepCopy(BulkRichardsonNumber, 0.0_Real);
   deepCopy(BulkRichardsonShear, 0.0_Real);
   deepCopy(UnresolvedShear, 0.0_Real);
   deepCopy(BuoyancyJump, 0.0_Real);

   // Compile-time bound for the per-thread edge scratch arrays below
   constexpr I4 MaxEdgesBound = HorzMesh::MaxEdgesBound;

   parallelFor(
       "KPP-OBLDepth", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          using namespace KPP;

          const Real UStar =
              Kokkos::fmax(0.0_Real, SurfaceFrictionVelocity(ICell));
          const Real BuoyFlux = SurfaceBuoyancyFlux(ICell);

          const I4 KMin     = MinLayerCell(ICell);
          const I4 KMax     = MaxLayerCell(ICell);
          const I4 KIntTop  = Kokkos::min(KMin + 1, NVertLayers);
          const I4 KIntDeep = Kokkos::min(KMax + 1, NVertLayers);

          const Real IceFrac = LocIceFraction(ICell);

          // KPP depths are measured below the free surface, so geometric
          // heights must be offset by the sea surface height.
          const Real Ssh = LocSshCell(ICell);

          // Default to the full water column; overwritten if Ri crosses.
          Real OBLDepth         = Ssh - ZInterface(ICell, KIntDeep);
          I4 KCross             = -1;
          const Real RiCritical = LocCriticalRichardson;
          // Ri is evaluated at cell centers while the reference average spans
          // the top epsilon*d; this factor corrects for that offset.
          const Real RiScaling   = 1.0_Real - 0.5_Real * LocSurfaceLayerExtent;
          const Real BuoyFluxEff = BuoyFlux * LocLangmuirFactor(ICell);

          // Unresolved shear coefficient, Large et al. (1994) Eq. (23):
          // Vt^2 = Cv * sqrt(-beta_T/(c_s*eps)) / (kappa^2 * Ri_crit) * d*N*w_s
          // CSUnres is c_s for the strongly-unstable scalar branch and VtCoef
          // collects the constant prefactor.
          const Real CSUnres = 24.0_Real * Kokkos::sqrt(17.0_Real);
          const Real VtCoef =
              Kokkos::sqrt(
                  0.2_Real /
                  Kokkos::max(1.0e-12_Real, CSUnres * LocSurfaceLayerExtent)) /
              (VonKar * VonKar);

          // -------------------------------------------------------------------
          // Velocities live on edges (C-grid), so the shear entering Ri is a
          // weighted average over the cell's edges. Weights are the MPAS kite
          // areas (0.25*dc*dv) normalized by the cell area.
          // -------------------------------------------------------------------
          const I4 NEdges                 = NEdgesOnCell(ICell);
          const I4 NEdgesEff              = Kokkos::min(NEdges, MaxEdgesBound);
          bool EdgeValid[MaxEdgesBound]   = {};
          Real EdgeWeights[MaxEdgesBound] = {};
          const Real InvAreaCell =
              1.0_Real / Kokkos::max(AreaCell(ICell), 1.0e-20_Real);
          for (I4 J = 0; J < NEdgesEff; ++J) {
             const I4 IEdge = EdgesOnCell(ICell, J);
             const I4 KEMin = MinLayerEdgeBot(IEdge);
             const I4 KEMax = MaxLayerEdgeTop(IEdge);
             EdgeValid[J] =
                 (KEMax >= KEMin && KEMin >= 0 && KEMin < NVertLayers);
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
             if (SumW < 1.0e-20_Real) {
                // Degenerate kite areas: fall back to equal weighting.
                I4 NEdgesValid = 0;
                for (I4 J = 0; J < NEdgesEff; ++J) {
                   if (EdgeValid[J]) {
                      ++NEdgesValid;
                   }
                }
                if (NEdgesValid > 0) {
                   const Real EqualW =
                       1.0_Real / static_cast<Real>(NEdgesValid);
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

          // -------------------------------------------------------------------
          // Bulk Richardson search, Large et al. (1994) Eq. (21):
          //   Ri_b(d) = (B_r - B(d)) d / (|V_r - V(d)|^2 + Vt^2(d))
          // where the reference values B_r, V_r are averaged over the top
          // epsilon*d of the column. Since epsilon*d grows monotonically with
          // the trial depth, the averaging window only ever extends downward,
          // so the running sums are carried across trial depths rather than
          // rebuilt from the surface at each one.
          // The OBL base is the first d at which Ri_b reaches RiCritical, and
          // the search stops there unless the full profile is wanted for the
          // Ri diagnostics.
          // -------------------------------------------------------------------

          // Cell surface-layer density average
          I4 KSurfaceAvg      = KMin;
          const Real ThickTop = Kokkos::abs(ZInterface(ICell, KMin + 1) -
                                            ZInterface(ICell, KMin));
          Real SumThickness   = Kokkos::max(ThickTop, 1.0e-12_Real);
          Real SumRho         = LocPotentialDensity(ICell, KMin) * SumThickness;

          // Per-edge surface-layer velocity averages
          I4 KSurfE[MaxEdgesBound]      = {};
          Real SumThickE[MaxEdgesBound] = {};
          Real SumUnE[MaxEdgesBound]    = {};
          Real SumVtE[MaxEdgesBound]    = {};

          for (I4 J = 0; J < NEdgesEff; ++J) {
             if (!EdgeValid[J]) {
                continue;
             }
             const I4 IEdge    = EdgesOnCell(ICell, J);
             const I4 KEMin    = MinLayerEdgeBot(IEdge);
             KSurfE[J]         = KEMin;
             const I4 KIntE0   = Kokkos::min(KEMin + 1, NVertLayers);
             const Real Thick0 = Kokkos::abs(ZInterface(ICell, KIntE0) -
                                             ZInterface(ICell, KEMin));
             SumThickE[J]      = Kokkos::max(Thick0, 1.0e-12_Real);
             const I4 KE0      = Kokkos::min(KEMin, NVertLayers - 1);
             SumUnE[J]         = LocNormalVelocity(IEdge, KE0) * SumThickE[J];
             SumVtE[J] = LocTangentialVelocity(IEdge, KE0) * SumThickE[J];
          }

          for (I4 K = KMin; K <= KMax; ++K) {
             const I4 KCell     = Kokkos::min(K, NVertLayers - 1);
             const I4 KInt      = Kokkos::min(K + 1, NVertLayers);
             const Real ZDepth  = Ssh - ZInterface(ICell, KInt);
             const Real ZCenter = Ssh - ZMid(ICell, KCell);
             if (ZDepth < 1.0e-12)
                continue;

             const Real SurfLayerDepth = LocSurfaceLayerExtent * ZDepth;

             // Advance cell surface average for density
             while (KSurfaceAvg < K &&
                    (Ssh - ZInterface(ICell, KSurfaceAvg + 1)) <
                        SurfLayerDepth) {
                ++KSurfaceAvg;
                const I4 KSA  = Kokkos::min(KSurfaceAvg, NVertLayers - 1);
                const Real DZ = Kokkos::abs(ZInterface(ICell, KSurfaceAvg + 1) -
                                            ZInterface(ICell, KSurfaceAvg));
                const Real ThickK = Kokkos::max(DZ, 1.0e-12_Real);
                SumThickness += ThickK;
                SumRho += LocPotentialDensity(ICell, KSA) * ThickK;
             }

             // Advance per-edge surface averages for velocity
             for (I4 J = 0; J < NEdgesEff; ++J) {
                if (!EdgeValid[J]) {
                   continue;
                }
                const I4 IEdge = EdgesOnCell(ICell, J);
                const I4 KEMax = MaxLayerEdgeTop(IEdge);
                while (KSurfE[J] < K &&
                       (Ssh - ZInterface(ICell, KSurfE[J] + 1)) <
                           SurfLayerDepth) {
                   ++KSurfE[J];
                   const I4 KE = Kokkos::min(
                       Kokkos::max(KSurfE[J], MinLayerEdgeBot(IEdge)), KEMax);
                   const Real DZ =
                       Kokkos::abs(ZInterface(ICell, KSurfE[J] + 1) -
                                   ZInterface(ICell, KSurfE[J]));
                   const Real ThickK = Kokkos::max(DZ, 1.0e-12_Real);
                   SumThickE[J] += ThickK;
                   SumUnE[J] += LocNormalVelocity(IEdge, KE) * ThickK;
                   SumVtE[J] += LocTangentialVelocity(IEdge, KE) * ThickK;
                }
             }

             const Real InvSumThickness =
                 1.0_Real / Kokkos::max(SumThickness, 1.0e-12_Real);
             const Real RhoAvgSurf = SumRho * InvSumThickness;

             // Buoyancy jump B_r - B(d), positive for stable stratification
             const Real RhoK              = LocPotentialDensity(ICell, KCell);
             const Real DeltaRho          = RhoK - RhoAvgSurf;
             const Real DeltaB            = Gravity * DeltaRho / RhoSw;
             LocBuoyancyJump(ICell, KInt) = DeltaB;

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
                       1.0_Real / Kokkos::max(SumThickE[J], 1.0e-12_Real);
                   const Real UnAvg = SumUnE[J] * InvThickE;
                   const Real VtAvg = SumVtE[J] * InvThickE;
                   const Real UnK   = LocNormalVelocity(IEdge, KE);
                   const Real VtK   = LocTangentialVelocity(IEdge, KE);
                   const Real DUn   = UnK - UnAvg;
                   const Real DVt   = VtK - VtAvg;
                   DeltaVSq += EdgeWeights[J] * (DUn * DUn + DVt * DVt);
                }
             }
             LocBulkRichardsonShear(ICell, KInt) =
                 Kokkos::max(DeltaVSq, 1.0e-15_Real);

             const Real SigmaLoc = Kokkos::fmin(
                 1.0_Real, Kokkos::fmax(0.0_Real, LocSurfaceLayerExtent));

             // Turbulent scalar velocity scale w_s at the surface-layer depth
             Real WTurb = 0.0_Real;
             if (UStar > 1.0e-12_Real) {
                const Real U3      = UStar * UStar * UStar;
                const Real Zeta    = SigmaLoc * ZDepth * VonKar * BuoyFluxEff /
                                     Kokkos::max(U3, 1.0e-20_Real);
                const Real PhiInvS = KPP::kppPhiInvScalar(Zeta);
                WTurb = VonKar * UStar * Kokkos::max(PhiInvS, 0.0_Real);
             } else if (BuoyFluxEff < 0.0_Real) {
                // Free convection limit: u* drops out and w_s ~ (c_s d B_0)^1/3
                const Real CS  = KPP::CMoS;
                const Real WS3 = -CS * SigmaLoc * ZDepth * VonKar * BuoyFluxEff;
                WTurb = VonKar * Kokkos::pow(Kokkos::max(WS3, 0.0_Real),
                                             1.0_Real / 3.0_Real);
             }

             // Unresolved turbulent shear Vt^2 (m^2/s^2), Large et al. Eq.
             // (23). Cv ramps from 2.1 to 1.7 as stratification strengthens.
             const Real NCntr = Kokkos::sqrt(
                 Kokkos::max(0.0_Real, LocBruntVaisalaFreqSq(ICell, KInt)));
             const Real Cv  = (NCntr < 0.002_Real)
                                  ? (2.1_Real - 200.0_Real * NCntr)
                                  : 1.7_Real;
             const Real Vt2 = Kokkos::max(
                 1.0e-10_Real, Cv * VtCoef * ZCenter * NCntr * WTurb /
                                   Kokkos::max(RiCritical, 1.0e-12_Real));
             LocUnresolvedShear(ICell, KInt) = Vt2;

             const Real VelScaleSq = DeltaVSq + Vt2;

             const Real RiBulk = RiScaling * DeltaB * ZCenter /
                                 Kokkos::max(VelScaleSq, 1.0e-12_Real);
             LocBulkRichardson(ICell, KInt) = RiBulk;

             if (KCross < 0 && RiBulk > RiCritical) {
                KCross = K;
                // Levels below the crossing only feed the Ri diagnostics
                if (!LocFullRiProfile)
                   break;
             }
          }

          if (KCross >= KMin) {
             if (KCross > KMin) {
                // Ri values are defined at cell centers, so interpolate on
                // center depths to keep the abscissa consistent.
                const I4 KAbove    = Kokkos::max(KMin, KCross - 1);
                const I4 KBelow    = Kokkos::min(KCross, NVertLayers - 1);
                const I4 KAboveRi  = Kokkos::min(KAbove + 1, NVertLayers);
                const I4 KBelowRi  = Kokkos::min(KBelow + 1, NVertLayers);
                const Real ZAbove  = Ssh - ZMid(ICell, KAbove);
                const Real ZBelow  = Ssh - ZMid(ICell, KBelow);
                const Real RiAbove = LocBulkRichardson(ICell, KAboveRi);
                const Real RiBelow = LocBulkRichardson(ICell, KBelowRi);

                const Real H = ZBelow - ZAbove;
                if (H > 1.0e-12_Real) {
                   // CVMix-style QUAD interpolation for OBL crossing:
                   // - first interior crossing uses zero slope at top point
                   // - deeper crossings use upstream slope
                   Real SlopeAbove = 0.0_Real;
                   if (KCross > KMin + 1) {
                      const I4 KPrev    = Kokkos::max(KMin, KAbove - 1);
                      const I4 KPrevRi  = Kokkos::min(KPrev + 1, NVertLayers);
                      const Real ZPrev  = Ssh - ZMid(ICell, KPrev);
                      const Real RiPrev = LocBulkRichardson(ICell, KPrevRi);
                      const Real DZPrev = ZAbove - ZPrev;
                      if (Kokkos::abs(DZPrev) > 1.0e-12_Real) {
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
                      if (Kokkos::abs(DRi) > 1.0e-12_Real) {
                         const Real Frac = Kokkos::fmax(
                             0.0_Real,
                             Kokkos::fmin(1.0_Real,
                                          (RiCritical - RiAbove) / DRi));
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
          OBLDepth               = Kokkos::fmax(OBLDepth, MinOBLDepth);
          if (IceFrac > LocIceFracThresholdForMinOBL) {
             OBLDepth = Kokkos::fmax(OBLDepth, LocMinimumOBLUnderSeaIce);
          }
          OBLDepth = Kokkos::fmin(OBLDepth, MaxOBLDepth);

          I4 KFinal = KMax;
          for (I4 K = KMin; K < KMax; ++K) {
             const Real ZAbove = Ssh - ZInterface(ICell, K);
             const Real ZBelow = Ssh - ZInterface(ICell, K + 1);
             if (OBLDepth >= ZAbove && OBLDepth <= ZBelow) {
                KFinal = K;
                break;
             }
          }

          LocBoundaryLayerDepth(ICell)      = OBLDepth;
          LocIndexBoundaryLayerDepth(ICell) = KFinal;
       });

   if (LocUseBLDSmoothing) {
      Array1DReal BoundaryLayerDepthSmooth("BoundaryLayerDepthSmooth",
                                           Mesh->NCellsSize);
      OMEGA_SCOPE(LocBoundaryLayerDepthSmooth, BoundaryLayerDepthSmooth);
      OMEGA_SCOPE(LocNCellsAll, Mesh->NCellsAll);

      // Area-weighted smoothing of the BLD over each cell and its neighbors
      // (MPAS-Ocean cvmix convention). This suppresses the grid-scale noise
      // that the discrete Ri crossing search can introduce.
      parallelFor(
          "KPP-OBLDepth-Smooth", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
             const I4 KMin = MinLayerCell(ICell);
             if (KMin < 0 || KMin >= NVertLayers) {
                LocBoundaryLayerDepthSmooth(ICell) =
                    LocBoundaryLayerDepth(ICell);
                return;
             }

             const I4 NEdges = NEdgesOnCell(ICell);
             Real AreaSum    = 0.0_Real;
             Real BLDSum     = 0.0_Real;
             I4 EdgeCount    = 0;

             for (I4 J = 0; J < NEdges; ++J) {
                const I4 INeighbor = CellsOnCell(ICell, J);
                if (INeighbor == LocNCellsAll) {
                   continue;
                }

                const I4 KMinNbr = MinLayerCell(INeighbor);
                if (KMinNbr < 0 || KMinNbr >= NVertLayers) {
                   continue;
                }

                const Real NbrArea = AreaCell(INeighbor);
                BLDSum += 2.0_Real * NbrArea * LocBoundaryLayerDepth(INeighbor);
                AreaSum += 2.0_Real * NbrArea;
                ++EdgeCount;
             }

             if (EdgeCount > 0) {
                const Real SelfArea = AreaCell(ICell);
                BLDSum += LocBoundaryLayerDepth(ICell) *
                          static_cast<Real>(EdgeCount) * SelfArea;
                AreaSum += static_cast<Real>(EdgeCount) * SelfArea;
             }

             if (AreaSum > 0.0_Real) {
                LocBoundaryLayerDepthSmooth(ICell) = BLDSum / AreaSum;
             } else {
                LocBoundaryLayerDepthSmooth(ICell) =
                    LocBoundaryLayerDepth(ICell);
             }
          });

      parallelFor(
          "KPP-OBLDepth-CommitSmooth", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell) {
             const I4 KMin = MinLayerCell(ICell);
             const I4 KMax = MaxLayerCell(ICell);
             if (KMin < 0 || KMax < KMin || KMin >= NVertLayers) {
                return;
             }

             const Real Ssh = LocSshCell(ICell);

             const I4 KIntTop             = Kokkos::min(KMin + 1, NVertLayers);
             const Real TopLayerThickness = Kokkos::abs(
                 ZInterface(ICell, KIntTop) - ZInterface(ICell, KMin));
             const Real MinOBLDepth = 0.5_Real * TopLayerThickness;
             const Real MaxOBLDepth = Ssh - ZMid(ICell, KMax);

             Real OBLDepth = LocBoundaryLayerDepthSmooth(ICell);
             OBLDepth      = Kokkos::fmax(OBLDepth, MinOBLDepth);
             OBLDepth      = Kokkos::fmin(OBLDepth, MaxOBLDepth);

             I4 KFinal = KMax;
             for (I4 K = KMin; K < KMax; ++K) {
                const Real ZAbove = Ssh - ZInterface(ICell, K);
                const Real ZBelow = Ssh - ZInterface(ICell, K + 1);
                if (OBLDepth >= ZAbove && OBLDepth <= ZBelow) {
                   KFinal = K;
                   break;
                }
             }

             LocBoundaryLayerDepth(ICell)      = OBLDepth;
             LocIndexBoundaryLayerDepth(ICell) = KFinal;
          });
   }

   LOG_INFO("KPPMix::computeOBLDepth: OBL depth computed");
}

/// Stage 2: Compute KPP mixing contribution or matched coefficients
void KPPMix::computeMixingCoefficients(
    const Array2DReal &PotentialDensity,
    const Array1DReal &SurfaceFrictionVelocity,
    const Array1DReal &SurfaceBuoyancyFlux, const Array2DReal &InteriorVertDiff,
    const Array2DReal &InteriorVertVisc) {

   using namespace KPP;

   (void)PotentialDensity;

   I4 NVertLayers = VCoord->NVertLayers;

   // =======================================================================
   // Capture data for use in lambda
   // =======================================================================
   OMEGA_SCOPE(LocBoundaryLayerDepth, BoundaryLayerDepth);
   OMEGA_SCOPE(LocIndexBoundaryLayerDepth, IndexBoundaryLayerDepth);
   OMEGA_SCOPE(LocVertDiff, VertDiff);
   OMEGA_SCOPE(LocVertVisc, VertVisc);
   OMEGA_SCOPE(LocVertNonLocalFlux, VertNonLocalFlux);
   OMEGA_SCOPE(LocTurbulentVelocityScale, TurbulentVelocityScale);
   OMEGA_SCOPE(LocSurfaceFrictionVelocity, SurfaceFrictionVelocity);
   OMEGA_SCOPE(LocSurfaceBuoyancyFlux, SurfaceBuoyancyFlux);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(ZInterface, VCoord->GeomZInterface);
   OMEGA_SCOPE(ZMid, VCoord->GeomZMid);
   OMEGA_SCOPE(LocSshCell, VCoord->SshCell);
   OMEGA_SCOPE(LocInteriorVertDiff, InteriorVertDiff);
   OMEGA_SCOPE(LocInteriorVertVisc, InteriorVertVisc);

   // Capture member variables for use in lambda
   const Real LocSurfaceLayerExtent = SurfaceLayerExtent;
   const KPPMatchType LocMatch      = MatchTechnique;
   // Non-local flux normalization constant from Large et al. (1994) eq. 20:
   // C_s = C* * kappa * (c_s * kappa * epsilon)^(1/3)
   // where C* = 10, c_s = CMoS = 98.9545, kappa = VonKar, epsilon =
   // SurfaceLayerExtent
   const Real LocNonLocalCs =
       10.0_Real * VonKar *
       Kokkos::pow(KPP::CMoS * VonKar * LocSurfaceLayerExtent,
                   1.0_Real / 3.0_Real);
   bool LocUseEnhancedDiffusion = UseEnhancedDiffusion;
   const Real LocKappa          = VonKar;
   const bool LocUseInteriorMix =
       InteriorVertDiff.data() != nullptr && InteriorVertVisc.data() != nullptr;

   // =======================================================================
   // Initialize with zero KPP contribution, or precomputed interior mixing for
   // matched-coefficient construction.
   // =======================================================================
   parallelFor(
       "KPP-Coeffs-Init", {Mesh->NCellsAll, NVertLayers + 1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          LocVertDiff(ICell, K) =
              LocUseInteriorMix ? LocInteriorVertDiff(ICell, K) : 0.0_Real;
          LocVertVisc(ICell, K) =
              LocUseInteriorMix ? LocInteriorVertVisc(ICell, K) : 0.0_Real;
          LocVertNonLocalFlux(ICell, K)       = 0.0;
          LocTurbulentVelocityScale(ICell, K) = 0.0;
       });

   // =======================================================================
   // Stage 2: Compute KPP profile-based mixing coefficients
   // =======================================================================

   parallelFor(
       "KPP-MixingCoeffs", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          // Get OBL information for this cell
          Real HOBL = LocBoundaryLayerDepth(ICell);

          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell);
          if (KMin < 0 || KMin >= NVertLayers || KMax < KMin) {
             return;
          }
          const I4 KMatch =
              Kokkos::min(KMax + 1, LocIndexBoundaryLayerDepth(ICell) + 1);

          // KPP depths are measured below the free surface, so geometric
          // heights must be offset by the sea surface height.
          const Real Ssh = LocSshCell(ICell);

          // =============================================================
          // Compute turbulent velocity scales
          // =============================================================
          Real UStar    = LocSurfaceFrictionVelocity(ICell);
          Real BuoyFlux = LocSurfaceBuoyancyFlux(ICell);

          // =============================================================
          // Compute mixing coefficients at each interface
          // =============================================================
          for (I4 K = KMin; K <= KMax + 1; ++K) {
             const I4 KIface   = Kokkos::min(Kokkos::max(K, 0), NVertLayers);
             const Real ZDepth = Ssh - ZInterface(ICell, KIface);

             // Check if within OBL using depth below the free surface.
             if (ZDepth <= HOBL && HOBL > 0.0_Real) {
                // Normalized depth in Omega sign convention: sigma in [-1,0].
                Real Sigma = -ZDepth / HOBL;
                Sigma = Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, Sigma));

                // CVMix-style turbulent scales: w = kappa*u*/phi in general,
                // with explicit free-convection limits when u*=0. The scales
                // are frozen at the surface-layer depth below the surface
                // layer, so SigmaLoc is capped at SurfaceLayerExtent.
                const Real SigmaCoord = -Sigma; // [0,1]
                const Real SigmaLoc   = Kokkos::fmin(
                    LocSurfaceLayerExtent, Kokkos::fmax(0.0_Real, SigmaCoord));

                Real Zeta   = 0.0_Real;
                Real WMTurb = 0.0_Real;
                Real WSTurb = 0.0_Real;

                if (UStar > 0.0_Real) {
                   const Real U3 = UStar * UStar * UStar;
                   Zeta          = SigmaLoc * HOBL * BuoyFlux * LocKappa /
                                   Kokkos::max(U3, 1.0e-20_Real);

                   // These return phi^{-1}; do not invert again.
                   const Real PhiInvM = KPP::kppPhiInvMomentum(Zeta);
                   const Real PhiInvS = KPP::kppPhiInvScalar(Zeta);

                   WMTurb = LocKappa * UStar * Kokkos::max(PhiInvM, 0.0_Real);
                   WSTurb = LocKappa * UStar * Kokkos::max(PhiInvS, 0.0_Real);
                } else if (BuoyFlux < 0.0_Real) {
                   // Free-convection edge case (u*=0, unstable forcing).
                   const Real CM  = KPP::CMoM;
                   const Real CS  = KPP::CMoS;
                   const Real WM3 = -CM * SigmaLoc * HOBL * LocKappa * BuoyFlux;
                   const Real WS3 = -CS * SigmaLoc * HOBL * LocKappa * BuoyFlux;
                   WMTurb = LocKappa * Kokkos::pow(Kokkos::max(0.0_Real, WM3),
                                                   1.0_Real / 3.0_Real);
                   WSTurb = LocKappa * Kokkos::pow(Kokkos::max(0.0_Real, WS3),
                                                   1.0_Real / 3.0_Real);
                }

                // For MatchBoth, the shape value the KPP profile must reach at
                // the OBL base so that it joins the interior coefficient there.
                const Real MatchViscShape =
                    (LocUseInteriorMix && LocMatch == KPPMatchType::MatchBoth &&
                     HOBL > 0.0_Real && WMTurb > 0.0_Real)
                        ? LocInteriorVertVisc(ICell, KMatch) /
                              Kokkos::max(HOBL * WMTurb, 1.0e-20_Real)
                        : 0.0_Real;
                const Real MatchDiffShape =
                    (LocUseInteriorMix && LocMatch == KPPMatchType::MatchBoth &&
                     HOBL > 0.0_Real && WSTurb > 0.0_Real)
                        ? LocInteriorVertDiff(ICell, KMatch) /
                              Kokkos::max(HOBL * WSTurb, 1.0e-20_Real)
                        : 0.0_Real;

                // ========================================================
                // Momentum mixing contribution.
                // ========================================================
                Real ShapeM =
                    (LocUseInteriorMix && LocMatch == KPPMatchType::MatchBoth)
                        ? KPP::kppShapeMatched(Sigma, MatchViscShape)
                        : KPP::kppShapeMomentum(Sigma);
                LocVertVisc(ICell, K) = HOBL * WMTurb * ShapeM;

                // ========================================================
                // Tracer mixing contribution.
                // ========================================================
                Real ShapeS =
                    (LocUseInteriorMix && LocMatch == KPPMatchType::MatchBoth)
                        ? KPP::kppShapeMatched(Sigma, MatchDiffShape)
                        : KPP::kppShapeScalar(Sigma);
                LocVertDiff(ICell, K)               = HOBL * WSTurb * ShapeS;
                LocTurbulentVelocityScale(ICell, K) = WSTurb;

                // ========================================================
                // Non-local flux: C_s * G(sigma), reusing the scalar
                // diffusivity shape so gamma and K share one profile.
                // C_s = C* * kappa * (c_s * kappa * epsilon)^(1/3)
                // per Large et al. (1994) eq. 20 (~6.33 with default constants)
                // ========================================================
                // Match CVMix behavior: apply non-local term only when
                // surface buoyancy forcing is unstable/neutral.
                if (BuoyFlux <= 0.0_Real) {
                   LocVertNonLocalFlux(ICell, K) = LocNonLocalCs * ShapeS;
                } else {
                   LocVertNonLocalFlux(ICell, K) = 0.0;
                }

             } else {
                // Below OBL: preserve interior values for MatchBoth, otherwise
                // no KPP contribution.
                LocVertDiff(ICell, K) = LocUseInteriorMix
                                            ? LocInteriorVertDiff(ICell, K)
                                            : 0.0_Real;
                LocVertVisc(ICell, K) = LocUseInteriorMix
                                            ? LocInteriorVertVisc(ICell, K)
                                            : 0.0_Real;
                LocVertNonLocalFlux(ICell, K)       = 0.0;
                LocTurbulentVelocityScale(ICell, K) = 0.0;
             }
          }

          // Optional enhanced diffusion/viscosity treatment at OBL base.
          // Match CVMix Appendix D weighting at the interface nearest HOBL:
          // the OBL base rarely lands on an interface, so the coefficient at
          // the neighboring interface KTarget is replaced by a quadratic blend
          // of the KPP value extrapolated from KKtup and the value already
          // there, weighted by where HOBL falls between the two cell centers.
          if (LocUseEnhancedDiffusion && HOBL > 0.0_Real) {
             const I4 KOBL = Kokkos::max(
                 KMin, Kokkos::min(LocIndexBoundaryLayerDepth(ICell), KMax));
             const Real ZMidOBL = Ssh - ZMid(ICell, KOBL);

             const bool TargetOutsideOBL = HOBL >= ZMidOBL;
             const I4 KKtup =
                 TargetOutsideOBL ? KOBL : Kokkos::max(KMin, KOBL - 1);
             const I4 KTarget = TargetOutsideOBL
                                    ? Kokkos::min(KOBL + 1, KMax + 1)
                                    : Kokkos::max(KMin + 1, KOBL);

             const Real ZKtup = Ssh - ZMid(ICell, KKtup);
             const Real ZNext = (KKtup < KMax)
                                    ? (Ssh - ZMid(ICell, KKtup + 1))
                                    : (Ssh - ZInterface(ICell, KKtup + 1));
             const Real Delta = Kokkos::fmax(
                 0.0_Real,
                 Kokkos::fmin(1.0_Real,
                              (HOBL - ZKtup) /
                                  Kokkos::max(ZNext - ZKtup, 1.0e-12_Real)));
             const Real OneMinusDelta = 1.0_Real - Delta;

             Real SigmaKtup = -ZKtup / HOBL;
             SigmaKtup =
                 Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, SigmaKtup));
             const Real SigmaCoord = -SigmaKtup;
             const Real SigmaLoc   = Kokkos::fmin(
                 LocSurfaceLayerExtent, Kokkos::fmax(0.0_Real, SigmaCoord));

             Real WMKtup = 0.0_Real;
             Real WSKtup = 0.0_Real;
             if (UStar > 0.0_Real) {
                const Real U3   = UStar * UStar * UStar;
                const Real Zeta = SigmaLoc * HOBL * BuoyFlux * LocKappa /
                                  Kokkos::max(U3, 1.0e-20_Real);
                WMKtup = LocKappa * UStar *
                         Kokkos::max(KPP::kppPhiInvMomentum(Zeta), 0.0_Real);
                WSKtup = LocKappa * UStar *
                         Kokkos::max(KPP::kppPhiInvScalar(Zeta), 0.0_Real);
             } else if (BuoyFlux < 0.0_Real) {
                const Real WM3 =
                    -KPP::CMoM * SigmaLoc * HOBL * LocKappa * BuoyFlux;
                const Real WS3 =
                    -KPP::CMoS * SigmaLoc * HOBL * LocKappa * BuoyFlux;
                WMKtup = LocKappa * Kokkos::pow(Kokkos::max(0.0_Real, WM3),
                                                1.0_Real / 3.0_Real);
                WSKtup = LocKappa * Kokkos::pow(Kokkos::max(0.0_Real, WS3),
                                                1.0_Real / 3.0_Real);
             }

             const Real MatchViscShape =
                 (LocUseInteriorMix && LocMatch == KPPMatchType::MatchBoth &&
                  HOBL > 0.0_Real && WMKtup > 0.0_Real)
                     ? LocInteriorVertVisc(ICell, KMatch) /
                           Kokkos::max(HOBL * WMKtup, 1.0e-20_Real)
                     : 0.0_Real;
             const Real MatchDiffShape =
                 (LocUseInteriorMix && LocMatch == KPPMatchType::MatchBoth &&
                  HOBL > 0.0_Real && WSKtup > 0.0_Real)
                     ? LocInteriorVertDiff(ICell, KMatch) /
                           Kokkos::max(HOBL * WSKtup, 1.0e-20_Real)
                     : 0.0_Real;

             const Real ViscKtup =
                 HOBL * WMKtup *
                 ((LocUseInteriorMix && LocMatch == KPPMatchType::MatchBoth)
                      ? KPP::kppShapeMatched(SigmaKtup, MatchViscShape)
                      : KPP::kppShapeMomentum(SigmaKtup));
             const Real DiffKtup =
                 HOBL * WSKtup *
                 ((LocUseInteriorMix && LocMatch == KPPMatchType::MatchBoth)
                      ? KPP::kppShapeMatched(SigmaKtup, MatchDiffShape)
                      : KPP::kppShapeScalar(SigmaKtup));

             const Real ViscProfile = LocVertVisc(ICell, KTarget);
             const Real DiffProfile = LocVertDiff(ICell, KTarget);

             const Real EnhVisc = OneMinusDelta * OneMinusDelta * ViscKtup +
                                  Delta * Delta * ViscProfile;
             const Real EnhDiff = OneMinusDelta * OneMinusDelta * DiffKtup +
                                  Delta * Delta * DiffProfile;

             const Real OldVisc = LocUseInteriorMix
                                      ? LocInteriorVertVisc(ICell, KTarget)
                                      : 0.0_Real;
             const Real OldDiff = LocUseInteriorMix
                                      ? LocInteriorVertDiff(ICell, KTarget)
                                      : 0.0_Real;
             const Real NewVisc = OneMinusDelta * OldVisc + Delta * EnhVisc;
             const Real NewDiff = OneMinusDelta * OldDiff + Delta * EnhDiff;

             LocVertVisc(ICell, KTarget) = NewVisc;
             LocVertDiff(ICell, KTarget) = NewDiff;

             // Keep the non-local term consistent with the rescaled diffusivity
             if (!TargetOutsideOBL && DiffProfile != 0.0_Real) {
                LocVertNonLocalFlux(ICell, KTarget) =
                    LocVertNonLocalFlux(ICell, KTarget) * NewDiff / DiffProfile;
             } else if (!TargetOutsideOBL) {
                LocVertNonLocalFlux(ICell, KTarget) = 0.0_Real;
             }
          }
       });

   LOG_INFO("KPPMix::computeMixingCoefficients: Phase 2 mixing coefficients "
            "computed");
}

/// Register fields with I/O system
void KPPMix::defineFields() {
   // BoundaryLayerDepth on cells
   std::vector<std::string> CellDims(1);
   CellDims[0] = "NCells";
   auto OBLDepthField =
       Field::create(OBLDepthFldName,                  // field name
                     "ocean boundary layer depth",     // long name
                     "m",                              // units
                     "",                               // CF standard name
                     0.0,                              // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     1,                                // number of dims
                     CellDims);

   auto OBLDepthIndexField =
       Field::create(OBLDepthIndexFldName,               // field name
                     "ocean boundary layer depth index", // long name
                     "",                                 // units
                     "",                                 // CF standard name
                     -1,                                 // min valid value
                     std::numeric_limits<I4>::max(),     // max valid value
                     1,                                  // number of dims
                     CellDims);

   // KPP non-local tracer flux profile on cell-layer interfaces
   std::vector<std::string> FluxDims(2);
   FluxDims[0] = "NCells";
   FluxDims[1] = "NVertLayersP1";
   auto NonLocalFluxField =
       Field::create(NonLocalFluxFldName,                 // field name
                     "KPP non-local tracer flux profile", // long name
                     "1",                                 // units
                     "",                                  // CF standard name
                     std::numeric_limits<Real>::lowest(), // min valid value
                     std::numeric_limits<Real>::max(),    // max valid value
                     2,                                   // number of dims
                     FluxDims);

   auto BulkRichardsonField =
       Field::create(BulkRichardsonFldName,               // field name
                     "bulk Richardson number",            // long name
                     "1",                                 // units
                     "",                                  // CF standard name
                     std::numeric_limits<Real>::lowest(), // min valid value
                     std::numeric_limits<Real>::max(),    // max valid value
                     2,                                   // number of dims
                     FluxDims);

   auto BulkRichardsonShearField =
       Field::create(BulkRichardsonShearFldName,       // field name
                     "bulk Richardson shear term",     // long name
                     "m2 s-2",                         // units
                     "",                               // CF standard name
                     0.0,                              // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     2,                                // number of dims
                     FluxDims);

   auto UnresolvedShearField =
       Field::create(UnresolvedShearFldName,           // field name
                     "KPP unresolved shear term Vt2",  // long name
                     "m2 s-2",                         // units
                     "",                               // CF standard name
                     0.0,                              // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     2,                                // number of dims
                     FluxDims);

   auto BuoyancyJumpField =
       Field::create(BuoyancyJumpFldName,                   // field name
                     "KPP buoyancy jump (density anomaly)", // long name
                     "m s-2",                               // units
                     "",                                    // CF standard name
                     std::numeric_limits<Real>::lowest(),   // min valid value
                     std::numeric_limits<Real>::max(),      // max valid value
                     2,                                     // number of dims
                     FluxDims);

   auto TurbulentVelScaleField =
       Field::create(TurbulentVelScaleFldName,         // field name
                     "KPP turbulent velocity scale",   // long name
                     "m s-1",                          // units
                     "",                               // CF standard name
                     0.0,                              // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     2,                                // number of dims
                     FluxDims);

   std::vector<std::string> LayerDims(2);
   LayerDims[0] = "NCells";
   LayerDims[1] = "NVertLayers";
   auto PotentialDensityField =
       Field::create(PotentialDensityFldName,             // field name
                     "KPP potential density",             // long name
                     "kg m-3",                            // units
                     "",                                  // CF standard name
                     std::numeric_limits<Real>::lowest(), // min valid value
                     std::numeric_limits<Real>::max(),    // max valid value
                     2,                                   // number of dims
                     LayerDims);

   // Group KPP-specific outputs for convenient stream selection.
   auto KPPGroup = FieldGroup::create("KPPMix");
   KPPGroup->addField(OBLDepthFldName);
   KPPGroup->addField(OBLDepthIndexFldName);
   KPPGroup->addField(NonLocalFluxFldName);
   KPPGroup->addField(BulkRichardsonFldName);
   KPPGroup->addField(BulkRichardsonShearFldName);
   KPPGroup->addField(UnresolvedShearFldName);
   KPPGroup->addField(BuoyancyJumpFldName);
   KPPGroup->addField(TurbulentVelScaleFldName);
   KPPGroup->addField(PotentialDensityFldName);
   KPPGroup->addField(SurfFricVelFldName);
   KPPGroup->addField(SurfBuoyFluxFldName);

   OBLDepthField->attachData<Array1DReal>(BoundaryLayerDepth, false);
   OBLDepthIndexField->attachData<Array1DI4>(IndexBoundaryLayerDepth, false);
   NonLocalFluxField->attachData<Array2DReal>(VertNonLocalFlux, false);
   BulkRichardsonField->attachData<Array2DReal>(BulkRichardsonNumber, false);
   BulkRichardsonShearField->attachData<Array2DReal>(BulkRichardsonShear,
                                                     false);
   UnresolvedShearField->attachData<Array2DReal>(UnresolvedShear, false);
   BuoyancyJumpField->attachData<Array2DReal>(BuoyancyJump, false);
   TurbulentVelScaleField->attachData<Array2DReal>(TurbulentVelocityScale,
                                                   false);
   PotentialDensityField->attachData<Array2DReal>(PotentialDensity, false);

   // Surface forcing fields on cells
   auto SurfFricVelField =
       Field::create(SurfFricVelFldName,                 // field name
                     "KPP surface friction velocity u*", // long name
                     "m s-1",                            // units
                     "",                                 // CF standard name
                     0.0,                                // min valid value
                     std::numeric_limits<Real>::max(),   // max valid value
                     1,                                  // number of dims
                     CellDims);
   SurfFricVelField->attachData<Array1DReal>(SurfaceFrictionVelocity, false);

   auto SurfBuoyFluxField =
       Field::create(SurfBuoyFluxFldName,                 // field name
                     "KPP surface buoyancy flux",         // long name
                     "m2 s-3",                            // units
                     "",                                  // CF standard name
                     std::numeric_limits<Real>::lowest(), // min valid value
                     std::numeric_limits<Real>::max(),    // max valid value
                     1,                                   // number of dims
                     CellDims);
   SurfBuoyFluxField->attachData<Array1DReal>(SurfaceBuoyancyFlux, false);

   OBLDepthField->addMetadata("_FillValue", FillValueReal);
   OBLDepthIndexField->addMetadata("_FillValue", FillValueI4);
   NonLocalFluxField->addMetadata("_FillValue", FillValueReal);
   BulkRichardsonField->addMetadata("_FillValue", FillValueReal);
   BulkRichardsonShearField->addMetadata("_FillValue", FillValueReal);
   UnresolvedShearField->addMetadata("_FillValue", FillValueReal);
   BuoyancyJumpField->addMetadata("_FillValue", FillValueReal);
   TurbulentVelScaleField->addMetadata("_FillValue", FillValueReal);
   PotentialDensityField->addMetadata("_FillValue", FillValueReal);
   SurfFricVelField->addMetadata("_FillValue", FillValueReal);
   SurfBuoyFluxField->addMetadata("_FillValue", FillValueReal);

   LOG_INFO("KPPMix::defineFields: registered {}, {}, {}, {}, {}, {}, {}, {}, "
            "{}, {}, {}",
            OBLDepthFldName, OBLDepthIndexFldName, NonLocalFluxFldName,
            BulkRichardsonFldName, BulkRichardsonShearFldName,
            UnresolvedShearFldName, BuoyancyJumpFldName,
            TurbulentVelScaleFldName, PotentialDensityFldName,
            SurfFricVelFldName, SurfBuoyFluxFldName);
}

} // namespace OMEGA
