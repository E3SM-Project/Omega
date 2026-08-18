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

/// Constructor for KPPMix
KPPMix::KPPMix(const std::string &Name_in, const HorzMesh *Mesh_in,
               const VertCoord *VCoord_in)
    : Name(Name_in), Mesh(Mesh_in), VCoord(VCoord_in) {

   // Allocate output arrays
   VertDiff = Array2DReal("VertDiff", Mesh->NCellsAll, VCoord->NVertLayers + 1);
   VertVisc = Array2DReal("VertVisc", Mesh->NCellsAll, VCoord->NVertLayers + 1);
   BoundaryLayerDepth = Array1DReal("BoundaryLayerDepth", Mesh->NCellsAll);
   IndexBoundaryLayerDepth =
       Array1DI4("IndexBoundaryLayerDepth", Mesh->NCellsAll);
   VertNonLocalFlux     = Array2DReal("VertNonLocalFlux", Mesh->NCellsAll,
                                      VCoord->NVertLayers + 1);
   BulkRichardsonNumber = Array2DReal("BulkRichardsonNumber", Mesh->NCellsAll,
                                      VCoord->NVertLayers + 1);
   BulkRichardsonShear  = Array2DReal("BulkRichardsonShear", Mesh->NCellsAll,
                                      VCoord->NVertLayers + 1);
   UnresolvedShear =
       Array2DReal("UnresolvedShear", Mesh->NCellsAll, VCoord->NVertLayers + 1);
   BuoyancyJump =
       Array2DReal("BuoyancyJump", Mesh->NCellsAll, VCoord->NVertLayers + 1);
   TurbulentVelocityScale = Array2DReal(
       "TurbulentVelocityScale", Mesh->NCellsAll, VCoord->NVertLayers + 1);
   PotentialDensity =
       Array2DReal("PotentialDensity", Mesh->NCellsAll, VCoord->NVertLayers);
   SurfaceFrictionVelocity =
       Array1DReal("SurfaceFrictionVelocity", Mesh->NCellsAll);
   SurfaceBuoyancyFlux = Array1DReal("SurfaceBuoyancyFlux", Mesh->NCellsAll);

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
   bool enable = true;
   Err += KPPConfig.get("Enable", enable);
   DefKPPMix->Enabled = enable;

   Err += KPPConfig.get("CriticalBulkRichardsonNumber",
                        DefKPPMix->CriticalRichardson);
   Err += KPPConfig.get("StopOBLSearch", DefKPPMix->StopOBLSearchMult);
   Err += KPPConfig.get("SurfaceLayerExtent", DefKPPMix->SurfaceLayerExtent);

   // KPP matching/profile semantics.
   Error MatchErr =
       KPPConfig.get("MatchTechnique", DefKPPMix->MatchTechniqueStr);
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

   // Keep active options focused on what is used in OMEGA.
   if (DefKPPMix->MatchTechniqueStr == "MatchGradient") {
      LOG_INFO("KPPMix::init: MatchGradient is deprecated/unused in OMEGA; "
               "mapping to SimpleShapes");
      DefKPPMix->MatchTechniqueStr = "SimpleShapes";
   }
   if (DefKPPMix->MatchTechniqueStr != "SimpleShapes" &&
       DefKPPMix->MatchTechniqueStr != "MatchBoth" &&
       DefKPPMix->MatchTechniqueStr != "ParabolicNonLocal") {
      LOG_INFO(
          "KPPMix::init: Unsupported MatchTechnique='{}', using SimpleShapes",
          DefKPPMix->MatchTechniqueStr);
      DefKPPMix->MatchTechniqueStr = "SimpleShapes";
   }

   // Wave and flux options
   Err += KPPConfig.get("UseLangmuirCirculation",
                        DefKPPMix->UseLangmuirCirculation);
   Err += KPPConfig.get("UseNonLocalFlux", DefKPPMix->UseNonLocalFlux);
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
            DefKPPMix->MatchTechniqueStr);
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
   Real maxAbsB0       = 0.0_Real;
   Real maxAbsUStar    = 0.0_Real;
   Real maxVertDiff    = 0.0_Real;
   Real maxVertVisc    = 0.0_Real;
   int maxAbsB0Cell    = -1;
   int maxAbsUStarCell = -1;
   int maxVertDiffCell = -1;
   int maxVertDiffK    = -1;
   int maxVertViscCell = -1;
   int maxVertViscK    = -1;
   for (int C = 0; C < NCellsAll; ++C) {
      const Real b0c = B0H(C);
      const Real usc = UStarH(C);
      const Real ab0 = Kokkos::abs(b0c);
      const Real aus = Kokkos::abs(usc);
      if (ab0 > maxAbsB0) {
         maxAbsB0     = ab0;
         maxAbsB0Cell = C;
      }
      if (aus > maxAbsUStar) {
         maxAbsUStar     = aus;
         maxAbsUStarCell = C;
      }
      const int KCMin = MinLayerCellH(C);
      const int KCMax = MaxLayerCellH(C) + 1;
      for (int K = KCMin; K <= KCMax; ++K) {
         const Real diff = VertDiffH(C, K);
         const Real visc = VertViscH(C, K);
         if (diff > maxVertDiff) {
            maxVertDiff     = diff;
            maxVertDiffCell = C;
            maxVertDiffK    = K;
         }
         if (visc > maxVertVisc) {
            maxVertVisc     = visc;
            maxVertViscCell = C;
            maxVertViscK    = K;
         }
      }
   }
   LOG_WARN("KPP debug domain post-coeff: max|b0|={} at cell={} max|u*|={} "
            "at cell={} maxKPPDiff={} at cell={},k={} maxKPPVisc={} at "
            "cell={},k={}",
            maxAbsB0, maxAbsB0Cell, maxAbsUStar, maxAbsUStarCell, maxVertDiff,
            maxVertDiffCell, maxVertDiffK, maxVertVisc, maxVertViscCell,
            maxVertViscK);

   const int ICell       = 0;
   const int KMin        = MinLayerCellH(ICell);
   const int KMax        = MaxLayerCellH(ICell);
   const int NVertLayers = VCoord->NVertLayers;

   if (KMin > KMax) {
      return;
   }

   const int KSurf       = Kokkos::min(KMin, NVertLayers - 1);
   const Real rho_surf   = DensityH(ICell, KSurf);
   const Real u_star     = UStarH(ICell);
   const Real u_star_eff = Kokkos::fmax(KPP::MIN_USTAR, u_star);
   const Real b0         = B0H(ICell);
   Real u10              = 0.0_Real;
   if (WindSpeed10m.extent(0) > 0) {
      const auto Wind10mH = createHostMirrorCopy(WindSpeed10m);
      u10                 = Wind10mH(ICell);
   }
   const Real langmuir_factor =
       UseLangmuirCirculation ? ComputeEnhancementFactor(u10, u_star_eff, 50.0)
                              : 1.0_Real;
   const Real b0_eff = b0 * langmuir_factor;

   LOG_WARN("KPP debug: cell={} h_obl={} m k_obl={} u*={} b0={} b0_eff={} "
            "langmuir={}",
            ICell, OBLDepthH(ICell), OBLIndexH(ICell), u_star, b0, b0_eff,
            langmuir_factor);

   const int KOblIface = Kokkos::min(
       NVertLayers, Kokkos::max(KMin, static_cast<int>(OBLIndexH(ICell)) + 1));
   LOG_WARN("KPP debug coeff target: cell={} k_obl={} iface={} diff={} "
            "visc={}",
            ICell, OBLIndexH(ICell), KOblIface, VertDiffH(ICell, KOblIface),
            VertViscH(ICell, KOblIface));

   const int KTop   = Kokkos::min(KMax, KMin + 3);
   const int k_obl  = OBLIndexH(ICell);
   const Real h_obl = OBLDepthH(ICell);

   for (int K = KMin; K <= KTop; ++K) {
      const int kCell    = Kokkos::min(K, NVertLayers - 1);
      const int kInt     = Kokkos::min(K + 1, NVertLayers);
      const Real z_depth = SshCellH(ICell) - ZInterfaceH(ICell, kInt);

      const Real rho_k     = DensityH(ICell, kCell);
      const Real delta_rho = rho_k - rho_surf;
      const Real delta_b   = Gravity * delta_rho / RhoSw;
      const Real w_turb =
          ComputeTurbulentVelocityScale(u_star_eff, b0_eff, z_depth);
      const Real ri_b = delta_b * z_depth / (w_turb * w_turb + 1.0e-12_Real);

      Real sigma = 0.0_Real;
      if (K <= k_obl) {
         sigma = -1.0_Real * static_cast<Real>(K - KMin) /
                 static_cast<Real>(k_obl - KMin + 1);
         sigma = Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, sigma));
      }

      const Real z_local = -sigma * h_obl;
      Real zeta          = 0.0_Real;
      const Real denom   = VonKar * b0;
      if (Kokkos::abs(denom) > 1.0e-16_Real) {
         const Real l_mo = (u_star_eff * u_star_eff * u_star_eff) / denom;
         if (Kokkos::abs(l_mo) > 1.0e-16_Real) {
            zeta = z_local / l_mo;
         }
      }

      const Real phi_m = KPP::KPPProfileM2(zeta);
      const Real phi_s = KPP::KPPProfileS2(zeta);

      LOG_WARN(
          "KPP debug top: cell={} k={} z={} ri_b={} zeta={} phi_m={} phi_s={}",
          ICell, K, z_depth, ri_b, zeta, phi_m, phi_s);
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
   Array1DReal LangmuirFactor("LangmuirFactor", Mesh->NCellsAll);
   const bool LocUseLangmuirCirculation      = UseLangmuirCirculation;
   const Real LocSurfaceLayerExtent          = SurfaceLayerExtent;
   const Real LocCriticalRichardson          = CriticalRichardson;
   const Real LocIceFracThresholdForLangmuir = IceFractionThresholdForLangmuir;
   parallelFor(
       "KPP-Langmuir", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          const Real iceFrac = IceFraction(ICell);
          if (LocUseLangmuirCirculation &&
              iceFrac < LocIceFracThresholdForLangmuir) {
             const Real uStar = SurfaceFrictionVelocity(ICell);
             const Real u10 =
                 (WindSpeed10m.extent(0) > 0) ? WindSpeed10m(ICell) : 0.0_Real;
             LangmuirFactor(ICell) = ComputeEnhancementFactor(u10, uStar, 50.0);
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
   const Real LocIceFracThresholdForMinOBL = IceFractionThresholdForMinimumOBL;
   const Real LocMinimumOBLUnderSeaIce     = MinimumOBLUnderSeaIce;
   const Real LocStopOBLSearchMult         = StopOBLSearchMult;

   deepCopy(BulkRichardsonNumber, 0.0_Real);
   deepCopy(BulkRichardsonShear, 0.0_Real);
   deepCopy(UnresolvedShear, 0.0_Real);
   deepCopy(BuoyancyJump, 0.0_Real);

   // Maximum edges around a cell (matches MaxMaxEdges in HorzOperators.h)
   constexpr I4 MAX_EDGES_ON_CELL = 10;

   parallelFor(
       "KPP-OBLDepth", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          using namespace KPP;

          const Real u_star =
              Kokkos::fmax(0.0_Real, SurfaceFrictionVelocity(ICell));
          const Real b0 = SurfaceBuoyancyFlux(ICell);

          const I4 KMin     = MinLayerCell(ICell);
          const I4 KMax     = MaxLayerCell(ICell);
          const I4 KIntTop  = Kokkos::min(KMin + 1, NVertLayers);
          const I4 KIntDeep = Kokkos::min(KMax + 1, NVertLayers);

          const Real iceFrac = LocIceFraction(ICell);

          // KPP depths are measured below the free surface, so geometric
          // heights must be offset by the sea surface height.
          const Real Ssh = LocSshCell(ICell);

          Real obl_depth     = Ssh - ZInterface(ICell, KIntDeep);
          I4 k_cross         = -1;
          const Real ri_crit = LocCriticalRichardson;
          const Real ri_stop_crit =
              Kokkos::max(1.0e-6_Real, LocStopOBLSearchMult) * ri_crit;
          const Real ri_scaling = 1.0_Real - 0.5_Real * LocSurfaceLayerExtent;
          const Real b0_eff     = b0 * LocLangmuirFactor(ICell);

          // CVMix default unresolved-shear constants.
          const Real c_s_unres = 24.0_Real * Kokkos::sqrt(17.0_Real);
          const Real vtc =
              Kokkos::sqrt(0.2_Real /
                           Kokkos::max(1.0e-12_Real,
                                       c_s_unres * LocSurfaceLayerExtent)) /
              (VonKar * VonKar);

          // -------------------------------------------------------------------
          // Initialize per-edge running sums for surface-layer velocity
          // averages
          // -------------------------------------------------------------------
          const I4 nEdges = NEdgesOnCell(ICell);
          // MPAS-style area fractions for edge averaging.
          // Use edge kite area divided by cell area.
          const I4 nEdgesEff = Kokkos::min(nEdges, MAX_EDGES_ON_CELL);
          bool edge_valid[MAX_EDGES_ON_CELL]   = {};
          Real edge_weights[MAX_EDGES_ON_CELL] = {};
          const Real inv_area_cell =
              1.0_Real / Kokkos::max(AreaCell(ICell), 1.0e-20_Real);
          for (I4 J = 0; J < nEdgesEff; ++J) {
             const I4 IEdge = EdgesOnCell(ICell, J);
             const I4 KEMin = MinLayerEdgeBot(IEdge);
             const I4 KEMax = MaxLayerEdgeTop(IEdge);
             edge_valid[J] =
                 (KEMax >= KEMin && KEMin >= 0 && KEMin < NVertLayers);
             if (edge_valid[J]) {
                edge_weights[J] =
                    0.25_Real * DcEdge(IEdge) * DvEdge(IEdge) * inv_area_cell;
             }
          }
          if (nEdgesEff > 0) {
             Real sum_w = 0.0_Real;
             for (I4 J = 0; J < nEdgesEff; ++J) {
                if (edge_valid[J]) {
                   sum_w += edge_weights[J];
                }
             }
             if (sum_w < 1.0e-20_Real) {
                I4 n_edges_valid = 0;
                for (I4 J = 0; J < nEdgesEff; ++J) {
                   if (edge_valid[J]) {
                      ++n_edges_valid;
                   }
                }
                if (n_edges_valid > 0) {
                   const Real equal_w =
                       1.0_Real / static_cast<Real>(n_edges_valid);
                   for (I4 J = 0; J < nEdgesEff; ++J) {
                      edge_weights[J] = edge_valid[J] ? equal_w : 0.0_Real;
                   }
                }
             } else {
                const Real inv_sum_w = 1.0_Real / sum_w;
                for (I4 J = 0; J < nEdgesEff; ++J) {
                   if (edge_valid[J]) {
                      edge_weights[J] *= inv_sum_w;
                   }
                }
             }
          }

          // -------------------------------------------------------------------
          // Cell surface-layer running sums for density
          // MOVED INSIDE K-LOOP TO RESET EACH ITERATION (FIX FOR PROGRESSIVE
          // ACCUMULATION BUG)
          // -------------------------------------------------------------------

          for (I4 k = KMin; k <= KMax; ++k) {
             // Initialize fresh surface layer averages for this candidate OBL
             // depth
             I4 k_surface_avg     = KMin;
             const Real thick_top = Kokkos::abs(ZInterface(ICell, KMin + 1) -
                                                ZInterface(ICell, KMin));
             Real sum_thickness   = Kokkos::max(thick_top, 1.0e-12_Real);
             Real sum_rho = LocPotentialDensity(ICell, KMin) * sum_thickness;

             // Initialize fresh per-edge surface layer averages for this
             // candidate OBL depth
             I4 k_surf_e[MAX_EDGES_ON_CELL]      = {};
             Real sum_thick_e[MAX_EDGES_ON_CELL] = {};
             Real sum_un_e[MAX_EDGES_ON_CELL]    = {};
             Real sum_vt_e[MAX_EDGES_ON_CELL]    = {};

             for (I4 J = 0; J < nEdgesEff; ++J) {
                if (!edge_valid[J]) {
                   continue;
                }
                const I4 IEdge    = EdgesOnCell(ICell, J);
                const I4 KEMin    = MinLayerEdgeBot(IEdge);
                k_surf_e[J]       = KEMin;
                const I4 kInt0    = Kokkos::min(KEMin + 1, NVertLayers);
                const Real thick0 = Kokkos::abs(ZInterface(ICell, kInt0) -
                                                ZInterface(ICell, KEMin));
                sum_thick_e[J]    = Kokkos::max(thick0, 1.0e-12_Real);
                const I4 ke0      = Kokkos::min(KEMin, NVertLayers - 1);
                sum_un_e[J] = LocNormalVelocity(IEdge, ke0) * sum_thick_e[J];
                sum_vt_e[J] =
                    LocTangentialVelocity(IEdge, ke0) * sum_thick_e[J];
             }
             const I4 kCell      = Kokkos::min(k, NVertLayers - 1);
             const I4 kInt       = Kokkos::min(k + 1, NVertLayers);
             const Real z_depth  = Ssh - ZInterface(ICell, kInt);
             const Real z_center = Ssh - ZMid(ICell, kCell);
             if (z_depth < 1.0e-12)
                continue;

             const Real surf_layer_depth = LocSurfaceLayerExtent * z_depth;

             // Advance cell surface average for density
             while (k_surface_avg < k &&
                    (Ssh - ZInterface(ICell, k_surface_avg + 1)) <
                        surf_layer_depth) {
                ++k_surface_avg;
                const I4 ksa = Kokkos::min(k_surface_avg, NVertLayers - 1);
                const Real dk =
                    Kokkos::abs(ZInterface(ICell, k_surface_avg + 1) -
                                ZInterface(ICell, k_surface_avg));
                const Real thick_k = Kokkos::max(dk, 1.0e-12_Real);
                sum_thickness += thick_k;
                sum_rho += LocPotentialDensity(ICell, ksa) * thick_k;
             }

             // Advance per-edge surface averages for velocity
             for (I4 J = 0; J < nEdgesEff; ++J) {
                if (!edge_valid[J]) {
                   continue;
                }
                const I4 IEdge = EdgesOnCell(ICell, J);
                const I4 KEMax = MaxLayerEdgeTop(IEdge);
                while (k_surf_e[J] < k &&
                       (Ssh - ZInterface(ICell, k_surf_e[J] + 1)) <
                           surf_layer_depth) {
                   ++k_surf_e[J];
                   const I4 ke = Kokkos::min(
                       Kokkos::max(k_surf_e[J], MinLayerEdgeBot(IEdge)), KEMax);
                   const Real dk =
                       Kokkos::abs(ZInterface(ICell, k_surf_e[J] + 1) -
                                   ZInterface(ICell, k_surf_e[J]));
                   const Real thick_k = Kokkos::max(dk, 1.0e-12_Real);
                   sum_thick_e[J] += thick_k;
                   sum_un_e[J] += LocNormalVelocity(IEdge, ke) * thick_k;
                   sum_vt_e[J] += LocTangentialVelocity(IEdge, ke) * thick_k;
                }
             }

             const Real inv_sum_thickness =
                 1.0_Real / Kokkos::max(sum_thickness, 1.0e-12_Real);
             const Real rho_avg_surf = sum_rho * inv_sum_thickness;

             const Real rho_k             = LocPotentialDensity(ICell, kCell);
             const Real delta_rho         = rho_k - rho_avg_surf;
             const Real delta_b           = Gravity * delta_rho / RhoSw;
             LocBuoyancyJump(ICell, kInt) = delta_b;

             // Edge-based velocity shear: average deltaV^2 over cell edges
             Real deltaVsq = 0.0_Real;
             if (nEdges > 0) {
                for (I4 J = 0; J < nEdgesEff; ++J) {
                   if (!edge_valid[J]) {
                      continue;
                   }
                   const I4 IEdge = EdgesOnCell(ICell, J);
                   const I4 KEMin = MinLayerEdgeBot(IEdge);
                   const I4 KEMax = MaxLayerEdgeTop(IEdge);
                   const I4 k_e   = Kokkos::min(Kokkos::max(k, KEMin), KEMax);
                   const Real inv_thick_e =
                       1.0_Real / Kokkos::max(sum_thick_e[J], 1.0e-12_Real);
                   const Real un_avg = sum_un_e[J] * inv_thick_e;
                   const Real vt_avg = sum_vt_e[J] * inv_thick_e;
                   const Real un_k   = LocNormalVelocity(IEdge, k_e);
                   const Real vt_k   = LocTangentialVelocity(IEdge, k_e);
                   const Real dun    = un_k - un_avg;
                   const Real dvt    = vt_k - vt_avg;
                   deltaVsq += edge_weights[J] * (dun * dun + dvt * dvt);
                }
             }
             LocBulkRichardsonShear(ICell, kInt) =
                 Kokkos::max(deltaVsq, 1.0e-15_Real);

             const Real sigma_loc = Kokkos::fmin(
                 1.0_Real, Kokkos::fmax(0.0_Real, LocSurfaceLayerExtent));

             Real w_turb = 0.0_Real;
             if (u_star > 1.0e-12_Real) {
                const Real u3   = u_star * u_star * u_star;
                const Real zeta = sigma_loc * z_depth * VonKar * b0_eff /
                                  Kokkos::max(u3, 1.0e-20_Real);
                const Real phi_inv_s = KPP::KPPProfileS2(zeta);
                w_turb = VonKar * u_star * Kokkos::max(phi_inv_s, 0.0_Real);
             } else if (b0_eff < 0.0_Real) {
                const Real c_s = KPP::C_MO_S;
                const Real ws3 = -c_s * sigma_loc * z_depth * VonKar * b0_eff;
                w_turb = VonKar * Kokkos::pow(Kokkos::max(ws3, 0.0_Real),
                                              1.0_Real / 3.0_Real);
             }
             const Real n_cntr = Kokkos::sqrt(
                 Kokkos::max(0.0_Real, LocBruntVaisalaFreqSq(ICell, kInt)));
             const Real cv  = (n_cntr < 0.002_Real)
                                  ? (2.1_Real - 200.0_Real * n_cntr)
                                  : 1.7_Real;
             const Real vt2 = Kokkos::max(
                 1.0e-10_Real, cv * vtc * z_center * n_cntr * w_turb /
                                   Kokkos::max(ri_crit, 1.0e-12_Real));
             LocUnresolvedShear(ICell, kInt) = vt2;

             const Real vel_scale2 = deltaVsq + vt2;

             const Real ri_b = ri_scaling * delta_b * z_center /
                               Kokkos::max(vel_scale2, 1.0e-12_Real);
             LocBulkRichardson(ICell, kInt) = ri_b;

             if (k_cross < 0 && ri_b > ri_stop_crit) {
                k_cross = k;
             }
          }

          if (k_cross >= KMin) {
             if (k_cross > KMin) {
                // Ri values are defined at cell centers, so interpolate on
                // center depths to keep the abscissa consistent.
                const I4 kAbove     = Kokkos::max(KMin, k_cross - 1);
                const I4 kBelow     = Kokkos::min(k_cross, NVertLayers - 1);
                const I4 kAboveRi   = Kokkos::min(kAbove + 1, NVertLayers);
                const I4 kBelowRi   = Kokkos::min(kBelow + 1, NVertLayers);
                const Real z_above  = Ssh - ZMid(ICell, kAbove);
                const Real z_below  = Ssh - ZMid(ICell, kBelow);
                const Real ri_above = LocBulkRichardson(ICell, kAboveRi);
                const Real ri_below = LocBulkRichardson(ICell, kBelowRi);

                const Real h = z_below - z_above;
                if (h > 1.0e-12_Real) {
                   // CVMix-style QUAD interpolation for OBL crossing:
                   // - first interior crossing uses zero slope at top point
                   // - deeper crossings use upstream slope
                   Real slope_above = 0.0_Real;
                   if (k_cross > KMin + 1) {
                      const I4 kPrev     = Kokkos::max(KMin, kAbove - 1);
                      const I4 kPrevRi   = Kokkos::min(kPrev + 1, NVertLayers);
                      const Real z_prev  = Ssh - ZMid(ICell, kPrev);
                      const Real ri_prev = LocBulkRichardson(ICell, kPrevRi);
                      const Real dz_prev = z_above - z_prev;
                      if (Kokkos::abs(dz_prev) > 1.0e-12_Real) {
                         slope_above = (ri_above - ri_prev) / dz_prev;
                      }
                   }

                   // In local coordinate t = z - z_above:
                   // Ri(t) = A t^2 + slope_above t + ri_above
                   const Real A =
                       (ri_below - ri_above - slope_above * h) / (h * h);
                   const Real C = ri_above - ri_stop_crit;

                   Real t_cross = h;
                   if (Kokkos::abs(A) < 1.0e-14_Real) {
                      // Degenerate quadratic -> linear fallback.
                      const Real d_ri = ri_below - ri_above;
                      if (Kokkos::abs(d_ri) > 1.0e-12_Real) {
                         const Real frac = Kokkos::fmax(
                             0.0_Real,
                             Kokkos::fmin(1.0_Real,
                                          (ri_stop_crit - ri_above) / d_ri));
                         t_cross = frac * h;
                      }
                   } else {
                      const Real disc =
                          slope_above * slope_above - 4.0_Real * A * C;
                      if (disc >= 0.0_Real) {
                         const Real sqrt_disc = Kokkos::sqrt(disc);
                         const Real t1 =
                             (-slope_above + sqrt_disc) / (2.0_Real * A);
                         const Real t2 =
                             (-slope_above - sqrt_disc) / (2.0_Real * A);

                         const bool t1_ok = (t1 >= 0.0_Real && t1 <= h);
                         const bool t2_ok = (t2 >= 0.0_Real && t2 <= h);
                         if (t1_ok && t2_ok) {
                            const Real mid = 0.5_Real * h;
                            t_cross =
                                (Kokkos::abs(t1 - mid) <= Kokkos::abs(t2 - mid))
                                    ? t1
                                    : t2;
                         } else if (t1_ok) {
                            t_cross = t1;
                         } else if (t2_ok) {
                            t_cross = t2;
                         } else {
                            t_cross = h;
                         }
                      }
                   }

                   t_cross   = Kokkos::fmax(0.0_Real, Kokkos::fmin(h, t_cross));
                   obl_depth = z_above + t_cross;
                } else {
                   obl_depth = z_below;
                }
             } else {
                // Match center-based OBL convention when crossing occurs in
                // the top interval.
                obl_depth = Ssh - ZMid(ICell, KMin);
             }
          } else {
             obl_depth = Ssh - ZInterface(ICell, KIntDeep);
          }

          const Real top_layer_thickness =
              Kokkos::abs(ZInterface(ICell, KIntTop) - ZInterface(ICell, KMin));
          const Real min_obl_depth = 0.5_Real * top_layer_thickness;
          const Real max_obl_depth = Ssh - ZMid(ICell, KMax);
          obl_depth                = Kokkos::fmax(obl_depth, min_obl_depth);
          if (iceFrac > LocIceFracThresholdForMinOBL) {
             obl_depth = Kokkos::fmax(obl_depth, LocMinimumOBLUnderSeaIce);
          }
          obl_depth = Kokkos::fmin(obl_depth, max_obl_depth);

          I4 k_final = KMax;
          for (I4 k = KMin; k < KMax; ++k) {
             const Real z_above = Ssh - ZInterface(ICell, k);
             const Real z_below = Ssh - ZInterface(ICell, k + 1);
             if (obl_depth >= z_above && obl_depth <= z_below) {
                k_final = k;
                break;
             }
          }

          LocBoundaryLayerDepth(ICell)      = obl_depth;
          LocIndexBoundaryLayerDepth(ICell) = k_final;
       });

   if (LocUseBLDSmoothing) {
      Array1DReal BoundaryLayerDepthSmooth("BoundaryLayerDepthSmooth",
                                           Mesh->NCellsAll);
      OMEGA_SCOPE(LocBoundaryLayerDepthSmooth, BoundaryLayerDepthSmooth);
      OMEGA_SCOPE(LocNCellsAll, Mesh->NCellsAll);

      parallelFor(
          "KPP-OBLDepth-Smooth", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
             const I4 KMin = MinLayerCell(ICell);
             if (KMin < 0 || KMin >= NVertLayers) {
                LocBoundaryLayerDepthSmooth(ICell) =
                    LocBoundaryLayerDepth(ICell);
                return;
             }

             const I4 nEdges = NEdgesOnCell(ICell);
             Real area_sum   = 0.0_Real;
             Real bld_sum    = 0.0_Real;
             I4 edge_count   = 0;

             for (I4 J = 0; J < nEdges; ++J) {
                const I4 INeighbor = CellsOnCell(ICell, J);
                if (INeighbor == LocNCellsAll) {
                   continue;
                }

                const I4 KMinNbr = MinLayerCell(INeighbor);
                if (KMinNbr < 0 || KMinNbr >= NVertLayers) {
                   continue;
                }

                const Real nbr_area = AreaCell(INeighbor);
                bld_sum +=
                    2.0_Real * nbr_area * LocBoundaryLayerDepth(INeighbor);
                area_sum += 2.0_Real * nbr_area;
                ++edge_count;
             }

             if (edge_count > 0) {
                const Real self_area = AreaCell(ICell);
                bld_sum += LocBoundaryLayerDepth(ICell) *
                           static_cast<Real>(edge_count) * self_area;
                area_sum += static_cast<Real>(edge_count) * self_area;
             }

             if (area_sum > 0.0_Real) {
                LocBoundaryLayerDepthSmooth(ICell) = bld_sum / area_sum;
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

             const I4 KIntTop = Kokkos::min(KMin + 1, NVertLayers);
             const Real top_layer_thickness = Kokkos::abs(
                 ZInterface(ICell, KIntTop) - ZInterface(ICell, KMin));
             const Real min_obl_depth = 0.5_Real * top_layer_thickness;
             const Real max_obl_depth = Ssh - ZMid(ICell, KMax);

             Real obl_depth = LocBoundaryLayerDepthSmooth(ICell);
             obl_depth      = Kokkos::fmax(obl_depth, min_obl_depth);
             obl_depth      = Kokkos::fmin(obl_depth, max_obl_depth);

             I4 k_final = KMax;
             for (I4 k = KMin; k < KMax; ++k) {
                const Real z_above = Ssh - ZInterface(ICell, k);
                const Real z_below = Ssh - ZInterface(ICell, k + 1);
                if (obl_depth >= z_above && obl_depth <= z_below) {
                   k_final = k;
                   break;
                }
             }

             LocBoundaryLayerDepth(ICell)      = obl_depth;
             LocIndexBoundaryLayerDepth(ICell) = k_final;
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
   bool LocUseNonLocalFlux          = UseNonLocalFlux;
   const Real LocSurfaceLayerExtent = SurfaceLayerExtent;
   I4 LocMatchTechnique = 0; // 0=SimpleShapes, 1=MatchBoth, 2=ParabolicNonLocal
   if (MatchTechniqueStr == "MatchBoth") {
      LocMatchTechnique = 1;
   } else if (MatchTechniqueStr == "ParabolicNonLocal") {
      LocMatchTechnique = 2;
   }
   // Non-local flux normalization constant from Large et al. (1994) eq. 20:
   // C_s = C* * kappa * (c_s * kappa * epsilon)^(1/3)
   // where C* = 10, c_s = C_MO_S = 98.9545, kappa = VonKar, epsilon =
   // SurfaceLayerExtent
   const Real LocNonLocalCs =
       10.0_Real * VonKar *
       Kokkos::pow(KPP::C_MO_S * VonKar * LocSurfaceLayerExtent,
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
          Real h_obl = LocBoundaryLayerDepth(ICell);

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
          Real u_star = LocSurfaceFrictionVelocity(ICell);
          Real b0     = LocSurfaceBuoyancyFlux(ICell);

          // =============================================================
          // Compute mixing coefficients at each interface
          // =============================================================
          for (I4 k = KMin; k <= KMax + 1; ++k) {
             const I4 k_iface   = Kokkos::min(Kokkos::max(k, 0), NVertLayers);
             const Real z_depth = Ssh - ZInterface(ICell, k_iface);

             // Check if within OBL using depth below the free surface.
             if (z_depth <= h_obl && h_obl > 0.0_Real) {
                // Normalized depth in Omega sign convention: sigma in [-1,0].
                Real sigma = -z_depth / h_obl;
                sigma = Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, sigma));

                // CVMix-style turbulent scales: w = kappa*u*/phi in general,
                // with explicit free-convection limits when u*=0.
                const Real sigma_coord = -sigma; // [0,1]
                const Real sigma_loc   = Kokkos::fmin(
                    LocSurfaceLayerExtent, Kokkos::fmax(0.0_Real, sigma_coord));

                Real zeta     = 0.0_Real;
                Real w_m_turb = 0.0_Real;
                Real w_s_turb = 0.0_Real;

                if (u_star > 0.0_Real) {
                   const Real u3 = u_star * u_star * u_star;
                   zeta          = sigma_loc * h_obl * b0 * LocKappa /
                          Kokkos::max(u3, 1.0e-20_Real);

                   // KPPProfileM2/S2 return phi^{-1}; do not invert again.
                   const Real phi_inv_m = KPP::KPPProfileM2(zeta);
                   const Real phi_inv_s = KPP::KPPProfileS2(zeta);

                   w_m_turb =
                       LocKappa * u_star * Kokkos::max(phi_inv_m, 0.0_Real);
                   w_s_turb =
                       LocKappa * u_star * Kokkos::max(phi_inv_s, 0.0_Real);
                } else if (b0 < 0.0_Real) {
                   // Free-convection edge case (u*=0, unstable forcing).
                   const Real c_m = KPP::C_MO_M;
                   const Real c_s = KPP::C_MO_S;
                   const Real wm3 = -c_m * sigma_loc * h_obl * LocKappa * b0;
                   const Real ws3 = -c_s * sigma_loc * h_obl * LocKappa * b0;
                   w_m_turb = LocKappa * Kokkos::pow(Kokkos::max(0.0_Real, wm3),
                                                     1.0_Real / 3.0_Real);
                   w_s_turb = LocKappa * Kokkos::pow(Kokkos::max(0.0_Real, ws3),
                                                     1.0_Real / 3.0_Real);
                }

                const Real match_visc_shape =
                    (LocUseInteriorMix && LocMatchTechnique == 1 &&
                     h_obl > 0.0_Real && w_m_turb > 0.0_Real)
                        ? LocInteriorVertVisc(ICell, KMatch) /
                              Kokkos::max(h_obl * w_m_turb, 1.0e-20_Real)
                        : 0.0_Real;
                const Real match_diff_shape =
                    (LocUseInteriorMix && LocMatchTechnique == 1 &&
                     h_obl > 0.0_Real && w_s_turb > 0.0_Real)
                        ? LocInteriorVertDiff(ICell, KMatch) /
                              Kokkos::max(h_obl * w_s_turb, 1.0e-20_Real)
                        : 0.0_Real;

                // ========================================================
                // Momentum mixing contribution.
                // ========================================================
                Real m1 = (LocUseInteriorMix && LocMatchTechnique == 1)
                              ? KPP::KPPProfileMatched(sigma, match_visc_shape)
                              : KPP::KPPProfileM1(sigma);
                LocVertVisc(ICell, k) = h_obl * w_m_turb * m1;

                // ========================================================
                // Tracer mixing contribution.
                // ========================================================
                Real s1 = (LocUseInteriorMix && LocMatchTechnique == 1)
                              ? KPP::KPPProfileMatched(sigma, match_diff_shape)
                              : KPP::KPPProfileS1(sigma);
                LocVertDiff(ICell, k)               = h_obl * w_s_turb * s1;
                LocTurbulentVelocityScale(ICell, k) = w_s_turb;

                // ========================================================
                // Non-local flux: C_s * G(σ)
                // C_s = C* * kappa * (c_s * kappa * epsilon)^(1/3)
                // per Large et al. (1994) eq. 20 (~6.33 with default constants)
                // ========================================================
                // Match CVMix behavior: apply non-local term only when
                // surface buoyancy forcing is unstable/neutral.
                if (LocUseNonLocalFlux && b0 <= 0.0_Real) {
                   Real g_sigma = 0.0_Real;
                   if (LocMatchTechnique == 2) {
                      g_sigma = KPP::KPPProfileGParabolicNonLocal(sigma);
                   } else if (LocMatchTechnique == 1) {
                      g_sigma = KPP::KPPProfileGMatchBoth(sigma);
                   } else {
                      g_sigma = KPP::KPPProfileG(sigma);
                   }
                   LocVertNonLocalFlux(ICell, k) = LocNonLocalCs * g_sigma;
                } else {
                   LocVertNonLocalFlux(ICell, k) = 0.0;
                }

             } else {
                // Below OBL: preserve interior values for MatchBoth, otherwise
                // no KPP contribution.
                LocVertDiff(ICell, k)               = LocUseInteriorMix
                                                          ? LocInteriorVertDiff(ICell, k)
                                                          : 0.0_Real;
                LocVertVisc(ICell, k)               = LocUseInteriorMix
                                                          ? LocInteriorVertVisc(ICell, k)
                                                          : 0.0_Real;
                LocVertNonLocalFlux(ICell, k)       = 0.0;
                LocTurbulentVelocityScale(ICell, k) = 0.0;
             }
          }

          // Optional enhanced diffusion/viscosity treatment at OBL base.
          // Match CVMix Appendix D weighting at the interface nearest h_obl.
          if (LocUseEnhancedDiffusion && h_obl > 0.0_Real) {
             const I4 k_obl = Kokkos::max(
                 KMin, Kokkos::min(LocIndexBoundaryLayerDepth(ICell), KMax));
             const Real z_mid_obl = Ssh - ZMid(ICell, k_obl);

             const bool target_outside_obl = h_obl >= z_mid_obl;
             const I4 k_ktup =
                 target_outside_obl ? k_obl : Kokkos::max(KMin, k_obl - 1);
             const I4 k_target = target_outside_obl
                                     ? Kokkos::min(k_obl + 1, KMax + 1)
                                     : Kokkos::max(KMin + 1, k_obl);

             const Real z_ktup = Ssh - ZMid(ICell, k_ktup);
             const Real z_next = (k_ktup < KMax)
                                     ? (Ssh - ZMid(ICell, k_ktup + 1))
                                     : (Ssh - ZInterface(ICell, k_ktup + 1));
             const Real delta  = Kokkos::fmax(
                 0.0_Real,
                 Kokkos::fmin(1.0_Real,
                               (h_obl - z_ktup) /
                                   Kokkos::max(z_next - z_ktup, 1.0e-12_Real)));
             const Real one_minus_delta = 1.0_Real - delta;

             Real sigma_ktup = -z_ktup / h_obl;
             sigma_ktup =
                 Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, sigma_ktup));
             const Real sigma_coord = -sigma_ktup;
             const Real sigma_loc   = Kokkos::fmin(
                 LocSurfaceLayerExtent, Kokkos::fmax(0.0_Real, sigma_coord));

             Real w_m_ktup = 0.0_Real;
             Real w_s_ktup = 0.0_Real;
             if (u_star > 0.0_Real) {
                const Real u3   = u_star * u_star * u_star;
                const Real zeta = sigma_loc * h_obl * b0 * LocKappa /
                                  Kokkos::max(u3, 1.0e-20_Real);
                w_m_ktup = LocKappa * u_star *
                           Kokkos::max(KPP::KPPProfileM2(zeta), 0.0_Real);
                w_s_ktup = LocKappa * u_star *
                           Kokkos::max(KPP::KPPProfileS2(zeta), 0.0_Real);
             } else if (b0 < 0.0_Real) {
                const Real wm3 =
                    -KPP::C_MO_M * sigma_loc * h_obl * LocKappa * b0;
                const Real ws3 =
                    -KPP::C_MO_S * sigma_loc * h_obl * LocKappa * b0;
                w_m_ktup = LocKappa * Kokkos::pow(Kokkos::max(0.0_Real, wm3),
                                                  1.0_Real / 3.0_Real);
                w_s_ktup = LocKappa * Kokkos::pow(Kokkos::max(0.0_Real, ws3),
                                                  1.0_Real / 3.0_Real);
             }

             const Real match_visc_shape =
                 (LocUseInteriorMix && LocMatchTechnique == 1 &&
                  h_obl > 0.0_Real && w_m_ktup > 0.0_Real)
                     ? LocInteriorVertVisc(ICell, KMatch) /
                           Kokkos::max(h_obl * w_m_ktup, 1.0e-20_Real)
                     : 0.0_Real;
             const Real match_diff_shape =
                 (LocUseInteriorMix && LocMatchTechnique == 1 &&
                  h_obl > 0.0_Real && w_s_ktup > 0.0_Real)
                     ? LocInteriorVertDiff(ICell, KMatch) /
                           Kokkos::max(h_obl * w_s_ktup, 1.0e-20_Real)
                     : 0.0_Real;

             const Real visc_ktup =
                 h_obl * w_m_ktup *
                 ((LocUseInteriorMix && LocMatchTechnique == 1)
                      ? KPP::KPPProfileMatched(sigma_ktup, match_visc_shape)
                      : KPP::KPPProfileM1(sigma_ktup));
             const Real diff_ktup =
                 h_obl * w_s_ktup *
                 ((LocUseInteriorMix && LocMatchTechnique == 1)
                      ? KPP::KPPProfileMatched(sigma_ktup, match_diff_shape)
                      : KPP::KPPProfileS1(sigma_ktup));

             const Real visc_profile = LocVertVisc(ICell, k_target);
             const Real diff_profile = LocVertDiff(ICell, k_target);

             const Real enh_visc =
                 one_minus_delta * one_minus_delta * visc_ktup +
                 delta * delta * visc_profile;
             const Real enh_diff =
                 one_minus_delta * one_minus_delta * diff_ktup +
                 delta * delta * diff_profile;

             const Real old_visc = LocUseInteriorMix
                                       ? LocInteriorVertVisc(ICell, k_target)
                                       : 0.0_Real;
             const Real old_diff = LocUseInteriorMix
                                       ? LocInteriorVertDiff(ICell, k_target)
                                       : 0.0_Real;
             const Real new_visc =
                 one_minus_delta * old_visc + delta * enh_visc;
             const Real new_diff =
                 one_minus_delta * old_diff + delta * enh_diff;

             LocVertVisc(ICell, k_target) = new_visc;
             LocVertDiff(ICell, k_target) = new_diff;

             if (!target_outside_obl && diff_profile != 0.0_Real) {
                LocVertNonLocalFlux(ICell, k_target) =
                    LocVertNonLocalFlux(ICell, k_target) * new_diff /
                    diff_profile;
             } else if (!target_outside_obl) {
                LocVertNonLocalFlux(ICell, k_target) = 0.0_Real;
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
