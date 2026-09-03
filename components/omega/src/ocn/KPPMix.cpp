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

} // namespace

/// Constructor for KPPMix
KPPMix::KPPMix(const std::string &InName, const HorzMesh *InMesh,
               const VertCoord *InVCoord)
    : Name(InName), Mesh(InMesh), VCoord(InVCoord) {

   // Allocate output arrays
   VertDiff  = Array2DReal("VertDiff", Mesh->NCellsSize, VCoord->NVertLayersP1);
   VertVisc  = Array2DReal("VertVisc", Mesh->NCellsSize, VCoord->NVertLayersP1);
   OSBLDepth = Array1DReal("OSBLDepth", Mesh->NCellsSize);
   OSBLDepthIndex = Array1DI4("OSBLDepthIndex", Mesh->NCellsSize);
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
   OSBLDepthFldName           = "OSBLDepth";
   OSBLDepthIndexFldName      = "OSBLDepthIndex";
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
      OSBLDepthFldName.append(Name);
      OSBLDepthIndexFldName.append(Name);
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
       KPPConfig.get("UseOSBLSmoothing", DefKPPMix->UseOSBLSmoothing);
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
   Err += KPPConfig.get("IceFractionThresholdForMinimumOSBL",
                        DefKPPMix->IceFractionThresholdForMinimumOSBL);
   Err += KPPConfig.get("MinimumOSBLUnderSeaIce",
                        DefKPPMix->MinimumOSBLUnderSeaIce);
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
   // Stage 1: Compute OSBL depth
   // =======================================================================
   computeOSBLDepth(PotentialDensity, NormalVelocity, TangentialVelocity,
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
   const auto OSBLDepthH    = createHostMirrorCopy(OSBLDepth);
   const auto OSBLIndexH    = createHostMirrorCopy(OSBLDepthIndex);
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

   const int KOblIface = Kokkos::min(
       NVertLayers, Kokkos::max(KMin, static_cast<int>(OBLIndexH(ICell)) + 1));

   const int KTop   = Kokkos::min(KMax, KMin + 3);
   const int KOSBL  = OSBLIndexH(ICell);
   const Real HOSBL = OSBLDepthH(ICell);

   for (int K = KMin; K <= KTop; ++K) {
      const int KCell   = Kokkos::min(K, NVertLayers - 1);
      const int KIface  = Kokkos::min(K + 1, NVertLayers);
      const Real ZDepth = SshCellH(ICell) - ZInterfaceH(ICell, KIface);

      const Real RhoK     = DensityH(ICell, KCell);
      const Real DeltaRho = RhoK - RhoSurf;
      const Real DeltaB   = Gravity * DeltaRho / RhoSw;
      const Real WTurb =
          computeTurbVelocityScale(UStarEff, BuoyFluxEff, ZDepth);
      const Real RiBulk =
          DeltaB * ZDepth / (WTurb * WTurb + KPP::NumericalTolerance);

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
   }
}

/// Stage 1: Compute OSBL depth using bulk Richardson search with edge-based
/// velocity shear following the MPAS CVMix reference implementation
/// (mpas_ocn_vmix_cvmix.F)
void KPPMix::computeOSBLDepth(const Array2DReal &PotentialDensity,
                              const Array2DReal &NormalVelocity,
                              const Array2DReal &TangentialVelocity,
                              const Array1DReal &SurfaceFrictionVelocity,
                              const Array1DReal &SurfaceBuoyancyFlux,
                              const Array2DReal &BruntVaisalaFreqSq,
                              const Array1DReal &IceFraction,
                              const Array1DReal &WindSpeed10m) {

   // =======================================================================
   // Compute Langmuir enhancement factors if wind speed is available
   // =======================================================================
   Array1DReal LangmuirFactor("LangmuirFactor", Mesh->NCellsSize);

   KPPLangmuirFactor LangmuirCalc;
   LangmuirCalc.UseLangmuirCirculation      = UseLangmuirCirculation;
   LangmuirCalc.IceFracThresholdForLangmuir = IceFractionThresholdForLangmuir;

   OMEGA_SCOPE(LocLangmuirFactor, LangmuirFactor);
   OMEGA_SCOPE(LocIceFraction, IceFraction);
   OMEGA_SCOPE(LocSurfFricVel, SurfaceFrictionVelocity);
   OMEGA_SCOPE(LocSurfBuoyFlux, SurfaceBuoyancyFlux);
   OMEGA_SCOPE(LocWindSpeed10m, WindSpeed10m);

   parallelFor(
       "KPP-Langmuir", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          LangmuirCalc(LocLangmuirFactor, ICell, LocIceFraction, LocSurfFricVel,
                       LocWindSpeed10m);
       });

   // =======================================================================
   // Stage 1: Compute OSBL depth using edge-based velocity shear
   // =======================================================================
   KPPOSBLDepthSearch OSBLDepthCalc(Mesh, VCoord);
   OSBLDepthCalc.CriticalRichardson = CriticalRichardson;
   OSBLDepthCalc.SurfaceLayerExtent = SurfaceLayerExtent;
   OSBLDepthCalc.IceFracThresholdForMinimumOSBL =
       IceFractionThresholdForMinimumOSBL;
   OSBLDepthCalc.MinimumOSBLUnderSeaIce = MinimumOSBLUnderSeaIce;
   OSBLDepthCalc.FullRiProfile          = DebugDiagnostics;

   OMEGA_SCOPE(LocPotentialDensity, PotentialDensity);
   OMEGA_SCOPE(LocNormalVelocity, NormalVelocity);
   OMEGA_SCOPE(LocTangentialVelocity, TangentialVelocity);
   OMEGA_SCOPE(LocBruntVaisalaFreqSq, BruntVaisalaFreqSq);
   OMEGA_SCOPE(LocOSBLDepth, OSBLDepth);
   OMEGA_SCOPE(LocOSBLDepthIndex, OSBLDepthIndex);
   OMEGA_SCOPE(LocBulkRichardson, BulkRichardsonNumber);
   OMEGA_SCOPE(LocBulkRichardsonShear, BulkRichardsonShear);
   OMEGA_SCOPE(LocUnresolvedShear, UnresolvedShear);
   OMEGA_SCOPE(LocBuoyancyJump, BuoyancyJump);

   deepCopy(BulkRichardsonNumber, 0.0_Real);
   deepCopy(BulkRichardsonShear, 0.0_Real);
   deepCopy(UnresolvedShear, 0.0_Real);
   deepCopy(BuoyancyJump, 0.0_Real);

   parallelFor(
       "KPP-OSBLDepth", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          OSBLDepthCalc(LocOSBLDepth, LocOSBLDepthIndex, LocBulkRichardson,
                        LocBulkRichardsonShear, LocUnresolvedShear,
                        LocBuoyancyJump, ICell, LocPotentialDensity,
                        LocNormalVelocity, LocTangentialVelocity,
                        LocSurfFricVel, LocSurfBuoyFlux, LocBruntVaisalaFreqSq,
                        LocIceFraction, LocLangmuirFactor);
       });

   if (UseOSBLSmoothing) {
      Array1DReal OSBLDepthSmooth("OSBLDepthSmooth", Mesh->NCellsSize);
      OMEGA_SCOPE(LocOSBLDepthSmooth, OSBLDepthSmooth);

      KPPOSBLSmooth BLDSmoothCalc(Mesh, VCoord);
      KPPOSBLCommit BLDCommitCalc(Mesh, VCoord);

      parallelFor(
          "KPP-OSBLDepth-Smooth", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
             BLDSmoothCalc(LocOSBLDepthSmooth, ICell, LocOSBLDepth);
          });

      parallelFor(
          "KPP-OSBLDepth-CommitSmooth", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell) {
             BLDCommitCalc(LocOSBLDepth, LocOSBLDepthIndex, ICell,
                           LocOSBLDepthSmooth);
          });
   }
}

/// Stage 2: Compute KPP mixing contribution or matched coefficients
void KPPMix::computeMixingCoefficients(
    const Array2DReal &PotentialDensity,
    const Array1DReal &SurfaceFrictionVelocity,
    const Array1DReal &SurfaceBuoyancyFlux, const Array2DReal &InteriorVertDiff,
    const Array2DReal &InteriorVertVisc) {

   (void)PotentialDensity;

   I4 NVertLayers = VCoord->NVertLayers;

   const bool LocUseInteriorMix =
       InteriorVertDiff.data() != nullptr && InteriorVertVisc.data() != nullptr;
   const bool LocUseMatchedShapes =
       LocUseInteriorMix && MatchTechnique == KPPMatchType::MatchBoth;

   OMEGA_SCOPE(LocOSBLDepth, OSBLDepth);
   OMEGA_SCOPE(LocOSBLDepthIndex, OSBLDepthIndex);
   OMEGA_SCOPE(LocVertDiff, VertDiff);
   OMEGA_SCOPE(LocVertVisc, VertVisc);
   OMEGA_SCOPE(LocVertNonLocalFlux, VertNonLocalFlux);
   OMEGA_SCOPE(LocTurbulentVelocityScale, TurbulentVelocityScale);
   OMEGA_SCOPE(LocSurfFricVel, SurfaceFrictionVelocity);
   OMEGA_SCOPE(LocSurfBuoyFlux, SurfaceBuoyancyFlux);
   OMEGA_SCOPE(LocInteriorVertDiff, InteriorVertDiff);
   OMEGA_SCOPE(LocInteriorVertVisc, InteriorVertVisc);

   // =======================================================================
   // Initialize with zero KPP contribution, or precomputed interior mixing for
   // matched-coefficient construction.
   // =======================================================================
   KPPCoeffsInit CoeffsInit;
   CoeffsInit.UseInteriorMix = LocUseInteriorMix;

   parallelFor(
       "KPP-Coeffs-Init", {Mesh->NCellsAll, NVertLayers + 1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          CoeffsInit(LocVertDiff, LocVertVisc, LocVertNonLocalFlux,
                     LocTurbulentVelocityScale, ICell, K, LocInteriorVertDiff,
                     LocInteriorVertVisc);
       });

   // =======================================================================
   // Stage 2: Compute KPP profile-based mixing coefficients
   // =======================================================================
   KPPMixingCoeffs MixingCoeffsCalc(Mesh, VCoord);
   MixingCoeffsCalc.SurfaceLayerExtent   = SurfaceLayerExtent;
   MixingCoeffsCalc.UseEnhancedDiffusion = UseEnhancedDiffusion;
   MixingCoeffsCalc.UseInteriorMix       = LocUseInteriorMix;
   MixingCoeffsCalc.UseMatchedShapes     = LocUseMatchedShapes;

   parallelFor(
       "KPP-MixingCoeffs", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          MixingCoeffsCalc(LocVertDiff, LocVertVisc, LocVertNonLocalFlux,
                           LocTurbulentVelocityScale, ICell, LocOSBLDepth,
                           LocOSBLDepthIndex, LocSurfFricVel, LocSurfBuoyFlux,
                           LocInteriorVertDiff, LocInteriorVertVisc);
       });
}

/// Constructor for KPPOSBLDepthSearch
KPPOSBLDepthSearch::KPPOSBLDepthSearch(const HorzMesh *Mesh,
                                       const VertCoord *VCoord)
    : NVertLayers(VCoord->NVertLayers), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop),
      ZInterface(VCoord->GeomZInterface), ZMid(VCoord->GeomZMid),
      SshCell(VCoord->SshCell), NEdgesOnCell(Mesh->NEdgesOnCell),
      EdgesOnCell(Mesh->EdgesOnCell), CellsOnCell(Mesh->CellsOnCell),
      AreaCell(Mesh->AreaCell), DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge) {}

/// Constructor for KPPOSBLSmooth
KPPOSBLSmooth::KPPOSBLSmooth(const HorzMesh *Mesh, const VertCoord *VCoord)
    : NVertLayers(VCoord->NVertLayers), NCellsAll(Mesh->NCellsAll),
      MinLayerCell(VCoord->MinLayerCell), NEdgesOnCell(Mesh->NEdgesOnCell),
      CellsOnCell(Mesh->CellsOnCell), AreaCell(Mesh->AreaCell) {}

/// Constructor for KPPOSBLCommit
KPPOSBLCommit::KPPOSBLCommit(const HorzMesh *Mesh, const VertCoord *VCoord)
    : NVertLayers(VCoord->NVertLayers), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell), ZInterface(VCoord->GeomZInterface),
      ZMid(VCoord->GeomZMid), SshCell(VCoord->SshCell) {}

/// Constructor for KPPMixingCoeffs
KPPMixingCoeffs::KPPMixingCoeffs(const HorzMesh *Mesh, const VertCoord *VCoord)
    : NVertLayers(VCoord->NVertLayers), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell), ZInterface(VCoord->GeomZInterface),
      ZMid(VCoord->GeomZMid), SshCell(VCoord->SshCell) {}

/// Register fields with I/O system
void KPPMix::defineFields() {
   // OSBLDepth on cells
   std::vector<std::string> CellDims(1);
   CellDims[0] = "NCells";
   auto OSBLDepthField =
       Field::create(OSBLDepthFldName,                     // field name
                     "ocean surface boundary layer depth", // long name
                     "m",                                  // units
                     "",                                   // CF standard name
                     0.0,                                  // min valid value
                     std::numeric_limits<Real>::max(),     // max valid value
                     1,                                    // number of dims
                     CellDims);

   auto OSBLDepthIndexField =
       Field::create(OSBLDepthIndexFldName,                      // field name
                     "ocean surface boundary layer depth index", // long name
                     "",                                         // units
                     "",                             // CF standard name
                     -1,                             // min valid value
                     std::numeric_limits<I4>::max(), // max valid value
                     1,                              // number of dims
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
   KPPGroup->addField(OSBLDepthFldName);
   KPPGroup->addField(OSBLDepthIndexFldName);
   KPPGroup->addField(NonLocalFluxFldName);
   KPPGroup->addField(BulkRichardsonFldName);
   KPPGroup->addField(BulkRichardsonShearFldName);
   KPPGroup->addField(UnresolvedShearFldName);
   KPPGroup->addField(BuoyancyJumpFldName);
   KPPGroup->addField(TurbulentVelScaleFldName);
   KPPGroup->addField(PotentialDensityFldName);
   KPPGroup->addField(SurfFricVelFldName);
   KPPGroup->addField(SurfBuoyFluxFldName);

   OSBLDepthField->attachData<Array1DReal>(OSBLDepth, false);
   OSBLDepthIndexField->attachData<Array1DI4>(OSBLDepthIndex, false);
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

   OSBLDepthField->addMetadata("_FillValue", FillValueReal);
   OSBLDepthIndexField->addMetadata("_FillValue", FillValueI4);
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
            OSBLDepthFldName, OSBLDepthIndexFldName, NonLocalFluxFldName,
            BulkRichardsonFldName, BulkRichardsonShearFldName,
            UnresolvedShearFldName, BuoyancyJumpFldName,
            TurbulentVelScaleFldName, PotentialDensityFldName,
            SurfFricVelFldName, SurfBuoyFluxFldName);
}

} // namespace OMEGA
