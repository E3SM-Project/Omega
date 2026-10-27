//===-- ocn/GLSMix.cpp - Generic Length Scale Mixing Implementation --*- C++
//-*-===//
//
/// \file
/// \brief Implementation of the Generic Length Scale (GLS) turbulence
/// closure
///
/// Solves prognostic equations for turbulent kinetic energy (TKE) and a
/// generic length-scale variable (psi) per water column (vertical
/// diffusion, shear/buoyancy production, and dissipation only -- horizontal
/// and vertical advection of TKE/psi are not included, see GLSMix.h and
/// doc/design/GLSMix.md), following Umlauf and Burchard (2003) and the
/// ROMS/REMORA gls_corstep.F formulation.
//
//===----------------------------------------------------------------------===//

#include "GLSMix.h"
#include "DataTypes.h"
#include "Error.h"
#include "Field.h"
#include "Forcing.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "KPPMix.h"
#include "Logging.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TriDiagSolvers.h"
#include "VertCoord.h"
#include <limits>

namespace OMEGA {

// Singleton instance
GLSMix *GLSMix::Instance = nullptr;

namespace {

bool stabilityTypeFromString(const std::string &Name, GLSStabilityType &Type) {
   if (Name == "Canuto_A") {
      Type = GLSStabilityType::CanutoA;
      return true;
   }
   if (Name == "Canuto_B") {
      Type = GLSStabilityType::CanutoB;
      return true;
   }
   if (Name == "Galperin") {
      Type = GLSStabilityType::Galperin;
      return true;
   }
   return false;
}

const char *stabilityTypeName(GLSStabilityType Type) {
   switch (Type) {
   case GLSStabilityType::CanutoB:
      return "Canuto_B";
   case GLSStabilityType::Galperin:
      return "Galperin";
   default:
      return "Canuto_A";
   }
}

} // namespace

GLSShearSquaredOnCell::GLSShearSquaredOnCell(const HorzMesh *Mesh,
                                             const VertCoord *VCoord)
    : GeomZMid(VCoord->GeomZMid), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnCell(Mesh->CellsOnCell), NEdgesOnCell(Mesh->NEdgesOnCell),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop), DcEdge(Mesh->DcEdge),
      DvEdge(Mesh->DvEdge) {}

/// Constructor for GLSMix
GLSMix::GLSMix(const std::string &InName, const HorzMesh *InMesh,
              const VertCoord *InVCoord)
    : Name(InName), Mesh(InMesh), VCoord(InVCoord) {

   VertVisc = Array2DReal("GLSVertVisc", Mesh->NCellsSize,
                          VCoord->NVertLayersP1);
   VertDiff = Array2DReal("GLSVertDiff", Mesh->NCellsSize,
                          VCoord->NVertLayersP1);
   TKE    = Array2DReal("GLSTke", Mesh->NCellsSize, VCoord->NVertLayersP1);
   Psi    = Array2DReal("GLSPsi", Mesh->NCellsSize, VCoord->NVertLayersP1);
   Lscale = Array2DReal("GLSLscale", Mesh->NCellsSize, VCoord->NVertLayersP1);
   Akk    = Array2DReal("GLSAkk", Mesh->NCellsSize, VCoord->NVertLayersP1);
   Akp    = Array2DReal("GLSAkp", Mesh->NCellsSize, VCoord->NVertLayersP1);

   // Initialize prognostic/diagnostic state to background values so a
   // restart read (if any) has sane defaults to fall back on.
   deepCopy(TKE, Kmin);
   deepCopy(Psi, Pmin);
   deepCopy(Lscale, 0.0_Real);
   deepCopy(Akk, AkkBak);
   deepCopy(Akp, AkpBak);
   deepCopy(VertVisc, AkvBak);
   deepCopy(VertDiff, AktBak);

   VertDiffFldName = "GLSVertDiff";
   VertViscFldName = "GLSVertVisc";
   TkeFldName      = "GLSTke";
   PsiFldName      = "GLSPsi";
   LscaleFldName   = "GLSLscale";
   if (Name != "Default") {
      VertDiffFldName.append(Name);
      VertViscFldName.append(Name);
      TkeFldName.append(Name);
      PsiFldName.append(Name);
      LscaleFldName.append(Name);
   }

   defineFields();
}

/// Destructor for GLSMix
GLSMix::~GLSMix() {}

/// Get singleton instance
GLSMix *GLSMix::getInstance() { return Instance; }

/// Destroy singleton instance
void GLSMix::destroyInstance() {
   delete Instance;
   Instance = nullptr;
}

/// Initialize GLSMix from configuration
void GLSMix::init() {
   if (!Instance) {
      Instance = new GLSMix("Default", HorzMesh::getDefault(),
                            VertCoord::getDefault());
   }

   GLSMix *DefGLSMix   = GLSMix::getInstance();
   Config *OmegaConfig = Config::getOmegaConfig();
   OMEGA_REQUIRE(OmegaConfig, "Null OmegaConfig pointer in GLSMix::init");

   Config VertMixConfig("VertMix");
   Error Err;
   Err += OmegaConfig->get(VertMixConfig);
   CHECK_ERROR_ABORT(Err, "GLSMix::init: VertMix group not found in Config");

   Config GLSConfig("GLS");
   Err += VertMixConfig.get(GLSConfig);
   if (Err.isFail()) {
      LOG_WARN("GLSMix::init: GLS subgroup not found, using defaults");
      return;
   }

   bool Enable = false;
   Err         = GLSConfig.get("Enable", Enable);
   if (Err.isSuccess())
      DefGLSMix->Enabled = Enable;

   auto getOpt = [&](const std::string &Key, Real &Val) {
      Error OptErr = GLSConfig.get(Key, Val);
      if (!OptErr.isSuccess())
         OptErr.reset();
   };

   getOpt("P", DefGLSMix->P);
   getOpt("M", DefGLSMix->M);
   getOpt("N", DefGLSMix->N);
   getOpt("Cmu0", DefGLSMix->Cmu0);
   getOpt("C1", DefGLSMix->C1);
   getOpt("C2", DefGLSMix->C2);
   getOpt("C3Minus", DefGLSMix->C3Minus);
   getOpt("C3Plus", DefGLSMix->C3Plus);
   getOpt("SigmaK", DefGLSMix->SigmaK);
   getOpt("SigmaPsi", DefGLSMix->SigmaPsi);
   getOpt("Kmin", DefGLSMix->Kmin);
   getOpt("Pmin", DefGLSMix->Pmin);
   getOpt("AkvBak", DefGLSMix->AkvBak);
   getOpt("AktBak", DefGLSMix->AktBak);
   getOpt("AkkBak", DefGLSMix->AkkBak);
   getOpt("AkpBak", DefGLSMix->AkpBak);
   getOpt("Zos", DefGLSMix->Zos);
   getOpt("Zob", DefGLSMix->Zob);
   getOpt("BottomDragCoeff", DefGLSMix->BottomDragCoeff);
   getOpt("LengthScaleCapFraction", DefGLSMix->LengthScaleCapFraction);

   bool CapEnable = DefGLSMix->LengthScaleCapEnable;
   Error CapErr   = GLSConfig.get("LengthScaleCapEnable", CapEnable);
   if (CapErr.isSuccess())
      DefGLSMix->LengthScaleCapEnable = CapEnable;

   // Clamp floors strictly positive so kernels can divide by them safely.
   DefGLSMix->Kmin = Kokkos::fmax(1.0e-12_Real, DefGLSMix->Kmin);
   DefGLSMix->Pmin = Kokkos::fmax(1.0e-16_Real, DefGLSMix->Pmin);

   std::string StabilityStr = "Canuto_A";
   Error StabErr            = GLSConfig.get("StabilityFunction", StabilityStr);
   if (!StabErr.isSuccess())
      StabErr.reset();
   if (!stabilityTypeFromString(StabilityStr, DefGLSMix->StabilityType)) {
      ABORT_ERROR("GLSMix::init: Invalid StabilityFunction='{}', must be "
                  "Canuto_A, Canuto_B, or Galperin",
                  StabilityStr);
   }

   DefGLSMix->Derived = computeGLSDerivedConstants(
       DefGLSMix->StabilityType, DefGLSMix->P, DefGLSMix->M, DefGLSMix->N,
       DefGLSMix->Cmu0, DefGLSMix->C1, DefGLSMix->C2, DefGLSMix->C3Minus,
       DefGLSMix->C3Plus, DefGLSMix->SigmaK, DefGLSMix->SigmaPsi);

   // GLS and KPP are mutually exclusive alternatives for boundary/interior
   // mixing; VertMix::computeVertMix() assumes exactly one of them merges
   // into the final VertDiff/VertVisc.
   KPPMix *KPPInstance = KPPMix::getInstance();
   if (DefGLSMix->Enabled && KPPInstance && KPPInstance->Enabled) {
      ABORT_ERROR("GLSMix::init: GLS and KPP cannot both be enabled; they "
                  "are alternative vertical mixing schemes. Disable one of "
                  "VertMix:GLS:Enable or VertMix:KPP:Enable");
   }

   LOG_INFO("GLSMix::init: GLS initialized enabled={} stability={} p={} "
            "m={} n={}",
            DefGLSMix->Enabled, stabilityTypeName(DefGLSMix->StabilityType),
            DefGLSMix->P, DefGLSMix->M, DefGLSMix->N);
} // end init

/// Convenience wrapper: fetch surface stress from the default Forcing
/// instance and call advance().
void GLSMix::computeGLSFields(const Array2DReal &NormalVelocity,
                              const Array2DReal &TangentialVelocity,
                              const Array2DReal &BruntVaisalaFreqSq,
                              const Array2DReal &PseudoThickCell,
                              const Array2DReal &KineticEnergyCell, Real DT) {
   if (!Enabled)
      return;

   const Forcing *ForcingState = Forcing::getDefault();
   if (!ForcingState) {
      LOG_WARN("GLSMix::computeGLSFields: Forcing has not been "
               "initialized, skipping GLS update");
      return;
   }

   advance(NormalVelocity, TangentialVelocity, BruntVaisalaFreqSq,
          PseudoThickCell, KineticEnergyCell,
          ForcingState->SfcStressForcing.ZonalStressCell,
          ForcingState->SfcStressForcing.MeridStressCell, DT);
}

/// Advance the TKE/psi equations by one time step and update VertDiff/
/// VertVisc from the resulting stability functions.
void GLSMix::advance(const Array2DReal &NormalVelocity,
                     const Array2DReal &TangentialVelocity,
                     const Array2DReal &BruntVaisalaFreqSq,
                     const Array2DReal &PseudoThickCell,
                     const Array2DReal &KineticEnergyCell,
                     const Array1DReal &ZonalStressCell,
                     const Array1DReal &MeridStressCell, Real DT) {
   if (!Enabled)
      return;

   Pacer::start("GLS:advance", 1);

   const I4 NCellsAll     = Mesh->NCellsAll;
   const I4 NVertLayersP1 = VCoord->NVertLayersP1;

   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(GeomZInterface, VCoord->GeomZInterface);
   OMEGA_SCOPE(ZonalStress, ZonalStressCell);
   OMEGA_SCOPE(MeridStress, MeridStressCell);

   OMEGA_SCOPE(LocTKE, TKE);
   OMEGA_SCOPE(LocPsi, Psi);
   OMEGA_SCOPE(LocLscale, Lscale);
   OMEGA_SCOPE(LocAkk, Akk);
   OMEGA_SCOPE(LocAkp, Akp);
   OMEGA_SCOPE(LocVertVisc, VertVisc);
   OMEGA_SCOPE(LocVertDiff, VertDiff);

   const GLSDerivedConstants LocDerived = Derived;
   const GLSStabilityType LocStabType   = StabilityType;
   const Real LocM = M, LocN = N;
   const Real LocCmu0 = Cmu0, LocC1 = C1, LocC2 = C2;
   const Real LocC3Minus = C3Minus, LocC3Plus = C3Plus;
   const Real LocSigmaK = SigmaK, LocSigmaPsi = SigmaPsi;
   const Real LocKmin = Kmin, LocPmin = Pmin;
   const Real LocAkvBak = AkvBak, LocAktBak = AktBak;
   const Real LocAkkBak = AkkBak, LocAkpBak = AkpBak;
   const Real LocZos = Zos, LocZob = Zob;
   const Real LocBottomDragCoeff = BottomDragCoeff;
   const bool LocCapEnable       = LengthScaleCapEnable;
   const Real LocCapFraction     = LengthScaleCapFraction;

   // ------------------------------------------------------------------
   // Step 1: shear-squared at interior layer interfaces.
   // ------------------------------------------------------------------
   Array2DReal ShearSq("GLS-ShearSq", Mesh->NCellsSize, NVertLayersP1);
   deepCopy(ShearSq, 0.0_Real);
   GLSShearSquaredOnCell ComputeShearSq(Mesh, VCoord);
   OMEGA_SCOPE(LocComputeShearSq, ComputeShearSq);

   parallelFor(
       "GLS-ShearSq", {NCellsAll, NVertLayersP1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell);
          if (K > KMin && K <= KMax) {
             LocComputeShearSq(ShearSq, ICell, K, NormalVelocity,
                              TangentialVelocity);
          }
       });

   // ------------------------------------------------------------------
   // Step 2: shear/buoyancy production and local capacity (thickness)
   // terms at interior interfaces, from the previous step's Akv/Akt.
   // ------------------------------------------------------------------
   Array2DReal CffHz("GLS-CffHz", Mesh->NCellsSize, NVertLayersP1);
   Array2DReal KProd("GLS-KProd", Mesh->NCellsSize, NVertLayersP1);
   Array2DReal PProd("GLS-PProd", Mesh->NCellsSize, NVertLayersP1);
   Array2DReal KDissExtra("GLS-KDissExtra", Mesh->NCellsSize, NVertLayersP1);
   Array2DReal PDissExtra("GLS-PDissExtra", Mesh->NCellsSize, NVertLayersP1);

   parallelFor(
       "GLS-Production", {NCellsAll, NVertLayersP1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell);
          if (K <= KMin || K > KMax)
             return;

          const Real ThickAbove = PseudoThickCell(ICell, K - 1);
          const Real ThickBelow = PseudoThickCell(ICell, K);
          CffHz(ICell, K)       = 0.5_Real * (ThickAbove + ThickBelow);

          const Real Strat2 = BruntVaisalaFreqSq(ICell, K);
          const Real AkvOld = LocVertVisc(ICell, K);
          const Real AktOld = LocVertDiff(ICell, K);
          const Real GlsC3  = (Strat2 > 0.0_Real) ? LocC3Minus : LocC3Plus;

          Real KprodVal = ShearSq(ICell, K) * (AkvOld - LocAkvBak) -
                         Strat2 * (AktOld - LocAktBak);
          Real PprodVal = LocC1 * ShearSq(ICell, K) * (AkvOld - LocAkvBak) -
                         GlsC3 * Strat2 * (AktOld - LocAktBak);

          Real Cff1 = 1.0_Real;
          Real Cff2 = 1.0_Real;
          if (KprodVal < 0.0_Real) {
             Cff1 = 0.0_Real;
             KprodVal += Strat2 * (AktOld - LocAktBak);
          }
          if (PprodVal < 0.0_Real) {
             Cff2 = 0.0_Real;
             PprodVal += GlsC3 * Strat2 * (AktOld - LocAktBak);
          }

          KProd(ICell, K) = KprodVal;
          PProd(ICell, K) = PprodVal;
          KDissExtra(ICell, K) =
              (1.0_Real - Cff1) * Strat2 * (AktOld - LocAktBak);
          PDissExtra(ICell, K) =
              (1.0_Real - Cff2) * GlsC3 * Strat2 * (AktOld - LocAktBak);
       });

   // ------------------------------------------------------------------
   // Step 3/4: implicit solves for TKE and psi. Boundary rows (surface and
   // bottom) are Dirichlet, derived from the surface wind stress / a
   // quadratic bottom-drag estimate; interior rows are diffusion (via Akk
   // or Akp from the previous step) plus a semi-implicit linear dissipation
   // term and the explicit production forcing from Step 2.
   // ------------------------------------------------------------------
   auto LConfig = TriDiagDiffSolver::makeLaunchConfig(NCellsAll, NVertLayersP1);
   const int LocVecLength = VecLength;

   // Solved values are written to separate arrays rather than back into
   // LocTKE/LocPsi in place, because Step 5 below needs both the old (pre-
   // solve) and new TKE/length-scale values together to time-average the
   // turbulent velocity scale (see the Ql computation).
   Array2DReal TkeSolved("GLS-TkeSolved", Mesh->NCellsSize, NVertLayersP1);
   Array2DReal PsiSolved("GLS-PsiSolved", Mesh->NCellsSize, NVertLayersP1);

   // --- TKE solve ---
   parallelForOuter(
       LConfig, KOKKOS_LAMBDA(int, const TeamMember &Team) {
          const int IStart = Team.league_rank() * LocVecLength;

          TriDiagDiffScratch Scratch(Team, NVertLayersP1);

          parallelForInner(Team, NVertLayersP1, [=](int K) {
             for (int IVec = 0; IVec < LocVecLength; ++IVec) {
                const int ICell = IStart + IVec;
                if (ICell >= NCellsAll) {
                   Scratch.G(K, IVec) = 0.0_Real;
                   Scratch.H(K, IVec) = 1.0_Real;
                   Scratch.X(K, IVec) = 0.0_Real;
                   continue;
                }

                const int KMin = MinLayerCell(ICell);
                const int KMax = MaxLayerCell(ICell) + 1;

                if (K < KMin || K > KMax) {
                   Scratch.G(K, IVec) = 0.0_Real;
                   Scratch.H(K, IVec) = 1.0_Real;
                   Scratch.X(K, IVec) = 0.0_Real;
                   continue;
                }

                // Surface/bottom Dirichlet target values for this column.
                // Needed both by the boundary rows themselves (below) and by
                // the adjacent interior rows, which receive the boundary's
                // diffusive influence as a known flux injected directly into
                // their own H/X (see note below on why G cannot carry it).
                const Real TauMag = Kokkos::sqrt(
                    ZonalStress(ICell) * ZonalStress(ICell) +
                    MeridStress(ICell) * MeridStress(ICell));
                const Real TkeSurfBC = Kokkos::max(
                    LocDerived.CmuFac3 * (TauMag / RhoSw), LocKmin);
                const Real BotSpeedSq =
                    2.0_Real * KineticEnergyCell(ICell, MaxLayerCell(ICell));
                const Real TkeBotBC = Kokkos::max(
                    LocDerived.CmuFac3 * LocBottomDragCoeff * BotSpeedSq,
                    LocKmin);

                if (K == KMin) {
                   // Dirichlet row: G must be 0 here, since the solver adds
                   // this row's own G into its own diagonal (H(0)+=G(0)),
                   // which would corrupt a Dirichlet assignment if nonzero.
                   Scratch.G(K, IVec) = 0.0_Real;
                   Scratch.H(K, IVec) = 1.0_Real;
                   Scratch.X(K, IVec) = TkeSurfBC;
                } else if (K == KMax) {
                   Scratch.G(K, IVec) = 0.0_Real;
                   Scratch.H(K, IVec) = 1.0_Real;
                   Scratch.X(K, IVec) = TkeBotBC;
                } else {
                   const Real TkeOldK =
                       Kokkos::max(LocTKE(ICell, K), LocKmin);
                   const Real PsiOldK =
                       Kokkos::max(LocPsi(ICell, K), LocPmin);
                   const Real GlsExpMExp1 =
                       1.0_Real / Kokkos::pow(PsiOldK, LocDerived.GlsExp1);
                   const Real TkeExpExp2 =
                       Kokkos::pow(TkeOldK, LocDerived.TkeExp2);
                   const Real DissCoeff =
                       LocDerived.CmuFac2 * GlsExpMExp1 * TkeExpExp2;

                   const Real CffHzK = CffHz(ICell, K);
                   Real HVal = CffHzK * (1.0_Real + DT * DissCoeff) +
                              DT * CffHzK * KDissExtra(ICell, K) / TkeOldK;
                   Real XVal = CffHzK * TkeOldK + DT * CffHzK * KProd(ICell, K);
                   Real GVal = 0.0_Real;

                   // Coupling to K+1: a normal matrix entry, unless K+1 is
                   // the bottom Dirichlet row, in which case that row's own
                   // G is forced to 0 above, so this row must inject the
                   // flux into its own equation directly instead.
                   if (K + 1 < KMax) {
                      const Real AkkLayer = 0.5_Real * (LocAkk(ICell, K) +
                                                       LocAkk(ICell, K + 1));
                      GVal = DT * AkkLayer / PseudoThickCell(ICell, K);
                   } else {
                      const Real AkkLayer =
                          0.5_Real * (LocAkk(ICell, K) + LocAkk(ICell, KMax));
                      const Real GBotFlux =
                          DT * AkkLayer / PseudoThickCell(ICell, K);
                      HVal += GBotFlux;
                      XVal += GBotFlux * TkeBotBC;
                   }

                   // Coupling to K-1: the solver derives this row's link to
                   // K-1 from G stored AT K-1 (via the Alpha recursion). If
                   // K-1 is the surface Dirichlet row, that G is 0 (as set
                   // above), so the flux must be injected here instead.
                   if (K - 1 == KMin) {
                      const Real AkkLayer =
                          0.5_Real * (LocAkk(ICell, KMin) + LocAkk(ICell, K));
                      const Real GSurfFlux =
                          DT * AkkLayer / PseudoThickCell(ICell, KMin);
                      HVal += GSurfFlux;
                      XVal += GSurfFlux * TkeSurfBC;
                   }

                   Scratch.H(K, IVec) = HVal;
                   Scratch.X(K, IVec) = XVal;
                   Scratch.G(K, IVec) = GVal;
                }
             }
          });

          Team.team_barrier();
          TriDiagDiffSolver::solve(Team, Scratch);
          Team.team_barrier();

          parallelForInner(Team, NVertLayersP1, [=](int K) {
             for (int IVec = 0; IVec < LocVecLength; ++IVec) {
                const int ICell = IStart + IVec;
                if (ICell >= NCellsAll)
                   continue;
                const int KMin = MinLayerCell(ICell);
                const int KMax = MaxLayerCell(ICell) + 1;
                if (K >= KMin && K <= KMax) {
                   TkeSolved(ICell, K) =
                       Kokkos::max(Scratch.X(K, IVec), LocKmin);
                }
             }
          });
       });

   // --- Psi solve ---
   parallelForOuter(
       LConfig, KOKKOS_LAMBDA(int, const TeamMember &Team) {
          const int IStart = Team.league_rank() * LocVecLength;

          TriDiagDiffScratch Scratch(Team, NVertLayersP1);

          parallelForInner(Team, NVertLayersP1, [=](int K) {
             for (int IVec = 0; IVec < LocVecLength; ++IVec) {
                const int ICell = IStart + IVec;
                if (ICell >= NCellsAll) {
                   Scratch.G(K, IVec) = 0.0_Real;
                   Scratch.H(K, IVec) = 1.0_Real;
                   Scratch.X(K, IVec) = 0.0_Real;
                   continue;
                }

                const int KMin = MinLayerCell(ICell);
                const int KMax = MaxLayerCell(ICell) + 1;

                if (K < KMin || K > KMax) {
                   Scratch.G(K, IVec) = 0.0_Real;
                   Scratch.H(K, IVec) = 1.0_Real;
                   Scratch.X(K, IVec) = 0.0_Real;
                   continue;
                }

                const Real ZosEff = Kokkos::max(LocZos, 1.0e-4_Real);
                const Real ZobEff = Kokkos::max(LocZob, 1.0e-4_Real);

                const Real TauMag = Kokkos::sqrt(
                    ZonalStress(ICell) * ZonalStress(ICell) +
                    MeridStress(ICell) * MeridStress(ICell));
                const Real TkeSurfBC = Kokkos::max(
                    LocDerived.CmuFac3 * (TauMag / RhoSw), LocKmin);
                const Real PsiSurfBC = Kokkos::max(
                    LocDerived.Cmu0ExpP * Kokkos::pow(TkeSurfBC, LocM) *
                        Kokkos::pow(VonKar * ZosEff, LocN),
                    LocPmin);
                const Real BotSpeedSq =
                    2.0_Real * KineticEnergyCell(ICell, MaxLayerCell(ICell));
                const Real TkeBotBC = Kokkos::max(
                    LocDerived.CmuFac3 * LocBottomDragCoeff * BotSpeedSq,
                    LocKmin);
                const Real PsiBotBC = Kokkos::max(
                    LocDerived.Cmu0ExpP * Kokkos::pow(TkeBotBC, LocM) *
                        Kokkos::pow(VonKar * ZobEff, LocN),
                    LocPmin);

                if (K == KMin) {
                   // See the TKE solve above for why G must be 0 on a
                   // Dirichlet row.
                   Scratch.G(K, IVec) = 0.0_Real;
                   Scratch.H(K, IVec) = 1.0_Real;
                   Scratch.X(K, IVec) = PsiSurfBC;
                } else if (K == KMax) {
                   Scratch.G(K, IVec) = 0.0_Real;
                   Scratch.H(K, IVec) = 1.0_Real;
                   Scratch.X(K, IVec) = PsiBotBC;
                } else {
                   const Real TkeOldK =
                       Kokkos::max(LocTKE(ICell, K), LocKmin);
                   const Real PsiOldK =
                       Kokkos::max(LocPsi(ICell, K), LocPmin);
                   const Real GlsExpMExp1 =
                       1.0_Real / Kokkos::pow(PsiOldK, LocDerived.GlsExp1);
                   const Real TkeExpExp2 =
                       Kokkos::pow(TkeOldK, LocDerived.TkeExp2);
                   const Real DissCoeff = LocC2 * LocDerived.CmuFac2 *
                                          GlsExpMExp1 * TkeExpExp2;

                   const Real CffHzK = CffHz(ICell, K);
                   Real HVal = CffHzK * (1.0_Real + DT * DissCoeff) +
                              DT * CffHzK * PDissExtra(ICell, K) / TkeOldK;
                   Real XVal = CffHzK * PsiOldK +
                              DT * CffHzK * PProd(ICell, K) * PsiOldK / TkeOldK;
                   Real GVal = 0.0_Real;

                   // See the TKE solve above: boundary-adjacent rows inject
                   // the Dirichlet neighbor's flux directly rather than
                   // relying on the solver's own coupling, since the
                   // boundary row's G is forced to 0.
                   if (K + 1 < KMax) {
                      const Real AkpLayer = 0.5_Real * (LocAkp(ICell, K) +
                                                       LocAkp(ICell, K + 1));
                      GVal = DT * AkpLayer / PseudoThickCell(ICell, K);
                   } else {
                      const Real AkpLayer =
                          0.5_Real * (LocAkp(ICell, K) + LocAkp(ICell, KMax));
                      const Real GBotFlux =
                          DT * AkpLayer / PseudoThickCell(ICell, K);
                      HVal += GBotFlux;
                      XVal += GBotFlux * PsiBotBC;
                   }

                   if (K - 1 == KMin) {
                      const Real AkpLayer =
                          0.5_Real * (LocAkp(ICell, KMin) + LocAkp(ICell, K));
                      const Real GSurfFlux =
                          DT * AkpLayer / PseudoThickCell(ICell, KMin);
                      HVal += GSurfFlux;
                      XVal += GSurfFlux * PsiSurfBC;
                   }

                   Scratch.H(K, IVec) = HVal;
                   Scratch.X(K, IVec) = XVal;
                   Scratch.G(K, IVec) = GVal;
                }
             }
          });

          Team.team_barrier();
          TriDiagDiffSolver::solve(Team, Scratch);
          Team.team_barrier();

          parallelForInner(Team, NVertLayersP1, [=](int K) {
             for (int IVec = 0; IVec < LocVecLength; ++IVec) {
                const int ICell = IStart + IVec;
                if (ICell >= NCellsAll)
                   continue;
                const int KMin = MinLayerCell(ICell);
                const int KMax = MaxLayerCell(ICell) + 1;
                if (K >= KMin && K <= KMax) {
                   PsiSolved(ICell, K) =
                       Kokkos::max(Scratch.X(K, IVec), LocPmin);
                }
             }
          });
       });

   // ------------------------------------------------------------------
   // Step 5: length-scale limiter, stability functions, and final
   // viscosity/diffusivity. The unconditional hard cap below is the fix
   // for the REMORA/ROMS length-scale-limiter bug: the Galperin bound only
   // applies under stable stratification (N^2 > 0), so without this cap
   // the length scale (and hence viscosity) is unbounded under
   // neutral/unstable stratification. See doc/design/GLSMix.md.
   // ------------------------------------------------------------------
   constexpr Real Eps = 1.0e-10_Real;

   parallelFor(
       "GLS-StabilityFunctions", {NCellsAll, NVertLayersP1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell) + 1;
          if (K < KMin || K > KMax)
             return;

          const Real ZosEff = Kokkos::max(LocZos, 1.0e-4_Real);
          const Real ZobEff = Kokkos::max(LocZob, 1.0e-4_Real);

          if (K == KMin) {
             const Real NewVisc =
                 LocAkvBak + VonKar * ZosEff * LocCmu0 *
                                 Kokkos::sqrt(TkeSolved(ICell, K));
             LocVertVisc(ICell, K) = 0.0_Real; // not used by momentum mixing
             LocVertDiff(ICell, K) = 0.0_Real; // (see VertMix boundary zeroing)
             LocAkk(ICell, K)      = LocAkkBak + NewVisc / LocSigmaK;
             LocAkp(ICell, K)      = LocAkpBak + NewVisc / LocSigmaPsi;
             LocLscale(ICell, K)   = 0.0_Real;
             LocTKE(ICell, K)      = TkeSolved(ICell, K);
             LocPsi(ICell, K)      = PsiSolved(ICell, K);
             return;
          }
          if (K == KMax) {
             const Real NewVisc =
                 LocAkvBak + VonKar * ZobEff * LocCmu0 *
                                 Kokkos::sqrt(TkeSolved(ICell, K));
             LocVertVisc(ICell, K) = 0.0_Real;
             LocVertDiff(ICell, K) = 0.0_Real;
             LocAkk(ICell, K)      = LocAkkBak + NewVisc / LocSigmaK;
             LocAkp(ICell, K)      = LocAkpBak + NewVisc / LocSigmaPsi;
             LocLscale(ICell, K)   = 0.0_Real;
             LocTKE(ICell, K)      = TkeSolved(ICell, K);
             LocPsi(ICell, K)      = PsiSolved(ICell, K);
             return;
          }

          // LocTKE/LocLscale here are still the pre-solve (old) values; the
          // just-solved new values live in TkeSolved/PsiSolved until the end
          // of this kernel, when LocTKE/LocPsi are finally updated below.
          Real TkeNew = Kokkos::max(TkeSolved(ICell, K), LocKmin);
          Real PsiNew = Kokkos::max(PsiSolved(ICell, K), LocPmin);

          // Bound psi by the dissipation-rate comparison (Warner et al.
          // 2005), before computing the length scale.
          const Real GlsComparison =
              LocDerived.GlsFac5 * Kokkos::pow(TkeNew, LocDerived.TkeExp4) *
              Kokkos::pow(Kokkos::sqrt(Kokkos::max(
                             0.0_Real, BruntVaisalaFreqSq(ICell, K))) +
                             Eps,
                         -LocN);
          PsiNew = (LocN >= 0.0_Real) ? Kokkos::min(PsiNew, GlsComparison)
                                      : Kokkos::max(PsiNew, GlsComparison);

          const Real LsUnlmt =
              Kokkos::max(Eps, Kokkos::pow(PsiNew, LocDerived.GlsExp1) *
                                   LocDerived.CmuFac1 *
                                   Kokkos::pow(TkeNew, -LocDerived.TkeExp1));

          const Real Strat2 = BruntVaisalaFreqSq(ICell, K);
          Real LsLmt =
              (Strat2 > 0.0_Real)
                  ? Kokkos::min(LsUnlmt,
                               Kokkos::sqrt(0.56_Real * TkeNew /
                                           (Kokkos::max(0.0_Real, Strat2) +
                                            Eps)))
                  : LsUnlmt;

          if (LocCapEnable) {
             const Real Depth = Kokkos::abs(
                 GeomZInterface(ICell, KMin) - GeomZInterface(ICell, KMax));
             LsLmt = Kokkos::min(LsLmt, LocCapFraction * Depth);
          }

          PsiNew = Kokkos::max(
              LocDerived.Cmu0ExpP * Kokkos::pow(TkeNew, LocM) *
                  Kokkos::pow(LsLmt, LocN),
              LocPmin);

          Real Gh = Kokkos::min(LocDerived.Gh0, -Strat2 * LsLmt * LsLmt /
                                                    (2.0_Real * TkeNew));
          Gh      = Kokkos::min(Gh, Gh - (Gh - LocDerived.Ghcri) *
                                            (Gh - LocDerived.Ghcri) /
                                            (Gh + LocDerived.Gh0 -
                                             2.0_Real * LocDerived.Ghcri));
          Gh = Kokkos::max(Gh, LocDerived.Ghmin);

          Real Gm = computeCanutoGmBound(LocStabType, LocDerived, Gh);
          Gm = Kokkos::min(
              Gm, ShearSq(ICell, K) * LsLmt * LsLmt / (2.0_Real * TkeNew));

          Real Sm = 0.0_Real, Sh = 0.0_Real;
          computeStabilityFunctions(LocStabType, LocDerived, Gh, Gm, Sm, Sh);

          const Real Ql =
              LocDerived.Sqrt2 * 0.5_Real *
              (LsLmt * Kokkos::sqrt(TkeNew) +
               LocLscale(ICell, K) * Kokkos::sqrt(LocTKE(ICell, K)));

          LocVertVisc(ICell, K) = LocAkvBak + Sm * Ql;
          LocVertDiff(ICell, K) = LocAktBak + Sh * Ql;
          LocAkk(ICell, K)      = LocAkkBak + Sm * Ql / LocSigmaK;
          LocAkp(ICell, K)      = LocAkpBak + Sm * Ql / LocSigmaPsi;
          LocLscale(ICell, K)   = LsLmt;
          LocTKE(ICell, K)      = TkeNew;
          LocPsi(ICell, K)      = PsiNew;
       });

   Pacer::stop("GLS:advance", 1);
} // advance

/// Define IO and restart fields
void GLSMix::defineFields() {

   int NDims = 2;
   std::vector<std::string> DimNames(NDims);
   DimNames[0] = "NCells";
   DimNames[1] = "NVertLayersP1";

   auto VertDiffField = Field::create(
       VertDiffFldName, "GLS vertical diffusivity at layer interfaces",
       "m2 s-1", "", 0.0, std::numeric_limits<Real>::max(), NDims, DimNames);
   auto VertViscField = Field::create(
       VertViscFldName, "GLS vertical viscosity at layer interfaces", "m2 s-1",
       "", 0.0, std::numeric_limits<Real>::max(), NDims, DimNames);
   auto TkeField = Field::create(
       TkeFldName, "GLS turbulent kinetic energy", "m2 s-2", "", 0.0,
       std::numeric_limits<Real>::max(), NDims, DimNames);
   auto PsiField = Field::create(
       PsiFldName, "GLS generic length-scale variable", "", "", 0.0,
       std::numeric_limits<Real>::max(), NDims, DimNames);
   auto LscaleField = Field::create(
       LscaleFldName, "GLS diagnostic turbulence length scale", "m", "", 0.0,
       std::numeric_limits<Real>::max(), NDims, DimNames);

   std::string GroupName = "GLS";
   if (Name != "Default")
      GroupName.append(Name);
   auto GLSGroup = FieldGroup::create(GroupName);
   GLSGroup->addField(VertDiffFldName);
   GLSGroup->addField(VertViscFldName);
   GLSGroup->addField(TkeFldName);
   GLSGroup->addField(PsiFldName);
   GLSGroup->addField(LscaleFldName);

   // TKE and psi are prognostic state and must survive restarts; the
   // remaining fields are recomputed every step and are diagnostic-only.
   if (!FieldGroup::exists("Restart"))
      FieldGroup::create("Restart");
   FieldGroup::addFieldToGroup(TkeFldName, "Restart");
   FieldGroup::addFieldToGroup(PsiFldName, "Restart");

   VertDiffField->attachData<Array2DReal>(VertDiff);
   VertViscField->attachData<Array2DReal>(VertVisc);
   TkeField->attachData<Array2DReal>(TKE);
   PsiField->attachData<Array2DReal>(Psi);
   LscaleField->attachData<Array2DReal>(Lscale);

} // defineFields

} // namespace OMEGA
