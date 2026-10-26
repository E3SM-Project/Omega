#ifndef OMEGA_GLS_MIX_H
#define OMEGA_GLS_MIX_H
//===-- ocn/GLSMix.h - Generic Length Scale Mixing -------*- C++ -*-===//
//
/// \file
/// \brief Generic Length Scale (GLS) two-equation turbulence closure
///
/// This header defines the GLSMix class, a config-selectable alternative to
/// KPP for computing vertical viscosity and diffusivity. GLS solves
/// prognostic equations for turbulent kinetic energy (TKE) and a generic
/// length-scale variable (psi) following Umlauf and Burchard (2003) and the
/// ROMS/REMORA gls_corstep.F formulation, and is intended for use in
/// shallow, strongly-stratified basins (e.g. the Great Lakes) where KPP's
/// open-ocean boundary-layer assumptions are a poor fit.
///
/// Scope of this implementation (see doc/design/GLSMix.md for the full
/// rationale): the TKE/psi equations are solved per water column (vertical
/// diffusion, production, and dissipation only); horizontal and vertical
/// advection of TKE/psi are not included, matching common practice in
/// column-based closures (e.g. GOTM, CVMix-TKE).
//
//===----------------------------------------------------------------------===//

#include "Config.h"
#include "DataTypes.h"
#include "GLSConstants.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"
#include <string>

namespace OMEGA {

class OceanState;
class AuxiliaryState;

/// @brief Shear-squared at layer interfaces, kite-area averaged over the
/// edges of each cell (same averaging as GradRichardsonNum in VertMix.h, but
/// returning the raw shear magnitude rather than a Richardson number).
class GLSShearSquaredOnCell {
 public:
   GLSShearSquaredOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void operator()(const Array2DReal &ShearSq, I4 ICell, I4 K,
                                   const Array2DReal &NormalVelocity,
                                   const Array2DReal &TangentialVelocity)
       const {
      const I4 K1 = K - 1;
      const I4 K2 = K;

      Real Weighted = 0.0_Real;
      Real WeightSum = 1.0e-12_Real;

      for (I4 J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);

         if (K > MaxLayerEdgeTop(JEdge))
            continue;

         const Real DNormVel =
             NormalVelocity(JEdge, K1) - NormalVelocity(JEdge, K2);
         const Real DTanVel =
             TangentialVelocity(JEdge, K1) - TangentialVelocity(JEdge, K2);
         const Real DzEdge =
             0.5_Real * (GeomZMid(ICell, K1) + GeomZMid(CellsOnCell(ICell, J), K1) -
                        (GeomZMid(ICell, K2) +
                         GeomZMid(CellsOnCell(ICell, J), K2)));
         const Real ShearEdge =
             (DNormVel * DNormVel + DTanVel * DTanVel) / (DzEdge * DzEdge);

         const Real Weight = 0.25_Real * DcEdge(JEdge) * DvEdge(JEdge);
         Weighted += Weight * ShearEdge;
         WeightSum += Weight;
      }

      ShearSq(ICell, K) = Weighted / WeightSum;
   }

 private:
   Array2DReal GeomZMid;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnCell;
   Array1DI4 NEdgesOnCell;
   Array1DI4 MaxLayerEdgeTop;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
};

/// @brief Generic Length Scale (GLS) two-equation turbulence closure
///
/// Solves prognostic equations for turbulent kinetic energy k and a generic
/// length-scale variable psi = (cmu0)^p k^m l^n, with p/m/n selecting the
/// closure family (k-epsilon, k-omega, k-kl/Mellor-Yamada 2.5, or generic),
/// and a choice of stability-function formulation (Canuto A/B, Galperin).
/// Produces vertical viscosity/diffusivity at layer interfaces, the same
/// grid as VertMix::VertDiff/VertVisc, so it plugs into the same merge seam
/// KPPMix uses in VertMix::computeVertMix (GLS and KPP are mutually
/// exclusive alternatives, validated at init()).
class GLSMix {

 public:
   /// @brief Singleton instance management
   static GLSMix *getInstance();
   static void init();
   static void destroyInstance();

   /// @brief Advance the TKE/psi equations by one full model time step and
   /// update VertDiff/VertVisc (and the internal Akk/Akp/Lscale diagnostics)
   /// from the resulting stability functions. Must be called once per step,
   /// after the current step's velocity and Brunt-Vaisala frequency have
   /// been updated, and before VertMix::computeVertMix() merges the result.
   /// Takes raw surface-stress arrays (rather than a Forcing pointer) so it
   /// can be unit tested without a full Forcing setup; computeGLSFields()
   /// below is the convenience wrapper that fetches them from Forcing.
   void advance(const Array2DReal &NormalVelocity,
               const Array2DReal &TangentialVelocity,
               const Array2DReal &BruntVaisalaFreqSq,
               const Array2DReal &PseudoThickCell,
               const Array2DReal &KineticEnergyCell,
               const Array1DReal &ZonalStressCell,
               const Array1DReal &MeridStressCell, Real DT);

   /// @brief Convenience wrapper: fetches the surface stress from the
   /// default Forcing instance and calls advance(). This is what
   /// VertMix::VertMixImplicit() calls each step.
   void computeGLSFields(const Array2DReal &NormalVelocity,
                         const Array2DReal &TangentialVelocity,
                         const Array2DReal &BruntVaisalaFreqSq,
                         const Array2DReal &PseudoThickCell,
                         const Array2DReal &KineticEnergyCell, Real DT);

   // =======================================================================
   // Output fields (same grid as VertMix::VertDiff/VertVisc)
   // =======================================================================

   /// @brief Vertical viscosity at layer interfaces (m^2/s)
   Array2DReal VertVisc;
   /// @brief Vertical diffusivity at layer interfaces (m^2/s)
   Array2DReal VertDiff;

   // =======================================================================
   // Prognostic state and diagnostics (persistent across time steps)
   // =======================================================================

   Array2DReal TKE;    ///< Turbulent kinetic energy (m^2/s^2)
   Array2DReal Psi;    ///< Generic length-scale variable
   Array2DReal Lscale; ///< Diagnostic turbulence length scale (m)
   Array2DReal Akk;    ///< TKE vertical-diffusion coefficient (m^2/s)
   Array2DReal Akp;    ///< Psi vertical-diffusion coefficient (m^2/s)

   // =======================================================================
   // Configuration parameters
   // =======================================================================

   bool Enabled = false; ///< Enable/disable GLS mixing

   /// Closure-family exponents (defaults: k-epsilon)
   Real P = 3.0_Real;
   Real M = 1.5_Real;
   Real N = -1.0_Real;

   GLSStabilityType StabilityType = GLSStabilityType::CanutoA;

   Real Cmu0     = 0.5477_Real;
   Real C1       = 1.44_Real;
   Real C2       = 1.92_Real;
   Real C3Minus  = -0.4_Real; ///< Used under stable stratification (N^2 > 0)
   Real C3Plus   = 1.0_Real;  ///< Used under unstable stratification
   Real SigmaK   = 1.0_Real;
   Real SigmaPsi = 1.3_Real;

   Real Kmin = 7.6e-6_Real;
   Real Pmin = 1.0e-12_Real;

   Real AkvBak = 5.0e-6_Real;
   Real AktBak = 1.0e-6_Real;
   Real AkkBak = 5.0e-6_Real;
   Real AkpBak = 5.0e-6_Real;

   Real Zos = 0.02_Real; ///< Surface roughness length (m)
   Real Zob = 0.02_Real; ///< Bottom roughness length (m)

   /// Quadratic bottom-drag coefficient used only to estimate a near-bottom
   /// friction velocity for the TKE/psi bottom boundary condition (GLS does
   /// not otherwise participate in the momentum bottom-drag tendency).
   Real BottomDragCoeff = 3.0e-3_Real;

   /// Cap the turbulence length scale at LengthScaleCapFraction * local
   /// water-column depth, independent of stratification. Ported from a
   /// REMORA/ROMS bug fix: the Galperin length-scale bound only applies
   /// under stable stratification (N^2 > 0), so an unconditional hard cap is
   /// required to keep the length scale (and hence viscosity) bounded under
   /// neutral/unstable conditions. See doc/design/GLSMix.md.
   bool LengthScaleCapEnable     = true;
   Real LengthScaleCapFraction   = 0.4_Real;

   std::string VertDiffFldName;
   std::string VertViscFldName;
   std::string TkeFldName;
   std::string PsiFldName;
   std::string LscaleFldName;
   std::string Name;

 private:
   /// @brief Private constructor for singleton pattern
   GLSMix(const std::string &InName, const HorzMesh *InMesh,
         const VertCoord *InVCoord);

   /// @brief Private destructor
   ~GLSMix();

   static GLSMix *Instance;

   const HorzMesh *Mesh;
   const VertCoord *VCoord;

   GLSDerivedConstants Derived;

   /// @brief Register fields with the I/O and restart systems
   void defineFields();

   // Delete copy and move constructors/assignment
   GLSMix(const GLSMix &)            = delete;
   GLSMix &operator=(const GLSMix &) = delete;
   GLSMix(GLSMix &&)                 = delete;
   GLSMix &operator=(GLSMix &&)      = delete;

}; // class GLSMix

} // namespace OMEGA

#endif // OMEGA_GLS_MIX_H
