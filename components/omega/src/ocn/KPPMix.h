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
#include "HorzMesh.h"
#include "HorzOperators.h"
#include "KPPConstants.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"
#include "VertCoord.h"
#include <string>

namespace OMEGA {

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

   Real CriticalRichardson = 0.3; ///< Ri_crit for OBL criterion
   Real StopOBLSearchMult  = 1.0; ///< Safety multiplier for search
   Real SurfaceLayerExtent = 0.1; ///< Surface layer fraction of OBL

   bool UseLangmuirCirculation = true;  ///< Apply wave enhancement
   bool UseNonLocalFlux        = true;  ///< Apply non-local tracer flux
   bool DebugDiagnostics       = false; ///< Print per-step KPP diagnostics

   // Ice/Langmuir controls (kept configurable to match reference semantics)
   Real IceFractionThresholdForLangmuir = 0.05; ///< Disable Langmuir above this
   Real IceFractionThresholdForMinimumOBL = 0.15; ///< Apply min OBL above this
   Real MinimumOBLUnderSeaIce = 5.0; ///< Min OBL depth under sea ice (m)

   Real BackgroundVisc = 1.0e-4; ///< Background viscosity below OBL (m²/s)
   Real BackgroundDiff = 1.0e-5; ///< Background diffusivity below OBL (m²/s)

   // KPP matching/profile controls (CVMix-style semantics)
   std::string MatchTechniqueStr =
       "SimpleShapes"; ///< SimpleShapes, MatchGradient, MatchBoth,
                       ///< ParabolicNonLocal
   std::string InterpType2Str = "LMD94"; ///< Linear, Quadratic, Cubic, LMD94
   bool UseEnhancedDiffusion  = true;    ///< Apply enhanced mixing at OBL base
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
   KPPMix(const std::string &Name_in, const HorzMesh *Mesh_in,
          const VertCoord *VCoord_in);

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
