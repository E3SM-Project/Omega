#ifndef OMEGA_TENDENCYTERMS_H
#define OMEGA_TENDENCYTERMS_H
//===-- ocn/TendencyTerms.h - Tendency Terms --------------------*- C++ -*-===//
//
/// \file
/// \brief Contains functors for calculating tendency terms
///
/// This header defines functors to be called by the time-stepping scheme
/// to calculate tendencies used to update state variables.
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Eos.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "VertCoord.h"

#include <cmath> // for std::copysign

namespace OMEGA {

//------------------------------------------------------------------------------
// Shared surface flux helpers
//
// Used by the thickness, tracer, and KPP surface forcing functors so the same
// flux definitions are not restated in each.
//------------------------------------------------------------------------------
/// Surface heat flux needed by the KPP buoyancy forcing: direct fluxes plus
/// the phase-change terms, but without the enthalpy carried by mass fluxes.
///
/// @param SaTop Absolute salinity of the top layer
/// @param PTopDb Pressure at the top layer (decibar)
KOKKOS_INLINE_FUNCTION Real sfcHeatFluxWithoutMassEnthalpy(
    I4 ICell, Real SaTop, Real PTopDb, EosType EosChoice,
    const Array1DReal &LongWaveHeatFluxUp,
    const Array1DReal &LongWaveHeatFluxDown,
    const Array1DReal &ShortWaveHeatFlux, const Array1DReal &SensibleHeatFlux,
    const Array1DReal &SeaIceHeatFlux, const Array1DReal &SeaIceFreshWaterFlux,
    const Array1DReal &SeaIceSaltFlux, const Array1DReal &LatentHeatFluxEvap,
    const Array1DReal &SnowFlux, const Array1DReal &IceRunoffFlux) {

   // Enthalpy of liq water is counted in SeaIceHeatFlux(ICell)
   // Here we estimate it to remove it
   // We approximate: melt*q_icepack(Smelt) + cong*q_icepack(SSS)
   // by (melt+cong)*q_teos10(SaSeaIce) where SaSeaIce either SaTop or Smelt
   // Melting ice adds water at its own bulk salinity; freezing takes
   // water at the local surface salinity
   const Real SaSeaIce =
       SeaIceFreshWaterFlux(ICell) > 0.0_Real
           ? SeaIceSaltFlux(ICell) / SeaIceFreshWaterFlux(ICell) * Salt2PPt
           : SaTop;
   const Real SeaIceLiqEnthalpyEstimate =
       SeaIceFreshWaterFlux(ICell) * Cp0Sw *
       Eos::calcCtFreezing(EosChoice, SaSeaIce, PTopDb, 0.0_Real);

   return LongWaveHeatFluxUp(ICell) + LongWaveHeatFluxDown(ICell) +
          ShortWaveHeatFlux(ICell) + SensibleHeatFlux(ICell) +
          SeaIceHeatFlux(ICell) - SeaIceLiqEnthalpyEstimate +
          LatentHeatFluxEvap(ICell) +
          (SnowFlux(ICell) + IceRunoffFlux(ICell)) * -LatIce;
}

/// Net surface freshwater mass flux (kg/m^2/s)
KOKKOS_INLINE_FUNCTION Real sfcFreshWaterFlux(
    I4 ICell, const Array1DReal &SnowFlux, const Array1DReal &RainFlux,
    const Array1DReal &EvaporationFlux, const Array1DReal &SeaIceFreshWaterFlux,
    const Array1DReal &IceRunoffFlux, const Array1DReal &RiverRunoffFlux) {
   return SnowFlux(ICell) + RainFlux(ICell) + EvaporationFlux(ICell) +
          SeaIceFreshWaterFlux(ICell) + IceRunoffFlux(ICell) +
          RiverRunoffFlux(ICell);
}

/// Divergence of pseudo-thickness flux at cell centers, for updating
/// pseudo-thickness arrays
class PseudoThicknessFluxDivOnCell {
 public:
   bool Enabled = false;

   /// constructor declaration
   PseudoThicknessFluxDivOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes cell index, vertical chunk index, and pseudo-thickness
   /// flux array as inputs, outputs the tendency array
   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array2DReal &Tend, I4 ICell,
                                   const Array2DReal &PseudoThicknessFlux,
                                   const Array2DReal &NormalVelEdge) const {

      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      ScratchArray1DReal DivTmp(teamScratch(Team), NVertLayers);
      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { DivTmp(K) = 0; });

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);

         const int MinLyrEdgeBot = MinLayerEdgeBot(JEdge);
         const int MaxLyrEdgeTop = MaxLayerEdgeTop(JEdge);

         parallelForInner(
             Team, Range{MinLyrEdgeBot, MaxLyrEdgeTop}, INNER_LAMBDA(int K) {
                DivTmp(K) -= DvEdge(JEdge) * EdgeSignOnCell(ICell, J) *
                             PseudoThicknessFlux(JEdge, K) *
                             NormalVelEdge(JEdge, K) * InvAreaCell;
             });
      }

      const int MinLyrCell = MinLayerCell(ICell);
      const int MaxLyrCell = MaxLayerCell(ICell);

      parallelForInner(
          Team, Range{MinLyrCell, MaxLyrCell},
          INNER_LAMBDA(int K) { Tend(ICell, K) -= DivTmp(K); });
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
   Array2DReal EdgeSignOnCell;
   I4 NVertLayers;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Horizontal advection of potential vorticity defined on edges, for
/// momentum equation
class PotentialVortHAdvOnEdge {
 public:
   bool Enabled = false;

   /// constructor declaration
   PotentialVortHAdvOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and arrays for
   /// normalized relative vorticity, normalized planetary vorticity, layer
   /// thickness on edges, and normal velocity on edges as inputs,
   /// outputs the tendency array
   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array2DReal &Tend, I4 IEdge,
                                   const Array2DReal &NormRVortEdge,
                                   const Array2DReal &NormFEdge,
                                   const Array2DReal &FluxPseudoThickEdge,
                                   const Array2DReal &NormVelEdge) const {

      ScratchArray1DReal VortTmp(teamScratch(Team), NVertLayers);

      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { VortTmp(K) = 0; });

      for (int J = 0; J < NEdgesOnEdge(IEdge); ++J) {
         I4 JEdge = EdgesOnEdge(IEdge, J);
         parallelForInner(
             Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                Real NormVort =
                    (NormRVortEdge(IEdge, K) + NormFEdge(IEdge, K) +
                     NormRVortEdge(JEdge, K) + NormFEdge(JEdge, K)) *
                    0.5_Real;

                VortTmp(K) += WeightsOnEdge(IEdge, J) *
                              FluxPseudoThickEdge(JEdge, K) *
                              NormVelEdge(JEdge, K) * NormVort;
             });
      }

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             Tend(IEdge, K) += EdgeMask(IEdge, K) * VortTmp(K);
          });
   }

 private:
   Array1DI4 NEdgesOnEdge;
   Array2DI4 EdgesOnEdge;
   Array2DReal WeightsOnEdge;
   I4 NVertLayers;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Gradient of kinetic energy defined on edges, for momentum equation
class KEGradOnEdge {
 public:
   bool Enabled = false;

   /// constructor declaration
   KEGradOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and kinetic energy
   /// array as inputs, outputs the tendency array
   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array2DReal &Tend, I4 IEdge,
                                   const Array2DReal &KECell) const {

      const I4 JCell0 = CellsOnEdge(IEdge, 0);
      const I4 JCell1 = CellsOnEdge(IEdge, 1);

      const Real InvDcEdge = 1._Real / DcEdge(IEdge);

      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             Tend(IEdge, K) -= EdgeMask(IEdge, K) *
                               (KECell(JCell1, K) - KECell(JCell0, K)) *
                               InvDcEdge;
          });
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Gradient of sea surface height defined on edges multipled by gravitational
/// acceleration, for momentum equation
/// NOTE: This term is only appropriate for shallow water (Omega v0) simulations
class SSHGradOnEdge {
 public:
   bool Enabled = false;

   /// constructor declaration
   SSHGradOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and array of
   /// pseudo-thickness/SSH, outputs tendency array
   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array2DReal &Tend, I4 IEdge,
                                   const Array1DReal &SshCell) const {

      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);

      const Real InvDcEdge = 1._Real / DcEdge(IEdge);

      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             Tend(IEdge, K) -= EdgeMask(IEdge, K) * Gravity *
                               (SshCell(ICell1) - SshCell(ICell0)) * InvDcEdge;
          });
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Laplacian horizontal mixing, for momentum equation
class VelocityDiffusionOnEdge {
 public:
   bool Enabled = false;

   Real ViscDel2;

   /// constructor declaration
   VelocityDiffusionOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and arrays for
   /// divergence of horizontal velocity (defined at cell centers) and relative
   /// vorticity (defined at vertices), outputs tendency array
   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array2DReal &Tend, I4 IEdge,
                                   const Array2DReal &DivCell,
                                   const Array2DReal &RVortVertex) const {

      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);

      const I4 IVertex0 = VerticesOnEdge(IEdge, 0);
      const I4 IVertex1 = VerticesOnEdge(IEdge, 1);

      const Real DcEdgeInv = 1._Real / DcEdge(IEdge);
      const Real DvEdgeInv = 1._Real / DvEdge(IEdge);

      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             const Real Del2U =
                 ((DivCell(ICell1, K) - DivCell(ICell0, K)) * DcEdgeInv -
                  (RVortVertex(IVertex1, K) - RVortVertex(IVertex0, K)) *
                      DvEdgeInv);

             Tend(IEdge, K) +=
                 EdgeMask(IEdge, K) * ViscDel2 * MeshScalingDel2(IEdge) * Del2U;
          });
   }

 private:
   Array2DI4 CellsOnEdge;
   Array2DI4 VerticesOnEdge;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   Array1DReal MeshScalingDel2;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Biharmonic horizontal mixing, for momentum equation
class VelocityHyperDiffOnEdge {
 public:
   bool Enabled = false;

   Real ViscDel4;
   Real DivFactor;

   /// Constructor declaration
   VelocityHyperDiffOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes the edge index, vertical chunk index, and arrays for
   /// the laplacian of divergence of horizontal velocity and the laplacian of
   /// the relative vorticity, outputs tendency array
   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array2DReal &Tend, I4 IEdge,
                                   const Array2DReal &Del2DivCell,
                                   const Array2DReal &Del2RVortVertex) const {

      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);

      const I4 IVertex0 = VerticesOnEdge(IEdge, 0);
      const I4 IVertex1 = VerticesOnEdge(IEdge, 1);

      const Real DcEdgeInv = 1._Real / DcEdge(IEdge);
      const Real DvEdgeInv = 1._Real / DvEdge(IEdge);

      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             const Real Del2U =
                 (DivFactor *
                      (Del2DivCell(ICell1, K) - Del2DivCell(ICell0, K)) *
                      DcEdgeInv -
                  (Del2RVortVertex(IVertex1, K) -
                   Del2RVortVertex(IVertex0, K)) *
                      DvEdgeInv);

             Tend(IEdge, K) -=
                 EdgeMask(IEdge, K) * ViscDel4 * MeshScalingDel4(IEdge) * Del2U;
          });
   }

 private:
   Array2DI4 CellsOnEdge;
   Array2DI4 VerticesOnEdge;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   Array1DReal MeshScalingDel4;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Surface stress forcing (eg. wind)
class SfcStressForcingOnEdge {
 public:
   bool Enabled = false;

   /// constructor declaration
   SfcStressForcingOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes the edge index, vertical chunk index, and arrays for
   /// normal surface stress and edge pseudo-thickness, outputs tendency array
   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array2DReal &Tend, I4 IEdge,
                                   const Array1DReal &NormalStressEdge,
                                   const Array2DReal &PseudoThickEdge) const {
      const I4 KMin = MinLayerEdgeBot(IEdge);
      const I4 KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             if (K == KMin) {
                const Real InvThickEdge = 1._Real / PseudoThickEdge(IEdge, K);
                Tend(IEdge, K) += EdgeMask(IEdge, K) * InvThickEdge *
                                  NormalStressEdge(IEdge) / RhoSw;
             }
          });
   }

 private:
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Bottom drag
class BottomDragOnEdge {
 public:
   bool Enabled = false;
   Real Coeff;

   /// constructor declaration
   BottomDragOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes the edge index and arrays for
   /// horizontal velocity, kinetic energy,
   /// and edge pseudo-thickness, outputs tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge,
                                   const Array2DReal &NormalVelEdge,
                                   const Array2DReal &KECell,
                                   const Array2DReal &PseudoThickEdge) const {
      const I4 KBot = MaxLayerEdgeTop(IEdge);

      // Land edges and the outermost edges of the halo have no active layer
      // on both sides, and MaxLayerEdgeTop is set to -1 for them. Unlike the
      // other tendency terms, which are chunked over [MinLayerEdgeBot,
      // MaxLayerEdgeTop] and so skip such edges automatically, bottom drag
      // indexes the bottom layer directly and must exclude them explicitly.
      if (KBot < 0)
         return;

      const I4 JCell0 = CellsOnEdge(IEdge, 0);
      const I4 JCell1 = CellsOnEdge(IEdge, 1);

      const Real VelNormEdge =
          Kokkos::sqrt(KECell(JCell0, KBot) + KECell(JCell1, KBot));

      const Real InvThickEdge = 1._Real / PseudoThickEdge(IEdge, KBot);
      Tend(IEdge, KBot) -= EdgeMask(IEdge, KBot) * Coeff * VelNormEdge *
                           InvThickEdge * NormalVelEdge(IEdge, KBot);
   }

 private:
   I4 NVertLayers;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeMask;
   Array1DI4 MaxLayerEdgeTop;
};

/// Coupled freshwater flux forcing for thickness equation.
class SfcThicknessForcingOnCell {
 public:
   bool Enabled = false;

   SfcThicknessForcingOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 ICell,
                                   const Array1DReal &SnowFlux,
                                   const Array1DReal &RainFlux,
                                   const Array1DReal &EvaporationFlux,
                                   const Array1DReal &SeaIceFreshWaterFlux,
                                   const Array1DReal &IceRunoffFlux,
                                   const Array1DReal &RiverRunoffFlux,
                                   const Array1DReal &SeaIceSaltFlux) const {

      const I4 KTop = MinLayerCell(ICell);
      if (KTop > MaxLayerCell(ICell)) {
         return;
      }

      const Real FreshWaterFlux = sfcFreshWaterFlux(
          ICell, SnowFlux, RainFlux, EvaporationFlux, SeaIceFreshWaterFlux,
          IceRunoffFlux, RiverRunoffFlux);

      Tend(ICell, KTop) += (FreshWaterFlux + SeaIceSaltFlux(ICell)) / RhoSw;
   }

 private:
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

/// Coupled surface flux forcing for active tracers.
class SfcTracerForcingOnCell {
 public:
   bool Enabled = false;

   SfcTracerForcingOnCell(const HorzMesh *Mesh, const VertCoord *VCoord,
                          I4 TempTracerIndex, I4 SaltTracerIndex,
                          const Eos *EosInst);

   KOKKOS_FUNCTION void
   operator()(const Array3DReal &Tend, I4 ICell, const Array3DReal &TracerCell,
              const Array1DReal &LatentHeatFluxEvap,
              const Array1DReal &SensibleHeatFlux,
              const Array1DReal &LongWaveHeatFluxUp,
              const Array1DReal &LongWaveHeatFluxDown,
              const Array1DReal &SeaIceHeatFlux,
              const Array1DReal &ShortWaveHeatFlux, const Array1DReal &SnowFlux,
              const Array1DReal &RainFlux, const Array1DReal &IceRunoffFlux,
              const Array1DReal &RiverRunoffFlux,
              const Array1DReal &EvaporationFlux,
              const Array1DReal &SeaIceSaltFlux) const {

      const I4 KTop = MinLayerCell(ICell);
      if (KTop > MaxLayerCell(ICell)) {
         return;
      }

      if (TempIndex >= 0) {
         const Real CtTop = TracerCell(TempIndex, ICell, KTop);

         // CT tendencies are due to direct heat fluxes + pot enthalpy fluxes
         // Each mass flux has an associated potential enthalpy flux.
         // Levels of simplification can be done here. For now:
         // - We approximate PotEnthalpyIce(Tinsitu, P=0) ~ -LatIce (constant);
         // Altrntively, we could use cnst PotEnthalpyIce(0, 0) ​= −333360
         // J/kg a 0.1% / 340 J/kg difference with the LatIce value from pcd.
         // The full expression is gsw_pot_enthalpy_ice(T, P).
         const Real PotEnthalpyIce = -LatIce;
         // - We assume dry snow and use the same PotEnthalpyIce for snow.
         // - We assume liquid water comes in at the specific enthalpy as the
         // top ocean layer, (i.e. mass flux is CT-neutral). The enthalpy is
         // capped by a lower bound of pot enthalpy of freshwater at 0.0 C.
         // Another choice would be to add it at the same in situ temp as ocean
         // i.e. CT(Sa=0, max(0, T)).
         //- We assume evaporation removes the same specific enthalpy as the
         // top ocean layer; not capped, to keep the mass flux CT-neutral.
         const Real CtLim =
             (EosChoice == EosType::Teos10Eos) ? Ct0Fw : 0.0_Real;
         const Real PotEnthalpyFwIn  = Cp0Sw * Kokkos::max(CtLim, CtTop);
         const Real PotEnthalpyFwOut = Cp0Sw * CtTop;

         const Real HeatFlux =
             LongWaveHeatFluxUp(ICell) + LongWaveHeatFluxDown(ICell) +
             ShortWaveHeatFlux(ICell) + SensibleHeatFlux(ICell) +
             SeaIceHeatFlux(ICell) + // includes enthalpy of meltwater already
             (RainFlux(ICell) + RiverRunoffFlux(ICell)) * PotEnthalpyFwIn +
             LatentHeatFluxEvap(ICell) +
             EvaporationFlux(ICell) * PotEnthalpyFwOut +
             (SnowFlux(ICell) + IceRunoffFlux(ICell)) * PotEnthalpyIce;

         Tend(TempIndex, ICell, KTop) += HeatFlux * HFluxFac;
      }

      if (SaltIndex >= 0) {
         Tend(SaltIndex, ICell, KTop) += SeaIceSaltFlux(ICell) * SFluxFac;
      }
   }

 private:
   I4 TempIndex;
   I4 SaltIndex;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   EosType EosChoice;
};

/// Potential density referenced to the surface, used by the KPP boundary
/// layer depth search.
class PotentialDensityOnCell {
 public:
   bool Enabled = false;

   PotentialDensityOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// Fills the reference pressure used for the displaced specific volume:
   /// every layer in a column is referenced to that column's surface pressure
   KOKKOS_FUNCTION void
   computeRefPressure(const Array2DReal &RefPressure, I4 ICell, I4 K,
                      const Array2DReal &PressureMid) const {
      RefPressure(ICell, K) = PressureMid(ICell, MinLayerCell(ICell));
   }

   /// Inverts the surface-referenced specific volume to give potential density
   KOKKOS_FUNCTION void operator()(const Array2DReal &PotentialDensity,
                                   I4 ICell, I4 K,
                                   const Array2DReal &SpecVolDisplaced) const {
      PotentialDensity(ICell, K) =
          1.0_Real / Kokkos::max(1.0e-12_Real, SpecVolDisplaced(ICell, K));
   }

 private:
   Array1DI4 MinLayerCell;
};

/// Surface forcing inputs consumed by the KPP boundary layer scheme: friction
/// velocity, buoyancy flux, and the surface tracer fluxes that scale the
/// non-local term.
class KPPSurfaceForcingOnCell {
 public:
   bool Enabled = false;

   /// Store surface tracer fluxes for the KPP non-local tracer tendency
   bool UpdateSurfaceTracerFlux = false;

   /// When false no coupled tracer forcing is active, so only the friction
   /// velocity is set and the buoyancy and tracer fluxes stay zero
   bool UseTracerForcing = false;

   KPPSurfaceForcingOnCell(const HorzMesh *Mesh, const VertCoord *VCoord,
                           I4 TempTracerIndex, I4 SaltTracerIndex,
                           const Eos *EosInst);

   KOKKOS_FUNCTION void operator()(
       const Array1DReal &FrictionVelocity, const Array1DReal &BuoyancyFlux,
       const Array2DReal &SurfaceTracerFlux, const Array1DReal &IceFraction,
       I4 ICell, const Array2DReal &ConservTemp, const Array2DReal &AbsSalinity,
       const Array2DReal &PressureMid, const Array2DReal &SpecVol,
       const Array1DReal &ZonalStress, const Array1DReal &MeridStress,
       const Array1DReal &LatentHeatFluxEvap,
       const Array1DReal &SensibleHeatFlux,
       const Array1DReal &LongWaveHeatFluxUp,
       const Array1DReal &LongWaveHeatFluxDown,
       const Array1DReal &SeaIceHeatFlux, const Array1DReal &ShortWaveHeatFlux,
       const Array1DReal &SnowFlux, const Array1DReal &RainFlux,
       const Array1DReal &EvaporationFlux,
       const Array1DReal &SeaIceFreshWaterFlux,
       const Array1DReal &IceRunoffFlux, const Array1DReal &RiverRunoffFlux,
       const Array1DReal &SeaIceSaltFlux) const {

      const Real TauX   = ZonalStress(ICell);
      const Real TauY   = MeridStress(ICell);
      const Real TauMag = Kokkos::sqrt(TauX * TauX + TauY * TauY);

      FrictionVelocity(ICell) =
          Kokkos::sqrt(Kokkos::max(0.0_Real, TauMag / RhoSw));
      BuoyancyFlux(ICell) = 0.0_Real;
      if (UpdateSurfaceTracerFlux) {
         SurfaceTracerFlux(TempIndex, ICell) = 0.0_Real;
         SurfaceTracerFlux(SaltIndex, ICell) = 0.0_Real;
      }
      // Sea ice coupling is not wired in yet, so KPP sees an ice-free ocean
      IceFraction(ICell) = 0.0_Real;

      if (!UseTracerForcing) {
         return;
      }

      const I4 KTop     = MinLayerCell(ICell);
      const Real SaTop  = AbsSalinity(ICell, KTop);
      const Real CtTop  = ConservTemp(ICell, KTop);
      const Real PTopDb = PressureMid(ICell, KTop) * Pa2Db;

      // KPP needs the surface heat flux without the mass-flux enthalpy that
      // SfcTracerForcingOnCell carries in the hT tendency
      const Real HeatFlux = sfcHeatFluxWithoutMassEnthalpy(
          ICell, SaTop, PTopDb, EosChoice, LongWaveHeatFluxUp,
          LongWaveHeatFluxDown, ShortWaveHeatFlux, SensibleHeatFlux,
          SeaIceHeatFlux, SeaIceFreshWaterFlux, SeaIceSaltFlux,
          LatentHeatFluxEvap, SnowFlux, IceRunoffFlux);
      const Real FreshWaterFlux = sfcFreshWaterFlux(
          ICell, SnowFlux, RainFlux, EvaporationFlux, SeaIceFreshWaterFlux,
          IceRunoffFlux, RiverRunoffFlux);

      const Real TempFlux = HeatFlux * HFluxFac;
      const Real SaltFlux =
          SeaIceSaltFlux(ICell) * SFluxFac - FreshWaterFlux * SaTop / RhoSw;

      const Real SpVol  = Kokkos::max(1.0e-12_Real, SpecVol(ICell, KTop));
      const Real RhoTop = 1.0_Real / SpVol;

      Real Alpha = 0.0_Real;
      Real Beta  = 0.0_Real;
      if (EosChoice == EosType::Teos10Eos) {
         Alpha = Teos10Coeff.calcAlpha(SaTop, CtTop, PTopDb, SpVol);
         Beta  = Teos10Coeff.calcBeta(SaTop, CtTop, PTopDb, SpVol);
      } else if (EosChoice == EosType::LinearEos) {
         Alpha = -LinearDRhodT / RhoTop;
         Beta  = LinearDRhodS / RhoTop;
      }

      BuoyancyFlux(ICell) = Gravity * (Alpha * TempFlux - Beta * SaltFlux);
      if (UpdateSurfaceTracerFlux) {
         SurfaceTracerFlux(TempIndex, ICell) = TempFlux;
         SurfaceTracerFlux(SaltIndex, ICell) = SaltFlux;
      }
   }

 private:
   I4 TempIndex;
   I4 SaltIndex;
   Real LinearDRhodT;
   Real LinearDRhodS;
   Array1DI4 MinLayerCell;
   EosType EosChoice;
   Teos10BruntVaisalaFreqSq Teos10Coeff;
};

// Tracer horizontal advection term
class TracerHorzAdvOnCell {
 public:
   bool Enabled       = false;
   bool ForceLowOrder = false;
   // coefficient for blending high-order terms
   Real Coef3rdOrder = 0.25;
   TracerHorzAdvOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);
   void init();
   KOKKOS_FUNCTION void operator()(const TeamMember &Team, const I4 L,
                                   const I4 IEdge,
                                   const Array3DReal &TracerCell,
                                   const Array2DReal &FluxPseudoThickEdge,
                                   const Array2DReal &NormVelEdge) const {

      ScratchArray1DReal FlxTmp(teamScratch(Team), NVertLayers);

      const auto LTracerCell =
          subviewUnmanaged(TracerCell, L, Kokkos::ALL, Kokkos::ALL);

      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { FlxTmp(K) = 0; });

      // Stay at low order at boundaries
      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             const I4 JCell0 = CellsOnEdge(IEdge, 0);
             const I4 JCell1 = CellsOnEdge(IEdge, 1);
             const Real NormalThicknessFlux =
                 FluxPseudoThickEdge(IEdge, K) * NormVelEdge(IEdge, K);
             const Real TracerWgt =
                 DvEdge(IEdge) * 0.5_Real * NormalThicknessFlux;
             FlxTmp(K) += TracerWgt * (1._Real - AdvMaskHighOrder(IEdge, K)) *
                          (LTracerCell(JCell1, K) + LTracerCell(JCell0, K));
          });

      // High order (3rd or 4th) fluxes elsewhere when requested
      //    - If HorzTracerFluxOrder = 2, NAdvCellsForEdge = 0 and
      //      this loop is skipped.
      for (int I = 0; I < NAdvCellsForEdge(IEdge); ++I) {
         const I4 ICell = AdvCellsForEdge(IEdge, I);
         parallelForInner(
             Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                const Real NormalThicknessFlux =
                    FluxPseudoThickEdge(IEdge, K) * NormVelEdge(IEdge, K);
                const Real TracerWgt =
                    (AdvCoefs(I, IEdge) +
                     Coef3rdOrder *
                         std::copysign(1._Real, NormalThicknessFlux) *
                         AdvCoefs3rd(I, IEdge)) *
                    NormalThicknessFlux;
                FlxTmp(K) += TracerWgt * LTracerCell(ICell, K) *
                             AdvMaskHighOrder(IEdge, K);
             });
      }

      const auto LHighOrderFlxHorz =
          subviewUnmanaged(HighOrderFlxHorz, L, Kokkos::ALL, Kokkos::ALL);

      parallelForInner(
          Team, Range{KMin, KMax},
          INNER_LAMBDA(int K) { LHighOrderFlxHorz(IEdge, K) = FlxTmp(K); });
   }

   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array3DReal &Tend, const I4 L,
                                   const I4 ICell) const {

      const auto LTend = subviewUnmanaged(Tend, L, Kokkos::ALL, Kokkos::ALL);
      const auto LHighOrderFlxHorz =
          subviewUnmanaged(HighOrderFlxHorz, L, Kokkos::ALL, Kokkos::ALL);

      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      ScratchArray1DReal TendTmp(teamScratch(Team), NVertLayers);
      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { TendTmp(K) = 0; });

      for (int I = 0; I < NEdgesOnCell(ICell); ++I) {
         const I4 IEdge = EdgesOnCell(ICell, I);

         const int MinLyrEdgeBot = MinLayerEdgeBot(IEdge);
         const int MaxLyrEdgeTop = MaxLayerEdgeTop(IEdge);

         parallelForInner(
             Team, Range{MinLyrEdgeBot, MaxLyrEdgeTop}, INNER_LAMBDA(int K) {
                TendTmp(K) += EdgeSignOnCell(ICell, I) *
                              LHighOrderFlxHorz(IEdge, K) * InvAreaCell;
             });
      }

      const int MinLyrCell = MinLayerCell(ICell);
      const int MaxLyrCell = MaxLayerCell(ICell);

      parallelForInner(
          Team, Range{MinLyrCell, MaxLyrCell},
          INNER_LAMBDA(int K) { LTend(ICell, K) += TendTmp(K); });
   }

 private:
   const HorzMesh *HorzontalMesh;
   const VertCoord *VerticalCoord;
   Array1DI4 NAdvCellsForEdge;
   Array2DI4 AdvCellsForEdge;
   Array2DI4 AdvMaskHighOrder;
   Array2DReal AdvCoefs;
   Array2DReal AdvCoefs3rd;
   Array3DReal HighOrderFlxHorz;

   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeSignOnCell;
   Array1DReal DvEdge;
   Array1DReal AreaCell;

   I4 NVertLayers;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

// Tracer horizontal diffusion term
class TracerDiffOnCell {
 public:
   bool Enabled = false;

   Real EddyDiff2;

   TracerDiffOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   operator()(const TeamMember &Team, const Array3DReal &Tend, I4 L, I4 ICell,
              const Array3DReal &TracerCell,
              const Array2DReal &MeanPseudoThickEdge) const {

      const auto LTend = subviewUnmanaged(Tend, L, Kokkos::ALL, Kokkos::ALL);
      const auto LTracerCell =
          subviewUnmanaged(TracerCell, L, Kokkos::ALL, Kokkos::ALL);

      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      ScratchArray1DReal DiffTmp(teamScratch(Team), NVertLayers);
      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { DiffTmp(K) = 0; });

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);

         const I4 JCell0 = CellsOnEdge(JEdge, 0);
         const I4 JCell1 = CellsOnEdge(JEdge, 1);

         const Real RTemp =
             MeshScalingDel2(JEdge) * DvEdge(JEdge) / DcEdge(JEdge);

         const int MinLyrEdgeBot = MinLayerEdgeBot(JEdge);
         const int MaxLyrEdgeTop = MaxLayerEdgeTop(JEdge);

         parallelForInner(
             Team, Range{MinLyrEdgeBot, MaxLyrEdgeTop}, INNER_LAMBDA(int K) {
                const Real TracerGrad =
                    (LTracerCell(JCell1, K) - LTracerCell(JCell0, K));

                DiffTmp(K) -= EdgeMask(JEdge, K) * EdgeSignOnCell(ICell, J) *
                              RTemp * MeanPseudoThickEdge(JEdge, K) *
                              TracerGrad;
             });
      }
      const int MinLyrCell = MinLayerCell(ICell);
      const int MaxLyrCell = MaxLayerCell(ICell);

      parallelForInner(
          Team, Range{MinLyrCell, MaxLyrCell}, INNER_LAMBDA(int K) {
             LTend(ICell, K) += EddyDiff2 * DiffTmp(K) * InvAreaCell;
          });
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeSignOnCell;
   Array1DReal DvEdge;
   Array1DReal DcEdge;
   Array1DReal AreaCell;
   Array1DReal MeshScalingDel2;
   I4 NVertLayers;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

// Tracer biharmonic horizontal mixing term
class TracerHyperDiffOnCell {
 public:
   bool Enabled = false;

   Real EddyDiff4;

   TracerHyperDiffOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array3DReal &Tend, I4 L, I4 ICell,
                                   const Array3DReal &TrDel2Cell) const {

      const auto LTend = subviewUnmanaged(Tend, L, Kokkos::ALL, Kokkos::ALL);
      const auto LTrDel2Cell =
          subviewUnmanaged(TrDel2Cell, L, Kokkos::ALL, Kokkos::ALL);

      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      ScratchArray1DReal HypTmp(teamScratch(Team), NVertLayers);
      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { HypTmp(K) = 0; });

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);

         const I4 JCell0 = CellsOnEdge(JEdge, 0);
         const I4 JCell1 = CellsOnEdge(JEdge, 1);

         const Real RTemp =
             MeshScalingDel4(JEdge) * DvEdge(JEdge) / DcEdge(JEdge);

         const int MinLyrEdgeBot = MinLayerEdgeBot(JEdge);
         const int MaxLyrEdgeTop = MaxLayerEdgeTop(JEdge);

         parallelForInner(
             Team, Range{MinLyrEdgeBot, MaxLyrEdgeTop}, INNER_LAMBDA(int K) {
                const Real Del2TrGrad =
                    (LTrDel2Cell(JCell1, K) - LTrDel2Cell(JCell0, K));

                HypTmp(K) -= EdgeMask(JEdge, K) * EdgeSignOnCell(ICell, J) *
                             RTemp * Del2TrGrad;
             });
      }
      const int MinLyrCell = MinLayerCell(ICell);
      const int MaxLyrCell = MaxLayerCell(ICell);

      parallelForInner(
          Team, Range{MinLyrCell, MaxLyrCell}, INNER_LAMBDA(int K) {
             LTend(ICell, K) -= EddyDiff4 * HypTmp(K) * InvAreaCell;
          });
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeSignOnCell;
   Array1DReal DvEdge;
   Array1DReal DcEdge;
   Array1DReal AreaCell;
   Array1DReal MeshScalingDel4;
   I4 NVertLayers;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Surface tracer restoring term
class SurfaceTracerRestoringOnCell {
 public:
   bool Enabled;
   Real PistonVelocity  = 1.585e-5; ///< piston velocity
   I4 NTracersToRestore = 0;        ///< number of tracers to restore
   Array1DI4 TracerIdsToRestore;    ///< tracer IDs to restore
   /// Need to add under sea ice restoring option when that is available

   /// constructor declaration
   SurfaceTracerRestoringOnCell(const HorzMesh *Mesh);

   /// The functor takes the cell index and the array for the tracer surface
   /// restoring values, outputs tendency array
   KOKKOS_FUNCTION void
   operator()(const Array3DReal &Tend, I4 L, I4 ICell, I4 KMin,
              const Array2DReal &TracersMonthlySurfClimoCell,
              const Array3DReal &TracerCell) const {

      Tend(L, ICell, KMin) +=
          PistonVelocity *
          (TracersMonthlySurfClimoCell(L, ICell) - TracerCell(L, ICell, KMin));
   }
};

} // namespace OMEGA
#endif
