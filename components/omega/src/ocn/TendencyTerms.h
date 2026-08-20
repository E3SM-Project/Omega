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
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "VertCoord.h"

#include <cmath> // for std::copysign
#include <iomanip>

namespace OMEGA {

/// Divergence of pseudo-thickness flux at cell centers, for updating
/// pseudo-thickness arrays
class PseudoThicknessFluxDivOnCell {
 public:
   bool Enabled = false;

   /// constructor declaration
   PseudoThicknessFluxDivOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes cell index, vertical chunk index, and pseudo-thickness
   /// flux array as inputs, outputs the tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 ICell, I4 KChunk,
                                   const Array2DReal &PseudoThicknessFlux,
                                   const Array2DReal &NormalVelEdge) const {

      const I4 KStartCell = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLenCell = chunkLength(KChunk, KStartCell, MaxLayerCell(ICell));
      const I4 KEndCell = KStartCell + KLenCell - 1;
      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real DivTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);

         const I4 KStartEdge = Kokkos::max(KStartCell, MinLayerEdgeBot(JEdge));
         const I4 KEndEdge   = Kokkos::min(KEndCell, MaxLayerEdgeTop(JEdge));

         for (int K = KStartEdge; K <= KEndEdge; ++K) {
            const I4 KVec = K - KStartCell;
            DivTmp[KVec] -= DvEdge(JEdge) * EdgeSignOnCell(ICell, J) *
                            PseudoThicknessFlux(JEdge, K) *
                            NormalVelEdge(JEdge, K) * InvAreaCell;
         }
      }

      for (int KVec = 0; KVec < KLenCell; ++KVec) {
         const I4 K = KStartCell + KVec;
         Tend(ICell, K) -= DivTmp[KVec];
      }
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
   Array2DReal EdgeSignOnCell;
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
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &NormRVortEdge,
                                   const Array2DReal &NormFEdge,
                                   const Array2DReal &FluxPseudoThickEdge,
                                   const Array2DReal &NormVelEdge) const {

      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));
      Real VortTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnEdge(IEdge); ++J) {
         I4 JEdge = EdgesOnEdge(IEdge, J);
         for (int KVec = 0; KVec < KLen; ++KVec) {
            const I4 K    = KStart + KVec;
            Real NormVort = (NormRVortEdge(IEdge, K) + NormFEdge(IEdge, K) +
                             NormRVortEdge(JEdge, K) + NormFEdge(JEdge, K)) *
                            0.5_Real;

            VortTmp[KVec] += WeightsOnEdge(IEdge, J) *
                             FluxPseudoThickEdge(JEdge, K) *
                             NormVelEdge(JEdge, K) * NormVort;
         }
      }

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(IEdge, K) += EdgeMask(IEdge, K) * VortTmp[KVec];
      }
   }

 private:
   Array1DI4 NEdgesOnEdge;
   Array2DI4 EdgesOnEdge;
   Array2DReal WeightsOnEdge;
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
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &KECell) const {

      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));
      const I4 JCell0 = CellsOnEdge(IEdge, 0);
      const I4 JCell1 = CellsOnEdge(IEdge, 1);
      const Real InvDcEdge = 1._Real / DcEdge(IEdge);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(IEdge, K) -= EdgeMask(IEdge, K) *
                           (KECell(JCell1, K) - KECell(JCell0, K)) * InvDcEdge;
      }
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
class SSHGradOnEdge {
 public:
   bool Enabled = false;

   /// constructor declaration
   SSHGradOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and array of
   /// pseudo-thickness/SSH, outputs tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array1DReal &SshCell) const {

      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));
      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);
      const Real InvDcEdge = 1._Real / DcEdge(IEdge);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(IEdge, K) -= EdgeMask(IEdge, K) * Gravity *
                           (SshCell(ICell1) - SshCell(ICell0)) * InvDcEdge;
      }
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
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &DivCell,
                                   const Array2DReal &RVortVertex) const {

      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));
      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);

      const I4 IVertex0 = VerticesOnEdge(IEdge, 0);
      const I4 IVertex1 = VerticesOnEdge(IEdge, 1);

      const Real DcEdgeInv = 1._Real / DcEdge(IEdge);
      const Real DvEdgeInv = 1._Real / DvEdge(IEdge);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         const Real Del2U =
             ((DivCell(ICell1, K) - DivCell(ICell0, K)) * DcEdgeInv -
              (RVortVertex(IVertex1, K) - RVortVertex(IVertex0, K)) *
                  DvEdgeInv);

         Tend(IEdge, K) +=
             EdgeMask(IEdge, K) * ViscDel2 * MeshScalingDel2(IEdge) * Del2U;
      }
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
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &Del2DivCell,
                                   const Array2DReal &Del2RVortVertex) const {

      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));
      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);

      const I4 IVertex0 = VerticesOnEdge(IEdge, 0);
      const I4 IVertex1 = VerticesOnEdge(IEdge, 1);

      const Real DcEdgeInv = 1._Real / DcEdge(IEdge);
      const Real DvEdgeInv = 1._Real / DvEdge(IEdge);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         const Real Del2U =
             (DivFactor * (Del2DivCell(ICell1, K) - Del2DivCell(ICell0, K)) *
                  DcEdgeInv -
              (Del2RVortVertex(IVertex1, K) - Del2RVortVertex(IVertex0, K)) *
                  DvEdgeInv);

         Tend(IEdge, K) -=
             EdgeMask(IEdge, K) * ViscDel4 * MeshScalingDel4(IEdge) * Del2U;
      }
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
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array1DReal &NormalStressEdge,
                                   const Array2DReal &PseudoThickEdge) const {
      if (KChunk == 0) {
         const I4 K = MinLayerEdgeBot(IEdge);

         const Real InvThickEdge = 1._Real / PseudoThickEdge(IEdge, K);
         Tend(IEdge, K) += EdgeMask(IEdge, K) * InvThickEdge *
                           NormalStressEdge(IEdge) / RhoSw;
      }
   }

 private:
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
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

// Tracer horizontal advection term
class TracerHorzAdvOnCell {
 public:
   bool Enabled           = true;
   bool ForceLowOrder     = false;
   bool FCT               = false;
   bool ComputeBudgets    = false;
   bool MonotonicityCheck = false;
   // coefficient for blending high-order terms
   Real Coef3rdOrder = 0.25;
   TracerHorzAdvOnCell(const HorzMesh *Mesh, const VertCoord *VCoord,
                       const VertAdv *VAdv);
   KOKKOS_FUNCTION ~TracerHorzAdvOnCell() {}
   void init();
   KOKKOS_FUNCTION void operator()(const I4 L, const I4 IEdge, const I4 K,
                                   const Array3DReal &TracerCell,
                                   const Array2DReal &FluxPseudoThickEdge,
                                   const Array2DReal &NormVelEdge) const {
      HighOrderFlxHorz(L, IEdge, K) = 0;

      // Stay at low order at boundaries
      const I4 JCell0 = CellsOnEdge(IEdge, 0);
      const I4 JCell1 = CellsOnEdge(IEdge, 1);
      const Real NormalThicknessFlux =
          FluxPseudoThickEdge(IEdge, K) * NormVelEdge(IEdge, K);
      const Real TracerWgt = DvEdge(IEdge) * 0.5_Real * NormalThicknessFlux;
      HighOrderFlxHorz(L, IEdge, K) +=
          TracerWgt * (1._Real - AdvMaskHighOrder(IEdge, K)) *
          (TracerCell(L, JCell1, K) + TracerCell(L, JCell0, K));

      // High order (3rd or 4th) fluxes elsewhere when requested
      //    - If HorzTracerFluxOrder = 2, NAdvCellsForEdge = 0 and
      //      this loop is skipped.
      for (int I = 0; I < NAdvCellsForEdge(IEdge); ++I) {
         const I4 ICell = AdvCellsForEdge(IEdge, I);
         const Real NormalThicknessFlux =
             FluxPseudoThickEdge(IEdge, K) * NormVelEdge(IEdge, K);
         const Real TracerWgt =
             (AdvCoefs(I, IEdge) +
              Coef3rdOrder * std::copysign(1._Real, NormalThicknessFlux) *
                  AdvCoefs3rd(I, IEdge)) *
             NormalThicknessFlux;
         HighOrderFlxHorz(L, IEdge, K) +=
             TracerWgt * TracerCell(L, ICell, K) * AdvMaskHighOrder(IEdge, K);
      }
   }

   KOKKOS_FUNCTION void operator()(const Array3DReal &Tend, const I4 L,
                                   const I4 ICell, const I4 K) const {
      const Real InvAreaCell = 1._Real / AreaCell(ICell);
      Tend(L, ICell, K)      = 0;

      for (int I = 0; I < NEdgesOnCell(ICell); ++I) {
         const I4 IEdge = EdgesOnCell(ICell, I);
         Tend(L, ICell, K) += EdgeSignOnCell(ICell, I) *
                              HighOrderFlxHorz(L, IEdge, K) * InvAreaCell;
      }
   }

   KOKKOS_FUNCTION void
   FCTProvisionaLayerThicknesses(const TeamMember &Team, const I4 ICell,
                                 const Real Dt,
                                 const Array2DReal &FluxPseudoThickEdge,
                                 const Array2DReal &LayerThickness,
                                 const Array2DReal &NormVelEdge) const {
      const int KMin = MinLayerCell(ICell);
      const int KMax = MaxLayerCell(ICell);
      parallelForInner(
          Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
             FCTProvisionaLayerThicknesses(ICell, K, Dt, FluxPseudoThickEdge,
                                           LayerThickness, NormVelEdge);
          });
   }

   KOKKOS_FUNCTION void
   FCTProvisionaLayerThicknesses(const I4 ICell, const I4 K, const Real Dt,
                                 const Array2DReal &FluxPseudoThickEdge,
                                 const Array2DReal &LayerThickness,
                                 const Array2DReal &NormVelEdge) const {

      const Real InvAreaCell = Dt / AreaCell(ICell);
      HProv(ICell, K)        = LayerThickness(ICell, K);
      for (I4 I = 0; I < NEdgesOnCell(ICell); ++I) {
         const I4 IEdge = EdgesOnCell(ICell, I);
         const Real SignedFactor =
             InvAreaCell * DvEdge(IEdge) * EdgeSignOnCell(ICell, I);
         // Provisional layer thickness is after horizontal
         // thickness flux only
         const Real NormalThicknessFlux =
             FluxPseudoThickEdge(IEdge, K) * NormVelEdge(IEdge, K);
         HProv(ICell, K) += SignedFactor * NormalThicknessFlux;
      }
      // New layer thickness is after horizontal and vertical
      // thickness flux
      HProvInv(ICell, K) = 1.0_Real / HProv(ICell, K);
      HNewInv(ICell, K) =
          1.0_Real / (HProv(ICell, K) - Dt * VerticalPseudoVelocity(ICell, K) +
                      Dt * VerticalPseudoVelocity(ICell, K + 1));
   }

   KOKKOS_FUNCTION void FCTTracerCurFill(const TeamMember &Team, const I4 L,
                                         const I4 ICell,
                                         const Array3DReal &TracerArray) const {
      const int KMin = MinLayerCell(ICell);
      const int KMax = MaxLayerCell(ICell);
      parallelForInner(
          Team, Range(KMin, KMax),
          INNER_LAMBDA(int K) { FCTTracerCurFill(L, ICell, K, TracerArray); });
   }

   KOKKOS_FUNCTION void FCTTracerCurFill(const I4 L, const I4 ICell, const I4 K,
                                         const Array3DReal &TracerArray) const {
      TracerCur(ICell, K) = TracerArray(L, ICell, K);
      TracerMin(ICell, K) = TracerArray(L, ICell, K);
      TracerMax(ICell, K) = TracerArray(L, ICell, K);
   }
   KOKKOS_FUNCTION void FCTTracerMinMax(const TeamMember &Team,
                                        const I4 ICell) const {
      const I4 KMin   = MinLayerCell(ICell);
      const I4 KMax   = MaxLayerCell(ICell);
      const I4 NEdges = NEdgesOnCell(ICell);
      for (I4 I = 0; I < NEdges; ++I) {
         const I4 ICell2 = CellsOnCell(ICell, I);
         const I4 KMin1  = Kokkos::max(KMin, MinLayerCell(ICell2));
         const I4 KMax1  = Kokkos::min(KMax, MaxLayerCell(ICell2));
         parallelForInner(
             Team, Range(KMin1, KMax1), INNER_LAMBDA(int K) {
                TracerMax(ICell, K) =
                    Kokkos::max(TracerMax(ICell, K), TracerCur(ICell2, K));
                TracerMin(ICell, K) =
                    Kokkos::min(TracerMin(ICell, K), TracerCur(ICell2, K));
             });
      }
   }

   KOKKOS_FUNCTION void
   FCTHighAndLowOrderFlux(const TeamMember &Team, const I4 IEdge,
                          const Array2DReal &FluxPseudoThickEdge,
                          const Array2DReal &NormVelEdge) const {

      const I4 ICell1 = CellsOnEdge(IEdge, 0);
      const I4 ICell2 = CellsOnEdge(IEdge, 1);

      const I4 KMin = MinLayerEdgeBot(IEdge);
      const I4 KMax = MaxLayerEdgeTop(IEdge);
      parallelForInner(
          Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
             HighOrderFlx(IEdge, K) = 0;
             LowOrderFlx(IEdge, K)  = 0;
          });
      // Compute 3rd or 4th fluxes where requested.
      for (int I = 0; I < NAdvCellsForEdge(IEdge); ++I) {
         const I4 ICell   = AdvCellsForEdge(IEdge, I);
         const Real Coef1 = AdvCoefs(I, IEdge);
         const Real Coef3 = AdvCoefs3rd(I, IEdge) * Coef3rdOrder;
         const I4 KMin    = MinLayerCell(ICell);
         const I4 KMax    = MaxLayerCell(ICell);
         parallelForInner(
             Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
                const Real NormalThicknessFlux =
                    FluxPseudoThickEdge(IEdge, K) * NormVelEdge(IEdge, K);
                const Real TracerWgt =
                    NormalThicknessFlux *
                    (Coef1 +
                     Coef3 * std::copysign(1.0_Real, NormalThicknessFlux));
                HighOrderFlx(IEdge, K) += TracerWgt * TracerCur(ICell, K) *
                                          AdvMaskHighOrder(IEdge, K);
             });
      }
      // Compute 2nd order fluxes where needed.
      // Also compute low order upwind horizontal flux (monotonic)
      // Remove low order flux from the high order flux
      // Store left over high order flux in highOrderFlx array
      parallelForInner(
          Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
             const Real NormalThicknessFlux =
                 FluxPseudoThickEdge(IEdge, K) * NormVelEdge(IEdge, K);
             const Real TracerWeight = (1.0_Real - AdvMaskHighOrder(IEdge, K)) *
                                       (DvEdge(IEdge) * 0.5_Real) *
                                       NormalThicknessFlux;
             LowOrderFlx(IEdge, K) =
                 DvEdge(IEdge) * (Kokkos::max(0.0_Real, NormalThicknessFlux) *
                                      TracerCur(ICell1, K) +
                                  Kokkos::min(0.0_Real, NormalThicknessFlux) *
                                      TracerCur(ICell2, K));
             HighOrderFlx(IEdge, K) +=
                 TracerWeight * (TracerCur(ICell1, K) + TracerCur(ICell2, K));
             HighOrderFlx(IEdge, K) -= LowOrderFlx(IEdge, K);
          });
   }

   KOKKOS_FUNCTION void FCTFluxInOut(const TeamMember &Team, const I4 ICell,
                                     const Real Dt,
                                     const Array2DReal &LayerThickness) const {
      const I4 KMin          = MinLayerCell(ICell);
      const I4 KMax          = MaxLayerCell(ICell);
      const Real InvAreaCell = 1._Real / AreaCell(ICell);
      parallelForInner(
          Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
             WorkTend(ICell, K) = 0;
             FlxIn(ICell, K)    = 0;
             FlxOut(ICell, K)   = 0;
          });
      // Finish computing the low order horizontal fluxes
      // Upwind fluxes are accumulated in workTend
      for (I4 I = 0; I < NEdgesOnCell(ICell); ++I) {
         const I4 IEdge          = EdgesOnCell(ICell, I);
         const Real SignedFactor = EdgeSignOnCell(ICell, I) * InvAreaCell;
         const I4 KMin           = MinLayerEdgeBot(IEdge);
         const I4 KMax           = MaxLayerEdgeTop(IEdge);
         parallelForInner(
             Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
                // Here workTend is the advection tendency due to the
                // upwind (low order) fluxes.
                WorkTend(ICell, K) += SignedFactor * LowOrderFlx(IEdge, K);
                // Accumulate remaining high order fluxes
                FlxOut(ICell, K) += Kokkos::min(
                    0.0_Real, SignedFactor * HighOrderFlx(IEdge, K));
                FlxIn(ICell, K) += Kokkos::max(
                    0.0_Real, SignedFactor * HighOrderFlx(IEdge, K));
             });
      }
      // Build the factors for the FCT
      // Computed using the bounds that were computed previously,
      // and the bounds on the newly updated value
      // Factors are placed in the flxIn and flxOut arrays
      parallelForInner(
          Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
             // Here workTend is the upwind tendency
             const Real TracerUpwindNew =
                 (TracerCur(ICell, K) * LayerThickness(ICell, K) +
                  Dt * WorkTend(ICell, K)) *
                 HProvInv(ICell, K);
             const Real TracerMinNew =
                 TracerUpwindNew + Dt * FlxOut(ICell, K) * HProvInv(ICell, K);
             const Real TracerMaxNew =
                 TracerUpwindNew + Dt * FlxIn(ICell, K) * HProvInv(ICell, K);
             const Real ScaleFactorIn =
                 (TracerMax(ICell, K) - TracerUpwindNew) /
                 (TracerMaxNew - TracerUpwindNew + Eps);
             FlxIn(ICell, K) =
                 Kokkos::min(1.0_Real, Kokkos::max(0.0_Real, ScaleFactorIn));
             const Real ScaleFactorOut =
                 (TracerUpwindNew - TracerMin(ICell, K)) /
                 (TracerUpwindNew - TracerMinNew + Eps);
             FlxOut(ICell, K) =
                 Kokkos::min(1.0_Real, Kokkos::max(0.0_Real, ScaleFactorOut));
          });
   }

   KOKKOS_FUNCTION void FCTRescaleHighOrderFlux(const TeamMember &Team,
                                                const I4 IEdge) const {
      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);
      parallelForInner(
          Team, Range(KMin, KMax),
          INNER_LAMBDA(int K) { FCTRescaleHighOrderFlux(IEdge, K); });
   }
   KOKKOS_FUNCTION void FCTRescaleHighOrderFlux(const I4 IEdge,
                                                const I4 K) const {
      const I4 ICell1 = CellsOnEdge(IEdge, 0);
      const I4 ICell2 = CellsOnEdge(IEdge, 1);
      HighOrderFlx(IEdge, K) =
          Kokkos::max(0.0_Real, HighOrderFlx(IEdge, K)) *
              Kokkos::min(FlxOut(ICell1, K), FlxIn(ICell2, K)) +
          Kokkos::min(0.0_Real, HighOrderFlx(IEdge, K)) *
              Kokkos::min(FlxIn(ICell1, K), FlxOut(ICell2, K));
   }

   KOKKOS_FUNCTION void
   FCTAccumulateHighOrderFlux(const TeamMember &Team, const I4 ICell,
                              const Real Dt, const Array2DReal &TracerArray,
                              const Array2DReal &LayerThickness) const {

      // Accumulate the scaled high order vertical tendencies
      // and the upwind tendencies
      const Real InvAreaCell1 = 1.0_Real / AreaCell(ICell);
      // Accumulate the scaled high order horizontal tendencies
      for (I4 I = 0; I < NEdgesOnCell(ICell); ++I) {
         const I4 IEdge          = EdgesOnCell(ICell, I);
         const Real SignedFactor = InvAreaCell1 * EdgeSignOnCell(ICell, I);
         const I4 KMin           = MinLayerEdgeBot(IEdge);
         const I4 KMax           = MaxLayerEdgeTop(IEdge);
         parallelForInner(
             Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
                // WorkTend on RHS is upwind tendency
                // WorkTend on LHS is total horiz advect tendency
                WorkTend(ICell, K) += SignedFactor * HighOrderFlx(IEdge, K);
             });
      }
      const int KMin = MinLayerCell(ICell);
      const int KMax = MaxLayerCell(ICell);
      parallelForInner(
          Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
             // workTend  on RHS is total horiz advection tendency
             // TracerCur on LHS is provisional tracer after
             //                     horizontal fluxes only.
             TracerCur(ICell, K) =
                 (TracerCur(ICell, K) * LayerThickness(ICell, K) +
                  Dt * WorkTend(ICell, K)) *
                 HProvInv(ICell, K);
             TracerArray(ICell, K) += WorkTend(ICell, K);
          });
   }

   KOKKOS_FUNCTION void
   FCTComputeBudgetAdvectionEdgeFlux(const TeamMember &Team, const I4 L,
                                     const I4 IEdge) const {
      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);
      parallelForInner(
          Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
             FCTComputeBudgetAdvectionEdgeFlux(L, IEdge, K);
          });
   }
   KOKKOS_FUNCTION void FCTComputeBudgetAdvectionEdgeFlux(const I4 L,
                                                          const I4 IEdge,
                                                          const I4 K) const {
      // Save u*h*T flux on edge for analysis. This variable will be
      // divided by h at the end of the time step.
      ActiveTracerHorizontalAdvectionEdgeFlux(L, IEdge, K) =
          (LowOrderFlx(IEdge, K) + HighOrderFlx(IEdge, K)) / DvEdge(IEdge);
   }

   KOKKOS_FUNCTION void
   FCTComputeBudgetAdvectionTendency(const TeamMember &Team, const I4 L,
                                     const I4 ICell) const {
      const int KMin = MinLayerCell(ICell);
      const int KMax = MaxLayerCell(ICell);
      parallelForInner(
          Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
             FCTComputeBudgetAdvectionTendency(L, ICell, K);
          });
   }
   KOKKOS_FUNCTION void FCTComputeBudgetAdvectionTendency(const I4 L,
                                                          const I4 ICell,
                                                          const I4 K) const {
      ActiveTracerHorizontalAdvectionTendency(L, ICell, K) = WorkTend(ICell, K);
   }

   KOKKOS_FUNCTION void FCTMonotonicityCheck(const TeamMember &Team,
                                             const I4 ICell) const {
      const int KMin = MinLayerCell(ICell);
      const int KMax = MaxLayerCell(ICell);
      parallelForInner(
          Team, Range(KMin, KMax),
          INNER_LAMBDA(int K) { FCTMonotonicityCheck(ICell, K); });
   }
   KOKKOS_FUNCTION void FCTMonotonicityCheck(const I4 ICell, const I4 K) const {
      // Check tracer values against local min,max to detect
      // non-monotone values and write warning if found
      if (TracerCur(ICell, K) < TracerMin(ICell, K) - Eps) {
         Kokkos::printf(
             "Horizontal minimum out of bounds on cell %d, level %d "
             "new tracer value %lg is smaller than previous tracer minimum "
             "%lg\n",
             ICell, K, TracerCur(ICell, K), TracerMin(ICell, K));
      }
      if (TracerCur(ICell, K) > TracerMax(ICell, K) + Eps) {
         Kokkos::printf(
             "Horizontal maximum out of bounds on cell %d, level %d "
             "new tracer value %lg is larger than previous tracer maximum "
             "%lg\n",
             ICell, K, TracerCur(ICell, K), TracerMax(ICell, K));
      }
   }

 protected:
   const Real Eps = 1.e-10_Real;
   const HorzMesh *HorzontalMesh;
   const VertCoord *VerticalCoord;
   const I4 NVertLayers;

   Array1DI4 NAdvCellsForEdge;
   Array2DI4 AdvCellsForEdge;
   Array2DI4 AdvMaskHighOrder;
   Array2DI4 CellsOnCell;
   Array1DI4 NEdgesOnCell;
   Array2DReal AdvCoefs;
   Array2DReal AdvCoefs3rd;
   Array3DReal HighOrderFlxHorz;
   Array2DReal TracerCur;

   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
   Array2DReal EdgeSignOnCell;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
   Array2DReal VerticalPseudoVelocity;
   Array2DReal HProvInv;
   Array2DReal HNewInv;
   Array2DReal HProv;
   Array2DReal TracerMax;
   Array2DReal TracerMin;
   Array2DReal HighOrderFlx;
   Array2DReal LowOrderFlx;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array2DReal WorkTend;
   Array2DReal FlxIn;
   Array2DReal FlxOut;
   Array3DReal ActiveTracerHorizontalAdvectionEdgeFlux;
   Array3DReal ActiveTracerHorizontalAdvectionTendency;
};

// Tracer horizontal diffusion term
class TracerDiffOnCell {
 public:
   bool Enabled = false;

   Real EddyDiff2;

   TracerDiffOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   operator()(const Array3DReal &Tend, I4 L, I4 ICell, I4 KChunk,
              const Array3DReal &TracerCell,
              const Array2DReal &MeanPseudoThickEdge) const {

      const I4 KStartCell = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLenCell = chunkLength(KChunk, KStartCell, MaxLayerCell(ICell));
      const I4 KEndCell = KStartCell + KLenCell - 1;
      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real DiffTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge      = EdgesOnCell(ICell, J);
         const I4 KStartEdge = Kokkos::max(KStartCell, MinLayerEdgeBot(JEdge));
         const I4 KEndEdge   = Kokkos::min(KEndCell, MaxLayerEdgeTop(JEdge));

         const I4 JCell0 = CellsOnEdge(JEdge, 0);
         const I4 JCell1 = CellsOnEdge(JEdge, 1);

         const Real RTemp =
             MeshScalingDel2(JEdge) * DvEdge(JEdge) / DcEdge(JEdge);

         for (int K = KStartEdge; K <= KEndEdge; ++K) {
            const I4 KVec = K - KStartCell;
            const Real TracerGrad =
                (TracerCell(L, JCell1, K) - TracerCell(L, JCell0, K));

            DiffTmp[KVec] -= EdgeMask(JEdge, K) * EdgeSignOnCell(ICell, J) *
                             RTemp * MeanPseudoThickEdge(JEdge, K) * TracerGrad;
         }
      }
      for (int KVec = 0; KVec < KLenCell; ++KVec) {
         const I4 K = KStartCell + KVec;
         Tend(L, ICell, K) += EddyDiff2 * DiffTmp[KVec] * InvAreaCell;
      }
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

   KOKKOS_FUNCTION void operator()(const Array3DReal &Tend, I4 L, I4 ICell,
                                   I4 KChunk,
                                   const Array3DReal &TrDel2Cell) const {

      const I4 KStartCell = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLenCell = chunkLength(KChunk, KStartCell, MaxLayerCell(ICell));
      const I4 KEndCell = KStartCell + KLenCell - 1;
      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real HypTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge      = EdgesOnCell(ICell, J);
         const I4 KStartEdge = Kokkos::max(KStartCell, MinLayerEdgeBot(JEdge));
         const I4 KEndEdge   = Kokkos::min(KEndCell, MaxLayerEdgeTop(JEdge));

         const I4 JCell0 = CellsOnEdge(JEdge, 0);
         const I4 JCell1 = CellsOnEdge(JEdge, 1);

         const Real RTemp =
             MeshScalingDel4(JEdge) * DvEdge(JEdge) / DcEdge(JEdge);

         for (int K = KStartEdge; K <= KEndEdge; ++K) {
            const I4 KVec = K - KStartCell;
            const Real Del2TrGrad =
                (TrDel2Cell(L, JCell1, K) - TrDel2Cell(L, JCell0, K));

            HypTmp[KVec] -= EdgeMask(JEdge, K) * EdgeSignOnCell(ICell, J) *
                            RTemp * Del2TrGrad;
         }
      }
      for (int KVec = 0; KVec < KLenCell; ++KVec) {
         const I4 K = KStartCell + KVec;
         Tend(L, ICell, K) -= EddyDiff4 * HypTmp[KVec] * InvAreaCell;
      }
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
