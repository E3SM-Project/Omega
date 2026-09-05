//===-- ocn/TendencyTerms.cpp - Tendency Terms ------------------*- C++ -*-===//
//
// The tendency terms that update state variables are implemented as functors,
// i.e. as classes that act like functions. This source defines the class
// constructors for these functors, which initialize the functor objects using
// the Mesh objects and info from the Config. The function call operators () are
// defined in the corresponding header file.
//
//===----------------------------------------------------------------------===//

#include "TendencyTerms.h"
#include "AuxiliaryState.h"
#include "DataTypes.h"
#include "Eos.h"
#include "Error.h"
#include "HorzMesh.h"
#include "HorzOperators.h"
#include "OceanState.h"
#include "Tracers.h"

namespace OMEGA {

PseudoThicknessFluxDivOnCell::PseudoThicknessFluxDivOnCell(
    const HorzMesh *Mesh, const VertCoord *VCoord)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell),
      EdgeSignOnCell(Mesh->EdgeSignOnCell), NVertLayers(VCoord->NVertLayers),
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

PotentialVortHAdvOnEdge::PotentialVortHAdvOnEdge(const HorzMesh *Mesh,
                                                 const VertCoord *VCoord)
    : NEdgesOnEdge(Mesh->NEdgesOnEdge), EdgesOnEdge(Mesh->EdgesOnEdge),
      WeightsOnEdge(Mesh->WeightsOnEdge), NVertLayers(VCoord->NVertLayers),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

KEGradOnEdge::KEGradOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord)
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

SSHGradOnEdge::SSHGradOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord)
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

VelocityDiffusionOnEdge::VelocityDiffusionOnEdge(const HorzMesh *Mesh,
                                                 const VertCoord *VCoord)
    : CellsOnEdge(Mesh->CellsOnEdge), VerticesOnEdge(Mesh->VerticesOnEdge),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge),
      MeshScalingDel2(Mesh->MeshScalingDel2), EdgeMask(VCoord->EdgeMask),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

VelocityHyperDiffOnEdge::VelocityHyperDiffOnEdge(const HorzMesh *Mesh,
                                                 const VertCoord *VCoord)
    : CellsOnEdge(Mesh->CellsOnEdge), VerticesOnEdge(Mesh->VerticesOnEdge),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge),
      MeshScalingDel4(Mesh->MeshScalingDel4), EdgeMask(VCoord->EdgeMask),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

SfcStressForcingOnEdge::SfcStressForcingOnEdge(const HorzMesh *Mesh,
                                               const VertCoord *VCoord)
    : Enabled(false), EdgeMask(VCoord->EdgeMask),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

BottomDragOnEdge::BottomDragOnEdge(const HorzMesh *Mesh,
                                   const VertCoord *VCoord)
    : Enabled(false), Coeff(0), CellsOnEdge(Mesh->CellsOnEdge),
      NVertLayers(VCoord->NVertLayers), EdgeMask(VCoord->EdgeMask),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

SfcThicknessForcingOnCell::SfcThicknessForcingOnCell(const HorzMesh *Mesh,
                                                     const VertCoord *VCoord)
    : MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {}

SfcTracerForcingOnCell::SfcTracerForcingOnCell(const HorzMesh *Mesh,
                                               const VertCoord *VCoord,
                                               I4 TempTracerIndex,
                                               I4 SaltTracerIndex,
                                               const Eos *EosInst)
    : TempIndex(TempTracerIndex), SaltIndex(SaltTracerIndex),
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell),
      EosChoice(EosInst->EosChoice) {}

TracerHorzAdvOnCell::TracerHorzAdvOnCell(const HorzMesh *Mesh,
                                         const VertCoord *VCoord)
    : HorzontalMesh(Mesh), VerticalCoord(VCoord),
      NAdvCellsForEdge("NumberOfCellsContribToAdvectionAtEdge",
                       Mesh->NEdgesAll),
      AdvCellsForEdge("IndexOfCellsContributingToAdvection", Mesh->NEdgesAll,
                      Mesh->MaxEdges2 + 2),
      AdvMaskHighOrder("MaskForHighOrderAdvectionTerms", Mesh->NEdgesAll,
                       VCoord->NVertLayers),
      AdvCoefs("CommonAdvectionCoefficients", Mesh->MaxEdges2 + 2,
               Mesh->NEdgesAll),
      AdvCoefs3rd("CommonAdvectionCoeffsForHighOrder", Mesh->MaxEdges2 + 2,
                  Mesh->NEdgesAll),
      HighOrderFlxHorz("HigherOrderHorizontalFlux", Tracers::getNumTracers(),
                       Mesh->NEdgesAll, VCoord->NVertLayers),
      NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell),
      NVertLayers(VCoord->NVertLayers), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

TracerDiffOnCell::TracerDiffOnCell(const HorzMesh *Mesh,
                                   const VertCoord *VCoord)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DvEdge(Mesh->DvEdge), DcEdge(Mesh->DcEdge), AreaCell(Mesh->AreaCell),
      MeshScalingDel2(Mesh->MeshScalingDel2), NVertLayers(VCoord->NVertLayers),
      EdgeMask(VCoord->EdgeMask), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

TracerHyperDiffOnCell::TracerHyperDiffOnCell(const HorzMesh *Mesh,
                                             const VertCoord *VCoord)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DvEdge(Mesh->DvEdge), DcEdge(Mesh->DcEdge), AreaCell(Mesh->AreaCell),
      MeshScalingDel4(Mesh->MeshScalingDel4), NVertLayers(VCoord->NVertLayers),
      EdgeMask(VCoord->EdgeMask), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

SurfaceTracerRestoringOnCell::SurfaceTracerRestoringOnCell(
    const HorzMesh *Mesh) {}

FrazilOnCell::FrazilOnCell(const HorzMesh *Mesh, const VertCoord *VCoord)
    : NCellsAll(Mesh->NCellsAll), TempTracerIndex(-1), SaltTracerIndex(-1),
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {
   Tracers::getIndex(TempTracerIndex, "Temperature");
   Tracers::getIndex(SaltTracerIndex, "Salinity");

   OMEGA_REQUIRE(TempTracerIndex >= 0,
                 "FrazilOnCell: Temperature tracer index is undefined");
   OMEGA_REQUIRE(SaltTracerIndex >= 0,
                 "FrazilOnCell: Salinity tracer index is undefined");
}

void FrazilOnCell::operator()(const Array2DReal &PseudoThicknessTend,
                              const Array3DReal &TracerTend,
                              const Array3DReal &TracerArray,
                              const Array2DReal &PressureMid,
                              const Array2DReal &PseudoThickness) const {
   auto *Frazil = Frazil::getDefault();
   if (!Enabled || !Frazil) {
      return;
   }

   deepCopy(Frazil->FrazilTTend, 0.0_Real);
   deepCopy(Frazil->FrazilSTend, 0.0_Real);
   deepCopy(Frazil->FrazilHTend, 0.0_Real);
   deepCopy(Frazil->AccMIce, 0.0_Real);
   deepCopy(Frazil->AccEIce, 0.0_Real);
   deepCopy(Frazil->AccMLiq, 0.0_Real);
   deepCopy(Frazil->AccELiq, 0.0_Real);
   deepCopy(Frazil->AccMSalt, 0.0_Real);

   const auto ConservTemp =
       Kokkos::subview(TracerArray, TempTracerIndex, Kokkos::ALL, Kokkos::ALL);
   const auto AbsSalinity =
       Kokkos::subview(TracerArray, SaltTracerIndex, Kokkos::ALL, Kokkos::ALL);

   Frazil->computeFrazil(ConservTemp, AbsSalinity, PressureMid,
                         PseudoThickness);

   const auto FrazilHTend = Frazil->FrazilHTend;
   const auto FrazilTTend = Frazil->FrazilTTend;
   const auto FrazilSTend = Frazil->FrazilSTend;
   const I4 TempIndex     = TempTracerIndex;
   const I4 SaltIndex     = SaltTracerIndex;

   OMEGA_SCOPE(LocPseudoThicknessTend, PseudoThicknessTend);
   OMEGA_SCOPE(LocTracerTend, TracerTend);
   OMEGA_SCOPE(LocFrazilHTend, FrazilHTend);
   OMEGA_SCOPE(LocFrazilTTend, FrazilTTend);
   OMEGA_SCOPE(LocFrazilSTend, FrazilSTend);
   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);

   parallelForOuter(
       {NCellsAll}, KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const int KMin = LocMinLayerCell(ICell);
          const int KMax = LocMaxLayerCell(ICell);

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 LocPseudoThicknessTend(ICell, K) += LocFrazilHTend(ICell, K);
                 LocTracerTend(TempIndex, ICell, K) += LocFrazilTTend(ICell, K);
                 LocTracerTend(SaltIndex, ICell, K) += LocFrazilSTend(ICell, K);
              });
       });
}

void TracerHorzAdvOnCell::init() {
   const HorzMesh *Mesh    = this->HorzontalMesh;
   const VertCoord *VCoord = this->VerticalCoord;
   const auto MaxEdges2    = Mesh->MaxEdges2;
   const auto NEdgesAll    = Mesh->NEdgesAll;
   const auto NCellsAll    = Mesh->NCellsAll;
   // Allocate Kokkos arrays in member data

   if (ForceLowOrder) {
      // Return when the 2nd-order tracer horz adv
      deepCopy(NAdvCellsForEdge, 0);
      deepCopy(AdvMaskHighOrder, 0);
      return;
   }

   SecondDerivativeOnCell secondDerivativeOnCell(Mesh);
   Array3DReal DerivTwo("DerivTwo", MaxEdges2 + 2, 2, NEdgesAll);
   parallelFor(
       {NCellsAll},
       KOKKOS_LAMBDA(int ICell) { secondDerivativeOnCell(DerivTwo, ICell); });
   // Compute masks and coefficients
   Kokkos::fence();
   MasksAndCoefficients masksAndCoefficients(
       Mesh, VCoord, DerivTwo, NAdvCellsForEdge, AdvCellsForEdge,
       AdvMaskHighOrder, AdvCoefs, AdvCoefs3rd);
   Kokkos::fence();
   parallelFor(
       {NEdgesAll}, KOKKOS_LAMBDA(int IEdge) { masksAndCoefficients(IEdge); });
   Kokkos::fence();
}
} // end namespace OMEGA

//===----------------------------------------------------------------------===//
