//===-- base/VertCoord.cpp - vertical coordinate ----------------*- C++ -*-===//
//
//===----------------------------------------------------------------------===//

#include "VertCoord.h"
#include "Dimension.h"
#include "IO.h"
#include "OmegaKokkos.h"

namespace OMEGA {

VertCoord *VertCoord::DefaultVertCoord = nullptr;
std::map<std::string, std::unique_ptr<VertCoord>> VertCoord::AllVertCoords;

//------------------------------------------------------------------------------
// Initialize default vertical coordinate, requires prior initialization of
// HorzMesh.
void VertCoord::init() {

   Decomp *DefDecomp = Decomp::getDefault();

   Config *OmegaConfig = Config::getOmegaConfig();

   VertCoord::DefaultVertCoord = create("Default", DefDecomp, OmegaConfig);

} // end init

VertCoord::VertCoord(const std::string &Name_, const Decomp *Decomp,
                     Config *Options) {

   Name = Name_;

   // Retrieve mesh filename from Decomp
   MeshFileName = Decomp->MeshFileName;

   // Open the mesh file for reading (assume IO has already been initialized)
   I4 Err;
   Err = IO::openFile(MeshFileID, MeshFileName, IO::ModeRead);

   // Set NVertLevels and NVertLevelsP1 and create the vertical dimension
   I4 NVertLevelsID;
   Err = IO::getDimFromFile(MeshFileID, "nVertLevels", NVertLevelsID,
                            NVertLevels);

   NVertLevelsP1 = NVertLevels + 1;

   auto VertDim = Dimension::create("NVertLevels", NVertLevels);

   // Retrieve mesh variables from Decomp
   NCellsOwned = Decomp->NCellsOwned;
   NCellsAll   = Decomp->NCellsAll;
   NCellsSize  = Decomp->NCellsSize;

   NEdgesOwned = Decomp->NEdgesOwned;
   NEdgesAll   = Decomp->NEdgesAll;
   NEdgesSize  = Decomp->NEdgesSize;

   NVerticesOwned = Decomp->NVerticesOwned;
   NVerticesAll   = Decomp->NVerticesAll;
   NVerticesSize  = Decomp->NVerticesSize;
   VertexDegree   = Decomp->VertexDegree;

   // Retrieve connectivity arrays from HorzMesh
   CellsOnEdge   = Decomp->CellsOnEdge;
   CellsOnVertex = Decomp->CellsOnVertex;

   // Allocate device arrays
   MaxLevelCell = Array1DI4("MaxLevelCell", NCellsSize);
   MinLevelCell = Array1DI4("MinLevelCell", NCellsSize);
   BottomDepth  = Array1DReal("BottomDepth", NCellsSize);

   PressureInterface =
       Array2DReal("PressureInterface", NCellsSize, NVertLevelsP1);
   PressureMid = Array2DReal("PressureMid", NCellsSize, NVertLevels);

   ZInterface = Array2DReal("ZInterface", NCellsSize, NVertLevelsP1);
   ZMid       = Array2DReal("ZMid", NCellsSize, NVertLevels);

   GeopotentialMid = Array2DReal("GeopotentialMid", NCellsSize, NVertLevels);
   LayerThicknessTarget =
       Array2DReal("LayerThicknessTarget", NCellsSize, NVertLevels);

   MaxLevelCellH         = createHostMirrorCopy(MaxLevelCell);
   MinLevelCellH         = createHostMirrorCopy(MinLevelCell);
   BottomDepthH          = createHostMirrorCopy(BottomDepth);
   PressureInterfaceH    = createHostMirrorCopy(PressureInterface);
   PressureMidH          = createHostMirrorCopy(PressureMid);
   ZInterfaceH           = createHostMirrorCopy(ZInterface);
   ZMidH                 = createHostMirrorCopy(ZMid);
   GeopotentialMidH      = createHostMirrorCopy(GeopotentialMid);
   LayerThicknessTargetH = createHostMirrorCopy(LayerThicknessTarget);

   readArrays(Decomp);

   minMaxLevelEdge();
   minMaxLevelVertex();

} // end constructor

VertCoord *VertCoord::create(const std::string &Name, const Decomp *Decomp,
                             Config *Options) {
   // Check to see if a VertCoord of the same name already exists and, if so,
   // exit with an error
   if (AllVertCoords.find(Name) != AllVertCoords.end()) {
      LOG_ERROR("Attempted to create a VertCoord with name {} but a VertCoord "
                "of that name already exists",
                Name);
      return nullptr;
   }

   // create a new VertCoord on the heap and put it in a map of unique_ptrs,
   // which will manage its lifetime
   auto *NewVertCoord = new VertCoord(Name, Decomp, Options);
   AllVertCoords.emplace(Name, NewVertCoord);

   return NewVertCoord;
} // end create

void VertCoord::readArrays(const Decomp *Decomp) {

   I4 NDims             = 1;
   IO::Rearranger Rearr = IO::RearrBox;

   I4 CellDecompI4;
   std::vector<I4> CellDims{Decomp->NCellsGlobal};
   std::vector<I4> CellID(NCellsAll);
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      CellID[Cell] = Decomp->CellIDH(Cell) - 1;
   }
   I4 DecErr = IO::createDecomp(CellDecompI4, IO::IOTypeI4, NDims, CellDims,
                                NCellsAll, CellID, Rearr);
   if (DecErr != 0) {
      LOG_CRITICAL("VertCoord: error creating cell I4 IO decomposition");
   }

   HostArray1DI4 TmpArrayI4H("TmpCellArray", NCellsSize);
   I4 ArrayID;
   const std::string MaxNameMPAS = "maxLevelCell";
   I4 MaxReadErr = IO::readArray(TmpArrayI4H.data(), NCellsAll, MaxNameMPAS,
                                 MeshFileID, CellDecompI4, ArrayID);

   if (MaxReadErr != 0) {
      LOG_WARN("VertCoord: error reading maxLevelCell from mesh file, "
               "using MaxLevelCell = NVertLevels - 1");
      deepCopy(TmpArrayI4H, NVertLevels - 1);
   } else {
      for (int ICell = 0; ICell < NCellsAll; ++ICell) {
         TmpArrayI4H(ICell) = TmpArrayI4H(ICell) - 1;
      }
   }
   TmpArrayI4H(NCellsAll) = -1;

   deepCopy(MaxLevelCellH, TmpArrayI4H);
   deepCopy(MaxLevelCell, TmpArrayI4H);

   const std::string MinNameMPAS = "minLevelCell";
   I4 MinReadErr = IO::readArray(TmpArrayI4H.data(), NCellsAll, MinNameMPAS,
                                 MeshFileID, CellDecompI4, ArrayID);

   if (MinReadErr != 0) {
      LOG_WARN("VertCoord: error reading minLevelCell from mesh file, "
               "using MinLevelCell = 0");
      deepCopy(TmpArrayI4H, 0);
   } else {
      for (int ICell = 0; ICell < NCellsAll; ++ICell) {
         TmpArrayI4H(ICell) = TmpArrayI4H(ICell) - 1;
      }
   }
   TmpArrayI4H(NCellsAll) = NVertLevelsP1;

   deepCopy(MinLevelCellH, TmpArrayI4H);
   deepCopy(MinLevelCell, TmpArrayI4H);

   I4 CellDecompR8;
   DecErr = IO::createDecomp(CellDecompR8, IO::IOTypeR8, NDims, CellDims,
                             NCellsAll, CellID, Rearr);
   if (DecErr != 0) {
      LOG_CRITICAL("VertCoord: error creating cell R8 IO decomposition");
   }

   HostArray1DR8 TmpArrayR8H("TmpCellArray", NCellsSize);
   const std::string BotNameMPAS = "bottomDepth";
   I4 BotReadErr = IO::readArray(TmpArrayR8H.data(), NCellsAll, BotNameMPAS,
                                 MeshFileID, CellDecompR8, ArrayID);
   if (BotReadErr != 0) {
      LOG_CRITICAL("VertCoord: error reading bottomDepth");
   }

   deepCopy(BottomDepthH, TmpArrayR8H);
   deepCopy(BottomDepth, TmpArrayR8H);

   DecErr = IO::destroyDecomp(CellDecompI4);
   if (DecErr != 0) {
      LOG_CRITICAL("VertCoord: error destroying cell I4 IO decomposition");
   }
   DecErr = IO::destroyDecomp(CellDecompR8);
   if (DecErr != 0) {
      LOG_CRITICAL("VertCoord: error destroying cell R8 IO decomposition");
   }
}

//------------------------------------------------------------------------------
// Destroys a local VertCoord and deallocates all arrays
VertCoord::~VertCoord() {} // end deconstructor

//------------------------------------------------------------------------------
// Removes a VertCoord from list by name
void VertCoord::erase(std::string InName) {
   AllVertCoords.erase(InName); // removes the VertCoord from the list and in
                                // the process, calls the destructor
} // end erase

//------------------------------------------------------------------------------
// Removes all VertCoords to clean up before exit
void VertCoord::clear() {

   AllVertCoords.clear(); // removes all VertCoords from the list and in the
                          // process, calls the destructors for each

} // end clear

//------------------------------------------------------------------------------
void VertCoord::minMaxLevelEdge() {

   MinLevelEdgeTop = Array1DI4("MinLevelEdgeTop", NEdgesSize);
   MinLevelEdgeBot = Array1DI4("MinLevelEdgeBot", NEdgesSize);
   MaxLevelEdgeTop = Array1DI4("MaxLevelEdgeTop", NEdgesSize);
   MaxLevelEdgeBot = Array1DI4("MaxLevelEdgeBot", NEdgesSize);

   OMEGA_SCOPE(LocNVertLevelsP1, NVertLevelsP1);
   OMEGA_SCOPE(LocCellsOnEdge, CellsOnEdge);
   OMEGA_SCOPE(LocMinLevelCell, MinLevelCell);
   OMEGA_SCOPE(LocMaxLevelCell, MaxLevelCell);
   OMEGA_SCOPE(LocMinLevelEdgeTop, MinLevelEdgeTop);
   OMEGA_SCOPE(LocMinLevelEdgeBot, MinLevelEdgeBot);
   OMEGA_SCOPE(LocMaxLevelEdgeTop, MaxLevelEdgeTop);
   OMEGA_SCOPE(LocMaxLevelEdgeBot, MaxLevelEdgeBot);
   parallelFor(
       {NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
          I4 Lvl1;
          I4 Lvl2;
          const I4 ICell1 = LocCellsOnEdge(IEdge, 0);
          const I4 ICell2 = LocCellsOnEdge(IEdge, 1);
          Lvl1            = LocMaxLevelCell(ICell1) == -1 ? LocNVertLevelsP1
                                                          : LocMinLevelCell(ICell1);
          Lvl2            = LocMaxLevelCell(ICell2) == -1 ? LocNVertLevelsP1
                                                          : LocMinLevelCell(ICell2);
          LocMinLevelEdgeTop(IEdge) = Kokkos::min(Lvl1, Lvl2);

          Lvl1 = LocMaxLevelCell(ICell1) == -1 ? 0 : LocMinLevelCell(ICell1);
          Lvl2 = LocMaxLevelCell(ICell2) == -1 ? 0 : LocMinLevelCell(ICell2);
          LocMinLevelEdgeBot(IEdge) = Kokkos::max(Lvl1, Lvl2);

          LocMaxLevelEdgeTop(IEdge) =
              Kokkos::min(LocMaxLevelCell(ICell1), LocMaxLevelCell(ICell2));
          LocMaxLevelEdgeBot(IEdge) =
              Kokkos::max(LocMaxLevelCell(ICell1), LocMaxLevelCell(ICell2));
       });

   OMEGA_SCOPE(LocNEdgesAll, NEdgesAll);
   parallelFor(
       {1}, KOKKOS_LAMBDA(const int &) {
          LocMinLevelEdgeTop(LocNEdgesAll) = LocNVertLevelsP1;
          LocMinLevelEdgeBot(LocNEdgesAll) = LocNVertLevelsP1;
          LocMaxLevelEdgeTop(LocNEdgesAll) = -1;
          LocMaxLevelEdgeBot(LocNEdgesAll) = -1;
       });

   MinLevelEdgeTopH = createHostMirrorCopy(MinLevelEdgeTop);
   MinLevelEdgeBotH = createHostMirrorCopy(MinLevelEdgeBot);
   MaxLevelEdgeTopH = createHostMirrorCopy(MaxLevelEdgeTop);
   MaxLevelEdgeBotH = createHostMirrorCopy(MaxLevelEdgeBot);
}

//------------------------------------------------------------------------------
void VertCoord::minMaxLevelVertex() {

   MinLevelVertexTop = Array1DI4("MinLevelVertexTop", NVerticesSize);
   MinLevelVertexBot = Array1DI4("MinLevelVertexBot", NVerticesSize);
   MaxLevelVertexTop = Array1DI4("MaxLevelVertexTop", NVerticesSize);
   MaxLevelVertexBot = Array1DI4("MaxLevelVertexBot", NVerticesSize);

   OMEGA_SCOPE(LocNVertLevelsP1, NVertLevelsP1);
   OMEGA_SCOPE(LocVertexDegree, VertexDegree);
   OMEGA_SCOPE(LocCellsOnVertex, CellsOnVertex);
   OMEGA_SCOPE(LocMinLevelCell, MinLevelCell);
   OMEGA_SCOPE(LocMaxLevelCell, MaxLevelCell);
   OMEGA_SCOPE(LocMinLevelVertexTop, MinLevelVertexTop);
   OMEGA_SCOPE(LocMinLevelVertexBot, MinLevelVertexBot);
   OMEGA_SCOPE(LocMaxLevelVertexTop, MaxLevelVertexTop);
   OMEGA_SCOPE(LocMaxLevelVertexBot, MaxLevelVertexBot);

   parallelFor(
       {NVerticesAll}, KOKKOS_LAMBDA(int IVertex) {
          I4 Lvl;
          I4 ICell = LocCellsOnVertex(IVertex, 0);
          Lvl      = LocMaxLevelCell(ICell) == -1 ? 0 : LocMinLevelCell(ICell);
          LocMinLevelVertexBot(IVertex) = Lvl;
          for (int I = 1; I < LocVertexDegree; ++I) {
             ICell = LocCellsOnVertex(IVertex, I);
             Lvl   = LocMaxLevelCell(ICell) == -1 ? 0 : LocMinLevelCell(ICell);
             LocMinLevelVertexBot(IVertex) =
                 Kokkos::max(LocMinLevelVertexBot(IVertex), Lvl);
          }

          ICell = LocCellsOnVertex(IVertex, 0);
          Lvl   = LocMaxLevelCell(ICell) == -1 ? LocNVertLevelsP1
                                               : LocMinLevelCell(ICell);
          LocMinLevelVertexTop(IVertex) = Lvl;
          for (int I = 1; I < LocVertexDegree; ++I) {
             ICell = LocCellsOnVertex(IVertex, I);
             Lvl   = LocMaxLevelCell(ICell) == -1 ? LocNVertLevelsP1
                                                  : LocMinLevelCell(ICell);
             LocMinLevelVertexTop(IVertex) =
                 Kokkos::min(LocMinLevelVertexTop(IVertex), Lvl);
          }

          ICell                         = LocCellsOnVertex(IVertex, 0);
          LocMaxLevelVertexBot(IVertex) = LocMaxLevelCell(ICell);
          for (int I = 1; I < LocVertexDegree; ++I) {
             ICell                         = LocCellsOnVertex(IVertex, I);
             LocMaxLevelVertexBot(IVertex) = Kokkos::max(
                 LocMaxLevelVertexBot(IVertex), LocMaxLevelCell(ICell));
          }

          ICell                         = LocCellsOnVertex(IVertex, 0);
          LocMaxLevelVertexTop(IVertex) = LocMaxLevelCell(ICell);
          for (int I = 1; I < LocVertexDegree; ++I) {
             ICell                         = LocCellsOnVertex(IVertex, I);
             LocMaxLevelVertexTop(IVertex) = Kokkos::min(
                 LocMaxLevelVertexTop(IVertex), LocMaxLevelCell(ICell));
          }
       });
   OMEGA_SCOPE(LocNVerticesAll, NVerticesAll);
   parallelFor(
       {1}, KOKKOS_LAMBDA(const int &) {
          LocMinLevelVertexTop(LocNVerticesAll) = LocNVertLevelsP1;
          LocMinLevelVertexBot(LocNVerticesAll) = LocNVertLevelsP1;
          LocMaxLevelVertexTop(LocNVerticesAll) = -1;
          LocMaxLevelVertexBot(LocNVerticesAll) = -1;
       });

   MinLevelVertexTopH = createHostMirrorCopy(MinLevelVertexTop);
   MinLevelVertexBotH = createHostMirrorCopy(MinLevelVertexBot);
   MaxLevelVertexTopH = createHostMirrorCopy(MaxLevelVertexTop);
   MaxLevelVertexBotH = createHostMirrorCopy(MaxLevelVertexBot);
}

//------------------------------------------------------------------------------
void VertCoord::initMovementWeights(Config *Options) {

   Error Err; // default successful error code

   Config VCoordConfig("VertCoord");
   Err += Options->get(VCoordConfig);
   CHECK_ERROR_ABORT(Err, "VertCoord: VertCoord group not found in Config");

   std::string MovementWeightType;
   Err += VCoordConfig.get("MovementWeightType", MovementWeightType);
   CHECK_ERROR_ABORT(Err,
                     "VertCoord: MovementWeightType not found in VertCoord");

   VertCoordMovementWeights =
       Array1DReal("VertCoordMovementWeights", NVertLevels);

   OMEGA_SCOPE(LocVertCoordMovementWeights, VertCoordMovementWeights);
   if (MovementWeightType == "Fixed") {
      deepCopy(VertCoordMovementWeights, 0._Real);
      parallelFor(
          {1}, KOKKOS_LAMBDA(const int &) {
             LocVertCoordMovementWeights(0) = 1._Real;
          });
   } else if (MovementWeightType == "Uniform") {
      deepCopy(VertCoordMovementWeights, 1._Real);
   } else {
      ABORT_ERROR("VertCoord: Unknown MovementWeightType requested");
   }

   VertCoordMovementWeightsH = createHostMirrorCopy(VertCoordMovementWeights);
}

//------------------------------------------------------------------------------
void VertCoord::computePressure(const Array2DReal &PressureInterface,
                                const Array2DReal &PressureMid,
                                const Array2DReal &LayerThickness,
                                const Array1DReal &SurfacePressure) {

   Real Gravity = 9.80616_Real;
   Real Rho0    = 1035._Real;

   OMEGA_SCOPE(LocMinLevelCell, MinLevelCell);
   OMEGA_SCOPE(LocMaxLevelCell, MaxLevelCell);

   const auto Policy = TeamPolicy(NCellsAll, OMEGA_TEAMSIZE, 1);
   Kokkos::parallel_for(
       "computePressure", Policy, KOKKOS_LAMBDA(const TeamMember &Member) {
          const I4 ICell = Member.league_rank();
          const I4 KMin  = LocMinLevelCell(ICell);
          const I4 KMax  = LocMaxLevelCell(ICell);
          const I4 Range = KMax - KMin + 1;

          PressureInterface(ICell, KMin) = SurfacePressure(ICell);
          Kokkos::parallel_scan(
              TeamThreadRange(Member, Range),
              [&](int K, Real &Accum, bool IsFinal) {
                 const I4 KLvl  = K + KMin;
                 Real Increment = Gravity * Rho0 * LayerThickness(ICell, KLvl);
                 Accum += Increment;

                 if (IsFinal) {
                    PressureInterface(ICell, KLvl + 1) =
                        SurfacePressure(ICell) + Accum;
                    PressureMid(ICell, KLvl) =
                        SurfacePressure(ICell) + Accum - 0.5 * Increment;
                 }
              });
       });
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
