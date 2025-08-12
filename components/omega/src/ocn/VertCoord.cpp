//===-- base/VertCoord.cpp - vertical coordinate ----------------*- C++ -*-===//
//
// The VertCoord class contains member variables related to the vertical mesh
// and methods for computing these variables in Omega. The class will also serve
// as the container for information defining the extent of the vertical
// dimension and the number of active vertical layers in each column.
//
//===----------------------------------------------------------------------===//

#include "VertCoord.h"
#include "Dimension.h"
#include "Field.h"
#include "IO.h"
#include "OmegaKokkos.h"

#include <limits>

namespace OMEGA {

// create static class members
VertCoord *VertCoord::DefaultVertCoord = nullptr;
std::map<std::string, std::unique_ptr<VertCoord>> VertCoord::AllVertCoords;

//------------------------------------------------------------------------------
// Initialize default vertical coordinate, requires prior initialization of
// Decomp.
void VertCoord::init() {

   Decomp *DefDecomp = Decomp::getDefault();

   Config *OmegaConfig = Config::getOmegaConfig();

   VertCoord::DefaultVertCoord = create("Default", DefDecomp, OmegaConfig);

} // end init

//------------------------------------------------------------------------------
// Construct a new VertCoord instance given a Decomp
VertCoord::VertCoord(const std::string &Name_, //< [in] Name for new VertCoord
                     const Decomp *Decomp,     //< [in] associated Decomp
                     Config *Options           //< [in] configuration options
) {

   Name = Name_;

   // Retrieve mesh filename from Decomp
   MeshFileName = Decomp->MeshFileName;

   // Open the mesh file for reading (assume IO has already been initialized)
   I4 Err;
   Err = IO::openFile(MeshFileID, MeshFileName, IO::ModeRead);

   // Set NVertLayers and NVertLayersP1 and create the vertical dimension
   I4 NVertLayersID;
   Err = IO::getDimFromFile(MeshFileID, "nVertLevels", NVertLayersID,
                            NVertLayers);
   if (Err != 0) {
      LOG_WARN("VertCoord: error reading nVertLevels from mesh file, "
               "using NVertLayers = 1");
      NVertLayers = 1;
   }
   NVertLayersP1 = NVertLayers + 1;

   auto VertDim = Dimension::create("NVertLayers", NVertLayers);

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
   MaxLayerCell = Array1DI4("MaxLayerCell", NCellsSize);
   MinLayerCell = Array1DI4("MinLayerCell", NCellsSize);
   BottomDepth  = Array1DReal("BottomDepth", NCellsSize);

   PressureInterface =
       Array2DReal("PressureInterface", NCellsSize, NVertLayersP1);
   PressureMid = Array2DReal("PressureMid", NCellsSize, NVertLayers);

   ZInterface = Array2DReal("ZInterface", NCellsSize, NVertLayersP1);
   ZMid       = Array2DReal("ZMid", NCellsSize, NVertLayers);

   GeopotentialMid = Array2DReal("GeopotentialMid", NCellsSize, NVertLayers);
   LayerThicknessTarget =
       Array2DReal("LayerThicknessTarget", NCellsSize, NVertLayers);

   MaxLayerCellH         = createHostMirrorCopy(MaxLayerCell);
   MinLayerCellH         = createHostMirrorCopy(MinLayerCell);
   BottomDepthH          = createHostMirrorCopy(BottomDepth);
   PressureInterfaceH    = createHostMirrorCopy(PressureInterface);
   PressureMidH          = createHostMirrorCopy(PressureMid);
   ZInterfaceH           = createHostMirrorCopy(ZInterface);
   ZMidH                 = createHostMirrorCopy(ZMid);
   GeopotentialMidH      = createHostMirrorCopy(GeopotentialMid);
   LayerThicknessTargetH = createHostMirrorCopy(LayerThicknessTarget);

   readArrays(Decomp);

   minMaxLayerEdge();
   minMaxLayerVertex();

   initMovementWeights(Options);

} // end constructor

//------------------------------------------------------------------------------
// Calls the VertCoord constructor and places it in the AllVertCoords map
VertCoord *
VertCoord::create(const std::string &Name, // [in] name for new VertCoord
                  const Decomp *Decomp,    // [in] associated Decomp
                  Config *Options          // [in] configuration options
) {
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

//------------------------------------------------------------------------------
// Define IO fields and metadata
void VertCoord::defineFields() {

   I4 Err = 0; // default error code

   // Set field names (append Name if not default)
   PressInterfFldName    = "PressureInterface";
   PressMidFldName       = "PressureMid";
   ZInterfFldName        = "ZInterface";
   ZMidFldName           = "ZMid";
   GeopotFldName         = "GeopotentialMid";
   LyrThickTargetFldName = "LayerThicknessTarget";

   if (Name != "Default") {
      PressInterfFldName.append(Name);
      PressMidFldName.append(Name);
      ZInterfFldName.append(Name);
      ZMidFldName.append(Name);
      GeopotFldName.append(Name);
      LyrThickTargetFldName.append(Name);
   }

   // Create fields for VertCoord variables
   const Real FillValue = -9.99e30;
   int NDims            = 2;
   std::vector<std::string> DimNames(NDims);
   DimNames[0] = "NCells";
   DimNames[1] = "NVertLayers";

   auto PressureInterfaceField = Field::create(
       PressInterfFldName,                      // field name
       "Pressure at vertical layer interfaces", // long name or description
       "Pa",                                    // units
       "sea_water_pressure",                    // CF standard Name
       0.0,                                     // min valid value
       std::numeric_limits<Real>::max(),        // max valid value
       FillValue,                               // scalar for undefined entries
       NDims,                                   // number of dimensions
       DimNames                                 // dimension names
   );

   auto PressureMidField = Field::create(
       PressMidFldName,                        // field name
       "Pressure at vertical layer midpoints", // long name or description
       "Pa",                                   // units
       "sea_water_pressure",                   // CF standard Name
       0.0,                                    // min valid value
       std::numeric_limits<Real>::max(),       // max valid value
       FillValue,                              // scalar for undefined entries
       NDims,                                  // number of dimensions
       DimNames                                // dimension names
   );

   auto ZInterfaceField = Field::create(
       ZInterfFldName,                               // field name
       "Cartesian Z coordinate at layer interfaces", // long name or description
       "m",                                          // units
       "height",                                     // CF standard Name
       std::numeric_limits<Real>::min(),             // min valid value
       std::numeric_limits<Real>::max(),             // max valid value
       FillValue, // scalar for undefined entries
       NDims,     // number of dimensions
       DimNames   // dimension names
   );

   auto ZMidField = Field::create(
       ZMidFldName,                                 // field name
       "Cartesian Z coordinate at layer midpoints", // long name or description
       "m",                                         // units
       "height",                                    // CF standard Name
       std::numeric_limits<Real>::min(),            // min valid value
       std::numeric_limits<Real>::max(),            // max valid value
       FillValue, // scalar for undefined entries
       NDims,     // number of dimensions
       DimNames   // dimension names
   );

   auto GeopotentialMidField = Field::create(
       GeopotFldName,                     // field name
       "Geopotential at layer midpoints", // long name or description
       "m^2 s^-2",                        // units
       "geopotential",                    // CF standard Name
       std::numeric_limits<Real>::min(),  // min valid value
       std::numeric_limits<Real>::max(),  // max valid value
       FillValue,                         // scalar for undefined entries
       NDims,                             // number of dimensions
       DimNames                           // dimension names
   );

   auto LayerThicknessTargetField =
       Field::create(LyrThickTargetFldName, // field name
                     "desired layer thickness based on total perturbation from "
                     "the reference thickness", // long name or description
                     "m",                       // units
                     "",                        // CF standard Name
                     0.0,                       // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValue, // scalar for undefined entries
                     NDims,     // number of dimensions
                     DimNames   // dimension names
       );

   // Create a field group for VertCoord fields
   GroupName = "VertCoord";
   if (Name != "Default") {
      GroupName.append(Name);
   }
   auto VCoordGroup = FieldGroup::create(GroupName);

   Err = VCoordGroup->addField(PressInterfFldName);
   if (Err != 0)
      LOG_ERROR("Error adding {} to field group {}", PressInterfFldName,
                GroupName);
   Err = VCoordGroup->addField(PressMidFldName);
   if (Err != 0)
      LOG_ERROR("Error adding {} to field group {}", PressMidFldName,
                GroupName);
   Err = VCoordGroup->addField(ZInterfFldName);
   if (Err != 0)
      LOG_ERROR("Error adding {} to field group {}", ZInterfFldName, GroupName);
   Err = VCoordGroup->addField(ZMidFldName);
   if (Err != 0)
      LOG_ERROR("Error adding {} to field group {}", ZMidFldName, GroupName);
   Err = VCoordGroup->addField(GeopotFldName);
   if (Err != 0)
      LOG_ERROR("Error adding {} to field group {}", GeopotFldName, GroupName);
   Err = VCoordGroup->addField(LyrThickTargetFldName);
   if (Err != 0)
      LOG_ERROR("Error adding {} to field group {}", LyrThickTargetFldName,
                GroupName);

   // Associate Field with data
   Err = PressureInterfaceField->attachData<Array2DReal>(PressureInterface);
   if (Err != 0)
      LOG_ERROR("Error attaching data array to field {}", PressInterfFldName);
   Err = PressureMidField->attachData<Array2DReal>(PressureMid);
   if (Err != 0)
      LOG_ERROR("Error attaching data array to field {}", PressMidFldName);
   Err = ZInterfaceField->attachData<Array2DReal>(ZInterface);
   if (Err != 0)
      LOG_ERROR("Error attaching data array to field {}", ZInterfFldName);
   Err = ZMidField->attachData<Array2DReal>(ZMid);
   if (Err != 0)
      LOG_ERROR("Error attaching data array to field {}", ZMidFldName);
   Err = GeopotentialMidField->attachData<Array2DReal>(GeopotentialMid);
   if (Err != 0)
      LOG_ERROR("Error attaching data array to field {}", GeopotFldName);
   Err =
       LayerThicknessTargetField->attachData<Array2DReal>(LayerThicknessTarget);
   if (Err != 0)
      LOG_ERROR("Error attaching data array to field {}",
                LyrThickTargetFldName);

} // end defineFields

//------------------------------------------------------------------------------
// Read desired quantities from the mesh file
void VertCoord::readArrays(const Decomp *Decomp //< [in] Decomp for mesh
) {

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
               "using MaxLayerCell = NVertLayers - 1");
      deepCopy(TmpArrayI4H, NVertLayers - 1);
   } else {
      for (int ICell = 0; ICell < NCellsAll; ++ICell) {
         TmpArrayI4H(ICell) = TmpArrayI4H(ICell) - 1;
      }
   }
   TmpArrayI4H(NCellsAll) = -1;

   deepCopy(MaxLayerCellH, TmpArrayI4H);
   deepCopy(MaxLayerCell, TmpArrayI4H);

   const std::string MinNameMPAS = "minLevelCell";
   I4 MinReadErr = IO::readArray(TmpArrayI4H.data(), NCellsAll, MinNameMPAS,
                                 MeshFileID, CellDecompI4, ArrayID);

   if (MinReadErr != 0) {
      LOG_WARN("VertCoord: error reading minLevelCell from mesh file, "
               "using MinLayerCell = 0");
      deepCopy(TmpArrayI4H, 0);
   } else {
      for (int ICell = 0; ICell < NCellsAll; ++ICell) {
         TmpArrayI4H(ICell) = TmpArrayI4H(ICell) - 1;
      }
   }
   TmpArrayI4H(NCellsAll) = NVertLayersP1;

   deepCopy(MinLayerCellH, TmpArrayI4H);
   deepCopy(MinLayerCell, TmpArrayI4H);

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
   deepCopy(BottomDepth, BottomDepthH);

   DecErr = IO::destroyDecomp(CellDecompI4);
   if (DecErr != 0) {
      LOG_CRITICAL("VertCoord: error destroying cell I4 IO decomposition");
   }
   DecErr = IO::destroyDecomp(CellDecompR8);
   if (DecErr != 0) {
      LOG_CRITICAL("VertCoord: error destroying cell R8 IO decomposition");
   }
} // end readArrays

//------------------------------------------------------------------------------
// Destroys a local VertCoord and deallocates all arrays
VertCoord::~VertCoord() {

   int Err;
   Err = FieldGroup::destroy(GroupName);
   if (Err != 0)
      LOG_ERROR("Error removing FieldGrup {}", GroupName);
   Err = Field::destroy(PressInterfFldName);
   if (Err != 0)
      LOG_ERROR("Error removing Field {}", PressInterfFldName);
   Err = Field::destroy(PressMidFldName);
   if (Err != 0)
      LOG_ERROR("Error removing Field {}", PressMidFldName);
   Err = Field::destroy(ZInterfFldName);
   if (Err != 0)
      LOG_ERROR("Error removing Field {}", ZInterfFldName);
   Err = Field::destroy(ZMidFldName);
   if (Err != 0)
      LOG_ERROR("Error removing Field {}", ZMidFldName);
   Err = Field::destroy(GeopotFldName);
   if (Err != 0)
      LOG_ERROR("Error removing Field {}", GeopotFldName);
   Err = Field::destroy(LyrThickTargetFldName);
   if (Err != 0)
      LOG_ERROR("Error removing Field {}", LyrThickTargetFldName);

} // end destructor

//------------------------------------------------------------------------------
// Removes a VertCoord from map by name
void VertCoord::erase(std::string Name) {
   AllVertCoords.erase(Name); // removes the VertCoord from the list and in the
                              // process, calls the destructor
} // end erase

//------------------------------------------------------------------------------
// Removes all VertCoords to clean up before exit
void VertCoord::clear() {

   AllVertCoords.clear(); // removes all VertCoords from the list and in the
                          // process, calls the destructors for each

} // end clear

//------------------------------------------------------------------------------
// Compute min and max layer indices for edges based on MinLayerCell and
// MaxLayerCell
void VertCoord::minMaxLayerEdge() {

   MinLayerEdgeTop = Array1DI4("MinLayerEdgeTop", NEdgesSize);
   MinLayerEdgeBot = Array1DI4("MinLayerEdgeBot", NEdgesSize);
   MaxLayerEdgeTop = Array1DI4("MaxLayerEdgeTop", NEdgesSize);
   MaxLayerEdgeBot = Array1DI4("MaxLayerEdgeBot", NEdgesSize);

   OMEGA_SCOPE(LocNVertLayersP1, NVertLayersP1);
   OMEGA_SCOPE(LocCellsOnEdge, CellsOnEdge);
   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);
   OMEGA_SCOPE(LocMinLayerEdgeTop, MinLayerEdgeTop);
   OMEGA_SCOPE(LocMinLayerEdgeBot, MinLayerEdgeBot);
   OMEGA_SCOPE(LocMaxLayerEdgeTop, MaxLayerEdgeTop);
   OMEGA_SCOPE(LocMaxLayerEdgeBot, MaxLayerEdgeBot);
   parallelFor(
       {NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
          I4 Lyr1;
          I4 Lyr2;
          const I4 ICell1 = LocCellsOnEdge(IEdge, 0);
          const I4 ICell2 = LocCellsOnEdge(IEdge, 1);
          Lyr1            = LocMaxLayerCell(ICell1) == -1 ? LocNVertLayersP1
                                                          : LocMinLayerCell(ICell1);
          Lyr2            = LocMaxLayerCell(ICell2) == -1 ? LocNVertLayersP1
                                                          : LocMinLayerCell(ICell2);
          LocMinLayerEdgeTop(IEdge) = Kokkos::min(Lyr1, Lyr2);

          Lyr1 = LocMaxLayerCell(ICell1) == -1 ? 0 : LocMinLayerCell(ICell1);
          Lyr2 = LocMaxLayerCell(ICell2) == -1 ? 0 : LocMinLayerCell(ICell2);
          LocMinLayerEdgeBot(IEdge) = Kokkos::max(Lyr1, Lyr2);

          LocMaxLayerEdgeTop(IEdge) =
              Kokkos::min(LocMaxLayerCell(ICell1), LocMaxLayerCell(ICell2));
          LocMaxLayerEdgeBot(IEdge) =
              Kokkos::max(LocMaxLayerCell(ICell1), LocMaxLayerCell(ICell2));
       });

   OMEGA_SCOPE(LocNEdgesAll, NEdgesAll);
   parallelFor(
       {1}, KOKKOS_LAMBDA(const int &) {
          LocMinLayerEdgeTop(LocNEdgesAll) = LocNVertLayersP1;
          LocMinLayerEdgeBot(LocNEdgesAll) = LocNVertLayersP1;
          LocMaxLayerEdgeTop(LocNEdgesAll) = -1;
          LocMaxLayerEdgeBot(LocNEdgesAll) = -1;
       });

   MinLayerEdgeTopH = createHostMirrorCopy(MinLayerEdgeTop);
   MinLayerEdgeBotH = createHostMirrorCopy(MinLayerEdgeBot);
   MaxLayerEdgeTopH = createHostMirrorCopy(MaxLayerEdgeTop);
   MaxLayerEdgeBotH = createHostMirrorCopy(MaxLayerEdgeBot);
} // end MinMaxLayerEdge

//------------------------------------------------------------------------------
// Compute min and max layer indices for vertices based on MinLayerCell and
// MaxLayerCell
void VertCoord::minMaxLayerVertex() {

   MinLayerVertexTop = Array1DI4("MinLayerVertexTop", NVerticesSize);
   MinLayerVertexBot = Array1DI4("MinLayerVertexBot", NVerticesSize);
   MaxLayerVertexTop = Array1DI4("MaxLayerVertexTop", NVerticesSize);
   MaxLayerVertexBot = Array1DI4("MaxLayerVertexBot", NVerticesSize);

   OMEGA_SCOPE(LocNVertLayersP1, NVertLayersP1);
   OMEGA_SCOPE(LocVertexDegree, VertexDegree);
   OMEGA_SCOPE(LocCellsOnVertex, CellsOnVertex);
   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);
   OMEGA_SCOPE(LocMinLayerVertexTop, MinLayerVertexTop);
   OMEGA_SCOPE(LocMinLayerVertexBot, MinLayerVertexBot);
   OMEGA_SCOPE(LocMaxLayerVertexTop, MaxLayerVertexTop);
   OMEGA_SCOPE(LocMaxLayerVertexBot, MaxLayerVertexBot);

   parallelFor(
       {NVerticesAll}, KOKKOS_LAMBDA(int IVertex) {
          I4 Lyr;
          I4 ICell = LocCellsOnVertex(IVertex, 0);
          Lyr      = LocMaxLayerCell(ICell) == -1 ? 0 : LocMinLayerCell(ICell);
          LocMinLayerVertexBot(IVertex) = Lyr;
          for (int I = 1; I < LocVertexDegree; ++I) {
             ICell = LocCellsOnVertex(IVertex, I);
             Lyr   = LocMaxLayerCell(ICell) == -1 ? 0 : LocMinLayerCell(ICell);
             LocMinLayerVertexBot(IVertex) =
                 Kokkos::max(LocMinLayerVertexBot(IVertex), Lyr);
          }

          ICell = LocCellsOnVertex(IVertex, 0);
          Lyr   = LocMaxLayerCell(ICell) == -1 ? LocNVertLayersP1
                                               : LocMinLayerCell(ICell);
          LocMinLayerVertexTop(IVertex) = Lyr;
          for (int I = 1; I < LocVertexDegree; ++I) {
             ICell = LocCellsOnVertex(IVertex, I);
             Lyr   = LocMaxLayerCell(ICell) == -1 ? LocNVertLayersP1
                                                  : LocMinLayerCell(ICell);
             LocMinLayerVertexTop(IVertex) =
                 Kokkos::min(LocMinLayerVertexTop(IVertex), Lyr);
          }

          ICell                         = LocCellsOnVertex(IVertex, 0);
          LocMaxLayerVertexBot(IVertex) = LocMaxLayerCell(ICell);
          for (int I = 1; I < LocVertexDegree; ++I) {
             ICell                         = LocCellsOnVertex(IVertex, I);
             LocMaxLayerVertexBot(IVertex) = Kokkos::max(
                 LocMaxLayerVertexBot(IVertex), LocMaxLayerCell(ICell));
          }

          ICell                         = LocCellsOnVertex(IVertex, 0);
          LocMaxLayerVertexTop(IVertex) = LocMaxLayerCell(ICell);
          for (int I = 1; I < LocVertexDegree; ++I) {
             ICell                         = LocCellsOnVertex(IVertex, I);
             LocMaxLayerVertexTop(IVertex) = Kokkos::min(
                 LocMaxLayerVertexTop(IVertex), LocMaxLayerCell(ICell));
          }
       });
   OMEGA_SCOPE(LocNVerticesAll, NVerticesAll);
   parallelFor(
       {1}, KOKKOS_LAMBDA(const int &) {
          LocMinLayerVertexTop(LocNVerticesAll) = LocNVertLayersP1;
          LocMinLayerVertexBot(LocNVerticesAll) = LocNVertLayersP1;
          LocMaxLayerVertexTop(LocNVerticesAll) = -1;
          LocMaxLayerVertexBot(LocNVerticesAll) = -1;
       });

   MinLayerVertexTopH = createHostMirrorCopy(MinLayerVertexTop);
   MinLayerVertexBotH = createHostMirrorCopy(MinLayerVertexBot);
   MaxLayerVertexTopH = createHostMirrorCopy(MaxLayerVertexTop);
   MaxLayerVertexBotH = createHostMirrorCopy(MaxLayerVertexBot);
} // end MinMaxLayerVertex

//------------------------------------------------------------------------------
// Store VertCoordMovementWeights based on config choice
void VertCoord::initMovementWeights(
    Config *Options // [in] configuration options
) {

   Error Err; // default successful error code

   Config VCoordConfig("VertCoord");
   Err += Options->get(VCoordConfig);
   CHECK_ERROR_ABORT(Err, "VertCoord: VertCoord group not found in Config");

   std::string MovementWeightType;
   Err += VCoordConfig.get("MovementWeightType", MovementWeightType);
   CHECK_ERROR_ABORT(Err,
                     "VertCoord: MovementWeightType not found in VertCoord");

   VertCoordMovementWeights =
       Array1DReal("VertCoordMovementWeights", NVertLayers);

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
// Compute the pressure at each layer interface and midpoint given the
// LayerThickness and SurfacePressure. Hierarchical parallelism is used with a
// parallel_for loop over all cells and a parallel_scan performing a prefix sum
// in each column to compute pressure from the top-most active layer to the
// bottom-most active layer.
void VertCoord::computePressure(
    const Array2DReal &PressureInterface, // [out] pressure at layer interfaces
    const Array2DReal &PressureMid,       // [out] pressure at layer midpoints
    const Array2DReal &LayerThickness,    // [in] pseudo thickness
    const Array1DReal &SurfacePressure    // [in] surface pressure
) {

   Real Gravity = 9.80616_Real; // gravitationl acceleration
   Real Rho0    = 1035._Real;   // reference density

   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);

   const auto Policy = TeamPolicy(NCellsAll, OMEGA_TEAMSIZE, 1);
   Kokkos::parallel_for(
       "computePressure", Policy, KOKKOS_LAMBDA(const TeamMember &Member) {
          const I4 ICell = Member.league_rank();
          const I4 KMin  = LocMinLayerCell(ICell);
          const I4 KMax  = LocMaxLayerCell(ICell);
          const I4 Range = KMax - KMin + 1;

          PressureInterface(ICell, KMin) = SurfacePressure(ICell);
          Kokkos::parallel_scan(
              TeamThreadRange(Member, Range),
              [=](int K, Real &Accum, bool IsFinal) {
                 const I4 KLyr  = K + KMin;
                 Real Increment = Gravity * Rho0 * LayerThickness(ICell, KLyr);
                 Accum += Increment;

                 if (IsFinal) {
                    PressureInterface(ICell, KLyr + 1) =
                        SurfacePressure(ICell) + Accum;
                    PressureMid(ICell, KLyr) =
                        SurfacePressure(ICell) + Accum - 0.5 * Increment;
                 }
              });
       });
} // end computePressure

//------------------------------------------------------------------------------
// Compute geometric height z at layer interfaces and midpoints given the
// LayerThickness, SpecVol, and BottomDepth. Hierarchical parallelism is used
// with a parallel_for loop over cells and a parallel_scan performing a prefix
// sum in each column to compute z from the bottom-most active layer to the
// top-most active layer
void VertCoord::computeZHeight(
    const Array2DReal &ZInterface,     // [out] Z coord at layer interfaces
    const Array2DReal &ZMid,           // [out] Z coord at layer midpoints
    const Array2DReal &LayerThickness, // [in] pseudo thickness
    const Array2DReal &SpecVol,        // [in] specific volume
    const Array1DReal &BottomDepth     // [in] bottom depth
) {

   Real Rho0 = 1035._Real; // reference density

   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);

   const auto Policy = TeamPolicy(NCellsAll, OMEGA_TEAMSIZE, 1);
   Kokkos::parallel_for(
       "computeZHeight", Policy, KOKKOS_LAMBDA(const TeamMember &Member) {
          const I4 ICell = Member.league_rank();
          const I4 KMin  = LocMinLayerCell(ICell);
          const I4 KMax  = LocMaxLayerCell(ICell);
          const I4 Range = KMax - KMin + 1;

          ZInterface(ICell, KMax + 1) = -BottomDepth(ICell);
          Kokkos::parallel_scan(
              TeamThreadRange(Member, Range),
              [=](int K, Real &Accum, bool IsFinal) {
                 const I4 KLyr = KMax - K;
                 Real DZ =
                     Rho0 * SpecVol(ICell, KLyr) * LayerThickness(ICell, KLyr);
                 Accum += DZ;
                 if (IsFinal) {
                    ZInterface(ICell, KLyr) = -BottomDepth(ICell) + Accum;
                    ZMid(ICell, KLyr) = -BottomDepth(ICell) + Accum - 0.5 * DZ;
                 }
              });
       });
} // end computeZHeight

//------------------------------------------------------------------------------
// Compute geopotential given Zmid, TidalPotential, and SelfAttractionLoading.
// Nested parallel_fors loop over all cells and all active layers in a column to
// compute the geopotential at the midpoint of each layer. The tidal potential
// and SAL are configurable, default-off features. When off these arrays will
// just be zeroes.
void VertCoord::computeGeopotential(
    const Array2DReal &GeopotentialMid,      // [out] geopotential
    const Array2DReal &ZMid,                 // [in] Z coord at layer midpoints
    const Array1DReal &TidalPotential,       // [in] tidal potential
    const Array1DReal &SelfAttractionLoading // [in] self attraction and loading
) {

   Real Gravity = 9.80616_Real; // gravitationl acceleration

   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);

   Kokkos::parallel_for(
       "computeGeopotential", TeamPolicy(NCellsAll, OMEGA_TEAMSIZE),
       KOKKOS_LAMBDA(const TeamMember &Member) {
          const I4 ICell   = Member.league_rank();
          const I4 KMin    = LocMinLayerCell(ICell);
          const I4 KMax    = LocMaxLayerCell(ICell);
          const I4 KRange  = KMax - KMin + 1;
          const I4 NChunks = (KRange + VecLength - 1) / VecLength;
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(Member, NChunks), [=](const int KChunk) {
                 const I4 KStart = KMin + KChunk * VecLength;
                 const I4 KEnd   = KStart + VecLength;

                 const I4 KLen =
                     KEnd > KMax + 1 ? KMax + 1 - KStart : VecLength;
                 for (int KVec = 0; KVec < KLen; ++KVec) {
                    const I4 K                = KStart + KVec;
                    GeopotentialMid(ICell, K) = Gravity * ZMid(ICell, K) +
                                                TidalPotential(ICell) +
                                                SelfAttractionLoading(ICell);
                 }
              });
       });
} // end compute Geopotential

//------------------------------------------------------------------------------
// Compute the desired target thickness, given PressureInterface,
// RefLayerThickness, and VertCoordMovementWeights. Hierarchical parallelsim is
// used with an outer parallel_for loop over cells, and 2 paralel_reduce
// reductions and a parallel_for over the active layers within a column.
void VertCoord::computeTargetThickness(
    const Array2DReal &LayerThicknessTarget, // [out] desired target thickness
    const Array2DReal &PressureInterface, // [in] pressure at layer interfaces
    const Array2DReal &RefLayerThickness, // [in] reference pseudo thickness
    const Array1DReal &VertCoordMovementWeights // [in] movement weights
) {

   Real Gravity = 9.80616_Real; // gravitationl acceleration
   Real Rho0    = 1035._Real;   // reference density

   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);

   Kokkos::parallel_for(
       "computeTargetThickness", TeamPolicy(NCellsAll, OMEGA_TEAMSIZE),
       KOKKOS_LAMBDA(const TeamMember &Member) {
          const I4 ICell = Member.league_rank();
          const I4 KMin  = LocMinLayerCell(ICell);
          const I4 KMax  = LocMaxLayerCell(ICell);

          Real Coeff = (PressureInterface(ICell, KMax + 1) -
                        PressureInterface(ICell, KMin)) /
                       (Gravity * Rho0);

          Real SumWh = 0;
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange(Member, KMin, KMax + 1),
              [=](const int K, Real &LocalWh) {
                 LocalWh +=
                     VertCoordMovementWeights(K) * RefLayerThickness(ICell, K);
              },
              SumWh);

          Real SumRefH = 0;
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange(Member, KMin, KMax + 1),
              [=](const int K, Real &LocalSum) {
                 LocalSum += RefLayerThickness(ICell, K);
              },
              SumRefH);
          Coeff -= SumRefH;

          const I4 KRange  = KMax - KMin + 1;
          const I4 NChunks = (KRange + VecLength - 1) / VecLength;

          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(Member, NChunks), [=](const int KChunk) {
                 const I4 KStart = KMin + KChunk * VecLength;
                 const I4 KEnd   = KStart + VecLength;

                 const I4 KLen =
                     KEnd > KMax + 1 ? KMax + 1 - KStart : VecLength;
                 for (int KVec = 0; KVec < KLen; ++KVec) {
                    const I4 K = KStart + KVec;
                    LayerThicknessTarget(ICell, K) =
                        RefLayerThickness(ICell, K) *
                        (1._Real + Coeff * VertCoordMovementWeights(K) / SumWh);
                 }
              });
       });
}

//------------------------------------------------------------------------------
// Get default VertCoord
VertCoord *VertCoord::getDefault() { return VertCoord::DefaultVertCoord; }

//------------------------------------------------------------------------------
// Get VertCoord by name
VertCoord *VertCoord::get(const std::string Name ///< [in] Name of VertCoord
) {

   // look for an instance of this name
   auto it = AllVertCoords.find(Name);

   // if found, return the VertCoord pointer
   if (it != AllVertCoords.end()) {
      return it->second.get();

      // otherwise print error and return null pointer
   } else {
      LOG_ERROR("VertCoord::get: Attempt to retrieve non-existant VertCoord:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }

} // end get VertCoord

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
