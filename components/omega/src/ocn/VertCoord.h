#ifndef OMEGA_VERTCOORD_H
#define OMEGA_VERTCOORD_H
//===-- base/VertCoord.h - vertical coordinate --------------*- C++ -*-===//
//
/// \file
/// \brief
///
///
//
//===----------------------------------------------------------------------===//

#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Error.h"
#include "HorzMesh.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"

#include <memory>
#include <string>

namespace OMEGA {

class VertCoord {

 private:
   // Variables from Decomp
   I4 NCellsOwned;
   I4 NCellsAll;
   I4 NCellsSize;
   I4 NEdgesOwned;
   I4 NEdgesAll;
   I4 NEdgesSize;
   I4 NVerticesOwned;
   I4 NVerticesAll;
   I4 NVerticesSize;
   I4 VertexDegree;
   Array2DI4 CellsOnEdge;
   Array2DI4 CellsOnVertex;

   std::string MeshFileName;
   int MeshFileID;

   static VertCoord *DefaultVertCoord;
   static std::map<std::string, std::unique_ptr<VertCoord>> AllVertCoords;

   // methods

   void readArrays(const Decomp *MeshDecomp);

   /// construct a new vertical coordinate object
   VertCoord(const std::string &Name, const Decomp *MeshDecomp,
             Config *Options);

   // Forbid copy and move construction
   VertCoord(const VertCoord &) = delete;
   VertCoord(VertCoord &&)      = delete;

 public:
   I4 NVertLayers;
   I4 NVertLayersP1;

   // Variables computed
   Array2DReal PressureInterface;
   Array2DReal PressureMid;
   Array2DReal ZInterface;
   Array2DReal ZMid;
   Array2DReal GeopotentialMid;
   Array2DReal LayerThicknessTarget;

   HostArray2DReal PressureInterfaceH;
   HostArray2DReal PressureMidH;
   HostArray2DReal ZInterfaceH;
   HostArray2DReal ZMidH;
   HostArray2DReal GeopotentialMidH;
   HostArray2DReal LayerThicknessTargetH;

   // Vertical loop bounds
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeTop;
   Array1DI4 MaxLayerEdgeTop;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeBot;
   Array1DI4 MinLayerVertexTop;
   Array1DI4 MaxLayerVertexTop;
   Array1DI4 MinLayerVertexBot;
   Array1DI4 MaxLayerVertexBot;

   HostArray1DI4 MinLayerCellH;
   HostArray1DI4 MaxLayerCellH;
   HostArray1DI4 MinLayerEdgeTopH;
   HostArray1DI4 MaxLayerEdgeTopH;
   HostArray1DI4 MinLayerEdgeBotH;
   HostArray1DI4 MaxLayerEdgeBotH;
   HostArray1DI4 MinLayerVertexTopH;
   HostArray1DI4 MaxLayerVertexTopH;
   HostArray1DI4 MinLayerVertexBotH;
   HostArray1DI4 MaxLayerVertexBotH;

   void minMaxLayerEdge();
   void minMaxLayerVertex();

   void initMovementWeights(Config *Options);

   // p star coordinate variables
   Array1DReal VertCoordMovementWeights;
   Array2DReal RefLayerThickness;

   HostArray1DReal VertCoordMovementWeightsH;
   HostArray2DReal RefLayerThicknessH;

   // BottomDepth read from mesh file
   Array1DReal BottomDepth;

   HostArray1DReal BottomDepthH;

   std::string Name;

   // methods

   /// Initialize Omega vertical coordinate
   static void init();

   /// Creates a new vertical coordinate object by calling the constructor and
   /// puts it in the AllVertCoords map
   static VertCoord *create(const std::string &Name, const Decomp *MeshDecomp,
                            Config *Options);

   /// Destructor - deallocates all memory and deletes a VertCoord
   ~VertCoord();

   static void clear();

   /// Remove a VertCoord by name
   static void erase(std::string InName);

   /// Retrieve the default VertCoord
   static VertCoord *getDefault();

   /// Retreive a VertCoord by name
   static VertCoord *get(std::string name);

   /// Sums the mass thickness times g from the top layer down, starting with
   /// the surface pressure
   void computePressure(const Array2DReal &PressureInterface,
                        const Array2DReal &PressureMid,
                        const Array2DReal &LayerThickness,
                        const Array1DReal &SurfacePressure);

   /// Sum the mass thickness time specific volume from the bottom layer up,
   /// starting with the bottom elevation
   void computeZHeight(const Array2DReal &ZInterface, const Array2DReal &ZMid,
                       const Array2DReal &LayerThickness,
                       const Array2DReal &SpecVol,
                       const Array1DReal &BottomDepth);

   /// Sum the z height times g, the tidal potential, and self attraction and
   /// loading
   void computeGeopotential(const Array2DReal &GeopotentialMid,
                            const Array2DReal &ZMid,
                            const Array1DReal &TidalPotential,
                            const Array1DReal &SelfAttractionLoading);

   /// Determine mass thickness used for the target vertical coordinate
   void computeTargetThickness(const Array2DReal &LayerThicknessTarget,
                               const Array2DReal &PressureInterface,
                               const Array2DReal &RefLayerThickness,
                               const Array1DReal &VertCoordMovementWeights);

}; // end class VertCoord

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_VERTCOORD_H
