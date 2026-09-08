#include "VelocityReconAuxVars.h"
#include "DataTypes.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

VelocityReconAuxVars::VelocityReconAuxVars(const std::string &AuxStateSuffix,
                                           const HorzMesh *Mesh,
                                           const VertCoord *VCoord)
    : VelocityZonalCell("VelocityZonalCell" + AuxStateSuffix, Mesh->NCellsSize,
                        VCoord->NVertLayers),
      VelocityMeridionalCell("VelocityMeridionalCell" + AuxStateSuffix,
                             Mesh->NCellsSize, VCoord->NVertLayers),
      OnSphere(Mesh->OnSphere), NEdgesReconOnCell(Mesh->NEdgesReconOnCell),
      ReconStencilCell(Mesh->ReconStencilCell),
      ReconWeightsCell(Mesh->ReconWeightsCell), LatCell(Mesh->LatCell),
      LonCell(Mesh->LonCell), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell) {}

void VelocityReconAuxVars::registerFields(
    const std::string &AuxGroupName, // name of Auxiliary field group
    const std::string &MeshName      // name of horizontal mesh
) const {

   // Create fields
   int NDims = 2;
   std::vector<std::string> DimNames(NDims);
   std::string DimSuffix;
   if (MeshName == "Default") {
      DimSuffix = "";
   } else {
      DimSuffix = MeshName;
   }

   DimNames[0] = "NCells" + DimSuffix;
   DimNames[1] = "NVertLayers";

   // Zonal velocity on cells
   auto VelocityZonalCellField = Field::create(
       VelocityZonalCell.label(),                      // field name
       "zonal velocity reconstructed at cell centers", // long name/describe
       "m s^-1",                                       // units
       "eastward_sea_water_velocity",                  // CF standard Name
       std::numeric_limits<Real>::lowest(),            // min valid value
       std::numeric_limits<Real>::max(),               // max valid value
       NDims,                                          // number of dimensions
       DimNames                                        // dimension names
   );

   // Meridional velocity on cells
   auto VelocityMeridionalCellField = Field::create(
       VelocityMeridionalCell.label(),                      // field name
       "meridional velocity reconstructed at cell centers", // long name
       "m s^-1",                                            // units
       "northward_sea_water_velocity",                      // CF standard Name
       std::numeric_limits<Real>::lowest(),                 // min valid value
       std::numeric_limits<Real>::max(),                    // max valid value
       NDims,   // number of dimensions
       DimNames // dimension names
   );

   // Add fields to FieldGroup
   FieldGroup::addFieldToGroup(VelocityZonalCell.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(VelocityMeridionalCell.label(), AuxGroupName);

   // Attach data
   VelocityZonalCellField->attachData<Array2DReal>(VelocityZonalCell);
   VelocityMeridionalCellField->attachData<Array2DReal>(VelocityMeridionalCell);
}

void VelocityReconAuxVars::unregisterFields() const {
   Field::destroy(VelocityZonalCell.label());
   Field::destroy(VelocityMeridionalCell.label());
}

} // namespace OMEGA
