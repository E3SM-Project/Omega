#include "WindForcingAuxVars.h"
#include "DataTypes.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

WindForcingAuxVars::WindForcingAuxVars(const std::string &AuxStateSuffix,
                                       const HorzMesh *Mesh, int NVertLevels)
    : NormalWindEdge("NormalWindEdge" + AuxStateSuffix, Mesh->NEdgesSize),
      ZonalWindCell("ZonalWindCell" + AuxStateSuffix, Mesh->NCellsSize),
      MeridWindCell("MeridWindCell" + AuxStateSuffix, Mesh->NCellsSize),
      WindRelNormCell("WindRelNormCell" + AuxStateSuffix, Mesh->NCellsSize),
      NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      DvEdge(Mesh->DvEdge), AngleEdge(Mesh->AngleEdge),
      AreaCell(Mesh->AreaCell) {}

void WindForcingAuxVars::registerFields(
    const std::string &AuxGroupName, // name of Auxiliary field group
    const std::string &MeshName      // name of horizontal mesh
) const {

   int Err = 0; // error flag for some calls

   // Create fields
   const Real FillValue = -9.99e30;
   int NDims            = 1;
   std::vector<std::string> DimNames(NDims);
   std::string DimSuffix;
   if (MeshName == "Default") {
      DimSuffix = "";
   } else {
      DimSuffix = MeshName;
   }

   // Zonal wind
   DimNames[0] = "NCells" + DimSuffix;
   auto ZonalWindCellField =
       Field::create(ZonalWindCell.label(),            // field name
                     "zonal wind",                     // long name/describe
                     "m s^-1",                         // units
                     "",                               // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValue, // scalar for undefined entries
                     1,         // number of dimensions
                     DimNames   // dim names
       );

   // Meridional wind
   auto MeridWindCellField =
       Field::create(MeridWindCell.label(), // field name
                     "meridional wind",     // long Name or description
                     "m s^-1",              // units
                     "",                    // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValue, // scalar used for undefined entries
                     NDims,     // number of dimensions
                     DimNames   // dimension names
       );

   // Add fields to FieldGroup
   Err = FieldGroup::addFieldToGroup(ZonalWindCell.label(), AuxGroupName);
   if (Err != 0)
      LOG_ERROR("Error adding field {} to group {}", ZonalWindCell.label(),
                AuxGroupName);
   Err = FieldGroup::addFieldToGroup(MeridWindCell.label(), AuxGroupName);
   if (Err != 0)
      LOG_ERROR("Error adding field {} to group {}", MeridWindCell.label(),
                AuxGroupName);

   // Attach data
   Err = ZonalWindCellField->attachData<Array1DReal>(ZonalWindCell);
   if (Err != 0)
      LOG_ERROR("Error attaching data to field {}", ZonalWindCell.label());

   Err = MeridWindCellField->attachData<Array1DReal>(MeridWindCell);
   if (Err != 0)
      LOG_ERROR("Error attaching data to field {}", MeridWindCell.label());
}

void WindForcingAuxVars::unregisterFields() const {
   int Err = 0;
   Err     = Field::destroy(ZonalWindCell.label());
   if (Err != 0)
      LOG_ERROR("Error destroying field {}", ZonalWindCell.label());
   Err = Field::destroy(MeridWindCell.label());
   if (Err != 0)
      LOG_ERROR("Error destroying field {}", MeridWindCell.label());
}

} // namespace OMEGA
