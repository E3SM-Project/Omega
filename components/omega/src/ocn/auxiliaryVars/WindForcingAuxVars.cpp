#include "WindForcingAuxVars.h"
#include "DataTypes.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

WindForcingAuxVars::WindForcingAuxVars(const std::string &AuxStateSuffix,
                                       const HorzMesh *Mesh, int NVertLevels)
    : NormalStressEdge("NormalStressEdge" + AuxStateSuffix, Mesh->NEdgesSize),
      ZonalStressCell("WindStressZonal" + AuxStateSuffix, Mesh->NCellsSize),
      MeridStressCell("WindStressMeridional" + AuxStateSuffix,
                      Mesh->NCellsSize),
      CellsOnEdge(Mesh->CellsOnEdge), AngleEdge(Mesh->AngleEdge), Interp(Mesh) {
}

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

   // Zonal wind stress
   DimNames[0] = "NCells" + DimSuffix;
   auto ZonalStressCellField =
       Field::create(ZonalStressCell.label(),          // field name
                     "zonal wind stress",              // long name/describe
                     "N m^{-2}",                       // units
                     "",                               // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValue, // scalar for undefined entries
                     NDims,     // number of dimensions
                     DimNames   // dim names
       );

   // Meridional wind stress
   auto MeridStressCellField =
       Field::create(MeridStressCell.label(),  // field name
                     "meridional wind stress", // long Name or description
                     "N m^{-2}",               // units
                     "",                       // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValue, // scalar used for undefined entries
                     NDims,     // number of dimensions
                     DimNames   // dimension names
       );

   // Add fields to FieldGroup
   Err = FieldGroup::addFieldToGroup(ZonalStressCell.label(), AuxGroupName);
   if (Err != 0)
      LOG_ERROR("Error adding field {} to group {}", ZonalStressCell.label(),
                AuxGroupName);
   Err = FieldGroup::addFieldToGroup(MeridStressCell.label(), AuxGroupName);
   if (Err != 0)
      LOG_ERROR("Error adding field {} to group {}", MeridStressCell.label(),
                AuxGroupName);

   // Attach data
   Err = ZonalStressCellField->attachData<Array1DReal>(ZonalStressCell);
   if (Err != 0)
      LOG_ERROR("Error attaching data to field {}", ZonalStressCell.label());

   Err = MeridStressCellField->attachData<Array1DReal>(MeridStressCell);
   if (Err != 0)
      LOG_ERROR("Error attaching data to field {}", MeridStressCell.label());
}

void WindForcingAuxVars::unregisterFields() const {
   int Err = 0;
   Err     = Field::destroy(ZonalStressCell.label());
   if (Err != 0)
      LOG_ERROR("Error destroying field {}", ZonalStressCell.label());
   Err = Field::destroy(MeridStressCell.label());
   if (Err != 0)
      LOG_ERROR("Error destroying field {}", MeridStressCell.label());
}

} // namespace OMEGA
