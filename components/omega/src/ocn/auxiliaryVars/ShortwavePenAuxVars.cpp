#include "ShortwavePenAuxVars.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

ShortwavePenAuxVars::ShortwavePenAuxVars(const std::string &AuxStateSuffix,
                                         const HorzMesh *Mesh)
    : ExtinctionCoeffRedCell("ExtinctionCoeffRedCell" + AuxStateSuffix,
                             Mesh->NCellsSize),
      ExtinctionCoeffBlueCell("ExtinctionCoeffBlueCell" + AuxStateSuffix,
                              Mesh->NCellsSize) {}

void ShortwavePenAuxVars::registerFields(const std::string &GroupName,
                                         const std::string &MeshName) const {
   const int NDims = 1;
   std::string DimSuffix;
   if (MeshName != "Default") {
      DimSuffix = MeshName;
   }
   const std::vector<std::string> DimNames = {"NCells" + DimSuffix};

   auto ExtinctionCoeffRedCellField = Field::create(
       ExtinctionCoeffRedCell.label(), "red-band extinction coefficient",
       "m^-1", "", std::numeric_limits<Real>::min(),
       std::numeric_limits<Real>::max(), NDims, DimNames);
   auto ExtinctionCoeffBlueCellField = Field::create(
       ExtinctionCoeffBlueCell.label(), "blue-band extinction coefficient",
       "m^-1", "", std::numeric_limits<Real>::min(),
       std::numeric_limits<Real>::max(), NDims, DimNames);

   FieldGroup::addFieldToGroup(ExtinctionCoeffRedCell.label(), GroupName);
   FieldGroup::addFieldToGroup(ExtinctionCoeffBlueCell.label(), GroupName);

   ExtinctionCoeffRedCellField->attachData<Array1DReal>(ExtinctionCoeffRedCell);
   ExtinctionCoeffBlueCellField->attachData<Array1DReal>(
       ExtinctionCoeffBlueCell);
}

void ShortwavePenAuxVars::unregisterFields() const {
   Field::destroy(ExtinctionCoeffRedCell.label());
   Field::destroy(ExtinctionCoeffBlueCell.label());
}

} // namespace OMEGA
