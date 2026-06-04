#include "TracerForcingAuxVars.h"
#include "Eos.h"
#include "Field.h"
#include "Tracers.h"
#include "VertCoord.h"

#include <limits>

namespace OMEGA {

TracerForcingAuxVars::TracerForcingAuxVars(const std::string &AuxStateSuffix,
                                           const HorzMesh *Mesh)
    : SnowFluxCell("snowFlux" + AuxStateSuffix, Mesh->NCellsSize),
      RainFluxCell("rainFlux" + AuxStateSuffix, Mesh->NCellsSize),
      EvaporationFluxCell("evaporationFlux" + AuxStateSuffix, Mesh->NCellsSize),
      SeaIceFreshWaterFluxCell("seaIceFreshWaterFlux" + AuxStateSuffix,
                               Mesh->NCellsSize),
      IceRunoffFluxCell("iceRunoffFlux" + AuxStateSuffix, Mesh->NCellsSize),
      RiverRunoffFluxCell("riverRunoffFlux" + AuxStateSuffix, Mesh->NCellsSize),
      LatentHeatFluxCell("latentHeatFlux" + AuxStateSuffix, Mesh->NCellsSize),
      SensibleHeatFluxCell("sensibleHeatFlux" + AuxStateSuffix,
                           Mesh->NCellsSize),
      LongWaveHeatFluxUpCell("longWaveHeatFluxUp" + AuxStateSuffix,
                             Mesh->NCellsSize),
      LongWaveHeatFluxDownCell("longWaveHeatFluxDown" + AuxStateSuffix,
                               Mesh->NCellsSize),
      SeaIceHeatFluxCell("seaIceHeatFlux" + AuxStateSuffix, Mesh->NCellsSize),
      ShortWaveHeatFluxCell("shortWaveHeatFlux" + AuxStateSuffix,
                            Mesh->NCellsSize),
      SeaIceSaltFluxCell("seaIceSalinityFlux" + AuxStateSuffix,
                         Mesh->NCellsSize),
      SurfInsituTemperature("surfInsituTemperature" + AuxStateSuffix,
                            Mesh->NCellsSize) {
   deepCopy(SnowFluxCell, 0.0_Real);
   deepCopy(RainFluxCell, 0.0_Real);
   deepCopy(EvaporationFluxCell, 0.0_Real);
   deepCopy(SeaIceFreshWaterFluxCell, 0.0_Real);
   deepCopy(IceRunoffFluxCell, 0.0_Real);
   deepCopy(RiverRunoffFluxCell, 0.0_Real);
   deepCopy(LatentHeatFluxCell, 0.0_Real);
   deepCopy(SensibleHeatFluxCell, 0.0_Real);
   deepCopy(LongWaveHeatFluxUpCell, 0.0_Real);
   deepCopy(LongWaveHeatFluxDownCell, 0.0_Real);
   deepCopy(SeaIceHeatFluxCell, 0.0_Real);
   deepCopy(ShortWaveHeatFluxCell, 0.0_Real);
   deepCopy(SeaIceSaltFluxCell, 0.0_Real);
   deepCopy(SurfInsituTemperature, 0.0_Real);
}

void TracerForcingAuxVars::registerFields(const std::string &MeshName) const {
   const Real FillValue = -9.99e30;
   const int NDims      = 1;
   std::vector<std::string> DimNames(NDims);

   std::string DimSuffix;
   if (MeshName == "Default") {
      DimSuffix = "";
   } else {
      DimSuffix = MeshName;
   }

   DimNames[0] = "NCells" + DimSuffix;

   auto SnowFluxField = Field::create(
       SnowFluxCell.label(), "snow freshwater flux", "kg m^-2 s^-1", "",
       std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::max(),
       FillValue, NDims, DimNames);
   auto RainFluxField = Field::create(
       RainFluxCell.label(), "rain freshwater flux", "kg m^-2 s^-1", "",
       std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::max(),
       FillValue, NDims, DimNames);
   auto EvaporationFluxField = Field::create(
       EvaporationFluxCell.label(), "evaporation freshwater flux",
       "kg m^-2 s^-1", "", std::numeric_limits<Real>::lowest(),
       std::numeric_limits<Real>::max(), FillValue, NDims, DimNames);
   auto SeaIceFreshWaterFluxField = Field::create(
       SeaIceFreshWaterFluxCell.label(), "sea-ice freshwater flux",
       "kg m^-2 s^-1", "", std::numeric_limits<Real>::lowest(),
       std::numeric_limits<Real>::max(), FillValue, NDims, DimNames);
   auto IceRunoffFluxField = Field::create(
       IceRunoffFluxCell.label(), "ice runoff freshwater flux", "kg m^-2 s^-1",
       "", std::numeric_limits<Real>::lowest(),
       std::numeric_limits<Real>::max(), FillValue, NDims, DimNames);
   auto RiverRunoffFluxField = Field::create(
       RiverRunoffFluxCell.label(), "river runoff freshwater flux",
       "kg m^-2 s^-1", "", std::numeric_limits<Real>::lowest(),
       std::numeric_limits<Real>::max(), FillValue, NDims, DimNames);

   auto LatentHeatFluxField = Field::create(
       LatentHeatFluxCell.label(), "latent heat flux", "W m^-2", "",
       std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::max(),
       FillValue, NDims, DimNames);
   auto SensibleHeatFluxField = Field::create(
       SensibleHeatFluxCell.label(), "sensible heat flux", "W m^-2", "",
       std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::max(),
       FillValue, NDims, DimNames);
   auto LongWaveHeatFluxUpField = Field::create(
       LongWaveHeatFluxUpCell.label(), "upward longwave heat flux", "W m^-2",
       "", std::numeric_limits<Real>::lowest(),
       std::numeric_limits<Real>::max(), FillValue, NDims, DimNames);
   auto LongWaveHeatFluxDownField = Field::create(
       LongWaveHeatFluxDownCell.label(), "downward longwave heat flux",
       "W m^-2", "", std::numeric_limits<Real>::lowest(),
       std::numeric_limits<Real>::max(), FillValue, NDims, DimNames);
   auto SeaIceHeatFluxField = Field::create(
       SeaIceHeatFluxCell.label(), "sea-ice heat flux", "W m^-2", "",
       std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::max(),
       FillValue, NDims, DimNames);
   auto ShortWaveHeatFluxField = Field::create(
       ShortWaveHeatFluxCell.label(), "shortwave heat flux", "W m^-2", "",
       std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::max(),
       FillValue, NDims, DimNames);

   auto SeaIceSaltFluxField = Field::create(
       SeaIceSaltFluxCell.label(), "sea-ice salt flux", "kg m^-2 s^-1", "",
       std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::max(),
       FillValue, NDims, DimNames);

   auto SurfInsituTemperatureField = Field::create(
       SurfInsituTemperature.label(),
       "insitu (potential) temperature at surface layer", "degrees Celsius", "",
       std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::max(),
       FillValue, NDims, DimNames);

   FieldGroup::addFieldToGroup(SnowFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(RainFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(EvaporationFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(SeaIceFreshWaterFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(IceRunoffFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(RiverRunoffFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(LatentHeatFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(SensibleHeatFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(LongWaveHeatFluxUpCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(LongWaveHeatFluxDownCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(SeaIceHeatFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(ShortWaveHeatFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(SeaIceSaltFluxCell.label(), "Forcing");

   SnowFluxField->attachData<Array1DReal>(SnowFluxCell);
   RainFluxField->attachData<Array1DReal>(RainFluxCell);
   EvaporationFluxField->attachData<Array1DReal>(EvaporationFluxCell);
   SeaIceFreshWaterFluxField->attachData<Array1DReal>(SeaIceFreshWaterFluxCell);
   IceRunoffFluxField->attachData<Array1DReal>(IceRunoffFluxCell);
   RiverRunoffFluxField->attachData<Array1DReal>(RiverRunoffFluxCell);
   LatentHeatFluxField->attachData<Array1DReal>(LatentHeatFluxCell);
   SensibleHeatFluxField->attachData<Array1DReal>(SensibleHeatFluxCell);
   LongWaveHeatFluxUpField->attachData<Array1DReal>(LongWaveHeatFluxUpCell);
   LongWaveHeatFluxDownField->attachData<Array1DReal>(LongWaveHeatFluxDownCell);
   SeaIceHeatFluxField->attachData<Array1DReal>(SeaIceHeatFluxCell);
   ShortWaveHeatFluxField->attachData<Array1DReal>(ShortWaveHeatFluxCell);
   SurfInsituTemperatureField->attachData<Array1DReal>(SurfInsituTemperature);
   SeaIceSaltFluxField->attachData<Array1DReal>(SeaIceSaltFluxCell);
}

void TracerForcingAuxVars::unregisterFields() const {
   Field::destroy(SnowFluxCell.label());
   Field::destroy(RainFluxCell.label());
   Field::destroy(EvaporationFluxCell.label());
   Field::destroy(SeaIceFreshWaterFluxCell.label());
   Field::destroy(IceRunoffFluxCell.label());
   Field::destroy(RiverRunoffFluxCell.label());
   Field::destroy(LatentHeatFluxCell.label());
   Field::destroy(SensibleHeatFluxCell.label());
   Field::destroy(LongWaveHeatFluxUpCell.label());
   Field::destroy(LongWaveHeatFluxDownCell.label());
   Field::destroy(SeaIceHeatFluxCell.label());
   Field::destroy(ShortWaveHeatFluxCell.label());
   Field::destroy(SeaIceSaltFluxCell.label());
   Field::destroy(SurfInsituTemperature.label());
}

void TracerForcingAuxVars::computeSurfInsituTemp(const Array3DReal &TracerArray,
                                                 const VertCoord *VCoord,
                                                 const Eos *EosInst) const {
   const int IndxTemp = Tracers::IndxTemp;
   const int IndxSalt = Tracers::IndxSalt;

   // Skip computation if temperature or salinity tracers are not defined
   if (IndxTemp < 0 || IndxSalt < 0) {
      return;
   }

   OMEGA_SCOPE(LocMinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(LocSurfInsituTemp, SurfInsituTemperature);

   int NCellsOwned = SurfInsituTemperature.extent_int(0);

   parallelFor(
       "TracerForcingAux:computeSurfInsituTemp", {NCellsOwned},
       KOKKOS_LAMBDA(int ICell) {
          const int KMin = LocMinLayerCell(ICell);
          const int KMax = LocMaxLayerCell(ICell);

          // Only compute for valid ocean cells
          if (KMin <= KMax) {
             const Real ConservTemp = TracerArray(IndxTemp, ICell, KMin);
             const Real AbsSalinity = TracerArray(IndxSalt, ICell, KMin);

             // Call EOS function to compute potential temperature from
             // conservative temperature at surface (reference pressure = 0)
             LocSurfInsituTemp(ICell) =
                 EosInst->calcPtFromCt(AbsSalinity, ConservTemp);
          }
       });
}
} // namespace OMEGA
