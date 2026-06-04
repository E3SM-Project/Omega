#ifndef OMEGA_AUX_TRACER_FORCING_H
#define OMEGA_AUX_TRACER_FORCING_H

#include "DataTypes.h"
#include "HorzMesh.h"

#include <string>

namespace OMEGA {

// Forward declarations. Full definitions not needed in this header since only
// pointers are used.
class VertCoord;
class Eos;

class TracerForcingAuxVars {
 public:
   Array1DReal SnowFluxCell;
   Array1DReal RainFluxCell;
   Array1DReal EvaporationFluxCell;
   Array1DReal SeaIceFreshWaterFluxCell;
   Array1DReal IceRunoffFluxCell;
   Array1DReal RiverRunoffFluxCell;

   Array1DReal LatentHeatFluxCell;
   Array1DReal SensibleHeatFluxCell;
   Array1DReal LongWaveHeatFluxUpCell;
   Array1DReal LongWaveHeatFluxDownCell;
   Array1DReal SeaIceHeatFluxCell;
   Array1DReal ShortWaveHeatFluxCell;

   Array1DReal SeaIceSaltFluxCell;

   Array1DReal SurfInsituTemperature;

   TracerForcingAuxVars(const std::string &AuxStateSuffix,
                        const HorzMesh *Mesh);

   void registerFields(const std::string &MeshName) const;
   void unregisterFields() const;

   /// Compute surface insitu temperature from conservative temperature
   void computeSurfInsituTemp(const Array3DReal &TracerArray,
                              const VertCoord *VCoord,
                              const Eos *EosInst) const;
};

} // namespace OMEGA

#endif
