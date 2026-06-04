#ifndef OMEGA_FORCING_H
#define OMEGA_FORCING_H

#include "Config.h"
#include "DataTypes.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "VertCoord.h"
#include "auxiliaryVars/MomForcingAuxVars.h"
#include "auxiliaryVars/SurfTracerRestAuxVars.h"
#include "auxiliaryVars/TracerForcingAuxVars.h"

#include <memory>
#include <string>

namespace OMEGA {

class Forcing {
 public:
   std::string Name;

   TracerForcingAuxVars TracerForcingAux;
   SurfTracerRestAuxVars SurfTracerRestAux;
   MomForcingAuxVars MomForcingAux;

   ~Forcing();

   static void init();

   static Forcing *create(const std::string &Name, const HorzMesh *Mesh,
                          Halo *MeshHalo, const VertCoord *VCoord,
                          int NTracers);

   static Forcing *getDefault();

   static Forcing *get(const std::string &Name);

   static bool exists(const std::string &Name);

   static void erase(const std::string &Name);

   static void clear();

   void readConfigOptions(Config *OmegaConfig);

   void registerFields(const std::string &MeshName) const;

   void unregisterFields() const;

   void computeSrfStressForcingOnEdge() const;

   void computeSurfInsituTemp(const Array3DReal &TracerArray) const;

   I4 exchangeHalo() const;

 private:
   Forcing(const std::string &Name, const HorzMesh *Mesh, Halo *MeshHalo,
           const VertCoord *VCoord, int NTracers);

   Forcing(const Forcing &) = delete;
   Forcing(Forcing &&)      = delete;

   const HorzMesh *Mesh;
   Halo *MeshHalo;
   const VertCoord *VCoord;

   static Forcing *DefaultForcing;
   static std::map<std::string, std::unique_ptr<Forcing>> AllForcing;
};

} // namespace OMEGA

#endif
