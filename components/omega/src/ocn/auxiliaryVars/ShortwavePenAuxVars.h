#ifndef OMEGA_AUX_SHORTWAVE_PEN_H
#define OMEGA_AUX_SHORTWAVE_PEN_H

#include "DataTypes.h"
#include "HorzMesh.h"

#include <string>

namespace OMEGA {

class ShortwavePenAuxVars {
 public:
   Array1DReal ExtinctionCoeffRedCell;
   Array1DReal ExtinctionCoeffBlueCell;

   ShortwavePenAuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh);

   void registerFields(const std::string &GroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;
};

} // namespace OMEGA

#endif
