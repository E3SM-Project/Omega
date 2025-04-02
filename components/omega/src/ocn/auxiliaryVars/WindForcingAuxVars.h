#ifndef OMEGA_AUX_WIND_H
#define OMEGA_AUX_WIND_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"

#include <string>

namespace OMEGA {

class WindForcingAuxVars {
 public:
   Array1DReal NormalStressEdge;
   Array1DReal ZonalStressCell;
   Array1DReal MeridStressCell;

   WindForcingAuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh,
                      int NVertLevels);

   KOKKOS_FUNCTION void computeVarsOnEdge(int IEdge, int KChunk) const {
      if (KChunk == 0) {
         const int JCell0 = CellsOnEdge(IEdge, 0);
         const int JCell1 = CellsOnEdge(IEdge, 1);
         const Real ZonalStressEdge =
             0.5_Real * (ZonalStressCell(JCell0) + ZonalStressCell(JCell1));
         const Real MeridStressEdge =
             0.5_Real * (MeridStressCell(JCell0) + MeridStressCell(JCell1));

         NormalStressEdge(IEdge) =
             Kokkos::cos(AngleEdge(IEdge)) * ZonalStressEdge +
             Kokkos::sin(AngleEdge(IEdge)) * MeridStressEdge;
      }
   }

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal AngleEdge;
};

} // namespace OMEGA
#endif
