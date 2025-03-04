#ifndef OMEGA_AUX_WIND_H
#define OMEGA_AUX_WIND_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"

#include <string>

namespace OMEGA {

class WindForcingAuxVars {
 public:
   Array1DReal NormalWindEdge;
   Array1DReal WindRelNormCell;

   WindForcingAuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh,
                      int NVertLevels);

   KOKKOS_FUNCTION void
   computeVarsOnCell(int ICell, int KChunk,
                     const Array2DReal &NormalVelEdge) const {
      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real WindRelNormCellTmp = 0;

      if (KChunk == 0) {
         for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
            const int JEdge     = EdgesOnCell(ICell, J);
            const Real AreaEdge = 0.5_Real * DvEdge(JEdge) * DcEdge(JEdge);
            const Real WindRelVel =
                NormalVelEdge(JEdge, 0) - NormalWindEdge(JEdge);
            WindRelNormCellTmp +=
                AreaEdge * InvAreaCell * WindRelVel * WindRelVel;
         }
         WindRelNormCell(ICell) = Kokkos::sqrt(WindRelNormCellTmp);
      }
   }

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
};

} // namespace OMEGA
#endif
