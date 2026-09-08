#ifndef OMEGA_AUX_VELOCITY_RECON_H
#define OMEGA_AUX_VELOCITY_RECON_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <string>

namespace OMEGA {

/// Zonal and meridional velocity components reconstructed at cell centers
/// from the edge-normal velocity, using the least-squares stencil and
/// weights supplied by the mesh file. On a planar mesh the reconstructed
/// vector already lies in the plane of the mesh, so the two components are
/// the Cartesian x and y components instead.
///
/// These are diagnostic: nothing in the Omega equations uses them, so they
/// are computed once per time step rather than once per time stepper stage
/// (see AuxiliaryState::computeVelocityRecon).
///
/// This holds the same reconstruction as the VectorReconOnCell horizontal
/// operator, which is its reference implementation and is what the
/// operator unit test exercises, with the valid layer range of a column
/// taken into account.
class VelocityReconAuxVars {
 public:
   Array2DReal VelocityZonalCell;
   Array2DReal VelocityMeridionalCell;

   VelocityReconAuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh,
                        const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   computeVarsOnCell(int ICell, int K, const Array2DReal &NormalVelEdge) const {

      // Leave layers outside the valid range of this column at the fill
      // value that Field::attachData wrote
      if (K < MinLayerCell(ICell) || K > MaxLayerCell(ICell))
         return;

      // Accumulate the Cartesian components of the reconstructed vector
      Real Ux = 0._Real, Uy = 0._Real, Uz = 0._Real;

      for (int J = 0; J < NEdgesReconOnCell(ICell); ++J) {
         const I4 JEdge   = ReconStencilCell(ICell, J);
         const Real Field = NormalVelEdge(JEdge, K);

         Ux += ReconWeightsCell(ICell, 0, J) * Field;
         Uy += ReconWeightsCell(ICell, 1, J) * Field;
         Uz += ReconWeightsCell(ICell, 2, J) * Field;
      }

      if (OnSphere) {
         // cartesian to local geographic
         const Real CLat = Kokkos::cos(LatCell(ICell));
         const Real SLat = Kokkos::sin(LatCell(ICell));
         const Real CLon = Kokkos::cos(LonCell(ICell));
         const Real SLon = Kokkos::sin(LonCell(ICell));

         VelocityZonalCell(ICell, K) = -SLon * Ux + CLon * Uy;
         VelocityMeridionalCell(ICell, K) =
             -(CLon * Ux + SLon * Uy) * SLat + Uz * CLat;
      } else {
         VelocityZonalCell(ICell, K)      = Ux;
         VelocityMeridionalCell(ICell, K) = Uy;
      }
   }

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   bool OnSphere;
   Array1DI4 NEdgesReconOnCell;
   Array2DI4 ReconStencilCell;
   Array3DReal ReconWeightsCell;
   Array1DReal LatCell;
   Array1DReal LonCell;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

} // namespace OMEGA
#endif
