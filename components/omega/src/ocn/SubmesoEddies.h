#ifndef OMEGA_SUBMESOEDDIES_H
#define OMEGA_SUBMESOEDDIES_H

#include "HorzMesh.h"
#include "VertCoord.h"

namespace OMEGA {

// Generic linear interpolation routine. This should be in OmegaMath.h or
// something like that once it exists.
KOKKOS_INLINE_FUNCTION Real linearInterp(Real x, Real y1, Real x1, Real y2,
                                         Real x2) {
   const Real A = (y1 - y2) / (x1 - x2);
   const Real B = y1 - A * x1;
   return (x1 == x2) ? y1 : A * x + B;
}

// A class for the submesoscale eddy parametrization.
// It groups the variables related to this parametrization
// and provides methods to compute mixed layer depth, buoyancy gradient,
// and eddy velocity.
class SubmesoEddies {
 public:
   // Public methods

   // Initialize SubmesoEddies from config, default mesh, and default vertical
   // coordinate
   static void init();

   // Get instance of SubmesoEddies
   static SubmesoEddies *getInstance();

   // Destroy instance (frees Kokkos views)
   static void destroyInstance();

   // Compute time scale array from time scale constant and Coriolis parameter
   void computeTimeScale();

   // Compute mixed layer index and  depth based on the density difference
   // criterion
   void computeDenMixLayerDepth(const Array2DReal &SpecVol);

   // Compute buoyancy gradient
   void computeBuoyGrad(const Array2DReal &SpecVol,
                        const Array2DReal &MeanPseudoThickEdge,
                        const Array2DReal &GeomZMid,
                        const Array2DReal &BruntVaisalaFreqSq);

   // Compute eddy velocity
   void computeEddyVelocity(const Array2DReal &BruntVaisalaFreqSq,
                            const Array2DReal &MeanPseudoThickEdge);

   // Public member variables

   // Enable the parametrization
   bool Enable;

   // Parametrization constants

   // Minimum width of submesoscale fronts
   Real LfMin;

   // Efficiency coefficient
   Real Ce;

   // Maximum edge length
   Real DsMax;

   // Time scale constant
   Real Tau;

   // Reference depth
   Real ReferenceDepth = 10;

   // Density threshold for determining the mixed layer depth
   Real DenThreshold = 0.03;

   // Mixed layer index
   Array1DI4 DenMixLayerIndex;

   // Mixed layer depth
   Array1DReal DenMixLayerDepth;

   // Buoyancy gradient
   Array2DReal GradBuoyEdgeInterface;

   // Eddy Velocity
   Array2DReal EddyVelocity;

   // Time scale array
   Array1DReal TimeScale;

 private:
   // Private methods

   // Constructor
   SubmesoEddies(const HorzMesh *Mesh, const VertCoord *VCoord);

   // Destructor
   ~SubmesoEddies() = default;

   // Delete copy and move constructors and assignment operators
   SubmesoEddies(const SubmesoEddies &)            = delete;
   SubmesoEddies &operator=(const SubmesoEddies &) = delete;
   SubmesoEddies(SubmesoEddies &&)                 = delete;
   SubmesoEddies &operator=(SubmesoEddies &&)      = delete;

   // Define fields and metadata
   void defineFields();

   // Private member variables

   // Instance pointer
   static SubmesoEddies *Instance;

   // Pointer to horizontal mesh
   const HorzMesh *Mesh;

   // Pointer to vertical coordinate
   const VertCoord *VCoord;
};

} // namespace OMEGA
#endif
