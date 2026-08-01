#ifndef OMEGA_PGRAD_H
#define OMEGA_PGRAD_H
//===-- ocn/PGrad.h - Pressure Gradient -----------------*- C++ -*-===//
///
/// Implements the PressureGrad class which provides a centered and a
/// finite-volume pressure gradient option and dispatches computations to
/// functor objects. This follows the patterns used in Eos.h/Eos.cpp.
//
//===----------------------------------------------------------------------===//

#include "Config.h"
#include "Eos.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"
#include <memory>

namespace OMEGA {

enum class PressureGradType {
   Centered,    ///< existing 2nd-order Montgomery scheme
   FiniteVolume ///< layer-integrated finite-volume scheme
   // , <FutureVariant>  ///< e.g. a 6th-order option, added when implemented
};

/// Mean-preserving vertical reconstruction of ConservTemp and AbsSalinity in
/// pressure used by the FiniteVolume scheme. The degree of the reconstruction
/// sets the scheme's exact set: linear deviations make the scheme exact for
/// profiles that vary linearly with pressure.
enum class PressureGradVertRecon {
   Linear ///< linear deviations (Phase 1)
   // , PPM  ///< parabolic (PPM-style) deviations (Phase 2)
};

/// Largest number of quadrature points supported within an edge layer
inline constexpr I4 MaxPGradQuadPoints = 4;

// Centered pressure gradient functor
class PressureGradCentered {
 public:
   bool Enabled;

   // constructor declaration
   PressureGradCentered(const HorzMesh *Mesh,   ///< [in] Horizontal mesh
                        const VertCoord *VCoord ///< [in] Vertical coordinate
   );

   // Compute centered pressure gradient contribution for given edge and
   // vertical chunk. This appends results into the Tend array (in-place).
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &PressureMid,
                                   const Array2DReal &PressureInterface,
                                   const Array2DReal &GeomZInterface,
                                   const Array1DReal &TidalPotential,
                                   const Array1DReal &SelfAttractionLoading,
                                   const Array2DReal &SpecVol) const {

      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));

      const I4 ICell0      = CellsOnEdge(IEdge, 0);
      const I4 ICell1      = CellsOnEdge(IEdge, 1);
      const Real InvDcEdge = 1.0_Real / DcEdge(IEdge);

      Real GradGeoPot =
          (TidalPotential(ICell1) - TidalPotential(ICell0)) * InvDcEdge +
          (SelfAttractionLoading(ICell1) - SelfAttractionLoading(ICell0)) *
              InvDcEdge;

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         Real MontPotCell0K =
             PressureInterface(ICell0, K) * SpecVol(ICell0, K) +
             Gravity * GeomZInterface(ICell0, K);
         Real MontPotCell1K =
             PressureInterface(ICell1, K) * SpecVol(ICell1, K) +
             Gravity * GeomZInterface(ICell1, K);
         Real GradMontPotK = (MontPotCell1K - MontPotCell0K) * InvDcEdge;

         Real MontPotCell0Kp1 =
             PressureInterface(ICell0, K + 1) * SpecVol(ICell0, K) +
             Gravity * GeomZInterface(ICell0, K + 1);
         Real MontPotCell1Kp1 =
             PressureInterface(ICell1, K + 1) * SpecVol(ICell1, K) +
             Gravity * GeomZInterface(ICell1, K + 1);
         Real GradMontPotKp1 = (MontPotCell1Kp1 - MontPotCell0Kp1) * InvDcEdge;
         Real GradMontPot    = 0.5_Real * (GradMontPotK + GradMontPotKp1);

         Real PGradAlpha =
             0.5_Real * (PressureMid(ICell1, K) + PressureMid(ICell0, K)) *
             (SpecVol(ICell1, K) - SpecVol(ICell0, K)) * InvDcEdge;
         Tend(IEdge, K) +=
             EdgeMask(IEdge, K) * (-GradMontPot + PGradAlpha - GradGeoPot);
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

// Finite-volume pressure gradient functor (placeholder)
class PressureGradFiniteVolume {
 public:
   bool Enabled;

   // Options cached from the PressureGrad config group by the PressureGrad
   // constructor, which is also where they are validated. Phase 1 implements
   // HorzOrder 2 and VertRecon Linear only.
   //
   // QuadraturePoints is a pure accuracy knob: the matched-pressure integrand
   // is zero pointwise for any profile the reconstruction resolves exactly, so
   // no quadrature rule can break the robustness property.
   I4 HorzOrder                    = 2;
   PressureGradVertRecon VertRecon = PressureGradVertRecon::Linear;
   I4 QuadraturePoints             = 2;

   // constructor declaration
   PressureGradFiniteVolume(
       const HorzMesh *Mesh,   ///< [in] Horizontal mesh
       const VertCoord *VCoord ///< [in] Vertical coordinate
   );

   KOKKOS_FUNCTION void operator()(
       const Array2DReal &Tend, I4 IEdge, I4 KChunk,
       const Array2DReal &PressureMid, const Array2DReal &PressureInterface,
       const Array2DReal &GeomZInterface, const Array1DReal &TidalPotential,
       const Array1DReal &SelfAttractionLoading, const Array2DReal &SpecVol,
       const Array2DReal &ConservTemp, const Array2DReal &AbsSalinity,
       const Array2DReal &SpecVolDCt, const Array2DReal &SpecVolDSa,
       const Array2DReal &SpecVolDP) const {

      // Placeholder: for now, no-op (future finite-volume implementation)
      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(IEdge, K) += 0.0_Real;
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

// Pressure gradient manager class
class PressureGrad {
 public:
   // Flag to indicate if pressure gradient term is enabled
   bool Enabled;

   // Initialize the default instance
   static void init();

   // Create a new pressure gradient object and add to map
   static PressureGrad *create(const std::string &Name, const HorzMesh *Mesh,
                               const VertCoord *VCoord, Config *Options);

   // Get the default instance
   static PressureGrad *getDefault();

   // Get instance by name
   static PressureGrad *
   get(const std::string &Name ///< [in] Name of PressureGrad
   );

   // Deallocates arrays and deletes instance
   static void clear();

   // Remove pressure gradient object by name
   static void erase(const std::string &Name ///< [in] Name of PressureGrad
   );

   // Destructor
   ~PressureGrad();

   // Accessors for the configured scheme and its options
   PressureGradType getType() const { return PressureGradChoice; }
   I4 getHorzOrder() const { return FiniteVolumePGrad.HorzOrder; }
   PressureGradVertRecon getVertRecon() const {
      return FiniteVolumePGrad.VertRecon;
   }
   I4 getQuadraturePoints() const { return FiniteVolumePGrad.QuadraturePoints; }

   // Compute pressure gradient tendencies and add into Tend array. The
   // FiniteVolume scheme additionally needs the layer-mean conservative
   // temperature and absolute salinity, which it reconstructs in pressure,
   // and the specific volume derivatives held by Eos. The Centered scheme
   // ignores them.
   void computePressureGrad(Array2DReal &Tend, const Array2DReal &PressureMid,
                            const Array2DReal &PressureInterface,
                            const Array2DReal &SpecVol,
                            const Array2DReal &GeomZInterface,
                            const Array2DReal &PseudoThick,
                            const Array2DReal &ConservTemp,
                            const Array2DReal &AbsSalinity,
                            const Eos *EqState) const;

 private:
   // Construct a new pressure gradient object
   PressureGrad(const HorzMesh *Mesh, const VertCoord *VCoord, Config *Options);

   // forbid copy and move construction
   PressureGrad(const PressureGrad &) = delete;
   PressureGrad(PressureGrad &&)      = delete;

   // Pointer to default pressure gradient object
   static PressureGrad *DefaultPGrad;

   // Mesh-related sizes
   I4 NEdgesAll     = 0;
   I4 NEdgesOwned   = 0;
   I4 NVertLayers   = 0;
   I4 NVertLayersP1 = 0;

   // Data required for computation (stored copies of VCoord arrays)
   Array1DI4 MinLayerEdgeBot; ///< min vertical layer on each edge
   Array1DI4 MaxLayerEdgeTop; ///< max vertical layer on each edge

   // Temporary: to be moveed to tidal forcing module in future
   Array1DReal TidalPotential; ///< Tidal potential for tidal forcing
   Array1DReal
       SelfAttractionLoading; ///< Self attraction and loading for tidal forcing

   // Instances of functors
   PressureGradCentered CenteredPGrad;
   PressureGradFiniteVolume FiniteVolumePGrad;

   // Choice from config
   PressureGradType PressureGradChoice = PressureGradType::Centered;

   // Map of all pressure gradient objects by name
   static std::map<std::string, std::unique_ptr<PressureGrad>> AllPGrads;

}; // end class PressureGrad

} // namespace OMEGA
#endif
