//===-- ocn/PGrad.cpp - Pressure Gradient Term -----------------*- C++ -*-===//
//
// Implements the PGrad manager and two discretizations: Centered and
// FiniteVolume.
//
//===----------------------------------------------------------------------===//

#include "PGrad.h"
#include "Eos.h"
#include "Error.h"
#include "Field.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

namespace OMEGA {

PressureGrad *PressureGrad::DefaultPGrad = nullptr;
std::map<std::string, std::unique_ptr<PressureGrad>> PressureGrad::AllPGrads;

//------------------------------------------------------------------------------
// Initialize the PressureGrad. Assumes that HorzMesh and VertCoord have already
// been initialized.
void PressureGrad::init() {

   // Retrieve default mesh and vertical coordinate
   HorzMesh *DefMesh    = HorzMesh::getDefault();
   VertCoord *DefVCoord = VertCoord::getDefault();

   // Retrieve omega config
   Config *OmegaConfig = Config::getOmegaConfig();

   // Create the default PressureGrad and set pointer to it
   PressureGrad::DefaultPGrad =
       PressureGrad::create("Default", DefMesh, DefVCoord, OmegaConfig);

} // end init

//------------------------------------------------------------------------------
// Create a new PressureGrad object and add it to the map
PressureGrad *
PressureGrad::create(const std::string &Name, /// [in] Name for PressureGrad
                     const HorzMesh *Mesh,    ///< [in] Horizontal mesh
                     const VertCoord *VCoord, ///< [in] Vertical coordinate
                     Config *Options) {       ///< [in] Configuration options

   // Check to see if a PressureGrad of the same name already exists and
   // if so, exit with an error
   if (AllPGrads.find(Name) != AllPGrads.end()) {
      LOG_ERROR("Attempted to create a PressureGrad with name {} but a "
                "PressureGrad of that name already exists",
                Name);
      return nullptr;
   }

   // create a new PressureGrad on the heap and put it in a map of
   // unique_ptrs, which will manage its lifetime
   auto *NewPGrad = new PressureGrad(Mesh, VCoord, Options);
   AllPGrads.emplace(Name, NewPGrad);

   return NewPGrad;

} // end create

//------------------------------------------------------------------------------
// Get the default pressure gradient instance
PressureGrad *PressureGrad::getDefault() {

   return DefaultPGrad;

} // end get default

//------------------------------------------------------------------------------
// Constructor for PressureGrad
PressureGrad::PressureGrad(
    const HorzMesh *Mesh,    ///< [in] Horizontal mesh
    const VertCoord *VCoord, ///< [in] Vertical coordinate
    Config *Options)         ///< [in] Configuration options
    : MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop), CenteredPGrad(Mesh, VCoord),
      FiniteVolumePGrad(Mesh, VCoord) {

   // store mesh sizes
   NEdgesAll     = Mesh->NEdgesAll;
   NEdgesOwned   = Mesh->NEdgesOwned;
   NVertLayers   = VCoord->NVertLayers;
   NVertLayersP1 = NVertLayers + 1;

   // Read config options for PressureGrad type
   // and enable the appropriate functor
   Config PGradConfig("PressureGrad");
   Error Err;
   Err += Options->get(PGradConfig);
   CHECK_ERROR_ABORT(Err,
                     "PressureGrad: PressureGrad group not found in Config");
   std::string PGradTypeStr;
   Err += PGradConfig.get("PressureGradType", PGradTypeStr);

   if (PGradTypeStr == "centered" || PGradTypeStr == "Centered") {
      PressureGradChoice          = PressureGradType::Centered;
      this->CenteredPGrad.Enabled = true;
   } else if (PGradTypeStr == "finiteVolume" ||
              PGradTypeStr == "FiniteVolume") {
      PressureGradChoice              = PressureGradType::FiniteVolume;
      this->FiniteVolumePGrad.Enabled = true;
   } else {
      LOG_INFO(
          "PGrad: Unknown PressureGradType in config, defaulting to centered");
   }

   // Read the FiniteVolume sub-options. All three are optional so that a
   // configuration written before they existed still parses; where a key is
   // absent the Phase 1 default is used. Phase 2 will add values to these
   // keys rather than new keys, so a Phase 1 configuration will continue to
   // work unchanged.
   I4 HorzOrder = FiniteVolumePGrad.HorzOrder;
   if (PGradConfig.existsVar("HorzOrder"))
      Err += PGradConfig.get("HorzOrder", HorzOrder);

   std::string VertReconStr = "linear";
   if (PGradConfig.existsVar("VerticalReconstruction"))
      Err += PGradConfig.get("VerticalReconstruction", VertReconStr);

   I4 QuadraturePoints = FiniteVolumePGrad.QuadraturePoints;
   if (PGradConfig.existsVar("QuadraturePoints"))
      Err += PGradConfig.get("QuadraturePoints", QuadraturePoints);

   CHECK_ERROR_ABORT(Err, "PressureGrad: error reading PressureGrad options");

   // Validate. Values reserved for Phase 2 are rejected outright rather than
   // silently falling back to the Phase 1 setting, which would make a Phase 2
   // run look like a Phase 1 pass.
   if (HorzOrder != 2)
      ABORT_ERROR("PressureGrad: HorzOrder {} is not implemented in Phase 1 "
                  "of the FiniteVolume pressure gradient; only HorzOrder 2 "
                  "(the two-cell stencil) is available",
                  HorzOrder);

   if (VertReconStr == "linear" || VertReconStr == "Linear") {
      FiniteVolumePGrad.VertRecon = PressureGradVertRecon::Linear;
   } else if (VertReconStr == "ppm" || VertReconStr == "PPM") {
      ABORT_ERROR("PressureGrad: VerticalReconstruction '{}' is not "
                  "implemented in Phase 1 of the FiniteVolume pressure "
                  "gradient; only 'linear' is available",
                  VertReconStr);
   } else {
      ABORT_ERROR("PressureGrad: unknown VerticalReconstruction '{}'; valid "
                  "values are 'linear' (Phase 1) and 'ppm' (Phase 2)",
                  VertReconStr);
   }

   if (QuadraturePoints < 1 || QuadraturePoints > MaxPGradQuadPoints)
      ABORT_ERROR("PressureGrad: QuadraturePoints {} is out of range; must be "
                  "between 1 and {}",
                  QuadraturePoints, MaxPGradQuadPoints);

   FiniteVolumePGrad.HorzOrder        = HorzOrder;
   FiniteVolumePGrad.QuadraturePoints = QuadraturePoints;

   // Temporary: initialization of tidal potential and SAL
   TidalPotential = Array1DReal("TidalPotential", Mesh->NCellsSize);
   SelfAttractionLoading =
       Array1DReal("SelfAttractionLoading", Mesh->NCellsSize);

} // end constructor

//------------------------------------------------------------------------------
// Destructor for PressureGrad
PressureGrad::~PressureGrad() {

   // No operations needed, Kokkos arrays removed when no longer in scope

} // end destructor

//------------------------------------------------------------------------------
// Remove PressureGrad instances before exit
void PressureGrad::clear() {

   AllPGrads.clear();

   DefaultPGrad = nullptr; // prevent dangling pointe'r

} // end clear

//------------------------------------------------------------------------------
// Remove PressureGrad from list by name
void PressureGrad::erase(const std::string &Name) {

   AllPGrads.erase(Name);

} // end erase

//------------------------------------------------------------------------------
// Get pressure gradient instance by name
PressureGrad *PressureGrad::get(const std::string &Name ///< [in] Name of
) {

   auto it = AllPGrads.find(Name);

   if (it != AllPGrads.end()) {
      return it->second.get();
   } else {
      LOG_ERROR("PressureGrad::get: Attempt to retrieve non-existent "
                "PressureGrad:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }

} // end get pressure gradient

//------------------------------------------------------------------------------
// Compute pressure gradient tendencies and add into Tend array
void PressureGrad::computePressureGrad(
    Array2DReal &Tend, const Array2DReal &PressureMid,
    const Array2DReal &PressureInterface, const Array2DReal &SpecVol,
    const Array2DReal &GeomZInterface, const Array2DReal &PseudoThick,
    const Array2DReal &ConservTemp, const Array2DReal &AbsSalinity,
    const Eos *EqState) const {

   OMEGA_SCOPE(LocCenteredPGrad, CenteredPGrad);
   OMEGA_SCOPE(LocFiniteVolumePGrad, FiniteVolumePGrad);
   OMEGA_SCOPE(LocMinLayerEdgeBot, MinLayerEdgeBot);
   OMEGA_SCOPE(LocMaxLayerEdgeTop, MaxLayerEdgeTop);
   OMEGA_SCOPE(LocTidalPotential, TidalPotential);
   OMEGA_SCOPE(LocSelfAttractionLoading, SelfAttractionLoading);

   if (PressureGradChoice == PressureGradType::Centered) {

      // computes centered geopotential and pressure gradient tendency
      parallelForOuter(
          "pgrad-centered", {NEdgesAll},
          KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
             const int KMin   = LocMinLayerEdgeBot(IEdge);
             const int KMax   = LocMaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocCenteredPGrad(Tend, IEdge, KChunk, PressureMid,
                                     PressureInterface, GeomZInterface,
                                     LocTidalPotential,
                                     LocSelfAttractionLoading, SpecVol);
                 });
          });

   } else {

      // The specific volume derivatives come from Eos, which computes them
      // from the same single equation-of-state evaluation as SpecVol
      const Array2DReal SpecVolDCt = EqState->SpecVolDCt;
      const Array2DReal SpecVolDSa = EqState->SpecVolDSa;
      const Array2DReal SpecVolDP  = EqState->SpecVolDP;

      // computes finite-volume geopotential and pressure gradient tendency
      parallelForOuter(
          "pgrad-finitevolume", {NEdgesAll},
          KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
             const int KMin   = LocMinLayerEdgeBot(IEdge);
             const int KMax   = LocMaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocFiniteVolumePGrad(
                        Tend, IEdge, KChunk, PressureMid, PressureInterface,
                        GeomZInterface, LocTidalPotential,
                        LocSelfAttractionLoading, SpecVol, ConservTemp,
                        AbsSalinity, SpecVolDCt, SpecVolDSa, SpecVolDP);
                 });
          });
   }
} // end compute pressure gradient

//------------------------------------------------------------------------------
// Constructor for centered pressure gradient functor
PressureGradCentered::PressureGradCentered(
    const HorzMesh *Mesh,   ///< [in] Horizontal mesh
    const VertCoord *VCoord ///< [in] Vertical coordinate
    )
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

//------------------------------------------------------------------------------
// Constructor for finite volume pressure gradient functor
PressureGradFiniteVolume::PressureGradFiniteVolume(
    const HorzMesh *Mesh,   ///< [in] Horizontal mesh
    const VertCoord *VCoord ///< [in] Vertical coordinate
    )
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

} // namespace OMEGA
