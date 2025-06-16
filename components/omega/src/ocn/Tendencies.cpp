//===-- ocn/Tendencies.cpp - Tendencies ------------------*- C++ -*-===//
//
// The Tendencies class is responsible for managing tendencies of state
// variables and tracers. It contains arrays that store the tendency data
// and provides methods for computing different tendency groups for use
// within the timestepping algorithm. At initialization, it determines which
// tendency terms are enabled.
//
//===----------------------------------------------------------------------===//

#include "Tendencies.h"
#include "CustomTendencyTerms.h"
#include "Error.h"
#include "Tracers.h"

namespace OMEGA {

Tendencies *Tendencies::DefaultTendencies = nullptr;
std::map<std::string, std::unique_ptr<Tendencies>> Tendencies::AllTendencies;

//------------------------------------------------------------------------------
// Initialize the tendencies. Assumes that HorzMesh as alread been initialized.
void Tendencies::init() {
   Error Err; // error code

   HorzMesh *DefHorzMesh = HorzMesh::getDefault();

   I4 NVertLevels = DefHorzMesh->NVertLevels;
   I4 NTracers    = Tracers::getNumTracers();

   // Get TendConfig group
   Config *OmegaConfig = Config::getOmegaConfig();
   Config TendConfig("Tendencies");
   Err += OmegaConfig->get(TendConfig);
   CHECK_ERROR_ABORT(Err, "Tendencies: Tendencies group not found in Config");

   // Check if use the customized tendencies. If it is not found in the
   // config, we assume it is not used (false)
   bool UseCustomTendency = false;
   Err += TendConfig.get("UseCustomTendency", UseCustomTendency);

   /// Instances of custom tendencies - empty by default
   CustomTendencyType CustomThickTend;
   CustomTendencyType CustomVelTend;

   if (UseCustomTendency) {
      // Check if use manufactured tendency terms if it is not found in
      // the config file, we will assume it is not used (false)
      bool ManufacturedTend = false;
      Error ManufacturedTendErr =
          TendConfig.get("ManufacturedSolutionTendency", ManufacturedTend);

      if (ManufacturedTend) {
         ManufacturedSolution ManufacturedSol;
         ManufacturedSol.init();

         CustomThickTend = ManufacturedSol.ManufacturedThickTend;
         CustomVelTend   = ManufacturedSol.ManufacturedVelTend;

      } // if ManufacturedTend

   } // end if UseCustomTendency

   // Ceate default tendencies
   Tendencies::DefaultTendencies =
       create("Default", DefHorzMesh, NVertLevels, NTracers, &TendConfig,
              CustomThickTend, CustomVelTend);

   DefaultTendencies->readTendConfig(&TendConfig);

} // end init

//------------------------------------------------------------------------------
// Destroys the tendencies
Tendencies::~Tendencies() {

   // No operations needed, Kokkos arrays removed when no longer in scope

} // end destructor

//------------------------------------------------------------------------------
// Removes all tendencies instances before exit
void Tendencies::clear() { AllTendencies.clear(); } // end clear

//------------------------------------------------------------------------------
// Removes tendencies from list by name
void Tendencies::erase(const std::string &Name) {

   AllTendencies.erase(Name);

} // end erase

//------------------------------------------------------------------------------
// Get default tendencies
Tendencies *Tendencies::getDefault() {

   return Tendencies::DefaultTendencies;

} // end get default

//------------------------------------------------------------------------------
// Get tendencies by name
Tendencies *Tendencies::get(const std::string &Name ///< [in] Name of tendencies
) {

   auto it = AllTendencies.find(Name);

   if (it != AllTendencies.end()) {
      return it->second.get();
   } else {
      LOG_ERROR(
          "Tendencies::get: Attempt to retrieve non-existent tendencies:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }

} // end get tendencies

//------------------------------------------------------------------------------
// read and set config options
void Tendencies::readTendConfig(
    Config *TendConfig ///< [in] Tendencies subconfig
) {
   Error Err; // error code

   Err += TendConfig->get("ThicknessFluxTendencyEnable",
                          this->ThicknessFluxDiv.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: ThicknessFluxTendencyEnable not found in TendConfig");

   Err += TendConfig->get("PVTendencyEnable", this->PotientialVortHAdv.Enabled);
   CHECK_ERROR_ABORT(Err,
                     "Tendencies: PVTendencyEnable not found in TendConfig");

   Err += TendConfig->get("KETendencyEnable", this->KEGrad.Enabled);
   CHECK_ERROR_ABORT(Err,
                     "Tendencies: KETendencyEnable not found in TendConfig");

   Err += TendConfig->get("SSHTendencyEnable", this->SSHGrad.Enabled);
   CHECK_ERROR_ABORT(Err,
                     "Tendencies: SSHTendencyEnable not found in TendConfig");

   Err += TendConfig->get("VelDiffTendencyEnable",
                          this->VelocityDiffusion.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: VelDiffTendencyEnable not found in TendConfig");

   Err += TendConfig->get("VelHyperDiffTendencyEnable",
                          this->VelocityHyperDiff.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: VelHyperDiffTendencyEnable not found in TendConfig");

   if (this->VelocityDiffusion.Enabled) {
      Err += TendConfig->get("ViscDel2", this->VelocityDiffusion.ViscDel2);
      CHECK_ERROR_ABORT(Err, "Tendencies: ViscDel2 not found in TendConfig");
   }

   if (this->VelocityHyperDiff.Enabled) {
      Err += TendConfig->get("ViscDel4", this->VelocityHyperDiff.ViscDel4);
      CHECK_ERROR_ABORT(Err, "Tendencies: ViscDel4 not found in TendConfig");
      Err += TendConfig->get("DivFactor", this->VelocityHyperDiff.DivFactor);
      CHECK_ERROR_ABORT(Err, "Tendencies: DivFactor not found in TendConfig");
   }

   Err += TendConfig->get("TracerHorzAdvTendencyEnable",
                          this->TracerHorzAdv.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: TracerHorzAdvTendencyEnable not found in TendConfig");

   Err += TendConfig->get("TracerDiffTendencyEnable",
                          this->TracerDiffusion.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: TracerDiffTendencyEnable not found in TendConfig");

   Err +=
       TendConfig->get("WindForcingTendencyEnable", this->WindForcing.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: WindForcingTendencyEnable not found in TendConfig");

   Err += TendConfig->get("Density0", this->WindForcing.SaltWaterDensity);
   CHECK_ERROR_ABORT(Err, "Tendencies: Density0 not found in TendConfig");

   Err += TendConfig->get("BottomDragTendencyEnable", this->BottomDrag.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: BottomDragTendencyEnable not found in TendConfig");

   Err += TendConfig->get("BottomDragCoeff", this->BottomDrag.Coeff);
   CHECK_ERROR_ABORT(Err,
                     "Tendencies: BottomDragCoeff not found in TendConfig");

   Err += TendConfig->get("TracerHorzAdvTendencyEnable",
                          this->TracerHorzAdv.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: TracerHorzAdvTendencyEnable not found in TendConfig");

   if (this->TracerDiffusion.Enabled) {
      Err += TendConfig->get("EddyDiff2", this->TracerDiffusion.EddyDiff2);
      CHECK_ERROR_ABORT(Err, "Tendencies: EddyDiff2 not found in TendConfig");
   }

   Err += TendConfig->get("TracerHyperDiffTendencyEnable",
                          this->TracerHyperDiff.Enabled);
   CHECK_ERROR_ABORT(
       Err,
       "Tendencies: TracerHyperDiffTendencyEnable not found in TendConfig");

   if (this->TracerHyperDiff.Enabled) {
      Err += TendConfig->get("EddyDiff4", this->TracerHyperDiff.EddyDiff4);
      CHECK_ERROR_ABORT(Err, "Tendencies: EddyDiff4 not found in TendConfig");
   }
}

//------------------------------------------------------------------------------
// Construct a new group of tendencies
Tendencies::Tendencies(const std::string &Name, ///< [in] Name for tendencies
                       const HorzMesh *Mesh,    ///< [in] Horizontal mesh
                       int NVertLevels, ///< [in] Number of vertical levels
                       int NTracersIn,  ///< [in] Number of tracers
                       Config *Options, ///< [in] Configuration options
                       CustomTendencyType InCustomThicknessTend,
                       CustomTendencyType InCustomVelocityTend)
    : ThicknessFluxDiv(Mesh), PotientialVortHAdv(Mesh), KEGrad(Mesh),
      SSHGrad(Mesh), VelocityDiffusion(Mesh), VelocityHyperDiff(Mesh),
      WindForcing(Mesh), BottomDrag(Mesh), TracerHorzAdv(Mesh),
      TracerDiffusion(Mesh), TracerHyperDiff(Mesh),
      CustomThicknessTend(InCustomThicknessTend),
      CustomVelocityTend(InCustomVelocityTend) {

   // Tendency arrays
   LayerThicknessTend =
       Array2DReal("LayerThicknessTend", Mesh->NCellsSize, NVertLevels);
   NormalVelocityTend =
       Array2DReal("NormalVelocityTend", Mesh->NEdgesSize, NVertLevels);
   TracerTend =
       Array3DReal("TracerTend", NTracersIn, Mesh->NCellsSize, NVertLevels);

   // Array dimension lengths
   NCellsAll = Mesh->NCellsAll;
   NEdgesAll = Mesh->NEdgesAll;
   NTracers  = NTracersIn;
   NChunks   = NVertLevels / VecLength;

} // end constructor

Tendencies::Tendencies(const std::string &Name, ///< [in] Name for tendencies
                       const HorzMesh *Mesh,    ///< [in] Horizontal mesh
                       int NVertLevels, ///< [in] Number of vertical levels
                       int NTracersIn,  ///< [in] Number of tracers
                       Config *Options) ///< [in] Configuration options
    : Tendencies(Name, Mesh, NVertLevels, NTracersIn, Options,
                 CustomTendencyType{}, CustomTendencyType{}) {}

//------------------------------------------------------------------------------
// Compute tendencies for layer thickness equation
void Tendencies::computeThicknessTendenciesOnly(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {

   OMEGA_SCOPE(LocLayerThicknessTend, LayerThicknessTend);
   OMEGA_SCOPE(LocThicknessFluxDiv, ThicknessFluxDiv);
   Array2DReal NormalVelEdge;
   State->getNormalVelocity(NormalVelEdge, VelTimeLevel);

   deepCopy(LocLayerThicknessTend, 0);

   // Compute thickness flux divergence
   const Array2DReal &ThickFluxEdge =
       AuxState->LayerThicknessAux.FluxLayerThickEdge;

   if (LocThicknessFluxDiv.Enabled) {
      parallelFor(
          {NCellsAll, NChunks}, KOKKOS_LAMBDA(int ICell, int KChunk) {
             LocThicknessFluxDiv(LocLayerThicknessTend, ICell, KChunk,
                                 ThickFluxEdge, NormalVelEdge);
          });
   }

   if (CustomThicknessTend) {
      CustomThicknessTend(LocLayerThicknessTend, State, AuxState,
                          ThickTimeLevel, VelTimeLevel, Time);
   }

} // end thickness tendency compute

//------------------------------------------------------------------------------
// Compute tendencies for normal velocity equation
void Tendencies::computeVelocityTendenciesOnly(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {

   OMEGA_SCOPE(LocNormalVelocityTend, NormalVelocityTend);
   OMEGA_SCOPE(LocPotientialVortHAdv, PotientialVortHAdv);
   OMEGA_SCOPE(LocKEGrad, KEGrad);
   OMEGA_SCOPE(LocSSHGrad, SSHGrad);
   OMEGA_SCOPE(LocVelocityDiffusion, VelocityDiffusion);
   OMEGA_SCOPE(LocVelocityHyperDiff, VelocityHyperDiff);
   OMEGA_SCOPE(LocWindForcing, WindForcing);
   OMEGA_SCOPE(LocBottomDrag, BottomDrag);

   deepCopy(LocNormalVelocityTend, 0);

   const Array2DReal &NormalVelEdge = State->NormalVelocity[VelTimeLevel];

   // Compute potential vorticity horizontal advection
   const Array2DReal &FluxLayerThickEdge =
       AuxState->LayerThicknessAux.FluxLayerThickEdge;
   const Array2DReal &NormRVortEdge = AuxState->VorticityAux.NormRelVortEdge;
   const Array2DReal &NormFEdge     = AuxState->VorticityAux.NormPlanetVortEdge;
   Array2DReal NormVelEdge;
   State->getNormalVelocity(NormVelEdge, VelTimeLevel);
   if (LocPotientialVortHAdv.Enabled) {
      parallelFor(
          {NEdgesAll, NChunks}, KOKKOS_LAMBDA(int IEdge, int KChunk) {
             LocPotientialVortHAdv(LocNormalVelocityTend, IEdge, KChunk,
                                   NormRVortEdge, NormFEdge, FluxLayerThickEdge,
                                   NormVelEdge);
          });
   }

   // Compute kinetic energy gradient
   const Array2DReal &KECell = AuxState->KineticAux.KineticEnergyCell;
   if (LocKEGrad.Enabled) {
      parallelFor(
          {NEdgesAll, NChunks}, KOKKOS_LAMBDA(int IEdge, int KChunk) {
             LocKEGrad(LocNormalVelocityTend, IEdge, KChunk, KECell);
          });
   }

   // Compute sea surface height gradient
   const Array2DReal &SSHCell = AuxState->LayerThicknessAux.SshCell;
   if (LocSSHGrad.Enabled) {
      parallelFor(
          {NEdgesAll, NChunks}, KOKKOS_LAMBDA(int IEdge, int KChunk) {
             LocSSHGrad(LocNormalVelocityTend, IEdge, KChunk, SSHCell);
          });
   }

   // Compute del2 horizontal diffusion
   const Array2DReal &DivCell     = AuxState->KineticAux.VelocityDivCell;
   const Array2DReal &RVortVertex = AuxState->VorticityAux.RelVortVertex;
   if (LocVelocityDiffusion.Enabled) {
      parallelFor(
          {NEdgesAll, NChunks}, KOKKOS_LAMBDA(int IEdge, int KChunk) {
             LocVelocityDiffusion(LocNormalVelocityTend, IEdge, KChunk, DivCell,
                                  RVortVertex);
          });
   }

   // Compute del4 horizontal diffusion
   const Array2DReal &Del2DivCell = AuxState->VelocityDel2Aux.Del2DivCell;
   const Array2DReal &Del2RVortVertex =
       AuxState->VelocityDel2Aux.Del2RelVortVertex;
   if (LocVelocityHyperDiff.Enabled) {
      parallelFor(
          {NEdgesAll, NChunks}, KOKKOS_LAMBDA(int IEdge, int KChunk) {
             LocVelocityHyperDiff(LocNormalVelocityTend, IEdge, KChunk,
                                  Del2DivCell, Del2RVortVertex);
          });
   }

   // Compute wind forcing
   const auto &NormalStressEdge = AuxState->WindForcingAux.NormalStressEdge;
   const auto &MeanLayerThickEdge =
       AuxState->LayerThicknessAux.MeanLayerThickEdge;
   if (LocWindForcing.Enabled) {
      parallelFor(
          {NEdgesAll, NChunks}, KOKKOS_LAMBDA(int IEdge, int KChunk) {
             LocWindForcing(LocNormalVelocityTend, IEdge, KChunk,
                            NormalStressEdge, MeanLayerThickEdge);
          });
   }

   // Compute bottom drag
   if (LocBottomDrag.Enabled) {
      parallelFor(
          {NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
             LocBottomDrag(LocNormalVelocityTend, IEdge, NormalVelEdge, KECell,
                           MeanLayerThickEdge);
          });
   }

   if (CustomVelocityTend) {
      CustomVelocityTend(LocNormalVelocityTend, State, AuxState, ThickTimeLevel,
                         VelTimeLevel, Time);
   }

} // end velocity tendency compute

void Tendencies::computeTracerTendenciesOnly(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {
   OMEGA_SCOPE(LocTracerTend, TracerTend);
   OMEGA_SCOPE(LocTracerHorzAdv, TracerHorzAdv);
   OMEGA_SCOPE(LocTracerDiffusion, TracerDiffusion);
   OMEGA_SCOPE(LocTracerHyperDiff, TracerHyperDiff);

   deepCopy(LocTracerTend, 0);

   // compute tracer horizotal advection
   const Array2DReal &NormalVelEdge = State->NormalVelocity[VelTimeLevel];
   const Array3DReal &HTracersEdge  = AuxState->TracerAux.HTracersEdge;
   if (LocTracerHorzAdv.Enabled) {
      parallelFor(
          {NTracers, NCellsAll, NChunks},
          KOKKOS_LAMBDA(int L, int ICell, int KChunk) {
             LocTracerHorzAdv(LocTracerTend, L, ICell, KChunk, NormalVelEdge,
                              HTracersEdge);
          });
   }

   // compute tracer diffusion
   const Array2DReal &MeanLayerThickEdge =
       AuxState->LayerThicknessAux.MeanLayerThickEdge;
   if (LocTracerDiffusion.Enabled) {
      parallelFor(
          {NTracers, NCellsAll, NChunks},
          KOKKOS_LAMBDA(int L, int ICell, int KChunk) {
             LocTracerDiffusion(LocTracerTend, L, ICell, KChunk, TracerArray,
                                MeanLayerThickEdge);
          });
   }

   // compute tracer hyperdiffusion
   const Array3DReal &Del2TracersCell = AuxState->TracerAux.Del2TracersCell;
   if (LocTracerHyperDiff.Enabled) {
      parallelFor(
          {NTracers, NCellsAll, NChunks},
          KOKKOS_LAMBDA(int L, int ICell, int KChunk) {
             LocTracerHyperDiff(LocTracerTend, L, ICell, KChunk,
                                Del2TracersCell);
          });
   }

} // end tracer tendency compute

void Tendencies::computeThicknessTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {
   // only need LayerThicknessAux on edge
   Array2DReal LayerThick;
   Array2DReal NormVel;
   State->getLayerThickness(LayerThick, ThickTimeLevel);
   State->getNormalVelocity(NormVel, VelTimeLevel);
   OMEGA_SCOPE(LayerThicknessAux, AuxState->LayerThicknessAux);
   OMEGA_SCOPE(LayerThickCell, LayerThick);
   OMEGA_SCOPE(NormalVelEdge, NormVel);

   parallelFor(
       "computeLayerThickAux", {NEdgesAll, NChunks},
       KOKKOS_LAMBDA(int IEdge, int KChunk) {
          LayerThicknessAux.computeVarsOnEdge(IEdge, KChunk, LayerThickCell,
                                              NormalVelEdge);
       });

   computeThicknessTendenciesOnly(State, AuxState, ThickTimeLevel, VelTimeLevel,
                                  Time);
}

void Tendencies::computeVelocityTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {
   AuxState->computeMomAux(State, ThickTimeLevel, VelTimeLevel);
   computeVelocityTendenciesOnly(State, AuxState, ThickTimeLevel, VelTimeLevel,
                                 Time);
}

void Tendencies::computeTracerTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {
   OMEGA_SCOPE(TracerAux, AuxState->TracerAux);
   OMEGA_SCOPE(LayerThickCell, State->LayerThickness[ThickTimeLevel]);
   OMEGA_SCOPE(NormalVelEdge, State->NormalVelocity[VelTimeLevel]);

   parallelFor(
       "computeTracerAuxEdge", {NTracers, NEdgesAll, NChunks},
       KOKKOS_LAMBDA(int LTracer, int IEdge, int KChunk) {
          TracerAux.computeVarsOnEdge(LTracer, IEdge, KChunk, NormalVelEdge,
                                      LayerThickCell, TracerArray);
       });

   const auto &MeanLayerThickEdge =
       AuxState->LayerThicknessAux.MeanLayerThickEdge;
   parallelFor(
       "computeTracerAuxCell", {NTracers, NCellsAll, NChunks},
       KOKKOS_LAMBDA(int LTracer, int ICell, int KChunk) {
          TracerAux.computeVarsOnCells(LTracer, ICell, KChunk,
                                       MeanLayerThickEdge, TracerArray);
       });

   computeTracerTendenciesOnly(State, AuxState, TracerArray, ThickTimeLevel,
                               VelTimeLevel, Time);
}

//------------------------------------------------------------------------------
// Compute both layer thickness and normal velocity tendencies
void Tendencies::computeAllTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {

   AuxState->computeAll(State, TracerArray, ThickTimeLevel, VelTimeLevel);
   computeThicknessTendenciesOnly(State, AuxState, ThickTimeLevel, VelTimeLevel,
                                  Time);
   computeVelocityTendenciesOnly(State, AuxState, ThickTimeLevel, VelTimeLevel,
                                 Time);
   computeTracerTendenciesOnly(State, AuxState, TracerArray, ThickTimeLevel,
                               VelTimeLevel, Time);

} // end all tendency compute

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
