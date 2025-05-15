#include "AuxiliaryState.h"
#include "Config.h"
#include "Field.h"
#include "Logging.h"

namespace OMEGA {

// create the static class members
AuxiliaryState *AuxiliaryState::DefaultAuxState = nullptr;
std::map<std::string, std::unique_ptr<AuxiliaryState>>
    AuxiliaryState::AllAuxStates;

static std::string stripDefault(const std::string &Name) {
   return Name != "Default" ? Name : "";
}

// Constructor. Constructs the member auxiliary variables and registers their
// fields with IOStreams
AuxiliaryState::AuxiliaryState(const std::string &Name, const HorzMesh *Mesh,
                               Halo *MeshHalo, int NVertLevels, int NTracers)
    : Mesh(Mesh), MeshHalo(MeshHalo), Name(stripDefault(Name)),
      KineticAux(stripDefault(Name), Mesh, NVertLevels),
      LayerThicknessAux(stripDefault(Name), Mesh, NVertLevels),
      VorticityAux(stripDefault(Name), Mesh, NVertLevels),
      VelocityDel2Aux(stripDefault(Name), Mesh, NVertLevels),
      WindForcingAux(stripDefault(Name), Mesh, NVertLevels),
      TracerAux(stripDefault(Name), Mesh, NVertLevels, NTracers) {

   GroupName = "AuxiliaryState";
   if (Name != "Default") {
      GroupName.append(Name);
   }
   std::string AuxMeshName = Mesh->MeshName;

   auto AuxGroup = FieldGroup::create(GroupName);

   KineticAux.registerFields(GroupName, AuxMeshName);
   LayerThicknessAux.registerFields(GroupName, AuxMeshName);
   VorticityAux.registerFields(GroupName, AuxMeshName);
   VelocityDel2Aux.registerFields(GroupName, AuxMeshName);
   WindForcingAux.registerFields(GroupName, AuxMeshName);
   TracerAux.registerFields(GroupName, AuxMeshName);
}

// Destructor. Unregisters the fields with IOStreams and destroys this auxiliary
// state field group.
AuxiliaryState::~AuxiliaryState() {
   KineticAux.unregisterFields();
   LayerThicknessAux.unregisterFields();
   VorticityAux.unregisterFields();
   VelocityDel2Aux.unregisterFields();
   WindForcingAux.unregisterFields();
   TracerAux.unregisterFields();

   int Err = FieldGroup::destroy(GroupName);
   if (Err != 0)
      LOG_ERROR("Error destroying FieldGroup {}", GroupName);
}

// Compute the auxiliary variables needed for momentum equation
void AuxiliaryState::computeMomAux(const OceanState *State, int ThickTimeLevel,
                                   int VelTimeLevel) const {
   Array2DReal LayerThickCell;
   Array2DReal NormalVelEdge;
   State->getLayerThickness(LayerThickCell, ThickTimeLevel);
   State->getNormalVelocity(NormalVelEdge, VelTimeLevel);

   const int NVertLevels = LayerThickCell.extent_int(1);
   const int NChunks     = NVertLevels / VecLength;

   OMEGA_SCOPE(LocKineticAux, KineticAux);
   OMEGA_SCOPE(LocLayerThicknessAux, LayerThicknessAux);
   OMEGA_SCOPE(LocVorticityAux, VorticityAux);
   OMEGA_SCOPE(LocVelocityDel2Aux, VelocityDel2Aux);
   OMEGA_SCOPE(LocWindForcingAux, WindForcingAux);

   parallelFor(
       "vertexAuxState1", {Mesh->NVerticesAll, NChunks},
       KOKKOS_LAMBDA(int IVertex, int KChunk) {
          LocVorticityAux.computeVarsOnVertex(IVertex, KChunk, LayerThickCell,
                                              NormalVelEdge);
       });

   parallelFor(
       "cellAuxState1", {Mesh->NCellsAll, NChunks},
       KOKKOS_LAMBDA(int ICell, int KChunk) {
          LocKineticAux.computeVarsOnCell(ICell, KChunk, NormalVelEdge);
       });

   const auto &VelocityDivCell = KineticAux.VelocityDivCell;
   const auto &RelVortVertex   = VorticityAux.RelVortVertex;

   parallelFor(
       "edgeAuxState1", {Mesh->NEdgesAll, NChunks},
       KOKKOS_LAMBDA(int IEdge, int KChunk) {
          LocVorticityAux.computeVarsOnEdge(IEdge, KChunk);
          LocLayerThicknessAux.computeVarsOnEdge(IEdge, KChunk, LayerThickCell,
                                                 NormalVelEdge);
          LocVelocityDel2Aux.computeVarsOnEdge(IEdge, KChunk, VelocityDivCell,
                                               RelVortVertex);
          LocWindForcingAux.computeVarsOnEdge(IEdge, KChunk);
       });

   parallelFor(
       "vertexAuxState2", {Mesh->NVerticesAll, NChunks},
       KOKKOS_LAMBDA(int IVertex, int KChunk) {
          LocVelocityDel2Aux.computeVarsOnVertex(IVertex, KChunk);
       });

   parallelFor(
       "cellAuxState2", {Mesh->NCellsAll, NChunks},
       KOKKOS_LAMBDA(int ICell, int KChunk) {
          LocVelocityDel2Aux.computeVarsOnCell(ICell, KChunk);
       });

   parallelFor(
       "cellAuxState3", {Mesh->NCellsAll, NChunks},
       KOKKOS_LAMBDA(int ICell, int KChunk) {
          LocLayerThicknessAux.computeVarsOnCells(ICell, KChunk,
                                                  LayerThickCell);
       });
}

// Compute the auxiliary variables
void AuxiliaryState::computeAll(const OceanState *State,
                                const Array3DReal &TracerArray,
                                int ThickTimeLevel, int VelTimeLevel) const {
   Array2DReal LayerThickCell;
   Array2DReal NormalVelEdge;
   State->getLayerThickness(LayerThickCell, ThickTimeLevel);
   State->getNormalVelocity(NormalVelEdge, VelTimeLevel);

   const int NVertLevels = LayerThickCell.extent_int(1);
   const int NChunks     = NVertLevels / VecLength;
   const int NTracers    = TracerArray.extent_int(0);

   OMEGA_SCOPE(LocTracerAux, TracerAux);

   computeMomAux(State, ThickTimeLevel, VelTimeLevel);

   parallelFor(
       "edgeAuxState4", {NTracers, Mesh->NEdgesAll, NChunks},
       KOKKOS_LAMBDA(int LTracer, int IEdge, int KChunk) {
          LocTracerAux.computeVarsOnEdge(LTracer, IEdge, KChunk, NormalVelEdge,
                                         LayerThickCell, TracerArray);
       });

   const auto &MeanLayerThickEdge = LayerThicknessAux.MeanLayerThickEdge;

   parallelFor(
       "cellAuxState4", {NTracers, Mesh->NCellsAll, NChunks},
       KOKKOS_LAMBDA(int LTracer, int ICell, int KChunk) {
          LocTracerAux.computeVarsOnCells(LTracer, ICell, KChunk,
                                          MeanLayerThickEdge, TracerArray);
       });
}

void AuxiliaryState::computeAll(const OceanState *State,
                                const Array3DReal &TracerArray,
                                int TimeLevel) const {
   computeAll(State, TracerArray, TimeLevel, TimeLevel);
}

// Create a non-default auxiliary state
AuxiliaryState *AuxiliaryState::create(const std::string &Name,
                                       const HorzMesh *Mesh, Halo *MeshHalo,
                                       int NVertLevels, const int NTracers) {
   if (AllAuxStates.find(Name) != AllAuxStates.end()) {
      LOG_ERROR("Attempted to create a new AuxiliaryState with name {} but it "
                "already exists",
                Name);
      return nullptr;
   }

   auto *NewAuxState =
       new AuxiliaryState(Name, Mesh, MeshHalo, NVertLevels, NTracers);
   AllAuxStates.emplace(Name, NewAuxState);

   return NewAuxState;
}

// Create the default auxiliary state. Assumes that HorzMesh has been
// initialized.
int AuxiliaryState::init() {
   int Err                 = 0;
   const HorzMesh *DefMesh = HorzMesh::getDefault();
   Halo *DefHalo           = Halo::getDefault();

   int NVertLevels = DefMesh->NVertLevels;
   int NTracers    = Tracers::getNumTracers();

   AuxiliaryState::DefaultAuxState = AuxiliaryState::create(
       "Default", DefMesh, DefHalo, NVertLevels, NTracers);

   Config *OmegaConfig = Config::getOmegaConfig();
   Err                 = DefaultAuxState->readConfigOptions(OmegaConfig);

   return Err;
}

// Get the default auxiliary state
AuxiliaryState *AuxiliaryState::getDefault() {
   return AuxiliaryState::DefaultAuxState;
}

// Get auxiliary state by name
AuxiliaryState *AuxiliaryState::get(const std::string &Name) {
   // look for an instance of this name
   auto it = AllAuxStates.find(Name);

   // if found, return the pointer
   if (it != AllAuxStates.end()) {
      return it->second.get();

      // otherwise print error and return null pointer
   } else {
      LOG_ERROR("AuxiliaryState::get: Attempt to retrieve non-existent "
                "auxiliary state:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }
}

// Remove auxiliary state by name
void AuxiliaryState::erase(const std::string &Name) {
   AllAuxStates.erase(Name);
}

// Remove all auxiliary states
void AuxiliaryState::clear() { AllAuxStates.clear(); }

// Read and set config options
int AuxiliaryState::readConfigOptions(Config *OmegaConfig) {

   int Err = 0;

   Config AdvectConfig("Advection");
   Err = OmegaConfig->get(AdvectConfig);
   if (Err != 0) {
      LOG_CRITICAL("AuxiliaryState: Advection group not found in Config");
      return Err;
   }
   std::string FluxThickTypeStr;
   Err = AdvectConfig.get("FluxThicknessType", FluxThickTypeStr);
   if (Err != 0) {
      LOG_CRITICAL("AuxiliaryState: FluxThicknessType not found in "
                   "AdvectConfig");
      return Err;
   }

   if (FluxThickTypeStr == "Center") {
      this->LayerThicknessAux.FluxThickEdgeChoice = FluxThickEdgeOption::Center;
   } else if (FluxThickTypeStr == "Upwind") {
      this->LayerThicknessAux.FluxThickEdgeChoice = FluxThickEdgeOption::Upwind;
   } else {
      LOG_CRITICAL("AuxiliaryState: Unknown FluxThicknessType requested");
      Err = -1;
      return Err;
   }

   std::string FluxTracerTypeStr;
   Err = AdvectConfig.get("FluxTracerType", FluxTracerTypeStr);
   if (Err != 0) {
      LOG_CRITICAL("AuxiliaryState: FluxTracerType not found in "
                   "AdvectConfig");
      return Err;
   }

   if (FluxTracerTypeStr == "Center") {
      this->TracerAux.TracersOnEdgeChoice = FluxTracerEdgeOption::Center;
   } else if (FluxTracerTypeStr == "Upwind") {
      this->TracerAux.TracersOnEdgeChoice = FluxTracerEdgeOption::Upwind;
   } else {
      LOG_CRITICAL("AuxiliaryState: Unknown FluxTracerType requested");
      Err = -1;
      return Err;
   }

   Config WindStressConfig("WindStress");
   Err = OmegaConfig->get(WindStressConfig);

   std::string WindStressInterpTypeStr;
   Err = WindStressConfig.get("InterpType", WindStressInterpTypeStr);
   if (Err != 0) {
      LOG_CRITICAL("AuxiliaryState: InterpType not found in "
                   "WindStressConfig");
      return Err;
   }

   if (WindStressInterpTypeStr == "Isotropic") {
      this->WindForcingAux.InterpChoice = InterpCellToEdgeOption::Isotropic;
   } else if (WindStressInterpTypeStr == "Anisotropic") {
      this->WindForcingAux.InterpChoice = InterpCellToEdgeOption::Anisotropic;
   } else {
      LOG_CRITICAL("AuxiliaryState: Unknown InterpType requested");
      Err = -1;
      return Err;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Perform auxiliary state halo exchange
// Note that only non-computed auxiliary variables needs to be exchanged
I4 AuxiliaryState::exchangeHalo() {
   I4 Err = 0;

   Err +=
       MeshHalo->exchangeFullArrayHalo(WindForcingAux.ZonalStressCell, OnCell);
   Err +=
       MeshHalo->exchangeFullArrayHalo(WindForcingAux.MeridStressCell, OnCell);

   return Err;

} // end exchangeHalo

} // namespace OMEGA
