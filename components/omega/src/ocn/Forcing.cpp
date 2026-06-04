#include "Forcing.h"
#include "Eos.h"
#include "Field.h"
#include "Logging.h"
#include "Pacer.h"
#include "Tracers.h"

namespace OMEGA {

Forcing *Forcing::DefaultForcing = nullptr;
std::map<std::string, std::unique_ptr<Forcing>> Forcing::AllForcing;

static std::string stripDefault(const std::string &Name) {
   return Name != "Default" ? Name : "";
}

Forcing::Forcing(const std::string &Name, const HorzMesh *Mesh, Halo *MeshHalo,
                 const VertCoord *VCoord, int NTracers)
    : Name(stripDefault(Name)), TracerForcingAux(stripDefault(Name), Mesh),
      SurfTracerRestAux(stripDefault(Name), Mesh, NTracers),
      MomForcingAux(stripDefault(Name), Mesh), Mesh(Mesh), MeshHalo(MeshHalo),
      VCoord(VCoord) {}

Forcing::~Forcing() { unregisterFields(); }

void Forcing::registerFields(const std::string &MeshName) const {
   MomForcingAux.registerFields(MeshName);
   TracerForcingAux.registerFields(MeshName);
   SurfTracerRestAux.registerFields(MeshName);
}

void Forcing::unregisterFields() const {
   MomForcingAux.unregisterFields();
   TracerForcingAux.unregisterFields();
   SurfTracerRestAux.unregisterFields();
}

Forcing *Forcing::create(const std::string &Name, const HorzMesh *Mesh,
                         Halo *MeshHalo, const VertCoord *VCoord,
                         int NTracers) {
   if (AllForcing.find(Name) != AllForcing.end()) {
      LOG_ERROR("Attempted to create new Forcing with name {} but it "
                "already exists",
                Name);
      return nullptr;
   }

   auto *NewForcing = new Forcing(Name, Mesh, MeshHalo, VCoord, NTracers);
   AllForcing.emplace(Name, NewForcing);

   return NewForcing;
}

void Forcing::init() {
   if (DefaultForcing != nullptr) {
      return;
   }

   FieldGroup::create("Forcing");

   const HorzMesh *DefMesh    = HorzMesh::getDefault();
   Halo *DefHalo              = Halo::getDefault();
   const VertCoord *DefVCoord = VertCoord::getDefault();
   int NTracers               = Tracers::getNumTracers();

   DefaultForcing =
       Forcing::create("Default", DefMesh, DefHalo, DefVCoord, NTracers);

   if (DefaultForcing == nullptr) {
      ABORT_ERROR("Forcing: failed to initialize default forcing state");
   }

   DefaultForcing->registerFields(DefMesh->MeshName);

   Config *OmegaConfig = Config::getOmegaConfig();
   DefaultForcing->readConfigOptions(OmegaConfig);
}

Forcing *Forcing::getDefault() { return DefaultForcing; }

Forcing *Forcing::get(const std::string &Name) {
   auto it = AllForcing.find(Name);
   if (it != AllForcing.end()) {
      return it->second.get();
   }

   LOG_ERROR("Forcing::get: Attempt to retrieve non-existent forcing state:");
   LOG_ERROR("{} has not been defined or has been removed", Name);
   return nullptr;
}

bool Forcing::exists(const std::string &Name) {
   return AllForcing.find(Name) != AllForcing.end();
}

void Forcing::erase(const std::string &Name) { AllForcing.erase(Name); }

void Forcing::clear() {
   AllForcing.clear();
   DefaultForcing = nullptr;
   if (FieldGroup::exists("Forcing"))
      FieldGroup::destroy("Forcing");
}

void Forcing::readConfigOptions(Config *OmegaConfig) {
   Error Err;

   Config SrfStressConfig("SrfStress");
   Err += OmegaConfig->get(SrfStressConfig);

   std::string SrfStressInterpTypeStr;
   Err += SrfStressConfig.get("InterpType", SrfStressInterpTypeStr);
   CHECK_ERROR_ABORT(Err, "Forcing: InterpType not found in SrfStressConfig");

   if (SrfStressInterpTypeStr == "Isotropic") {
      this->MomForcingAux.InterpChoice = InterpCellToEdgeOption::Isotropic;
   } else if (SrfStressInterpTypeStr == "Anisotropic") {
      this->MomForcingAux.InterpChoice = InterpCellToEdgeOption::Anisotropic;
   } else {
      ABORT_ERROR("Forcing: Unknown InterpType requested");
   }
}

void Forcing::computeSrfStressForcingOnEdge() const {
   OMEGA_SCOPE(LocMomForcingAux, MomForcingAux);

   Pacer::start("Forcing:edgeAuxState1", 2);
   parallelFor(
       "Forcing:edgeAuxState1", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge) { LocMomForcingAux.computeVarsOnEdge(IEdge); });
   Pacer::stop("Forcing:edgeAuxState1", 2);
}

void Forcing::computeSurfInsituTemp(const Array3DReal &TracerArray) const {
   TracerForcingAux.computeSurfInsituTemp(TracerArray, VCoord,
                                          Eos::getInstance());
}

I4 Forcing::exchangeHalo() const {
   I4 Err = 0;

   Err +=
       MeshHalo->exchangeFullArrayHalo(MomForcingAux.ZonalStressCell, OnCell);
   Err +=
       MeshHalo->exchangeFullArrayHalo(MomForcingAux.MeridStressCell, OnCell);

   const I4 NTracers =
       SurfTracerRestAux.TracersMonthlySurfClimoCell.extent_int(0);
   for (I4 LTracer = 0; LTracer < NTracers; ++LTracer) {
      auto TracerSurfClimoCell = Kokkos::subview(
          SurfTracerRestAux.TracersMonthlySurfClimoCell, LTracer, Kokkos::ALL);
      Err += MeshHalo->exchangeFullArrayHalo(TracerSurfClimoCell, OnCell);
   }

   Err +=
       MeshHalo->exchangeFullArrayHalo(TracerForcingAux.SnowFluxCell, OnCell);
   Err +=
       MeshHalo->exchangeFullArrayHalo(TracerForcingAux.RainFluxCell, OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(TracerForcingAux.EvaporationFluxCell,
                                          OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(
       TracerForcingAux.SeaIceFreshWaterFluxCell, OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(TracerForcingAux.IceRunoffFluxCell,
                                          OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(TracerForcingAux.RiverRunoffFluxCell,
                                          OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(TracerForcingAux.LatentHeatFluxCell,
                                          OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(TracerForcingAux.SensibleHeatFluxCell,
                                          OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(
       TracerForcingAux.LongWaveHeatFluxUpCell, OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(
       TracerForcingAux.LongWaveHeatFluxDownCell, OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(TracerForcingAux.SeaIceHeatFluxCell,
                                          OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(
       TracerForcingAux.ShortWaveHeatFluxCell, OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(TracerForcingAux.SeaIceSaltFluxCell,
                                          OnCell);

   return Err;
}

} // namespace OMEGA
