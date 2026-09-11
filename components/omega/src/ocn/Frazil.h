#ifndef OMEGA_FRAZIL_H
#define OMEGA_FRAZIL_H
//===-- ocn/Frazil.h - Frazil Ice Formation -------------------*- C++ -*-===//
//
// The Frazil class manages frazil tendencies and accumulators.
// This initial implementation only has a teos-10 configuration.
// but carries scaffolding for other implementations.
//
//===----------------------------------------------------------------------===//

#include "Config.h"
// #include "DataTypes.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <map>
#include <memory>
#include <string>

#include <gswteos-10.h>

namespace OMEGA {

enum class FrazilType {
   FixedPropertyFrazil, ///< fixed-property frazil option
   TeosFrazil           ///< TEOS frazil option
};

class FixedPropertyFrazilFormation {
 public:
   FixedPropertyFrazilFormation();

   Real massLimit         = 0.1_Real;  // to do:  remove default
   Real FrazilIceSalinity = IceRefSal; // Global constant
   Real LatFrazil         = LatIce;    // Global constant

   Real FrazilPorosity =
       1.0_Real; // Internal for now; can move to config later.

   KOKKOS_FUNCTION void operator()(const Real SA, const Real CT, const Real PDb,
                                   const Real H, Real &SumIceThickness,
                                   Real &SumSalt, Real &SumEnergy, Real &HTend,
                                   Real &TTend, Real &STend,
                                   const Real Tfrz) const {

      const Real potential      = H * Cp0Sw * RhoSw * (CT - Tfrz);
      const Real freezingEnergy = Kokkos::max(0.0_Real, -potential);

      HTend = 0.0_Real;
      TTend = 0.0_Real;
      STend = 0.0_Real;

      Real newFrzThickness =
          freezingEnergy /
          (LatFrazil * RhoSw); // frazil (ice) mass in pseudo-thickness terms

      newFrzThickness = Kokkos::min(newFrzThickness, H * massLimit);
      Real newFrzEnergy =
          newFrzThickness *
          (-LatFrazil + Cp0Sw * Tfrz); // (<0; enthalpy of frazil, i.e. phase
                                       // change and enthalpy of melted equ)

      // MANUAL TOGGLE: uncomment line below to use porosity
      // const Real FrazilIceSalinity = FrazilPorosity * SA;

      const Real frazilSalinity = Kokkos::min(FrazilIceSalinity, SA);
      const Real newSaltContent =
          newFrzThickness * frazilSalinity; // in m.(g/kg)

      // HTend does not include the salt mass contribution to align with mpas-o.
      HTend = -newFrzThickness;
      // TTend includes the enthalpy associated with the mass flux to be
      // conservative with the melt impl. (this differs from old mpas-o)
      TTend =
          -(newFrzEnergy) / (Cp0Sw); // (E< 0 so TTend>0) // scaled to h.CT tend
      STend = -newSaltContent;

      SumIceThickness += newFrzThickness;
      SumSalt += newSaltContent;
      SumEnergy += newFrzEnergy;
   }
};

class FixedPropertyFrazilMelt {
 public:
   FixedPropertyFrazilMelt();

   Real massLimit = 0.1_Real; // to do:  remove default
   Real LatFrazil = LatIce;   // Global constant

   KOKKOS_FUNCTION void operator()(const Real SA, const Real CT, const Real PDb,
                                   const Real H, Real &SumIceThickness,
                                   Real &SumSalt, Real &SumEnergy, Real &HTend,
                                   Real &TTend, Real &STend,
                                   const Real Tfrz) const {
      constexpr Real Eps = 1.0e-12_Real;

      if (SumIceThickness <= Eps) { // potential leak if we dont redistribute
         SumIceThickness = 0.0_Real;
         SumSalt         = 0.0_Real;
         HTend           = 0.0_Real;
         TTend           = 0.0_Real;
         STend           = 0.0_Real;
         return;
      }
      const Real potential       = H * Cp0Sw * RhoSw * (CT - Tfrz);
      const Real availableEnergy = Kokkos::max(0.0_Real, potential);

      HTend = 0.0_Real;
      TTend = 0.0_Real;
      STend = 0.0_Real;

      Real meltThickness =
          availableEnergy /
          (LatFrazil * RhoSw); // mass in pseudo-thickness units
      meltThickness = Kokkos::min(meltThickness, SumIceThickness);
      meltThickness = Kokkos::min(meltThickness,
                                  H * massLimit); // also 0.1h lim on added mass
      const Real frazilFractionMelted =
          meltThickness / SumIceThickness; // mass fraction melted
      // const Real meltAverageSalinity = SumSalt / SumIceThickness;
      const Real meltEnergy = frazilFractionMelted * SumEnergy;

      HTend = frazilFractionMelted * (SumIceThickness); // (>0 so HTend>0)
      TTend = frazilFractionMelted * SumEnergy /
              Cp0Sw; // (SumE <0 thus TTend < 0 when melting for phase change)
      STend = +frazilFractionMelted * SumSalt; // (STend > 0 when melting)

      const Real frazilFractionLeft =
          Kokkos::max(0.0_Real, 1.0_Real - frazilFractionMelted);
      SumIceThickness = frazilFractionLeft * SumIceThickness;
      SumSalt         = frazilFractionLeft * SumSalt;
      SumEnergy =
          frazilFractionLeft * SumEnergy; // conservative by construction
   }
};

class FrazilMelt {
 public:
   /// constructor declaration
   FrazilMelt();

   // masslimit parameter (set in config)
   Real massLimit;

   //   The functor for FrazilMelt takes as inputs:
   //   the local ocean layer state (SA, CT, P, h),
   //   the accumulated frazil solid and liquid mass, energy, and salt
   //   and outputs the frazil tendencies (HTend, TTend, STend) and updated
   //   accumulators.
   //   Host-only: relies on GSW TEOS-10 routines that are not device-callable.
   //   This is a temporary implementation until a device-callable solution is
   //   available.
   void operator()(const Real SA, const Real CT, const Real P, const Real h,
                   Real &AccMIce, Real &AccMLiq, Real &AccMSalt, Real &AccELiq,
                   Real &AccEIce, Real &HTend, Real &TTend, Real &STend) const {

      constexpr Real Eps = 1.0e-12_Real;

      // this check on AccMIce and E should be done in the calling function, but
      // is here for safety we can do a better implementation of the checks
      if (AccMIce <= Eps ||
          AccEIce >= Eps) { // potential leak if we dont redistribute
         AccMIce = 0.0_Real;
         AccEIce = 0.0_Real;
         HTend   = 0.0_Real;
         TTend   = 0.0_Real;
         STend   = 0.0_Real;
         return;
      }

      if (AccMLiq <= Eps ||
          AccELiq >= Eps) { // potential leak if we dont redistribute
         AccMLiq = 0.0_Real;
         AccELiq = 0.0_Real;
         HTend   = 0.0_Real;
         TTend   = 0.0_Real;
         STend   = 0.0_Real;
         return;
      }

      // 1. we start by adding the solid ice to the ocean layer (no brine yet)
      const Real potEnthalpyIce = AccEIce / AccMIce;
      const Real newLayerMass =
          Kokkos::max(h + AccMIce, Eps); // max unnecessary but for safety
      const Real newLayerIceFraction = AccMIce / newLayerMass;

      // 2. we calculate the (mass- and energy-conserving) ocean layer evolution
      // ... how much (pure) ice does this layer melt?

      // typecasting for now but will be simplified once gsw functions are
      // ported
      const double SA_d     = static_cast<double>(SA);
      const double CT_d     = static_cast<double>(CT);
      const double P_d      = static_cast<double>(P);
      const double wIhIn_d  = static_cast<double>(newLayerIceFraction);
      const double pt0Ice_d = gsw_pt_from_pot_enthalpy_ice_poly(
          static_cast<double>(potEnthalpyIce));
      const double tIce_d = gsw_t_from_pt0_ice(pt0Ice_d, P_d);
      double SAnew_d      = SA_d;
      double CTnew_d      = CT_d;
      double wIhOut_d     = wIhIn_d;

      gsw_melting_ice_into_seawater(SA_d, CT_d, P_d, wIhIn_d, tIce_d, &SAnew_d,
                                    &CTnew_d, &wIhOut_d);

      const Real wIhOut = static_cast<Real>(
          wIhOut_d); // by def 0<= wIhOut <= 1 ;  To-do: check function behavior

      // 3. we calculate the mass fraction of the frazil (pure) ice that was
      // melted - limited by a total mass limit of 0.1h
      const Real solidMassMelted = Kokkos::max(
          0.0_Real,
          AccMIce - wIhOut * newLayerMass); // original - left-over solid ice,
      const Real frazilFractionMelted = Kokkos::min(
          solidMassMelted / AccMIce,
          h * massLimit / (AccMIce + AccMLiq)); // total added mass < 0.1h

      // the frazil fraction based on the solid ice also sets the (proportional)
      // contributions from the frazil brine
      HTend = +(frazilFractionMelted * (AccMLiq + AccMIce));
      TTend = +(frazilFractionMelted * (AccELiq + AccEIce)) / Cp0Sw;
      STend = +(frazilFractionMelted * AccMSalt);

      const Real frazilFractionLeft =
          Kokkos::max(0.0_Real, 1.0_Real - frazilFractionMelted);

      AccMIce  = frazilFractionLeft * AccMIce;
      AccMLiq  = frazilFractionLeft * AccMLiq;
      AccMSalt = frazilFractionLeft * AccMSalt;
      AccELiq  = frazilFractionLeft * AccELiq;
      AccEIce  = frazilFractionLeft * AccEIce;
   }
};

class FrazilFormation {
 public:
   /// constructor declaration
   FrazilFormation();

   /// Parameters for FrazilFormation (set by yaml file)
   Real phi;       ///< liquid mass fraction of new frazil (0 < Phi < 1)
   Real massLimit; ///< layer mass fraction limit for thickness tendency

   //   The functor for FrazilFormation takes as inputs:
   //   the local ocean layer state (SA, CT, P, h),
   //   the accumulated frazil solid and liquid mass, energy, and salt
   //   and outputs the frazil tendencies (HTend, TTend, STend) and updated
   //   accumulators.
   //   Host-only: relies on GSW TEOS-10 routines that are not device-callable.
   //   This is a temporary implementation until a device-callable solution is
   //   available.
   void operator()(const Real SA, const Real CT, const Real P, const Real h,
                   Real &AccMIce, Real &AccMLiq, Real &AccMSalt, Real &AccELiq,
                   Real &AccEIce, Real &HTend, Real &TTend, Real &STend) const {

      Real CTnew;
      Real SAnew;
      Real wIh            = 0.0_Real;
      Real solidMass      = 0.0_Real;
      Real liquidMass     = 0.0_Real;
      Real solidEnthalpy  = 0.0_Real;
      Real liquidEnthalpy = 0.0_Real;

      // type casting for now because gsw expects double
      // and Real can be either float or double depending on build configuration
      double SAnew_d = 0.0;
      double CTnew_d = 0.0;
      double wIh_d   = 0.0;
      // double PTnew_d = 0.0;

      gsw_frazil_properties_potential_poly(
          static_cast<double>(SA), static_cast<double>(Cp0Sw * CT),
          static_cast<double>(P), &SAnew_d, &CTnew_d, &wIh_d);
      double PTnew_d =
          gsw_pt_from_ct(SAnew_d, CTnew_d); // convert to potential temperature

      SAnew = static_cast<Real>(SAnew_d);
      CTnew = static_cast<Real>(CTnew_d);
      wIh   = static_cast<Real>(wIh_d);

      const Real oneMinusPhi = Kokkos::max(1.0e-12_Real, 1.0_Real - phi);
      // anything called mass below is in pseudo-thickness units (m) and needs
      // to be scaled by RhoSw for coupling
      solidMass      = h * Kokkos::min(wIh, oneMinusPhi * massLimit);
      liquidMass     = (phi / oneMinusPhi) * solidMass;
      solidEnthalpy  = solidMass * gsw_pot_enthalpy_from_pt_ice_poly(PTnew_d);
      liquidEnthalpy = liquidMass * Cp0Sw * CTnew;
      // per timestep (not scaled by dt here)
      HTend = -(solidMass +
                liquidMass); // because phi is a *mass* fraction, liquidMass
                             // includes the salt contribution to mass.
      TTend = -(liquidEnthalpy + solidEnthalpy) / Cp0Sw;
      STend = -(liquidMass * SAnew);

      // Local unit of mass is pseudo thickness (m)
      // these all need a RhoSw factor before coupling
      AccMIce += solidMass;           // m
      AccMLiq += liquidMass;          // m
      AccMSalt += liquidMass * SAnew; // (m)(g/kg)
      AccELiq += liquidEnthalpy;      // (m)(J/kg)
      AccEIce += solidEnthalpy;       // (m)(J/kg)
   }
};

class Frazil {
 public:
   static void init();
   /// Creates a new frazil object and stores it in the AllFrazil map.
   static Frazil *create(const std::string &Name);

   /// Retrieve frazil object by name.
   static Frazil *get(const std::string &Name);

   /// Retrieve default frazil object.
   static Frazil *getDefault();

   /// Destructor
   ~Frazil();

   /// Deallocates arrays
   static void clear();

   /// Remove frazil object by name.
   static void erase(std::string InName); ///< [in] name to remove

   Array2DReal FrazilTTend;
   Array2DReal FrazilSTend;
   Array2DReal FrazilHTend;
   Array1DReal AccMIce;
   Array1DReal AccEIce;
   Array1DReal AccMLiq;
   Array1DReal AccELiq;
   Array1DReal AccMSalt;

   void computeFrazil(const Array2DReal &CT, const Array2DReal &SA,
                      const Array2DReal &P, const Array2DReal &H);
   void computeFrazilFixedPropertyImpl(const Array2DReal &CT,
                                       const Array2DReal &SA,
                                       const Array2DReal &P,
                                       const Array2DReal &LayerH);
   void computeFrazilTeosImpl(const Array2DReal &CT, const Array2DReal &SA,
                              const Array2DReal &P, const Array2DReal &LayerH);
   bool conservationCheck = false;
   Real depthLimit        = -1.0_Real;

 private:
   static Frazil *DefaultFrazil;
   static std::map<std::string, std::unique_ptr<Frazil>> AllFrazil;

   Frazil(const HorzMesh *Mesh, const VertCoord *VCoord);

   // Forbid copy and move construction/assignment.
   Frazil(const Frazil &)            = delete;
   Frazil &operator=(const Frazil &) = delete;
   Frazil(Frazil &&)                 = delete;
   Frazil &operator=(Frazil &&)      = delete;

   FrazilType frazilChoice;
   FixedPropertyFrazilFormation computeFixedPropertyFrazilFormation;
   FixedPropertyFrazilMelt computeFixedPropertyFrazilMelt;
   FrazilFormation computeFrazilFormation;
   FrazilMelt computeFrazilMelt;
   I4 NCellsAll;
   I4 NChunks;

   const HorzMesh *MeshPtr;
   const VertCoord *VCoordPtr;

   void checkColumnConservation() const;
};

} // namespace OMEGA

#endif
