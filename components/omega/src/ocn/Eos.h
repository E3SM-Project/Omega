#ifndef OMEGA_EOS_H
#define OMEGA_EOS_H
//===-- ocn/Eos.h - Equation of State --------------------*- C++ -*-===//
//
/// \file
/// \brief Contains functors for calculating specific volume
///
/// This header defines functors to be called by the time-stepping scheme
/// to calculate the specific volume based on the choice of EOS
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Config.h"
#include "EosConstants.h"
#include "HorzMesh.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"
#include <string>

namespace OMEGA {

enum class EosType {
   Linear,       /// Linear equation of state
   TEOS10Poly75t /// Roquet et al. 2015 75 term expansion
};

//
/// TEOS10 75-term polynomial
class TEOS10Poly75t {
 public:
   Array2DReal specVolPcoeffs;

   /// constructor declaration
   TEOS10Poly75t(int NVertLevels);

   KOKKOS_FUNCTION void operator()(const Array2DReal &SpecVol, I4 ICell,
                                   I4 KChunk,
                                   const Array2DReal &ConservativeTemperature,
                                   const Array2DReal &AbsoluteSalinity,
                                   const Array2DReal &Pressure) const {

      OMEGA_SCOPE(LocSpecVolPcoeffs, specVolPcoeffs);
      const I4 KStart = KChunk * VecLength;
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         calcPCoeffs(LocSpecVolPcoeffs, KVec, ConservativeTemperature(ICell, K),
                     AbsoluteSalinity(ICell, K));
         SpecVol(ICell, K) =
             calcRefProfile(Pressure(ICell, K)) +
             calcDelta(LocSpecVolPcoeffs, KVec, Pressure(ICell, K));
      }
   }

   // Note that it assumes that we have called calcPCoeffs already
   //
   KOKKOS_FUNCTION void
   calcDisplacedSpecVol(const Array2DReal &SpecVolDisplaced, I4 ICell,
                        I4 KChunk, const Array2DReal &Pressure) const {

      OMEGA_SCOPE(LocSpecVolPcoeffs, specVolPcoeffs);
      const I4 KStart = KChunk * VecLength;
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K    = KStart + KVec;
         const I4 KTmp = Kokkos::min(K + 1, NVertLevels);
         SpecVolDisplaced(ICell, K) =
             calcRefProfile(Pressure(ICell, KTmp)) +
             calcDelta(LocSpecVolPcoeffs, KVec, Pressure(ICell, KTmp));
      }
   }

   //   This member function takes point-wise conservative temperature, absolute
   //   salinity and calculate the relevant coefficients stored as data members
   KOKKOS_FUNCTION void calcPCoeffs(const Array2DReal &specVolPcoeffs,
                                    const I4 K, const Real Ct,
                                    const Real Sa) const {
      const Real SAu    = 40 * 35.16504 / 35;
      const Real CTu    = 40.;
      const Real DeltaS = 24.;
      GSW_SPECVOL_COEFFICIENTS;
      const Real ss        = Kokkos::sqrt((Sa + DeltaS) / SAu);
      Real tt              = Ct / CTu;
      specVolPcoeffs(5, K) = V005;

      specVolPcoeffs(4, K) = V014 * tt + V104 * ss + V004;
      specVolPcoeffs(3, K) =
          (V023 * tt + V113 * ss + V013) * tt + (V203 * ss + V103) * ss + V003;
      specVolPcoeffs(2, K) =
          (((V042 * tt + V132 * ss + V032) * tt + (V222 * ss + V122) * ss +
            V022) *
               tt +
           ((V312 * ss + V212) * ss + V112) * ss + V012) *
              tt +
          (((V402 * ss + V302) * ss + V202) * ss + V102) * ss + V002;
      specVolPcoeffs(1, K) =
          ((((V051 * tt + V141 * ss + V041) * tt + (V231 * ss + V131) * ss +
             V031) *
                tt +
            ((V321 * ss + V221) * ss + V121) * ss + V021) *
               tt +
           (((V411 * ss + V311) * ss + V211) * ss + V111) * ss + V011) *
              tt +
          ((((V501 * ss + V401) * ss + V301) * ss + V201) * ss + V101) * ss +
          V001;
      specVolPcoeffs(0, K) =
          (((((V060 * tt + V150 * ss + V050) * tt + (V240 * ss + V140) * ss +
              V040) *
                 tt +
             ((V330 * ss + V230) * ss + V130) * ss + V030) *
                tt +
            (((V420 * ss + V320) * ss + V220) * ss + V120) * ss + V020) *
               tt +
           ((((V510 * ss + V410) * ss + V310) * ss + V210) * ss + V110) * ss +
           V010) *
              tt +
          (((((V600 * ss + V500) * ss + V400) * ss + V300) * ss + V200) * ss +
           V100) *
              ss +
          V000;

      // could insert a check here (abs(value)> 0 value or <e+33)
   }

   KOKKOS_FUNCTION Real calcDelta(const Array2DReal &specVolPcoeffs, I4 K,
                                  const Real P) const {
      const Real Pu = 1e4;
      Real pp       = P / Pu;

      Real delta = ((((specVolPcoeffs(5, K) * pp + specVolPcoeffs(4, K)) * pp +
                      specVolPcoeffs(3, K)) *
                         pp +
                     specVolPcoeffs(2, K)) *
                        pp +
                    specVolPcoeffs(1, K)) *
                       pp +
                   specVolPcoeffs(0, K);
      return delta;
   }
   KOKKOS_FUNCTION Real calcRefProfile(const Real P) const {
      const Real Pu  = 1e4;
      const Real V00 = -4.4015007269e-05;
      const Real V01 = 6.9232335784e-06;
      const Real V02 = -7.5004675975e-07;
      const Real V03 = 1.7009109288e-08;
      const Real V04 = -1.6884162004e-08;
      const Real V05 = 1.9613503930e-09;
      Real pp        = P / Pu;

      Real v0 =
          (((((V05 * pp + V04) * pp + V03) * pp + V02) * pp + V01) * pp + V00) *
          pp;
      return v0;
   }

 private:
   const int NVertLevels;
};

/// Linear Equation of State
class LinearEOS {
 public:
   Real dRhodT  = {-0.2};   // alpha in kg.m-3 degC-1
   Real dRhodS  = {0.8};    // beta in kg m-3
   Real RhoT0S0 = {1000.0}; // density at (T,S)=(0,0) in kg.m-3

   /// constructor declaration
   LinearEOS();

   //   The functor takes the full arrays of specific volume (inout),
   //   the indices ICell and KChunk, and the ocean tracers (conservative)
   //   temperature, and (absolute) salinity as inputs, and outputs the
   //   linear specific volume.

   KOKKOS_FUNCTION void operator()(const Array2DReal &SpecVol, I4 ICell,
                                   I4 KChunk,
                                   const Array2DReal &ConservativeTemperature,
                                   const Array2DReal &AbsoluteSalinity) const {
      const I4 KStart = KChunk * VecLength;
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         SpecVol(ICell, K) =
             1.0 / (RhoT0S0 + (dRhodT * ConservativeTemperature(ICell, K) +
                               dRhodS * AbsoluteSalinity(ICell, K)));
      }
   }
};

// Eos class
class Eos {
 public:
   EosType EosChoice;
   Array2DReal SpecVol;
   Array2DReal SpecVolDisplaced;
   std::string SpecVolFldName;          ///< Field name for SpecVol
   std::string SpecVolDisplacedFldName; ///< Field name for SpecVolDisplaced
   std::string EosGroupName;
   std::string Name;

   void computeSpecVol(const Array2DReal &SpecVol,
                       const Array2DReal &ConservativeTemperature,
                       const Array2DReal &AbsoluteSalinity,
                       const Array2DReal &Pressure) const;
   static I4 init();
   static Eos *create(const std::string &Name, const HorzMesh *Mesh,
                      int NVertLevels) {
      // Check to see if eos of the same name already exist and
      // if so, exit with an error
      if (AllEos.find(Name) != AllEos.end()) {
         LOG_ERROR("Attempted to create Eos with name {} but Eos of "
                   "that name already exists",
                   Name);
         return nullptr;
      }

      // create new eos on the heap and put it in a map of
      // unique_ptrs, which will manage its lifetime
      auto *NewEos = new Eos(Name, Mesh, NVertLevels);
      AllEos.emplace(Name, NewEos);

      return NewEos; // get(Name);
   }
   ~Eos();
   // Deallocates arrays
   static void clear();

   // Remove Eos object by name
   static void erase(const std::string &Name ///< [in]
   );
   // get default eos object
   static Eos *getDefault();
   // get eos object by name
   static Eos *get(const std::string &Name ///< [in]
   );

 private:
   I4 NCellsAll;
   I4 NChunks;
   // void truncateTempSal();
   TEOS10Poly75t computeSpecVolTEOS10Poly75t;
   LinearEOS computeSpecVolLinear;

   // constructor declaration
   Eos(const std::string &Name, ///< [in] Name for eos object
       const HorzMesh *Mesh,    ///< [in] Horizontal mesh
       int NVertLevels          ///< [in] Number of vertical levels
   );
   // pointer to default eos
   static Eos *DefaultEos;
   // map with all eos objects
   static std::map<std::string, std::unique_ptr<Eos>> AllEos;
   void defineFields();

}; // end class Eos

} // namespace OMEGA
#endif
