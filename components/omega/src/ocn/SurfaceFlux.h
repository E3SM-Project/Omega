#ifndef OMEGA_SURFACE_FLUX_H
#define OMEGA_SURFACE_FLUX_H
//===-- ocn/SurfaceFlux.h - Surface Flux Helpers ----------*- C++ -*-===//
//
/// \file
/// \brief Shared surface heat and freshwater flux calculations
//
//===----------------------------------------------------------------------===//

#include "Eos.h"
#include "GlobalConstants.h"

namespace OMEGA {

/// Enthalpy (W/m^2) of the sea-ice meltwater included in SeaIceHeatFlux,
/// taken at the freezing CT of the meltwater (sea-ice salinity when known,
/// otherwise the top-layer salinity).
KOKKOS_INLINE_FUNCTION Real
sfcSeaIceLiqEnthalpy(I4 ICell, Real SaTop, Real PTopDb, EosType EosChoice,
                     const Array1DReal &SeaIceFreshWaterFlux,
                     const Array1DReal &SeaIceSaltFlux) {
   const Real SeaIceFreshWater = SeaIceFreshWaterFlux(ICell);
   const Real SeaIceSalt       = SeaIceSaltFlux(ICell);
   const Real SaSeaIce = SeaIceFreshWater > 0.0_Real && SeaIceSalt > 0.0_Real
                             ? SeaIceSalt / SeaIceFreshWater * Salt2PPt
                             : SaTop;
   return SeaIceFreshWater * Cp0Sw *
          Eos::calcCtFreezing(EosChoice, SaSeaIce, PTopDb, 0.0_Real);
}

/// Surface heat flux used by KPP buoyancy forcing, excluding mass-flux
/// enthalpy carried by the tracer tendency.
KOKKOS_INLINE_FUNCTION Real sfcHeatFluxWithoutMassEnthalpy(
    I4 ICell, Real SaTop, Real PTopDb, EosType EosChoice,
    const Array1DReal &LongWaveHeatFluxUp,
    const Array1DReal &LongWaveHeatFluxDown,
    const Array1DReal &ShortWaveHeatFlux, const Array1DReal &SensibleHeatFlux,
    const Array1DReal &SeaIceHeatFlux, const Array1DReal &SeaIceFreshWaterFlux,
    const Array1DReal &SeaIceSaltFlux, const Array1DReal &LatentHeatFluxEvap,
    const Array1DReal &SnowFlux, const Array1DReal &IceRunoffFlux) {

   const Real SeaIceLiqEnthalpyEstimate = sfcSeaIceLiqEnthalpy(
       ICell, SaTop, PTopDb, EosChoice, SeaIceFreshWaterFlux, SeaIceSaltFlux);

   return LongWaveHeatFluxUp(ICell) + LongWaveHeatFluxDown(ICell) +
          ShortWaveHeatFlux(ICell) + SensibleHeatFlux(ICell) +
          SeaIceHeatFlux(ICell) - SeaIceLiqEnthalpyEstimate +
          LatentHeatFluxEvap(ICell) +
          (SnowFlux(ICell) + IceRunoffFlux(ICell)) * -LatIce;
}

/// Net surface freshwater mass flux (kg/m^2/s).
KOKKOS_INLINE_FUNCTION Real sfcFreshWaterFlux(
    I4 ICell, const Array1DReal &SnowFlux, const Array1DReal &RainFlux,
    const Array1DReal &EvaporationFlux, const Array1DReal &SeaIceFreshWaterFlux,
    const Array1DReal &IceRunoffFlux, const Array1DReal &RiverRunoffFlux) {
   return SnowFlux(ICell) + RainFlux(ICell) + EvaporationFlux(ICell) +
          SeaIceFreshWaterFlux(ICell) + IceRunoffFlux(ICell) +
          RiverRunoffFlux(ICell);
}

/// Fraction of the surface shortwave flux that penetrates below Depth (m),
/// using the three-band (near-infrared, red, blue) exponential profile.
KOKKOS_INLINE_FUNCTION Real sfcShortWaveTransmission(
    Real Depth, Real NearIrFraction, Real NearIrCoeff, Real RedFraction,
    Real Kr, Real BlueFraction, Real Kb) {
   return NearIrFraction * Kokkos::exp(-NearIrCoeff * Depth) +
          RedFraction * Kokkos::exp(-Kr * Depth) +
          BlueFraction * Kokkos::exp(-Kb * Depth);
}

/// Fraction of the surface shortwave flux that penetrates below
/// SurfaceBuoyancyDepth (m) in a column of thickness ColumnDepth (m).  All of
/// the shortwave is absorbed in columns shallower than SurfaceBuoyancyDepth.
KOKKOS_INLINE_FUNCTION Real sfcShortWaveTransmissionBelowDepth(
    Real SurfaceBuoyancyDepth, Real ColumnDepth, Real NearIrFraction,
    Real NearIrCoeff, Real RedFraction, Real Kr, Real BlueFraction, Real Kb) {
   if (ColumnDepth <= SurfaceBuoyancyDepth)
      return 0.0_Real;
   return sfcShortWaveTransmission(SurfaceBuoyancyDepth, NearIrFraction,
                                   NearIrCoeff, RedFraction, Kr, BlueFraction,
                                   Kb);
}

} // namespace OMEGA

#endif // OMEGA_SURFACE_FLUX_H
