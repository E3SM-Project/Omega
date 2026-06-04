(omega-dev-forcing)=

# Forcing

This page describes design and implementation details for forcing-related
pathways in Omega, currently this includes:

- Surface momemtum forcing (e.g. stress from wind forcing)
- Surface thickness and tracer forcing
- Surface tracer restoring

## Surface momentum forcing design

### Surface stress forcing data flow

1. External fields provide:
   - `WindStressZonal`
   - `WindStressMeridional`
2. Auxiliary-state compute builds `NormalStressEdge` from those fields.
3. Tendency term applies stress forcing to edge-normal velocity tendency.

### Surface momemtum forcing key classes/components

- `MomForcingAuxVars`
  - Stores surface stress cell fields and computed `NormalStressEdge`
  - Applies configured interpolation choice (`InterpType`)
- `AuxiliaryState::computeMomAux`
  - Calls `MomForcingAuxVars::computeVarsOnEdge`
- `SrfStressForcingOnEdge` tendency term
  - Adds contribution proportional to normal stress and inverse layer
    thickness in the surface layer

### Surface stress forcing config coupling

- `Omega.SrfStress.InterpType`
  - mapped to `InterpCellToEdgeOption`
- `Omega.Tendencies.SrfStressForcingTendencyEnable`
  - gates execution of surface stress (e.g. wind forcing) tendency kernel

## Surface flux forcing design

### Surface flux forcing data flow

**Thickness equation pathway:**

1. External fields provide freshwater and salt flux components:
   - `SnowFlux`, `RainFlux`, `EvaporationFlux`
   - `SeaIceFreshWaterFlux`, `IceRunoffFlux`, `RiverRunoffFlux`
   - `SeaIceSaltFlux`
2. `Forcing` stores the flux fields in `TracerForcingAuxVars`
3. The tendency term `SrfThicknessForcing`sums both the freshwater and salt mass fluxes, converted to be applied to
the surface layer pseudo-thickness.

**Tracer equation pathway:**

1. External fields provide heat and salt flux components:
   - `LatentHeatFlux`, `SensibleHeatFlux`
   - `LongWaveHeatFluxUp`, `LongWaveHeatFluxDown`
   - `SeaIceHeatFlux`, `ShortWaveHeatFlux`
   - `SeaIceSaltFlux`, `SnowFlux`, `IceRunoffFlux`
2. `Forcing` stores the flux fields in `TracerForcingAuxVars`
3. The tendency terms for each tracer `SrfTracerForcing` converts external heat fluxes to tendencies in conservative temperature,
 and external mass salt flux to salinity (g/kg) in the surface layer.

### Coupled flux forcing key classes/components

- `TracerForcingAuxVars`
  - Stores 13 coupled flux cell-centered fields: 7 freshwater fluxes and 6 heat
    fluxes, plus 1 salt flux component
  - Fields initialized to zero and registered in `Forcing` field group
- `SrfThicknessForcingOnCell` tendency term
  - Computes freshwater flux contribution: $\sum (\text{SnowFlux} + \text{RainFlux} + \text{EvaporationFlux} + \text{SeaIceFreshWaterFlux} + \text{IceRunoffFlux} + \text{RiverRunoffFlux} + \text{SeaIceSaltFlux}) / \rho_{sw}$
  - Applied only at surface layer (top active layer) using `MinLayerCell`
- `SrfTracerForcingOnCell` tendency term
  - For temperature: computes heat flux minus latent heat of fusion contribution: $(\sum \text{HeatFluxes} - (\text{SnowFlux} + \text{IceRunoffFlux}) L_i) \times H_{\text{FluxFac}}$
  - For salinity: applies salt flux with unit conversion: $\text{SeaIceSaltFlux} \times S_{\text{FluxFac}}$
  - Applied only at surface layer using `MinLayerCell`
  - Uses tracer index validation to apply to specific tracers only
- `Forcing`
  - Manages `TracerForcingAuxVars` instance
- `Tendencies`
  - Calls `SrfThicknessForcingOnCell` in `computeThicknessTendenciesOnly`
  - Calls `SrfTracerForcingOnCell` in `computeTracerTendenciesOnly` after surface tracer restoring

### Coupled flux forcing config coupling

- `Omega.Tendencies.SrfThicknessForcingTendencyEnable`
  - gates execution of coupled flux thickness kernel
  - controls freshwater and salt flux forcing on sea surface height
- `Omega.Tendencies.SrfTracerForcingTendencyEnable`
  - gates execution of coupled flux tracer kernel
  - controls heat flux forcing on temperature and salt flux forcing on salinity

## Surface tracer restoring design

### Surface tracer restoring data flow

1. External fields provide target values: `TracersMonthlySurfClimoCell` (values and units should match the state variables)
2. `SurfTracerRestAuxVars` stores `TracersMonthlySurfClimoCell` for later restoring use.
3. Tendency term applies restoring only at surface layer and only for tracers selected from `SrfRestoring.TracersToRestore`.

### Surface tracer restoring key classes/components

- `SurfTracerRestAuxVars`
  - Stores `TracersMonthlySurfClimoCell`
- `SrfTracerRestoringOnCell` tendency term
  - Applies `PistonVelocity * (target - state)` at the surface layer
- `Tendencies`
  - Parses `SrfRestoring.TracersToRestore` and resolves tracer indices
  - Builds `TracerIdsToRestore` and `NTracersToRestore`
  - Applies tracer-selection logic at call site in
    `computeTracerTendenciesOnly`
  - Aborts if restoring is enabled but no tracer IDs are available

### Surface tracer restoring config coupling

- `Omega.SrfRestoring.PistonVelocity`
  - tendency scaling
- `Omega.SrfRestoring.TracersToRestore`
  - tracer-level enable list used to build `TracerIdsToRestore`
- `Omega.Tendencies.SrfTracerRestoringEnable`
  - gates restoring tendency execution

## Notes

- If a tracer is not listed in `TracersToRestore`, no restoring tendency is applied to that tracer.
- If restoring is enabled but no tracer IDs are available at tendency compute-time, Omega aborts with an error.
- It is assumed that the incoming `TracersMonthlySurfClimoCell` fields (values and units) match the Omega state variables (i.e. conservative temperature and absolute salinity for Teos-10). If not, a pre-processing conversion should be implemented.
- Surface tracer restoring is active everywhere if enabled. A flag to turn it off under sea ice will need to be added in later development if this feature is desired.
- Unlike MPAS-Ocean, a `MaxDiff` clamping is not applied here. This check should instead be implemented in Ocean Validate when that is available.
- A global scaling to ensure zero-sum has not been implemented for the surface tracer restoring, but should be added in later development.
- At this stage, there is no temporal interpolation applied to the restoring targets, the raw `TracersMonthlySurfClimoCell` snapshot is used.
