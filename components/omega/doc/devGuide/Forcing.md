(omega-dev-forcing)=

# Forcing

This page describes design and implementation details for forcing-related
pathways in Omega, currently this includes:

- Surface stress forcing (e.g. wind stress)
- Surface thickness and tracer flux forcing (actively coupled or data-forced)
- Surface tracer restoring (soon to be ported as a field originating from the coupler)

## Surface stress forcing design

### Surface stress forcing data flow

1. External fields provide:
   - `SfcStressZonal`
   - `SfcStressMeridional`
2. Forcing compute builds `NormalStressEdge` from those fields.
3. Tendency term applies surface stress forcing to edge-normal velocity tendency.

### Surface stress forcing key classes/components

- `SfcStressForcingVars`
  - Stores surface stress cell fields and computed `NormalStressEdge`
  - Applies configured interpolation choice (`InterpType`)
- `Forcing`
  - Calls `SfcStressForcingVars::computeVarsOnEdge` from
    `computeSfcStressForcingOnEdge`
- `SfcStressForcingOnEdge` tendency term
  - Adds contribution proportional to normal stress and inverse layer
    thickness in the surface layer

### Surface stress forcing config coupling

- `Omega.SfcStress.InterpType`
  - mapped to `InterpCellToEdgeOption`
- `Omega.Tendencies.SfcStressForcingTendencyEnable`
  - gates execution of surface stress forcing tendency kernel

## Surface thickness and tracer flux forcing design

### Surface thickness and tracer flux forcing data flow

**Thickness equation pathway:**

1. External fields provide freshwater and salt flux components:
   - `SnowFlux`, `RainFlux`, `EvaporationFlux`
   - `SeaIceFreshWaterFlux`, `IceRunoffFlux`, `RiverRunoffFlux`
   - `SeaIceSaltFlux`
2. `Forcing` stores the flux fields in `TracerForcingVars`
3. The tendency term `SfcThicknessForcingOnCell` sums the freshwater and salt mass fluxes and applies them to
the surface layer pseudo-thickness.

**Tracer equation pathway:**

1. External fields provide heat and salt flux components:
   - `LatentHeatFluxEvap`, `SensibleHeatFlux`
   - `LongWaveHeatFluxUp`, `LongWaveHeatFluxDown`
   - `SeaIceHeatFlux`, `ShortWaveHeatFlux`
   - mass fluxes which add energy changes (`SnowFlux`, `RainFlux`, `IceRunoffFlux`, `RiverRunoffFlux`)
   - `SeaIceSaltFlux`
2. `Forcing` stores the flux fields in `TracerForcingVars`
3. The tendency term `SfcTracerForcingOnCell` converts the summed external heat fluxes to a conservative-temperature tendency,
  and applies the external sea-ice salt flux to the top layer salt content thus impacting salinity.

### Surface thickness and tracer flux forcing key classes/components

- `TracerForcingVars`
  - Stores 13 coupled flux cell-centered fields: 6 freshwater fluxes, 6 heat
    fluxes, and 1 salt flux component
  - Fields initialized to zero and registered in `Forcing` field group
- `SfcThicknessForcingOnCell` tendency term
  - Computes the layer mass contribution (converted to pseudo-thickness): $\sum (\text{SnowFlux} + \text{RainFlux} + \text{EvaporationFlux} + \text{SeaIceFreshWaterFlux} + \text{IceRunoffFlux} + \text{RiverRunoffFlux} + \text{SeaIceSaltFlux}) / \rho_{sw}$
  - Applied only at surface layer (top active layer) using `MinLayerCell`
- `SfcTracerForcingOnCell` tendency term
  - For temperature: adds the direct heat fluxes
    Q_{\text{sensible}} + Q_{\text{lw,up}} + Q_{\text{lw,down}} + Q_{\text{sw}}$,
    the enthalpy flux from sea-ice interactions (conduction, phase-change and enthalpy of mass flux)  Q_{\text{ice}}, the enthalpy of added freshwater mass $(\text{RainFlux} + \text{RiverRunoffFlux}) c^0_{p,sw} max(C_T(0,0,0),C_T^{\text{top}})$ (freshwater $T >= 0$), the evaporation terms ($Q_{\text{latentEvap}} + EvapFlux c^0_{p,sw} C_T^{\text{top}})
, the enthalpy change for frozen mass $(\text{SnowFlux} + \text{IceRunoffFlux}) PotEnthalpyIce$, where PotEnthalpyIce is the potential enthalpy of solid ice (which includes phase change and enthalpy of the melted mass). To first order, $PotEnthalpyIce$ is approximated by $-LatIce$, the engineering handbook value, neglecting 0.1%.
    These potential enthalpy terms are then scaled by $H_{\text{FluxFac}}$ to provide the layer tracer tendency.
  - For salinity: applies the mass salt flux with the appropriate unit conversion: $\text{SeaIceSaltFlux} \times S_{\text{FluxFac}}$
  - Applied only at surface layer using `MinLayerCell`
  - Uses tracer index validation to apply to specific tracers only
- `Forcing`
  - Manages `TracerForcingVars` instance
- `Tendencies`
  - Calls `SfcThicknessForcingOnCell` in `computePseudoThicknessTendenciesOnly`
  - Calls `SfcTracerForcingOnCell` in `computeTracerTendenciesOnly` after surface tracer restoring

### Surface thickness and tracer flux forcing config coupling

- `Omega.Tendencies.SfcThicknessForcingTendencyEnable`
  - gates execution of coupled flux thickness kernel
  - controls freshwater and salt flux forcing on sea surface height
- `Omega.Tendencies.SfcTracerForcingTendencyEnable`
  - gates execution of coupled flux tracer kernel
  - controls direct heat flux forcing on temperature and salt flux forcing on salinity

## Notes

- Currently all forcing is applied to the surface layer only. In the future, vertical spreading of river runoff contributions will be needed.
- `SeaIceFreshWaterFlux` is the pure freshwater mass from sea ice. The full mass flux from sea ice is `SeaIceFreshWaterFlux + SeaIceSaltFlux`

## Surface tracer restoring design

### Surface tracer restoring data flow

1. External fields provide target values: `TracersMonthlySurfClimoCell` (values and units should match the state variables)
2. Auxiliary-state compute forms restoring differences: `SurfTracerRestoringDiffsCell = target - tracer_surface`
3. Tendency term applies restoring only at surface layer and only for tracers selected from `SurfaceRestoring.TracersToRestore`.

### Surface tracer restoring key classes/components

- `SurfTracerRestAuxVars`
  - Inputs: `TracersMonthlySurfClimoCell`, tracer state array
  - Output: `SurfTracerRestoringDiffsCell`
  - Uses `MinLayerCell` to select surface layer index
- `SurfaceTracerRestoringOnCell` tendency term
  - Applies `PistonVelocity * SurfTracerRestoringDiffsCell` at surface
- `Tendencies`
  - Parses `SurfaceRestoring.TracersToRestore` and resolves tracer indices
  - Builds `TracerIdsToRestore` and `NTracersToRestore`
  - Applies tracer-selection logic at call site in
    `computeTracerTendenciesOnly`
  - Aborts if restoring is enabled but no tracer IDs are available

### Surface tracer restoring config coupling

- `Omega.SurfaceRestoring.PistonVelocity`
  - tendency scaling
- `Omega.SurfaceRestoring.TracersToRestore`
  - tracer-level enable list used to build `TracerIdsToRestore`
- `Omega.Tendencies.SurfaceTracerRestoringEnable`
  - gates restoring tendency execution

## Notes

- If a tracer is not listed in `TracersToRestore`, no restoring tendency is applied to that tracer.
- If restoring is enabled but no tracer IDs are available at tendency compute-time, Omega aborts with an error.
- It is assumed that the incoming `TracersMonthlySurfClimoCell` fields (values and units) match the Omega state variables (i.e. conservative temperature and absolute salinity for Teos-10). If not, a pre-processing conversion should be implemented.
- Surface tracer restoring is active everywhere if enabled. A flag to turn it off under sea ice will need to be added in later development if this feature is desired.
- Unlike MPAS-Ocean, a `MaxDiff` clamping is not applied here. This check should instead be implemented in Ocean Validate when that is available.
- A global scaling to ensure zero-sum has not been implemented for the surface tracer restoring, but should be added in later development.
- At this stage, there is no temporal interpolation applied to the restoring targets, the raw `TracersMonthlySurfClimoCell` snapshot is used.
