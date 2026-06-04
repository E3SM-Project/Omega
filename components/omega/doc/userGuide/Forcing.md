(omega-user-forcing)=

# Forcing

This page documents the user-facing configuration and behavior for current forcing in Omega:

- Wind forcing
- Coupled flux forcing
- Surface tracer restoring

## Wind forcing

Wind forcing adds momentum tendency from surface wind stress.

### Wind forcing configuration

Wind forcing behavior is controlled by two configuration blocks:

```yaml
Omega:
  SrfStress:
    InterpType: Isotropic

  Tendencies:
    SrfStressForcingTendencyEnable: true
```

- `SrfStress.InterpType`
  - `Isotropic`: isotropic cell-to-edge interpolation for wind stress
  - `Anisotropic`: anisotropic interpolation option
- `Tendencies.SrfStressForcingTendencyEnable`: switch to enable wind forcing tendency

### Required input fields

Wind forcing uses auxiliary wind-stress fields:

- `WindStressZonal`
- `WindStressMeridional`

These are used to form edge-normal stress (`NormalStressEdge`) that enters
momentum tendencies.

## Surface flux forcing

Surface flux forcing applies ocean-atmosphere and ocean-sea ice fluxes from the other model
components (atmosphere, sea ice) to the thickness and tracer equations. This enables
the ocean to respond to heat, freshwater, and salt exchanges at the surface. These fluxes can be from data or (active) coupled components.

### Surface flux forcing configuration

Surface flux forcing is controlled by two configuration flags:

```yaml
Omega:
  Tendencies:
    SrfThicknessForcingTendencyEnable: false
    SrfTracerForcingTendencyEnable: false
```

- `Tendencies.SrfThicknessForcingTendencyEnable`: enables coupled freshwater and salt flux forcing on thickness
- `Tendencies.SrfTracerForcingTendencyEnable`: enables coupled heat and salt flux forcing on tracers

### Required input fields

Coupled flux forcing uses 13 auxiliary fields organized by type:

**Freshwater mass fluxes (kg m⁻² s⁻¹):**
- `SnowFlux`: precipitation from snow
- `RainFlux`: precipitation from rain
- `EvaporationFlux`: evaporative water loss
- `SeaIceFreshWaterFlux`: freshwater input from sea-ice melt or formation
- `IceRunoffFlux`: runoff from land ice
- `RiverRunoffFlux`: runoff from rivers

**Heat fluxes (W m⁻²):**
- `LatentHeatFlux`: latent heat transfer
- `SensibleHeatFlux`: sensible heat transfer
- `LongWaveHeatFluxUp`: upward longwave radiation
- `LongWaveHeatFluxDown`: downward longwave radiation
- `SeaIceHeatFlux`: heat from sea-ice interaction
- `ShortWaveHeatFlux`: shortwave (solar) radiation

**Salt mass flux (kg m⁻² s⁻¹):**
- `SeaIceSaltFlux`: salt flux from sea-ice formation/melt processes

These fields are populated by external coupling components (typically atmosphere
and ice models). Omega assumes the incoming values match the documented units.
For now, there are assumed to come from a `forcing.nc` file, but later will be provided
by the equivalent `ocn_comp_mct.F`.

### Notes

- Coupled fluxes are applied only at the surface layer (top active layer) for each cell.
- Temperature tendency is computed as the sum of all heat fluxes minus latent heat
  of fusion for snow and ice runoff, converted to temperature tendency via
  $H_{\text{FluxFac}} = 1.0 / (\rho_{sw} c^0_{p,sw})$ where $c^0_{p,sw}$ is the reference
  specific heat of seawater defined by TEOS-10. This allows the conversion to the
  conservative temperature variable.
- Salinity tendency from salt flux is scaled by $S_{\text{FluxFac}} = 1.0e3 / \rho_{sw}$
  to account for unit conversion from kg/(m²·s) to salinity units (g/kg).
- Fluxes are assumed to be in the documented units (i.e. net mass fluxes);
  any unit conversion should be performed by the coupling component before providing flux
  values to Omega.
- The reference density used here ($\rho_{sw}$) is not a Boussinesq density, it is the
  conversion factor from mass to pseudo-thickness.
- No iceberg fluxes are included for now.

## Surface tracer restoring

Surface tracer restoring applies a piston-velocity tendency, or damping, at the ocean
surface for selected tracers. This is implemented to mitigate drifts in chosen tracers
(most often salinity) by nudging the model's simulated tracer values towards observed climatological values.
This process prevents oceanic regimes from shifting away from reality due to errors in surface freshwater
forcing (in the case of salinity restoring). Currently, it is applied everywhere when enabled.

### Surface tracer restoring configuration

Surface tracer restoring is controlled by two configuration blocks:

```yaml
Omega:
  SrfRestoring:
    TracersToRestore: [Temperature, Salinity]
    PistonVelocity: 1.585e-5

  Tendencies:
    SrfTracerRestoringEnable: true
```

- `TracersToRestore`: list of tracer names that restoring is applied to
- `PistonVelocity`: restoring rate coefficient
- `SrfTracerRestoringEnable`: switch to enable surface tracer restoring

When restoring is enabled, Omega resolves `TracersToRestore` into an internal
list of tracer IDs and applies restoring only to tracers in that list.

### Restoring target fields

Surface restoring uses auxiliary fields:

- `TracersMonthlySurfClimoCell`: restoring target climatological values

The restoring tendency `SrfTracerRestoring` is computed at the surface layer only using the
configured `PistonVelocity` and the inline target-minus-state difference.

## Notes

- If a tracer is not listed in `TracersToRestore`, no restoring tendency is
  applied to that tracer.
- If surface restoring is enabled but no tracer IDs are available at tendency
  compute-time, Omega aborts with an error.
