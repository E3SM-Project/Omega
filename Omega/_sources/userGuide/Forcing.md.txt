(omega-user-forcing)=

# Forcing

This page documents the user-facing configuration and behavior for current forcing in Omega:

- Surface stress forcing (e.g. wind stress)
- Surface tracer restoring

## Surface stress forcing

Surface stress forcing adds momentum tendency from surface stress (e.g. wind).

### Surface stress forcing configuration

Surface stress forcing behavior is controlled by two configuration blocks:

```yaml
Omega:
  SfcStress:
    InterpType: Isotropic

  Tendencies:
    SfcStressForcingTendencyEnable: true
```

- `SfcStress.InterpType`
  - `Isotropic`: isotropic cell-to-edge interpolation for surface stress
  - `Anisotropic`: anisotropic interpolation option
- `Tendencies.SfcStressForcingTendencyEnable`: switch to enable surface stress forcing tendency

### Required input fields

Surface stress forcing uses surface stress input fields:

- `SfcStressZonal`
- `SfcStressMeridional`

These are stored in forcing variables and used to form edge-normal stress
(`NormalStressEdge`) that enters momentum tendencies.

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
  SurfaceRestoring:
    TracersToRestore: [Temperature, Salinity]
    PistonVelocity: 1.585e-5

  Tendencies:
    SurfaceTracerRestoringEnable: true
```

- `TracersToRestore`: list of tracer names that restoring is applied to
- `PistonVelocity`: restoring rate coefficient
- `SurfaceTracerRestoringEnable`: switch to enable surface tracer restoring

When restoring is enabled, Omega resolves `TracersToRestore` into an internal
list of tracer IDs and applies restoring only to tracers in that list.

### Restoring target fields

Surface restoring uses auxiliary fields:

- `TracersMonthlySurfClimoCell`: restoring target climatological values
- `SurfTracerRestoringDiffsCell`: computed target-minus-state differences

The restoring tendency is computed at the surface layer only and is limited by
the configured `PistonVelocity` and target-minus-state difference.

## Notes

- If a tracer is not listed in `TracersToRestore`, no restoring tendency is
  applied to that tracer.
- If surface restoring is enabled but no tracer IDs are available at tendency
  compute-time, Omega aborts with an error.
