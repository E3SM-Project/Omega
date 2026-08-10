(omega-dev-kppmix)=

# KPP Boundary Layer Mixing

This page maps OMEGA KPP implementation details to runtime behavior and code
locations. It complements the design page by focusing on concrete APIs,
call flow, and developer test strategy.

Related pages:
- Design and theory: [Design KPP document](../design/KPPMix.md)
- User configuration and workflow: [User KPP guide](../userGuide/KPPMix.md)
- Broader vertical mixing: [Developer Vertical Mixing Coefficients](./VerticalMixingCoeff.md)

## Implementation Overview

OMEGA KPP is implemented in `KPPMix` as a singleton with two major compute
phases:

1. OBL depth diagnosis (`computeOBLDepth`)
2. Coefficient/profile construction (`computeMixingCoefficients`)

Main class/API surface is in `src/ocn/KPPMix.h` and implementation is in
`src/ocn/KPPMix.cpp`.

## Runtime Call Flow

### Tendency coupling

KPP coupling into tendencies occurs through:

- `Tendencies::computeAllTendencies(...)`
- `Tendencies::computeStageVerticalMixing(...)`

`computeStageVerticalMixing(...)` assembles required inputs:

- potential density from EOS specific volume
- Brunt-Vaisala frequency squared
- edge normal and reconstructed tangential velocity
- surface friction velocity from wind stress
- surface buoyancy flux from heat/freshwater forcing

Then it calls:

```c++
KPPInstance->computeKPPMix(...)
```

### RK4 interaction

Current RK4 behavior is:

1. Disable stage KPP recompute while stepping RK sub-stages by setting
   `StageVerticalMixingEnabled = false`.
2. After RK4 stage accumulation and time-level update, recompute auxiliary
   state and call `computeStageVerticalMixing(...)` once on the fully updated
   state.
3. Restore previous stage-mixing flag.
4. Apply implicit vertical mixing.

This behavior is implemented in the RK4 stepper and is important for
consistency with current coupling expectations.

## Configuration Mapping

KPP reads configuration from the `VertMix: KPP` subgroup during `KPPMix::init`.
Important keys and class members:

- `Enable` -> `Enabled`
- `CriticalBulkRichardsonNumber` -> `CriticalRichardson`
- `StopOBLSearch` -> `StopOBLSearchMult`
- `SurfaceLayerExtent` -> `SurfaceLayerExtent`
- `MatchTechnique` -> `MatchTechniqueStr`
- `InterpType2` -> `InterpType2Str`
- `UseEnhancedDiffusion` -> `UseEnhancedDiffusion`
- `UseLangmuirCirculation` -> `UseLangmuirCirculation`
- `UseNonLocalFlux` -> `UseNonLocalFlux`
- `IceFractionThresholdForLangmuir` -> `IceFractionThresholdForLangmuir`
- `IceFractionThresholdForMinimumOBL` -> `IceFractionThresholdForMinimumOBL`
- `MinimumOBLUnderSeaIce` -> `MinimumOBLUnderSeaIce`
- `BackgroundViscosity` -> `BackgroundVisc`
- `BackgroundDiffusivity` -> `BackgroundDiff`
- `DebugDiagnostics` -> `DebugDiagnostics`

See [User KPP guide](../userGuide/KPPMix.md) for defaults and runnable examples.

## Output and Diagnostic Fields

`KPPMix::defineFields()` registers KPP outputs for I/O. Frequently used outputs:

- `BoundaryLayerDepth`
- `VertNonLocalFlux`
- `BulkRichardsonNumber`
- `BulkRichardsonShear`
- `UnresolvedShear`
- `BuoyancyJump`
- `TurbulentVelocityScale`
- `PotentialDensity`
- `SurfaceFrictionVelocity`
- `SurfaceBuoyancyFlux`

These can be enabled in output stream contents to diagnose OBL and profile
behavior in experiments.

## Developer Notes

- `MatchGradient` is currently treated as deprecated/unused and remapped to
  `SimpleShapes` at init.
- Unsupported `MatchTechnique` values are guarded and fall back to
  `SimpleShapes` with a log message.
- When `DebugDiagnostics` is enabled in debug builds, targeted diagnostic
  logging is available; behavior is compile/build-mode aware.

## Testing Strategy

### Code-level checks

1. Verify KPP initialization with explicit and default YAML keys.
2. Verify stage call path executes with KPP enabled and is skipped when
   disabled.
3. Verify RK4 sequencing: no stage recompute during sub-stages, one recompute
   before implicit vertical mixing.

### Diagnostics-based checks

1. Output `BoundaryLayerDepth`, `BulkRichardsonNumber`, and `VertDiff`.
2. Confirm OBL depth and coefficient evolution under changing forcing.
3. Validate optional outputs (`VertNonLocalFlux`) only when enabled.

### Regression checks

Run CTests and appropriate OMEGA regression workflows documented in:
[Testing Code](./Testing.md).
