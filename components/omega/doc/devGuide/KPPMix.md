(omega-dev-kppmix)=

# KPP Boundary Layer Mixing

This guide maps Omega KPP implementation details to runtime behavior and code
locations. It complements the design page by focusing on concrete APIs,
call flow, and developer test strategy.

## Implementation Overview

Omega KPP is implemented in `KPPMix` as a singleton with two major compute
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

### Time stepper interaction

KPP is hooked into all three OMEGA time steppers
(`src/timeStepping/RungeKutta4Stepper.cpp`,
`src/timeStepping/RungeKutta2Stepper.cpp`,
`src/timeStepping/ForwardBackwardStepper.cpp`) through two mechanisms:

1. **Stage recompute**: `Tendencies::StageVerticalMixingEnabled` (default
   `true`) gates a call to `computeStageVerticalMixing(...)` inside
   `computeAllTendencies`, `computeVelocityTendencies`, and
   `computeTracerTendencies`. Whichever of these tendency functions a stepper
   calls during its stages will trigger a KPP recompute on that stage's
   state.
2. **Post-step recompute**: `TimeStepper::applyPostStepVerticalMixing(...)`
   (in `src/timeStepping/TimeStepper.cpp`) is called by every stepper's
   `doStep()` immediately after `State->updateTimeLevels()`. It recomputes
   auxiliary state and calls `computeStageVerticalMixing(...)` once more on
   the fully updated state, then applies implicit vertical mixing via
   `VertMix::VertMixImplicit(...)` if enabled.

Per-stepper call flow:

- **RungeKutta4Stepper**: calls `computeAllTendencies(...)` once per stage
  (base stage plus 3 provisional stages), so KPP recomputes 4 times during
  stepping, followed by `applyPostStepVerticalMixing(..., "RK4")`.
- **RungeKutta2Stepper**: calls `computeAllTendencies(...)` twice (initial
  stage, midpoint stage), so KPP recomputes twice during stepping, followed
  by `applyPostStepVerticalMixing(..., "RK2")`.
- **ForwardBackwardStepper**: calls `computeVelocityTendencies(...)` and
  `computeTracerTendencies(...)` separately, each triggering a KPP recompute,
  followed by `applyPostStepVerticalMixing(..., "ForwardBackward")`.

In every stepper, the post-step recompute uses the fully updated state and is
what determines the KPP diagnostics written to output for that step.

Note: `RungeKutta4Stepper::doStep` still saves/restores
`StageVerticalMixingEnabled` around its stage loop and ANDs it with
`KPPMix::Enabled`. This is currently a no-op with respect to gating stage
recompute, since `computeStageVerticalMixing` already early-returns when KPP
is disabled; do not assume stage recompute is suppressed during RK4
sub-stages when reading that code.

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
- `UseNonLocalFlux` -> `UseNonLocalFlux` (expert/debugging override only --
  non-local flux is on by default and required for physically correct
  tracer transport; this key is intentionally omitted from the User Guide so
  it is not disabled by non-experts)
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
3. Verify per-stepper sequencing: stage recompute at each stage of the active
   stepper (4 for RK4, 2 for RK2, 2 for Forward-Backward), plus one final
   recompute on the updated state before implicit vertical mixing, for all
   three steppers.

### Diagnostics-based checks

1. Output `BoundaryLayerDepth`, `BulkRichardsonNumber`, and `VertDiff`.
2. Confirm OBL depth and coefficient evolution under changing forcing.
3. Validate optional outputs (`VertNonLocalFlux`) only when enabled.

### Regression checks

Run CTests and appropriate OMEGA regression workflows documented in:
[Testing Code](./Testing.md).
