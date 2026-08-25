(omega-dev-kppmix)=

# KPP Boundary Layer Mixing

This guide maps Omega KPP implementation details to runtime behavior and code
locations. It complements the design page by focusing on concrete APIs,
call flow, and developer test strategy.

## Implementation Overview

The Omega implementation of KPP follows directly from the MPAS-Ocean implementation.
Notably, by default it does not match diffusivity and viscosity from interior mixing
(below the ocean surface boundary layer) sources at the base of the boundary
layer. Instead these separate sources are added directly to the KPP diagnosed
diffusivity and viscosity.  Matching can be enabled via the `MatchTechnique` parameter to
`MatchBoth` in the Omega yaml file.  Boundary layer depth is computed as the depth
where the bulk Richardson number exceeds a critical value.  It is then smoothed horizontally.
KPP diffusivity, viscosity, and non local fluxes are computed based on this boundary layer.

Omega KPP is implemented in `KPPMix` as a singleton with two major compute
phases:

1. OBL depth diagnosis (`computeOBLDepth`)
2. Coefficient/profile construction (`computeMixingCoefficients`)

Main class/API surface is in `src/ocn/KPPMix.h` and implementation is in
`src/ocn/KPPMix.cpp`.

### Depth convention

All KPP depths are measured downward from the free surface, not from the geoid.
`VertCoord::GeomZInterface` and `GeomZMid` are geometric heights relative to
`z = 0`, with `GeomZInterface(ICell, MinLayerCell)` equal to
`VertCoord::SshCell`. KPP therefore forms depths as
`SshCell(ICell) - GeomZ...(ICell, K)`. Layer thicknesses are differences of
geometric heights and are unaffected by the sea surface height.

## Notation for Readers New to KPP

KPP splits the water column at the ocean boundary layer (OBL) depth `h`, also
called the boundary layer depth (BLD). Inside the OBL, diffusivity and
viscosity are prescribed as a depth profile scaled by `h` and a turbulent
velocity scale; below it, only interior mixing applies.

| Symbol | Code name | Units | Meaning |
| --- | --- | --- | --- |
| `h` | `HOBL`, `BoundaryLayerDepth` | m | OBL depth below the free surface |
| `d` | `ZDepth`, `ZCenter` | m | depth below the free surface |
| `sigma` | `Sigma` | - | normalized depth in `[-1,0]`, Omega sign convention: 0 at the surface, -1 at the OBL base |
| `sigma_mu` | `SigmaMu` | - | `-Sigma`, in `[0,1]`; the CVMix/Large et al. convention |
| `u*` | `UStar` | m/s | surface friction velocity from wind stress |
| `B_0` | `BuoyFlux` | m^2/s^3 | surface buoyancy flux; **negative is destabilizing** (convection) |
| `L` | `LMoninObukhov` | m | Monin-Obukhov length, `u*^3 / (kappa B_0)` |
| `zeta` | `Zeta` | - | stability coordinate `d/L`; negative is unstable |
| `kappa` | `VonKar` | - | von Karman constant |
| `epsilon` | `SurfaceLayerExtent` | - | surface layer as a fraction of `h` |
| `w_m`, `w_s` | `WMTurb`, `WSTurb` | m/s | turbulent velocity scales for momentum and scalars |
| `Vt^2` | `Vt2`, `UnresolvedShear` | m^2/s^2 | unresolved turbulent shear, Large et al. (1994) Eq. 23 |
| `Ri_b` | `RiBulk`, `BulkRichardsonNumber` | - | bulk Richardson number, Eq. 21 |
| `G(sigma)` | `kppShape*` | - | non-dimensional profile shapes |
| `gamma` | `VertNonLocalFlux` | - | non-local (counter-gradient) tracer flux coefficient |

The OBL depth is the shallowest `d` at which `Ri_b(d)` reaches the critical
value; the search is done per cell in `computeOBLDepth`, and the crossing depth
is refined by a quadratic fit through the three nearest cell-center `Ri_b`
values.

### Shape and stability functions

The non-dimensional functions live in `src/ocn/KPPConstants.h` in namespace
`OMEGA::KPP`. All shape functions take `Sigma` in `[-1,0]` and convert to
`SigmaMu` internally, so callers never flip signs:

- `kppShapeMomentum`, `kppShapeScalar` -- gradient shapes for viscosity and
  diffusivity (`SimpleShapes`)
- `kppShapeMatched` -- gradient shape that additionally reaches a prescribed
  value at the OBL base, used by `MatchBoth`
- the non-local flux has no shape function of its own; it reuses whichever
  scalar shape is in effect, scaled by `C_s` instead of `h*w_s`
- `kppPhiInvMomentum`, `kppPhiInvScalar` -- **inverse** Monin-Obukhov stability
  functions `phi^-1(zeta)`; they already return the reciprocal, so do not
  invert them again at the call site

### Constants and defaults

`src/ocn/KPPConstants.h` is the single authoritative source for KPP default
values. The runtime-configurable members of `KPPMix` (`CriticalRichardson`,
`StopOBLSearchMult`, `SurfaceLayerExtent`, the two ice-fraction thresholds and
`MinimumOBLUnderSeaIce`) are initialized from those constants rather than from
inline literals, so a default is changed in exactly one place.

Per-thread edge scratch arrays in `computeOBLDepth` are sized from
`HorzMesh::MaxEdgesBound`, the shared compile-time bound on edges per cell;
KPP does not define its own maximum.

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
- `MatchTechnique` -> `MatchTechnique` (a `KPPMatchType` enum, not a string)
- `InterpType2` -> `InterpType2Str`
- `UseEnhancedDiffusion` -> `UseEnhancedDiffusion`
- `UseLangmuirCirculation` -> `UseLangmuirCirculation`
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

- `MatchTechnique` accepts only `SimpleShapes` and `MatchBoth`. Any other value
  aborts at init rather than falling back, so a typo cannot silently change the
  scheme being run.
- `MatchBoth` matches the interior coefficient *value* at the OBL base. The
  shape derivative there is zero, so the gradient is not yet matched despite the
  name.
- `MatchBoth` needs interior coefficients to be passed in; without them
  `ShapeAtBase` is zero and it degenerates exactly to `SimpleShapes`.
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
