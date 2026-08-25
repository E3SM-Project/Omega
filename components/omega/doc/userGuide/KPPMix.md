(omega-kppmix)=

# KPP Boundary Layer Mixing

This page explains how to enable, configure, and use Omega K-Profile
Parameterization (KPP) boundary layer mixing in runs.  The implementation follows
directly from the MPAS-Ocean implementation.

Related pages:
- KPP design/theory: [Design KPP document](../design/KPPMix.md)
- KPP implementation details: [Developer KPP document](../devGuide/KPPMix.md)
- Broader vertical mixing options: [Vertical Mixing Coefficients](./VerticalMixingCoeff.md)

## What KPP Provides

KPP computes:

- Ocean boundary layer depth (`BoundaryLayerDepth`)
- Vertical viscosity (`VertVisc`)
- Vertical diffusivity (`VertDiff`)
- Non-local tracer flux profile (`VertNonLocalFlux`)

It uses a bulk Richardson depth search followed by profile-based coefficient
construction. The non-local flux is a standard part of the KPP formulation and
is included by default.

## How KPP Is Used in Time Stepping

KPP is connected to all three OMEGA time steppers: Forward-Backward (the
default, split-explicit-style stepper), RungeKutta2, and RungeKutta4. For
whichever stepper is active, KPP recomputes boundary-layer depth and
coefficients at each internal stage of that stepper, and then once more on
the fully updated state after time levels are advanced, immediately before
implicit vertical mixing is applied:

- **Forward-Backward** (default): KPP recomputes when velocity tendencies
  are evaluated and again when tracer tendencies are evaluated, then once
  more on the updated state.
- **RungeKutta2**: KPP recomputes at the initial stage and at the midpoint
  stage, then once more on the updated state.
- **RungeKutta4**: KPP recomputes at each of the four RK4 stages, then once
  more on the updated state.

In all cases, the final recompute on the fully updated state is what
determines the KPP diagnostics written to output for that step.

## Configuration

KPP settings are under `VertMix: KPP` in `omega.yml`.

### Example

```yaml
VertMix:
  KPP:
    Enable: true
    CriticalBulkRichardsonNumber: 0.25
    StopOBLSearch: 1.0
    SurfaceLayerExtent: 0.1
    MatchTechnique: SimpleShapes
    InterpType2: LMD94
    UseEnhancedDiffusion: true
    UseBLDSmoothing: true
    UseLangmuirCirculation: true
    IceFractionThresholdForLangmuir: 0.05
    IceFractionThresholdForMinimumOBL: 0.15
    MinimumOBLUnderSeaIce: 5.0
    BackgroundViscosity: 1.0e-4
    BackgroundDiffusivity: 1.0e-5
    DebugDiagnostics: false
```

### Key Options

| Key | Meaning | Typical default |
|---|---|---|
| `Enable` | Enable KPP mixing | `true` |
| `CriticalBulkRichardsonNumber` | OBL depth criterion threshold | `0.25` |
| `StopOBLSearch` | Multiple of the critical Richardson number at which the OBL search stops descending | `1.0` |
| `SurfaceLayerExtent` | Surface layer thickness as a fraction of the OBL depth ($\epsilon$ in Large et al. 1994) | `0.1` |
| `MatchTechnique` | How the K profile meets interior mixing at the OBL base: `SimpleShapes` or `MatchBoth` | `SimpleShapes` |
| `InterpType2` | Interpolation type used near OBL matching/base logic | `LMD94` |
| `UseEnhancedDiffusion` | Enable enhanced diffusion treatment near OBL base | `true` |
| `UseBLDSmoothing` | Apply horizontal smoothing to the boundary layer depth | `true` |
| `UseLangmuirCirculation` | Apply Langmuir enhancement to the turbulent velocity scale | `true` |
| `IceFractionThresholdForLangmuir` | Above this ice fraction, disable Langmuir enhancement | `0.05` |
| `IceFractionThresholdForMinimumOBL` | Above this ice fraction, enforce minimum OBL depth | `0.15` |
| `MinimumOBLUnderSeaIce` | Minimum OBL depth under sea ice (m) | `5.0` |
| `BackgroundViscosity` | Background viscosity below the OBL (m^2/s) | `1.0e-4` |
| `BackgroundDiffusivity` | Background diffusivity below the OBL (m^2/s) | `1.0e-5` |
| `DebugDiagnostics` | Enable additional KPP diagnostics/logging in debug workflows | `false` |

Note that KPP reads its own `BackgroundViscosity` and `BackgroundDiffusivity`
from the `VertMix: KPP` group; these are separate from the `VertMix: Background`
values used by the other vertical mixing schemes.

## Output and Diagnostics

To diagnose KPP, include KPP fields in output stream contents. Common fields:

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

## Typical Workflow

1. Enable KPP and set baseline options in `omega.yml`.
2. Run a short case.
3. Inspect `BoundaryLayerDepth` and coefficient fields.
4. If needed, tune `CriticalBulkRichardsonNumber`, `MatchTechnique`, and
   `InterpType2`.
5. Re-run and compare diagnostics.

## Practical Notes

- Non-local tracer transport is required for physically correct KPP boundary
  layer tracer fluxes. It is applied through the
  `Tendencies: TracerNonLocalFluxTendencyEnable` flag, which is on by default
  and should be left enabled; see
  [Tendency Terms](./TendencyTerms.md).
- Use `DebugDiagnostics` sparingly for troubleshooting targeted cases.
- When studying sea-ice regions, review minimum-OBL and ice-threshold options.
