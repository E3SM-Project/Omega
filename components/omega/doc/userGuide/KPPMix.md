(omega-kppmix)=

# KPP Boundary Layer Mixing

This page explains how to enable, configure, and use OMEGA K Profile
Parameterization (KPP) boundary layer mixing in runs.

Related pages:
- KPP design/theory: [Design KPP document](../design/KPPMix.md)
- KPP implementation details: [Developer KPP document](../devGuide/KPPMix.md)
- Broader vertical mixing options: [Vertical Mixing Coefficients](./VerticalMixingCoeff.md)

## What KPP Provides

KPP computes:

- Ocean boundary layer depth (`BoundaryLayerDepth`)
- Vertical viscosity (`VertVisc`)
- Vertical diffusivity (`VertDiff`)
- Optional non-local tracer flux profile (`VertNonLocalFlux`)

It uses a bulk Richardson depth search followed by profile-based coefficient
construction.

## How KPP Is Used in Time Stepping

In the current RK4 workflow, stage-level KPP recomputation is gated off during
RK4 sub-stages and KPP is recomputed once after all four stages, on the fully
updated state, before implicit vertical mixing is applied.

This means KPP diagnostics in output correspond to the post-stage state for
each RK4 step.

## Configuration

KPP settings are under `VertMix: KPP` in `omega.yml`.

### Example

```yaml
VertMix:
  KPP:
    Enable: true
    UseNonLocalFlux: true
    CriticalBulkRichardsonNumber: 0.25
    MatchTechnique: SimpleShapes
    InterpType2: LMD94
    UseEnhancedDiffusion: true
    IceFractionThresholdForLangmuir: 0.05
    IceFractionThresholdForMinimumOBL: 0.15
    MinimumOBLUnderSeaIce: 5.0
    DebugDiagnostics: false
```

### Key Options

| Key | Meaning | Typical default |
|---|---|---|
| `Enable` | Enable KPP mixing | `true` |
| `UseNonLocalFlux` | Enable non-local tracer flux profile | `true` |
| `CriticalBulkRichardsonNumber` | OBL depth criterion threshold | `0.25` |
| `MatchTechnique` | KPP profile matching mode | `SimpleShapes` |
| `InterpType2` | Interpolation type used near OBL matching/base logic | `LMD94` |
| `UseEnhancedDiffusion` | Enable enhanced diffusion treatment near OBL base | `true` |
| `IceFractionThresholdForLangmuir` | Above this ice fraction, disable Langmuir enhancement | `0.05` |
| `IceFractionThresholdForMinimumOBL` | Above this ice fraction, enforce minimum OBL depth | `0.15` |
| `MinimumOBLUnderSeaIce` | Minimum OBL depth under sea ice (m) | `5.0` |
| `DebugDiagnostics` | Enable additional KPP diagnostics/logging in debug workflows | `false` |

KPP also uses background coefficients from:

```yaml
VertMix:
  Background:
    Viscosity: 1.0e-4
    Diffusivity: 1.0e-5
```

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

- Keep `UseNonLocalFlux` enabled when you want tracer non-local transport.
- Use `DebugDiagnostics` sparingly for troubleshooting targeted cases.
- When studying sea-ice regions, review minimum-OBL and ice-threshold options.
