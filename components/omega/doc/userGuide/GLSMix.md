(omega-glsmix)=

# Generic Length Scale (GLS) Mixing

This page explains how to enable and configure the Generic Length Scale
(GLS) turbulence closure, a config-selectable alternative to
[KPP](./KPPMix.md) for computing vertical viscosity and diffusivity.

Related pages:
- GLS design/theory: [Design GLS document](../design/GLSMix.md)
- GLS implementation details: [Developer GLS document](../devGuide/GLSMix.md)
- Broader vertical mixing options: [Vertical Mixing Coefficients](./VerticalMixingCoeff.md)

## When to use GLS instead of KPP

KPP diagnoses a boundary-layer depth and prescribes a profile shape within
it; this works well for the open-ocean mixed layer but is a poor fit for
shallow, strongly stratified coastal basins and lakes (e.g. the Great
Lakes), where the water column may be only a few tens of meters deep. GLS
instead solves two prognostic equations (turbulent kinetic energy and a
generic length-scale variable) and diagnoses viscosity/diffusivity
everywhere in the column from them, making it better suited to these
settings.

GLS and KPP are mutually exclusive: enabling both aborts at startup.

## Enabling GLS

In the YAML config, under `VertMix:`:

```yaml
VertMix:
  KPP:
    Enable: false
  GLS:
    Enable: true
    P: 3.0
    M: 1.5
    N: -1.0
    StabilityFunction: Canuto_A
```

The default `P: 3.0, M: 1.5, N: -1.0` with `StabilityFunction: Canuto_A` is
the k-epsilon closure, the recommended starting point. Other closure
families (k-omega, k-kl/Mellor-Yamada 2.5) are reachable through the same
`P`/`M`/`N` exponents; see the [design document](../design/GLSMix.md)
Section 3.2 for the exponent table. `StabilityFunction` may also be
`Canuto_B` or `Galperin`.

## Configuration options

| Option | Default | Meaning |
| --- | --- | --- |
| `Enable` | `false` | Enable GLS (mutually exclusive with KPP) |
| `P`, `M`, `N` | `3.0, 1.5, -1.0` | Closure-family exponents |
| `StabilityFunction` | `Canuto_A` | `Canuto_A`, `Canuto_B`, or `Galperin` |
| `Cmu0`, `C1`, `C2`, `C3Minus`, `C3Plus` | ROMS/REMORA defaults | Closure constants |
| `SigmaK`, `SigmaPsi` | `1.0`, `1.3` | Turbulent Schmidt numbers for TKE/psi diffusion |
| `Kmin`, `Pmin` | `7.6e-6`, `1.0e-12` | Floors on TKE/psi |
| `AkvBak`, `AktBak`, `AkkBak`, `AkpBak` | `5e-6`, `1e-6`, `5e-6`, `5e-6` | Background diffusivities |
| `Zos`, `Zob` | `0.02` m | Surface/bottom roughness lengths |
| `BottomDragCoeff` | `3.0e-3` | Quadratic drag coefficient used only for the TKE/psi bottom boundary condition |
| `LengthScaleCapEnable` | `true` | Cap the turbulence length scale (recommended on; see design doc Section 3.4) |
| `LengthScaleCapFraction` | `0.4` | Cap fraction of local water-column depth |

## Scope notes

This implementation solves the TKE/psi equations per water column (no
horizontal or vertical advection of TKE/psi), and does not include REMORA's
optional Craig-Banner wave breaking or Charnok roughness extras. See the
[design document](../design/GLSMix.md) Section 4.3 for the full list of
deliberate simplifications relative to REMORA/ROMS.
