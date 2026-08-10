(omega-design-kppmix)=
# KPP Boundary Layer Mixing

**Table of Contents**
1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Algorithmic Formulation](#3-algorithmic-formulation)
4. [Design](#4-design)
5. [Verification and Testing](#5-verification-and-testing)

## 1 Overview

This document describes the OMEGA implementation of K Profile Parameterization
(KPP) ocean boundary layer mixing. KPP computes boundary-layer depth, vertical
viscosity, vertical diffusivity, and an optional non-local tracer flux shape
used by tracer tendencies.

The implementation is in `KPPMix` and is integrated with the OMEGA tendency and
RK4 stepping workflow. Relative to broad vertical mixing documentation, this
page focuses specifically on KPP theory, algorithm choices, and verification.

Related pages:
- User usage/configuration: [KPP in the User Guide](../userGuide/KPPMix.md)
- Developer implementation details: [KPP in the Developer Guide](../devGuide/KPPMix.md)
- Broader vertical mixing context: [Vertical Mixing Coefficients](./VerticalMixingCoeff.md)

## 2 Requirements

### 2.1 Requirement: Boundary-layer depth from bulk Richardson criterion

The OBL depth must be diagnosed from a bulk Richardson criterion so that
mixing depth responds to evolving stratification, shear, and surface forcing.

### 2.2 Requirement: Coefficients must be computable in parallel over columns

The KPP implementation must operate over many columns in parallel using OMEGA
array/kernels, rather than serial single-column calls, to match accelerator
performance goals.

### 2.3 Requirement: Compatible with additive vertical-mixing framework

KPP viscosity/diffusivity fields must be compatible with existing OMEGA vertical
mixing infrastructure so they can be merged with other configured contributions.

### 2.4 Desired: Optional non-local flux and profile matching controls

KPP should support optional non-local tracer flux profiles and configurable
matching/interpolation choices to support scientific tuning studies.

### 2.5 Desired: Stable RK4 interaction

For RK4, KPP should be computed in a way that avoids repeated stage re-evaluation
when configuration requires a single post-stage update on the fully updated
state.

## 3 Algorithmic Formulation

The implementation follows a two-stage KPP structure.

### 3.1 Stage 1: OBL depth search

For each water column, OBL depth $h$ is diagnosed by searching downward until
bulk Richardson number reaches a critical value:

$$
Ri_b(z) = \frac{\Delta b(z)\, z}{|\Delta \mathbf{U}(z)|^2 + V_t^2(z)}
$$

with threshold

$$
Ri_b(h) = Ri_{crit}.
$$

Here, $\Delta b$ is buoyancy jump relative to the near-surface reference,
$|\Delta \mathbf{U}|^2$ is shear contribution, and $V_t^2$ is unresolved shear.
The code supports interpolation/matching choices near the crossing and applies
configured constraints such as minimum OBL under sea ice and a maximum by water
column depth.

### 3.2 Stage 2: KPP coefficients and optional non-local flux

Given diagnosed $h$, KPP computes interface coefficients using shape functions
in normalized depth $\sigma = -z/h$:

$$
K_m(\sigma) = h\, w_m(\sigma)\, M_1(\sigma),
$$

$$
K_s(\sigma) = h\, w_s(\sigma)\, S_1(\sigma),
$$

where $w_m$ and $w_s$ are turbulent velocity scales from Monin-Obukhov style
stability functions. Optional non-local tracer flux shape $G(\sigma)$ is
computed when enabled.

Below OBL, coefficients revert to configured background values, with optional
enhanced diffusion handling near the OBL base.

## 4 Design

### 4.1 Data types and parameters

#### 4.1.1 Parameters

KPP is configured from the `VertMix: KPP` YAML group. Key parameters include:

- `Enable`
- `UseNonLocalFlux`
- `CriticalBulkRichardsonNumber`
- `MatchTechnique`
- `InterpType2`
- `UseEnhancedDiffusion`
- `IceFractionThresholdForLangmuir`
- `IceFractionThresholdForMinimumOBL`
- `MinimumOBLUnderSeaIce`
- `BackgroundViscosity`
- `BackgroundDiffusivity`
- `DebugDiagnostics`

Defaults and usage examples are documented in the user guide page:
[KPP in the User Guide](../userGuide/KPPMix.md).

#### 4.1.2 Class/data structure

`KPPMix` is a singleton that owns persistent output fields, including:

- `BoundaryLayerDepth`, `IndexBoundaryLayerDepth`
- `VertDiff`, `VertVisc`
- `VertNonLocalFlux`
- diagnostics such as `BulkRichardsonNumber`, `BulkRichardsonShear`,
  `UnresolvedShear`, `BuoyancyJump`, and `TurbulentVelocityScale`

### 4.2 Methods

Main interface:

```c++
void computeKPPMix(const Array2DReal &PotentialDensity,
                   const Array2DReal &NormalVelocity,
                   const Array2DReal &TangentialVelocity,
                   const Array1DReal &SurfaceFrictionVelocity,
                   const Array1DReal &SurfaceBuoyancyFlux,
                   const Array2DReal &BruntVaisalaFreqSq,
                   const Array1DReal &IceFraction,
                   const Array1DReal &WindSpeed10m = Array1DReal());
```

Internal stages:
- `computeOBLDepth(...)`
- `computeMixingCoefficients(...)`

### 4.3 RK4 coupling behavior

In OMEGA RK4 stepping, stage-level KPP recomputation can be gated off and KPP
is recomputed once after all RK4 stages on the fully updated state before
implicit vertical mixing is applied. This behavior is part of the current
coupling design and is described in detail for developers and users in:

- [Developer KPP workflow](../devGuide/KPPMix.md)
- [User runtime notes](../userGuide/KPPMix.md)

## 5 Verification and Testing

### 5.1 Unit-level checks

Use targeted tests and diagnostics to verify:

- OBL depth search monotonicity and threshold crossing behavior
- Positive bounded coefficients and expected background behavior below OBL
- Correct enable/disable behavior for non-local flux and enhanced diffusion

Tests cover requirements: 2.1, 2.2, 2.3, 2.4.

### 5.2 Coupled/regression checks

Run regression cases and compare key diagnostics over time:

- `BoundaryLayerDepth`
- `BulkRichardsonNumber`
- `VertDiff`, `VertVisc`
- `VertNonLocalFlux` (when enabled)

For full OMEGA testing workflow, see the developer testing guide:
[Testing Code](../devGuide/Testing.md).

### 5.3 Configuration sensitivity checks

Perform short experiments varying:

- `CriticalBulkRichardsonNumber`
- `MatchTechnique`
- `InterpType2`
- `UseEnhancedDiffusion`
- sea-ice thresholds

to ensure expected qualitative and quantitative responses in OBL depth and
mixing intensity.
