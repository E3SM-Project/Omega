(omega-design-glsmix)=
# Generic Length Scale (GLS) Turbulence Closure

**Table of Contents**
1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Algorithmic Formulation](#3-algorithmic-formulation)
4. [Design](#4-design)
5. [Verification and Testing](#5-verification-and-testing)

## 1 Overview

This document describes the Omega implementation of the Generic Length Scale
(GLS) two-equation turbulence closure (Umlauf and Burchard, 2003, *A generic
length-scale equation for geophysical turbulence models*, J. Marine Research
61, 235-265), as a config-selectable alternative to
[KPP](omega-design-kppmix). GLS solves prognostic equations for turbulent
kinetic energy (TKE, `k`) and a generic length-scale variable (`psi`), and
diagnoses vertical viscosity and diffusivity from them everywhere in the
water column, rather than diagnosing a boundary-layer depth and prescribing
a profile shape as KPP does. This makes it better suited than KPP to
shallow, strongly stratified basins such as the Great Lakes, where KPP's
open-ocean boundary-layer assumptions are a poor fit.

The implementation follows the ROMS/REMORA `gls_corstep.F`/`REMORA_gls.cpp`
formulation for its equations, stability functions, and constants, and
carries over a REMORA/ROMS bug fix for the turbulence length-scale limiter
(Section 3.4) from the start, rather than reproducing the original bug.

## 2 Requirements

### 2.1 Requirement: config-selectable alternative to KPP

GLS and KPP are both schemes for computing the full `VertDiff`/`VertVisc`
profile; a run selects one or the other via `VertMix:GLS:Enable` /
`VertMix:KPP:Enable`. Enabling both is a configuration error.

### 2.2 Requirement: generalized closure family

The closure exponents `p`, `m`, `n` (defining `psi = cmu0^p k^m l^n`) and the
stability-function choice (Canuto A, Canuto B, Galperin) must be
independently configurable, so that k-epsilon (the default, `p=3, m=1.5,
n=-1`), k-omega, k-kl/Mellor-Yamada 2.5, and other GLS-family variants are
all reachable without code changes.

### 2.3 Requirement: bounded length scale under all stratification

The diagnosed turbulence length scale (and hence viscosity/diffusivity) must
remain bounded regardless of the sign of stratification. See Section 3.4.

### 2.4 Desired: restart-compatible prognostic state

TKE and psi are prognostic; a restarted run must reproduce the same
evolution as an uninterrupted run, so they are registered in the `Restart`
field group.

## 3 Algorithmic Formulation

### 3.1 Governing equations

$$
\frac{Dk}{Dt} = \frac{\partial}{\partial z}\left(\frac{A_k}{\sigma_k}\frac{\partial k}{\partial z}\right) + P_s + P_b - \varepsilon
$$

$$
\frac{D\psi}{Dt} = \frac{\partial}{\partial z}\left(\frac{A_k}{\sigma_\psi}\frac{\partial \psi}{\partial z}\right) + \frac{\psi}{k}\left(c_1 P_s + c_3 P_b - c_2 \varepsilon F_{wall}\right)
$$

with shear production $P_s = A_v (\partial u/\partial z)^2$, buoyancy
production $P_b = -A_t N^2$, dissipation
$\varepsilon = c_{\mu 0}^{3+p/n} k^{3/2+m/n} \psi^{-1/n}$, and
$\psi = c_{\mu 0}^p k^m l^n$ relating `psi` to the turbulence length scale
`l`. $F_{wall}$ (the Mellor-Yamada wall-proximity correction) is not
implemented; see Section 4.4.

### 3.2 Closure family (`p`, `m`, `n`)

| Closure | p | m | n |
| --- | --- | --- | --- |
| k-epsilon (default) | 3.0 | 1.5 | -1.0 |
| k-omega | -2.0 | 0.5 | -1.0 |
| k-kl / Mellor-Yamada 2.5 | 0.0 | 1.0 | 1.0 |

### 3.3 Stability functions

`Sm` (momentum) and `Sh` (tracers) relate the non-dimensional buoyancy
number `Gh = -N^2 l^2 / 2k` and shear number `Gm = shear^2 l^2 / 2k` to the
mixing coefficients. Canuto A/B use the quartic rational-function form of
Canuto et al. (2001); Galperin uses the simpler Mellor-Yamada
quasi-equilibrium form. Coefficient tables are in `GLSConstants.h`.

### 3.4 Length-scale limiter (bug fix)

The unlimited length scale `l_unlmt = psi^{1/n} cmu0^{-p/n} k^{-m/n}` is
bounded two ways:

1. **Galperin bound**, applied only under stable stratification
   (`N^2 > 0`): `l <= sqrt(0.56 k / N^2)`.
2. **Unconditional hard cap**: `l <= LengthScaleCapFraction * D`, where `D`
   is the local water-column depth, applied regardless of the sign of
   `N^2`.

The original REMORA implementation applied the hard cap only inside the
`N^2 > 0` branch (i.e., only where the Galperin bound already applied),
leaving the length scale unbounded under neutral or unstable stratification.
In one documented case this drove the diagnosed viscosity to roughly
`4e6 m^2/s` in a 23 m water column, homogenizing a 12 m layer in a single
time step. ROMS fixed this upstream (`gls_corstep.F`, January 2025) by
applying the cap unconditionally; the Omega implementation applies it
unconditionally from the start (`LengthScaleCapEnable`, default on).

## 4 Design

### 4.1 Data types and parameters

`GLSMix` (`src/ocn/GLSMix.h`/`.cpp`) is a singleton, structured like
`KPPMix`, reading a `GLS:` subgroup under `VertMix:` in
`configs/Default.yml`. Key parameters: `Enable`, `P`/`M`/`N`,
`StabilityFunction` (`Canuto_A`/`Canuto_B`/`Galperin`), `Cmu0`, `C1`, `C2`,
`C3Minus`/`C3Plus`, `SigmaK`/`SigmaPsi`, `Kmin`/`Pmin` (floors), background
diffusivities `AkvBak`/`AktBak`/`AkkBak`/`AkpBak`, roughness lengths
`Zos`/`Zob`, `BottomDragCoeff` (see 4.3), and
`LengthScaleCapEnable`/`LengthScaleCapFraction`.

Prognostic/diagnostic state (`Array2DReal`, `[NCells][NVertLayersP1]`,
layer-interface-centered, same grid as `VertMix::VertDiff`): `TKE`, `Psi`,
`Lscale`, `Akk`, `Akp` (GLS's own vertical-diffusion coefficients for TKE
and psi), and the public-facing `VertVisc`/`VertDiff`.

### 4.2 Methods

`GLSMix::init()` reads configuration and aborts if `StabilityFunction` is
unrecognized, or if both GLS and KPP are enabled.
`GLSMix::advance(NormalVelocity, TangentialVelocity, BruntVaisalaFreqSq,
PseudoThickCell, KineticEnergyCell, ZonalStressCell, MeridStressCell, DT)`
is the per-timestep physics kernel (unit-testable without `Forcing`).
`GLSMix::computeGLSFields(...)` is the convenience wrapper that fetches
`ZonalStressCell`/`MeridStressCell` from the default `Forcing` instance and
calls `advance()`; this is what `VertMix::VertMixImplicit()` calls once per
step, immediately before `computeVertMix()` merges the result (see
`src/ocn/VertMix.cpp`).

### 4.3 Scope simplifications relative to REMORA/ROMS

To keep the initial port tractable, three simplifications were made
deliberately (documented here rather than silently):

- **No horizontal or vertical advection of TKE/psi.** Both equations are
  solved per water column (vertical diffusion, production, and dissipation
  only), matching common practice in other ocean models' TKE/GLS closures
  (e.g. GOTM, CVMix-TKE). REMORA/ROMS do advect TKE/psi; this is a candidate
  future enhancement once the column physics is validated.
- **Dirichlet boundary conditions couple diffusively into the interior.**
  REMORA's tridiagonal solve treats the surface/bottom values as fixed and
  injects their influence via separate explicit flux terms
  (`tke_fluxt`/`gls_fluxb`, etc.), decoupled from the matrix. Omega's
  `TriDiagDiffSolver` convention instead requires a zero coupling
  coefficient (`G`) at a true Dirichlet row (otherwise the solver's own
  elimination would corrupt the fixed value by blending it with the
  interior solution); the equivalent diffusive flux is injected directly
  into the first/last interior row's own `H`/`X` instead. This is a
  standard, well-posed numerical treatment, but is not a line-for-line match
  to REMORA's bookkeeping.
- **Approximate bottom boundary condition.** The surface TKE/psi boundary
  condition uses the real wind stress (`Forcing::SfcStressForcing`, as KPP
  does). The bottom boundary condition estimates a near-bottom friction
  velocity from a configurable quadratic drag coefficient
  (`BottomDragCoeff`) and the model's near-bottom kinetic energy, rather
  than the actual bottom-stress diagnostic used by the momentum solver
  (which in Omega is only populated when `VertMix`'s own implicit bottom
  drag is enabled).

Not ported at all (REMORA optional extras, off by default there too):
Craig-Banner surface wave-breaking TKE injection, Charnok surface roughness,
and the Mellor-Yamada wall-proximity correction to psi dissipation
($F_{wall}$ in Section 3.1).

### 4.4 Integration with VertMix

`VertMix::computeVertMix()` treats GLS as a full-column replacement for the
Background/Shear/KPP path (GLS's own production/dissipation balance already
folds in shear and buoyancy effects, and the config validation in
`GLSMix::init()` guarantees KPP is not also enabled), copying
`GLSMix::VertDiff`/`VertVisc` directly into the merged output. Convective
mixing may still be layered on top, as it is for KPP.

## 5 Verification and Testing

### 5.1 Stability-function sanity (`GLS_STABILITY_TEST`)

For each stability-function choice, `Sm`/`Sh` are finite and non-negative
over a representative range of the buoyancy number `Gh`.
Tests requirement 2.2.

### 5.2 Length-scale cap regression (`GLS_LENGTHCAP_TEST`)

A shallow (23-layer, 1 m/layer), strongly sheared, neutral/unstable column
(`N^2 <= 0`, where the Galperin bound does not engage) is run for several
steps; `Lscale` must stay within `LengthScaleCapFraction * Depth` and
`VertVisc`/`VertDiff` must stay well below a generous sanity bound, in
contrast to the unbounded growth documented for the pre-fix REMORA code.
Tests requirement 2.3.

### 5.3 Config validation (`GLS_CONFIG_STABILITY_TEST`, `GLS_CONFIG_EXCLUSIVE_TEST`)

`GLSMix::init()` must abort on an unrecognized `StabilityFunction` string,
and on a config with both `GLS:Enable` and `KPP:Enable` true.
Tests requirements 2.1 and 2.2.
