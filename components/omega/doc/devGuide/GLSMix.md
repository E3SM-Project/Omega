(omega-dev-glsmix)=

# Generic Length Scale (GLS) Mixing

This guide maps the Omega GLS implementation to runtime behavior and code
locations. It complements the [design page](../design/GLSMix.md) by
focusing on concrete APIs, call flow, and developer test strategy.

## Implementation Overview

Omega GLS is implemented in `GLSMix` (`src/ocn/GLSMix.h`/`.cpp`) as a
singleton, structured like `KPPMix`. Stability-function coefficient tables
and the exponent-derived constants live in `src/ocn/GLSConstants.h`
(`namespace OMEGA::GLS`-equivalent free functions
`computeGLSDerivedConstants`/`computeStabilityFunctions`/
`computeCanutoGmBound`).

Unlike KPP (two-stage: boundary-layer search, then profile construction),
GLS has a single per-timestep entry point, `GLSMix::advance()`, performing:

1. Shear-squared at layer interfaces (`GLSShearSquaredOnCell`, the same
   edge-kite-weighted average as `GradRichardsonNum` in `VertMix.h`, but
   returning the raw shear rather than a Richardson-number ratio).
2. Shear/buoyancy production terms, from the *previous* step's
   `VertVisc`/`VertDiff` (i.e. lagged by one step, same as REMORA).
3. Implicit tridiagonal solve for TKE, then for psi, each via
   `TriDiagDiffSolver` (`src/base/TriDiagSolvers.h`).
4. Length-scale limiter, stability functions, and the final
   `VertVisc`/`VertDiff`/`Akk`/`Akp`/`Lscale` update.

## Call flow

`GLSMix::init()` is called from `OceanInit.cpp`, alongside
`VertMix::init()`/`KPPMix::init()`. Each step,
`VertMix::VertMixImplicit()` (`src/ocn/VertMix.cpp`) calls
`GLSMix::computeGLSFields()` immediately after computing this step's
Brunt-Vaisala frequency and before `computeVertMix()` merges the result —
this is the single call site regardless of which time stepper is active
(`VertMixImplicit()` is itself the one function all time steppers call for
implicit vertical mixing). `computeGLSFields()` fetches
`ZonalStressCell`/`MeridStressCell` from the default `Forcing` instance and
calls `advance()`, which takes those as plain arrays so it can be
unit-tested without a `Forcing` setup (see `GLS_LENGTHCAP_TEST` in
`test/ocn/GLSMixTest.cpp`).

`VertMix::computeVertMix()` branches on `GLSMix::getInstance()->Enabled`
near the top (before the Background/Shear computation): when enabled, GLS's
`VertDiff`/`VertVisc` are copied directly into the output, convective mixing
is optionally layered on top, and the function returns early, skipping the
Background/Shear/KPP merge path entirely.

## Grid convention

TKE, psi, `Lscale`, `Akk`, `Akp`, `VertVisc`, `VertDiff` are all
`Array2DReal`, `[NCells][NVertLayersP1]`, layer-interface-centered — the
same grid and `K`-increases-downward convention as
`VertMix::VertDiff`/`VertVisc` (`K = MinLayerCell(ICell)` is the surface
interface, `K = MaxLayerCell(ICell) + 1` is the bottom interface). This is
the *opposite* vertical direction from REMORA/ROMS, where `k` increases
upward from the bottom (`k=0`) to the surface (`k=N+1`); boundary-condition
formulas were transposed accordingly when porting from `REMORA_gls.cpp`.

## Tridiagonal solver boundary-condition convention

`TriDiagDiffSolver`'s system is `-G(i)x_{i-1} + (H(i)+G(i)+G(i+1))x_i -
G(i+1)x_{i+1} = X(i)`, and the solver internally folds `G(i)` into row `i`'s
own diagonal (`H(0) += G(0)` and similar in the elimination). This means a
true Dirichlet row **must** have `G = 0` at its own index — any nonzero `G`
there would blend the fixed value with the interior solution rather than
preserving it exactly. Consequently, the diffusive coupling from a boundary
value into the adjacent interior row cannot be expressed through that row's
own stored `G`; `GLSMix::advance()` instead injects it directly into the
first/last interior row's own `H`/`X` (see the comments at the `K + 1 <
KMax` / `K - 1 == KMin` branches in the TKE and psi solve kernels). This is
different from REMORA, which uses separate explicit flux terms
(`tke_fluxt`/`gls_fluxb`, etc.) for the same purpose; both are valid
numerical treatments of the same boundary condition.

## Where to find constants

- Closure exponent defaults and per-stability-function `Gh0`/`Ghcri`/`L1-L8`
  (Canuto) or quasi-equilibrium (Galperin) coefficients:
  `GLSConstants.h::computeGLSDerivedConstants()`.
- Config defaults: `configs/Default.yml`, `VertMix:GLS:`.
- Runtime overrides and validation (including the GLS/KPP mutual-exclusion
  check): `GLSMix::init()`.

## Developer test strategy

`test/ocn/GLSMixTest.cpp`, registered in `test/CMakeLists.txt`, following
the same `add_omega_test`/`add_omega_ctest` shared-executable pattern as
`KPPMixTest.cpp`:
- `stability`: pure function checks on `computeStabilityFunctions()` for
  each stability-function choice (no mesh required).
- `lengthcap`: full-column regression test of the length-scale-cap bug fix,
  calling `GLSMix::advance()` directly with synthetic inputs.
- `config-stability`, `config-exclusive`: expected-abort checks on
  `GLSMix::init()` (`WILL_FAIL` ctest properties), mirroring KPP's
  `config-gradient`/`config-unsupported`/`config-parabolic` groups.
