(omega-design-pressure-grad-high-order)=
# Higher-Order Horizontal Pressure Gradient

**Table of Contents**
1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Algorithmic Formulation](#3-algorithmic-formulation)
4. [Design](#4-design)
5. [Verification and Testing](#5-verification-and-testing)

## 1 Overview

This document specifies a higher-order discretization of the horizontal pressure
gradient force (PGF) in Omega's momentum equation. It is the "future design document"
promised in {ref}`Pressure Gradient <omega-design-pressure-grad>` (§2.3, §3.2),
which commits Omega to a high-order option "similar to
[Adcroft et al. 2008](https://doi.org/10.1016/j.ocemod.2008.02.001)" but defers the
details. The existing centered, second-order Montgomery-potential scheme
(`PressureGradCentered`) remains the default and the reference implementation;
this design adds the high-order option that the `PressureGrad` class and its
`PressureGradType` enum already anticipate (the placeholder high-order entry, named
`FiniteVolume` here; see §4.1.1).

The PGF is the single most error-sensitive term in a layered ocean model. In the
non-Boussinesq, hydrostatic momentum equation
({ref}`layered momentum equation <omega-v1-momentum-eq>`) the horizontal acceleration
from pressure and gravity is

$$
-\left(\alpha \nabla p + \nabla \Phi \right)_{e,k},
$$

in which $\alpha\nabla p$ and $\nabla\Phi = g\nabla z + \nabla(\phi_{TP}+\phi_{SAL})$
are individually large and nearly cancel. The residual is the dynamically relevant
baroclinic signal. When layers are tilted in geometric height — as they always are in
a general ALE coordinate, and severely so for thin, steeply sloped layers near shelf
breaks and (future) ice-shelf cavities — truncation error in this cancellation appears
directly as spurious velocity (the classic "pressure gradient error"). Experience with
MPAS-Ocean and published results for other models show that an accurate PGF goes a long
way toward making thin, steeply sloped layers usable.

This design is shaped by three competing principles:

1. **Accuracy.** The PGF must be markedly more accurate than the centered scheme,
   especially in the presence of sloping layers.
2. **Cost.** Omega uses TEOS-10 (a 75-term polynomial) for specific volume $\alpha$.
   Naively evaluating $\alpha$ at every quadrature point would dominate run time. The
   design must bound TEOS-10 evaluations.
3. **General ALE coordinate.** The scheme must remain accurate and stable for thin,
   steeply sloped layers. The PGF alone cannot guarantee this, but it is a necessary
   ingredient.

A guiding choice that follows from principles 2 and 3 is that the relevant metric of
success is **low absolute PGF error at affordable (coarse) resolution**, not high
asymptotic order of convergence. Omega will be strongly resource-limited in how fine a
mesh it can run, so a scheme whose error is small at coarse resolution is preferred over
one with a steeper convergence slope but larger coarse-resolution error. Order of
accuracy is configurable and a fourth-order variant is targeted, but it is a means to the
end of coarse-resolution accuracy, not the objective itself.

## 2 Requirements

### 2.1 Requirement: Higher accuracy than the centered scheme at affordable resolution

The high-order PGF must produce substantially lower absolute error than
`PressureGradCentered` at the coarse-to-moderate resolutions Omega can afford to run,
for representative stratified columns with horizontal gradients of temperature, salinity,
surface pressure, and coordinate slope. Convergence rate is a secondary diagnostic, not
the primary acceptance criterion.

### 2.2 Requirement: Bounded TEOS-10 cost

The number of TEOS-10 specific-volume evaluations performed by the PGF must be
independent of the quadrature order and bounded at approximately one evaluation (plus its
first derivatives) per cell per layer per time step — comparable to the cost Omega already
pays to compute the `Eos::SpecVol` field. The scheme must not require evaluating the full
TEOS-10 polynomial at each quadrature point, nor integrating the TEOS-10 polynomial
itself semi-analytically.

### 2.3 Requirement: Robustness for thin, steeply sloped layers via discrete hydrostatic consistency

The scheme must suppress spurious pressure-gradient accelerations when layers are thin and
steeply sloped in geometric height. The primary guarantee is **discrete hydrostatic
consistency**: for a physically resting ocean with horizontally uniform conservative
temperature and absolute salinity but coordinate surfaces tilted arbitrarily, the scheme
must return **exactly zero** horizontal PGF, to machine precision, regardless of the vertical
profile or the steepness of the tilt. This property is *structural* — it depends only on the
construction of the discretization, not on any background reference state — and therefore
holds uniformly across the global ocean. The residual error for non-uniform $T$/$S$ is then
genuine truncation error, reduced by the high-order, mean-preserving reconstruction (§3.4)
and by evaluating the equation of state consistently, including its pressure dependence
(compressibility $\alpha_p$), so that thermobaric effects under tilted layers are not aliased
into spurious flow.

### 2.4 Requirement: Consistency with the rest of the model state

The reconstructions used by the PGF must be **mean-preserving** (reproduce each layer's
prognostic mean $\Theta_{i,k}$, $S_{i,k}$ exactly) and **anchored to `VertCoord`**
(use the model's diagnostic interface pressures and geometric heights). The PGF must
therefore see the same layer-mean state the tracer and thickness equations evolve, and the
same hydrostatic pressure the rest of the model uses. It may assume a smoother sub-layer
profile for integration, but it must not introduce a second, inconsistent representation
of $\Theta$, $S$, $p$, or $z$.

### 2.5 Requirement: Runtime-selectable, backward-compatible option

The high-order scheme must be selectable at runtime through the existing `PressureGrad`
configuration group and `PressureGradType` enum, leaving `Centered` as the default. Its
sub-options (horizontal reconstruction order, vertical-reconstruction mode, quadrature) must
be configurable, and the centered scheme must be recoverable as the lowest-order limit.

### 2.6 Desired: Extensibility to sixth order and to tidal/geoid geopotential terms

The framework should accommodate a future higher-order variant (e.g. sixth order, added to
the `PressureGradType` enum when implemented) and the
tidal-potential and self-attraction-and-loading contributions to the geopotential
(supplied by `VertCoord`; see {ref}`Pressure Gradient <omega-design-pressure-grad>`
§2.4), without restructuring.

## 3 Algorithmic Formulation

### 3.1 Continuous target form

The high-order PGF starts from the layer-integrated, finite-volume form of the pressure
and geopotential terms derived in
{ref}`the governing equations <omega-design-governing-eqns-omega1>`
(the {ref}`full volume-integral form <omega-v1-vh-momentum-reynolds2>`). Neglecting the
turbulent correlations and dropping the resolved-component notation, the layer tendency is
({ref}`Pressure Gradient <omega-design-pressure-grad>` §3.2):

$$
T^p &= - \int_A \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \rho_0 \, \left( \nabla \Phi \right) \, d\tilde{z} \, dA \\
& - \int_{\partial A} \left( \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \rho_0 \left(\alpha p \right) \, d\tilde{z} \right) dl \\
& - \int_A \rho_0 \left[ \alpha \, p \, \nabla \tilde{z}_k^{\text{top}} \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} \, dA \\
& + \int_A \rho_0 \left[ \alpha \, p \, \nabla \tilde{z}_k^{\text{bot}} \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \, dA.
$$ (ho-target)

The four terms are: the geopotential (gravity) body force integrated over the layer
volume; the side-wall integral of $\alpha p$ (the pressure traction on the cell faces); and
two metric terms accounting for the pressure traction on the sloping top and bottom layer
interfaces. This is the Adcroft et al. (2008) finite-volume route — the net pressure force
on each control volume is obtained by integrating in-situ pressure over the faces, rather
than by forming the pointwise product $\alpha\nabla p$. We adopt this form because it
matches Omega's non-Boussinesq layer-integral momentum equation exactly and is the form
MOM6 uses successfully for the same problem.

Adcroft et al. (2008) were targeting an isopycnal model (constant density within a layer)
and the analytically integrable Wright equation of state. Omega differs on both counts:
density varies within a layer in a general ALE coordinate, and TEOS-10 is far more
expensive and not practical to integrate in closed form. The remainder of this section
adapts the form [](#ho-target) to those two realities.

### 3.2 A key simplification: pressure is linear in pseudo-height

Omega's vertical coordinate is pseudo-height,
$\tilde z \equiv -p/(\rho_0 g)$ ({ref}`omega-design-governing-eqns-omega1` §5). Pressure is
therefore, by definition, an **exactly linear function of $\tilde z$**:

$$
p(\tilde z) = p_{k}^{\text{top}} - \rho_0 g \left( \tilde z - \tilde z_k^{\text{top}} \right),
\qquad \tilde z_k^{\text{bot}} \le \tilde z \le \tilde z_k^{\text{top}},
$$ (p-linear)

anchored at the interface pressures $p_k^{\text{top}}$, $p_k^{\text{bot}}$ that `VertCoord`
already provides (`PressureInterface`). No reconstruction of $p$ within a layer is needed
or permitted — using [](#p-linear) is exactly the `VertCoord` anchoring required by
Requirement 2.4. Consequently, the only within-layer nonlinearity in the integrand
$\alpha p$ comes from the variation of $\alpha$ through the reconstructed $\Theta$ and $S$.

### 3.3 Reference-state expansion of the equation of state

This is the core device that satisfies Requirement 2.2. We **never integrate TEOS-10**.
Instead, for each cell $i$ and layer $k$ we expand $\alpha$ in a Taylor series about the
layer's reference state $(\Theta_{i,k}, S_{i,k}, p_{i,k}^{\text{mid}})$, where
$\Theta_{i,k}$ and $S_{i,k}$ are the prognostic layer means and $p_{i,k}^{\text{mid}}$ is
the mid-layer pressure from `VertCoord` (`PressureMid`):

$$
\alpha(\Theta, S, p) \approx \alpha_{0} + \alpha_{\Theta}\,(\Theta - \Theta_{i,k}) + \alpha_{S}\,(S - S_{i,k}) + \alpha_{p}\,(p - p_{i,k}^{\text{mid}}),
$$ (alpha-taylor)

with

$$
\alpha_0 = \alpha(\Theta_{i,k}, S_{i,k}, p_{i,k}^{\text{mid}}), \quad
\alpha_\Theta = \left.\frac{\partial \alpha}{\partial \Theta}\right|_0, \quad
\alpha_S = \left.\frac{\partial \alpha}{\partial S}\right|_0, \quad
\alpha_p = \left.\frac{\partial \alpha}{\partial p}\right|_0.
$$ (alpha-derivs)

The four coefficients in [](#alpha-derivs) are obtained from a **single** evaluation of the
TEOS-10 polynomial, since the first derivatives reuse the same polynomial coefficients (the
same quantities computed by `gsw.specvol_first_derivatives`, which the Polaris reference
solution already uses). Indeed $\alpha_0$ is exactly the `Eos::SpecVol` field Omega already
computes, so only the three derivative fields are additional work.

The first-order expansion [](#alpha-taylor) is the default. A second-order expansion (adding
$\alpha_{\Theta\Theta}$, $\alpha_{SS}$, $\alpha_{\Theta S}$, $\alpha_{pp}$, …) is an option
where stratification is strong; it adds derivative evaluations but still no in-integral
TEOS-10 calls.

Because $\Theta$, $S$ are reconstructed as low-order polynomials in $\tilde z$ (and across
edges) and $p$ is linear in $\tilde z$ by [](#p-linear), the expansion [](#alpha-taylor)
makes $\alpha$ a **low-order polynomial** whose product with $p$ integrates **exactly** by a
low-order Gauss rule. The semi-analytic integration is performed on the cheap Taylor model,
not on TEOS-10 — directly resolving the feasibility and cost concerns with applying
Adcroft-style analytic integration to TEOS-10.

### 3.4 Mean-preserving vertical reconstruction

Within layer $k$ of column $i$ we reconstruct the conservative temperature and absolute
salinity as

$$
\Theta(\tilde z) = \Theta_{i,k} + \Theta'_{i,k}(\tilde z), \qquad
S(\tilde z) = S_{i,k} + S'_{i,k}(\tilde z),
$$ (vert-recon)

where $\Theta'_{i,k}$, $S'_{i,k}$ are parabolic (PPM-style) deviations built from the
neighboring layer means and constrained to **integrate to zero over the layer**
($\int_{\tilde z_k^{\text{bot}}}^{\tilde z_k^{\text{top}}} \Theta'_{i,k}\, d\tilde z = 0$).
This mean-preserving property (Requirement 2.4) guarantees the PGF uses the same layer-mean
state as the rest of the model; the reconstruction only supplies the smoother sub-layer
shape needed to integrate $\alpha p$ to high order. The "constant density within a layer"
assumption of an isopycnal model is recovered exactly by the degenerate choice
$\Theta'_{i,k} = S'_{i,k} = 0$, which is available as the cheapest configuration.

### 3.5 Analytic layer integral of the side-wall term

Combining [](#p-linear), [](#alpha-taylor), and [](#vert-recon), the side-wall integrand
$\alpha p$ in a single column/layer is a polynomial in $\tilde z$: with first-order $\alpha$
and parabolic $\Theta', S'$, $\alpha$ is parabolic and $p$ is linear, so $\alpha p$ is cubic
and is integrated **exactly** by a two-point Gauss-Legendre rule over the layer:

$$
\Pi_{i,k} \equiv \int_{\tilde z_k^{\text{bot}}}^{\tilde z_k^{\text{top}}} \alpha\, p \; d\tilde z
\;=\; \tilde h_{i,k} \sum_{q} w_q \, \alpha\!\left(\tilde z_q\right) p\!\left(\tilde z_q\right),
$$ (sidewall-int)

with Gauss nodes $\tilde z_q$ and weights $w_q$ on the layer, and $\alpha(\tilde z_q)$,
$p(\tilde z_q)$ from [](#alpha-taylor) and [](#p-linear). No TEOS-10 calls occur inside the
sum. The quadrature order is configurable and is matched to the reconstruction so the
integral is exact for the reconstructed polynomial.

### 3.6 Horizontal reconstruction and the edge gradient

The discrete tendency lives at edge $e$ and is the edge-normal projection of the PGF. The
side-wall line integral $\int_{\partial A}(\cdots)\,dl$ in [](#ho-target) becomes, on the
TRiSK C-grid, the difference of the two adjacent columns' face contributions across the
edge, divided by the cell-center distance $d_e$:

$$
\left[\nabla_n \Pi\right]_{e,k} = \frac{1}{d_e}\sum_{i \in CE(e)} -n_{e,i}\, \widehat{\Pi}_{i,k},
$$ (edge-grad)

where $\widehat\Pi_{i,k}$ is the column integral [](#sidewall-int) reconstructed to the edge.
To reach fourth order, $\Theta$, $S$ (and hence $\alpha$ and $\Pi$) are reconstructed from
cell centers to the edge with a **cubic** reconstruction consistent with the TRiSK edge
stencils used elsewhere in Omega (the same neighborhood used for high-order tracer
reconstruction; cf. the third-order interface reconstruction noted at
{ref}`omega-design-governing-eqns-omega1` §10, White & Adcroft 2008). The lowest-order limit
(two-cell, centered) reproduces the operator already used by `PressureGradCentered`. We do
**not** assume anything constant across the edge — the horizontal density contrast between
adjacent columns is fully reconstructed.

### 3.7 Discrete hydrostatic consistency (the steep-layer robustness property)

Robustness for thin, steeply sloped layers (Requirement 2.3) is delivered not by subtracting
a background profile but by a structural property of the discretization. We require **discrete
hydrostatic consistency**:

> When $\Theta$ and $S$ are horizontally uniform — a physically resting, flat-isopycnal ocean
> — the discrete horizontal PGF at every edge and layer must be **exactly zero, to machine
> precision**, for arbitrary tilt of the coordinate surfaces $\tilde z_k^{\text{top}}$,
> $\tilde z_k^{\text{bot}}$ between the two columns and for an arbitrary vertical profile.

The finite-volume form [](#ho-target) is constructed to satisfy this exactly: with $\Theta$,
$S$ horizontally uniform, $\alpha(\Theta,S,p)$ is the same function of pressure in both
columns, so the side-wall integral [](#sidewall-int) and the sloping-interface metric terms
(§3.8) combine to the net pressure force on a control volume with horizontally uniform density,
which is balanced exactly by the geopotential term. The cancellation must be enforced
*discretely*: the pressure term and the geopotential term must share the same interface
locations (from `VertCoord`), the same edge reconstruction [](#edge-grad), and the same
quadrature nodes, so that their difference vanishes term by term rather than only in the
continuum limit.

This property is the reason the scheme is robust: it is *reference-free* (it depends on no
background state, so it holds identically in the warm pool, at the poles, and over a seamount)
and it is directly *testable* (§5.1, §5.2). It strictly supersedes a Shchepetkin–McWilliams
(2003)-style reference-profile subtraction, which only approximates this cancellation to the
extent the local column resembles a global reference and therefore degrades precisely in the
strong-gradient, steep-layer regions of interest. For non-uniform $\Theta$/$S$ the remaining
error is genuine truncation error, controlled by the reconstruction order (§3.4, §3.6) and by
the consistent treatment of compressibility $\alpha_p$ (§3.3); no formulation removes it, and
a reference subtraction would not help. Because the consistency cancellation is exact in
exact arithmetic, the only residual in the resting state is round-off, which is at the
machine-precision floor in double precision; single-precision builds (`OMEGA_SINGLE_PRECISION`)
may expose a small floor, addressed if needed as a local implementation-level conditioning
choice rather than a physics option.

### 3.8 Geopotential and metric terms

The geopotential body-force term (first line of [](#ho-target)) is
$g\,\nabla z + \nabla(\phi_{TP}+\phi_{SAL})$. The geometric height $z$ and geopotential are
provided by `VertCoord` (`GeomZMid`/`GeomZInterface`, `GeopotentialMid`; computed via
{ref}`omega-design-governing-eqns-omega1` Eqs. discrete-z and the geopotential relation).
Its layer average and edge-normal gradient are evaluated with the same high-order edge
reconstruction [](#edge-grad). The tidal-potential and self-attraction-and-loading
contributions enter through `VertCoord` and are differenced identically (Requirement 2.6);
they are zero in early Omega versions.

The two metric terms (third and fourth lines of [](#ho-target)) use the **same** edge-normal
operator to evaluate the interface slopes $\nabla \tilde z_k^{\text{top}}$,
$\nabla \tilde z_k^{\text{bot}}$, with $\alpha\,p$ evaluated at the interfaces from
[](#alpha-taylor) and [](#p-linear). These terms vanish for flat interfaces and supply the
leading correction for sloping layers; together with the side-wall integral [](#sidewall-int)
they form a discretely consistent finite-volume pressure force on the layer control volume.

### 3.9 Reduction to the centered scheme

As a consistency check, the high-order scheme collapses to the implemented
`PressureGradCentered` form
($T^p_{e,k} = -\nabla M + \tfrac12(p_0+p_1)\nabla\alpha - \nabla\Phi$, with
$M = \alpha p + g z$) in the joint limit of: constant in-layer reconstruction
($\Theta' = S' = 0$, §3.4), midpoint quadrature in place of [](#sidewall-int), and two-cell
centered horizontal differencing (§3.6). This guarantees the new code reproduces the existing
scheme in its lowest-order configuration and provides a direct path for regression testing.

### 3.10 Per-step algorithm summary

1. From `VertCoord`: read `PressureInterface`, `PressureMid`, `GeomZInterface`/`GeomZMid`,
   geopotential, and interface pseudo-heights (already computed diagnostically each step).
2. Per cell-layer: obtain $\alpha_0$ (= existing `Eos::SpecVol`) and the derivatives
   $\alpha_\Theta, \alpha_S, \alpha_p$ from one TEOS-10 evaluation ([](#alpha-derivs)).
3. Per cell-layer: build mean-preserving PPM deviations $\Theta', S'$ ([](#vert-recon)).
4. Per edge-layer: reconstruct edge quantities (cubic, §3.6) using interface locations and
   quadrature nodes shared by the pressure and geopotential terms (the discrete hydrostatic
   consistency requirement, §3.7); evaluate the analytic layer integral [](#sidewall-int),
   the geopotential gradient, and the metric terms (§3.8); assemble $T^p_{e,k}$ and accumulate
   into the tendency with `EdgeMask`.

## 4 Design

The high-order scheme extends the existing `PressureGrad` class in
`components/omega/src/ocn/PGrad.{h,cpp}`. The class manager, creation/retrieval/removal
methods, named-instance map, and dispatch logic are unchanged; the work is to flesh out the
`PressureGradHighOrder` functor (currently a no-op placeholder) and to widen the inputs the
manager hands to it.

### 4.1 Data types and parameters

#### 4.1.1 Parameters

Configuration lives under the existing `PressureGrad` YAML group. The `PressureGradType`
enum (PGrad.h:22) is updated to give the high-order scheme a descriptive name and to drop
the unimplemented second placeholder, leaving a commented stub for a future variant:

```c++
enum class PressureGradType {
   Centered,      // existing 2nd-order Montgomery scheme
   FiniteVolume   // high-order finite-volume analytic-integration scheme (this design)
   // , <FutureVariant>   // e.g. a 6th-order option, added when implemented
};
```

New sub-options for the high-order scheme:

```yaml
    PressureGrad:
       PressureGradType: 'FiniteVolume'   # Centered | FiniteVolume
       ReconstructionOrder: 4             # horizontal cell->edge order (2 = centered limit)
       VerticalReconstruction: 'ppm'      # 'constant' (isopycnal limit) | 'ppm'
       QuadraturePoints: 2                # per-layer Gauss points for the side-wall integral
```

The fourth-order target is `ReconstructionOrder: 4` with parabolic (`ppm`) vertical
reconstruction. The centered scheme is recovered by `ReconstructionOrder: 2` and
`VerticalReconstruction: 'constant'` (§3.9).

#### 4.1.2 New EOS support

`PressureGradHighOrder` needs $\alpha$ together with its first derivatives. The `Eos` class
(Eos.h) currently exposes `computeSpecVol`, `computeSpecVolDisp`, and
`computeBruntVaisalaFreqSq`, but no specific-volume derivatives. This design adds one method
and three device fields:

```c++
// New on the Eos class
Array2DReal SpecVolDThetaCons; ///< d(alpha)/d(ConservTemp) at cell centers
Array2DReal SpecVolDSalt;      ///< d(alpha)/d(AbsSalinity)  at cell centers
Array2DReal SpecVolDPressure;  ///< d(alpha)/d(Pressure)     at cell centers

/// Compute specific volume AND its first derivatives in one TEOS-10 pass
void computeSpecVolAndDerivs(const Array2DReal &ConservTemp,
                             const Array2DReal &AbsSalinity,
                             const Array2DReal &Pressure);
```

The TEOS-10 derivatives reuse the polynomial coefficients already assembled inside the
`Teos10Eos` functor (`calcPCoeffs`/`calcDelta`), so the marginal cost over `computeSpecVol`
is the derivative arithmetic only (Requirement 2.2). The linear and constant EOS options
supply trivial analytic derivatives.

#### 4.1.3 `PressureGradHighOrder` functor

The functor mirrors `PressureGradCentered` (cached mesh/coordinate arrays, `Enabled` flag,
`chunkStart`/`chunkLength` vertical iteration) but takes the additional reconstruction and
EOS-derivative inputs. Its `operator()` signature is widened from the placeholder to:

```c++
KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                const Array2DReal &PressureMid,
                                const Array2DReal &PressureInterface,
                                const Array2DReal &GeomZInterface,
                                const Array1DReal &TidalPotential,
                                const Array1DReal &SelfAttractionLoading,
                                const Array2DReal &SpecVol,
                                const Array2DReal &ConservTemp,    // new
                                const Array2DReal &AbsSalinity,    // new
                                const Array2DReal &SpecVolDThetaCons, // new
                                const Array2DReal &SpecVolDSalt,      // new
                                const Array2DReal &SpecVolDPressure)  // new
                                const;
```

(The existing centered signature is unchanged.) Additional cached members hold the cubic
edge-reconstruction stencil and weights. The functor implements §3.3–§3.8 per edge and
vertical chunk — sharing interface locations and quadrature nodes between the pressure and
geopotential terms to satisfy discrete hydrostatic consistency (§3.7) — accumulating into
`Tend` with `EdgeMask`, exactly as the centered functor does.

### 4.2 Methods

`PressureGrad::computePressureGrad` keeps its role of selecting the configured option, but
its input list grows so the high-order branch can reach $\Theta$, $S$, and the EOS
derivatives. The signature becomes:

```c++
void computePressureGrad(Array2DReal &Tend, const Array2DReal &PressureMid,
                         const Array2DReal &PressureInterface,
                         const Array2DReal &SpecVol,
                         const Array2DReal &GeomZInterface,
                         const Array2DReal &PseudoThick,
                         const Array2DReal &ConservTemp,        // new
                         const Array2DReal &AbsSalinity,        // new
                         const Eos *EqState) const;             // new: derivative fields
```

The `Centered` branch ignores the new arguments and is byte-for-byte unchanged; the
`FiniteVolume` branch dispatches to `PressureGradHighOrder` through the same
`parallelForOuter`/`parallelForInner` team pattern used today (PGrad.cpp). The call site in
`Tendencies.cpp` is updated to pass $\Theta$, $S$ (from the tracer state) and the `Eos`
instance, which it already references for `SpecVol`.

Creation, retrieval, and removal (`init`, `create`, `get`, `getDefault`, `clear`, `erase`)
are unchanged. The high-order constructor additionally caches the reconstruction stencil and
quadrature weights from config.

### 4.3 Consistency and follow-up work

The high-order PGF is deliberately more accurate in the vertical than the rest of Omega's
layer-mean operators. This is sanctioned by {ref}`omega-design-governing-eqns-omega1` §9,
which singles out the PGF as **the** exception to the piecewise-constant assumption ("we will
ignore most of these [within-layer deviation] terms in Omega. The exception is the pressure
gradient force … appropriate for the *simple* pressure gradient force targeted for *early*
versions of Omega. This assumption will be revisited at a later date."). The mean-preserving,
`VertCoord`-anchored constraints (Requirement 2.4) keep the PGF consistent with the model
state at the layer-mean level. For completeness, the following terms would need higher-order
vertical treatment to make the **whole** model consistent with this PGF; they are explicitly
out of scope here and flagged as follow-up:

- **Vertical remapping and advection** ({ref}`omega-design-vert-adv` and the ALE
  remap): a mean-preserving PPM remap consistent with §3.4.
- **EOS evaluation point**: buoyancy and $N^2$ ({ref}`omega-design-eos`) are currently
  evaluated at layer mid; consistency with the reconstructed profile would evaluate them on
  the same reconstruction.
- **Dropped within-layer products** $\overline{\delta\varphi\,\delta u}$ in the tracer and
  momentum equations ({ref}`omega-design-governing-eqns-omega1` §9), which become non-negligible
  only at much higher resolution.

## 5 Verification and Testing

Testing reuses and extends the Polaris `horiz_press_grad` task family
(`polaris/tasks/ocean/horiz_press_grad`), which is already Omega-only and built around a
quasi-analytic, layer-mean TEOS-10 reference solution (`reference.py`, surface-anchored,
4-point Gauss). That analytic reference remains valid as "truth" for the high-order scheme,
so the principal changes are scheme selection and revised pass criteria, plus two full-model
acceptance tests.

Because the design optimizes for **coarse-resolution accuracy** rather than convergence rate
(§1), the primary pass/fail gate is an absolute error tolerance at a representative coarse
resolution, with the convergence slope demoted to a loose secondary diagnostic.

### 5.1 Test: Two-column HPGA convergence (extend existing)

Extend the four existing variants — `temperature_gradient`, `salinity_gradient`,
`surface_pressure_gradient`, `ztilde_gradient` — to run both the centered and high-order
schemes:

- **Scheme selection.** Add `PressureGrad: { PressureGradType: FiniteVolume, … }` to
  `forward.yaml` and parametrize each task over `centered` vs. `high_order`. The forward step
  still runs a single time step with only `PressureGradTendencyEnable: true`, reading the PGF
  acceleration from `NormalVelocityTend`.
- **Reference.** `reference.py`/`analysis.py` compare `NormalVelocityTend` against the
  layer-mean analytic HPGA (unchanged). For the high-order scheme the layer-mean comparison
  remains the correct target, since the scheme is a finite-volume, layer-averaged
  discretization.
- **Primary pass criterion (new):** at a representative coarse resolution (e.g. the coarsest
  in `horiz_resolutions`), the absolute RMS HPGA error vs. the reference must be below a
  tolerance, **and** the high-order RMS error must be below the centered RMS error at that
  same resolution (the scheme must demonstrably help where it matters).
- **Consistency check (retained):** `omega_vs_polaris_rms_threshold` (~1e-10 m/s²) — Omega's
  forward output must still match the Python-computed HPGA, confirming the implementation
  matches the intended discretization.
- **Secondary diagnostic (loosened):** the `omega_vs_reference_convergence_rate_*` band and
  `omega_vs_reference_high_res_rms_threshold` are retained but relaxed/retuned; the
  convergence slope is reported and required only to be no worse than the centered scheme,
  not to hit a specific order. New cfg keys mirror the existing ones
  (`horiz_press_grad.cfg`), e.g. a coarse-resolution tolerance and a
  `high_order_vs_centered` ratio gate.

Tests Requirements 2.1, 2.2 (via the bounded-EOS implementation exercised), 2.4, 2.5.

### 5.2 Unit test: Discrete hydrostatic consistency (exact resting state)

This is the primary test of the steep-layer robustness property (Requirement 2.3, §3.7).
Construct a two-column (or seamount) configuration with **horizontally uniform** conservative
temperature and absolute salinity but with the coordinate interfaces
$\tilde z_k^{\text{top}}$, $\tilde z_k^{\text{bot}}$ deliberately tilted between the columns
(arbitrary, including steep slopes and thin layers), over an arbitrary vertical $T$/$S$
profile. **Pass:** the high-order PGF is zero at every edge and layer to **machine precision**
(double-precision builds; the threshold tracks `Real`'s epsilon and the magnitude of the
hydrostatic terms, not a physical tolerance). This verifies that the pressure and geopotential
terms cancel discretely, term by term, independent of any background state. It is implemented
as a fast C++ unit test and also exercised as a configuration of the Polaris two-column task.

### 5.3 Test: Seamount resting state (steep-layer robustness)

Use the Polaris seamount task with tilted layers over a seamount and a horizontally uniform
$T$/$S$ stratification, integrated for a fixed period. **Pass:** the maximum spurious velocity
stays below a threshold and is substantially smaller for the high-order scheme than for the
centered scheme. This is the dynamical counterpart of §5.2 (the same resting state, now run
forward in the full solver) and the direct test of Requirement 2.3.

### 5.4 Test: Overflow (full non-Boussinesq dynamics)

Use the Polaris overflow task to exercise the PGF within the full non-Boussinesq equations
with strongly sloping layers and active dynamics. **Pass:** the solution remains stable and
the down-slope evolution agrees with the reference behavior; spurious mixing/velocity
attributable to PGF error is reduced relative to the centered scheme. This tests Requirements
2.1 and 2.3 under realistic, coupled conditions.

### 5.5 Test: Reduction to the centered scheme (regression)

Configure the high-order option in its lowest-order limit (§3.9:
`ReconstructionOrder: 2`, `VerticalReconstruction: constant`, midpoint quadrature) and confirm
it reproduces `PressureGradCentered` to round-off on the two-column test. This guards
Requirement 2.5 and protects the existing default during refactoring.
