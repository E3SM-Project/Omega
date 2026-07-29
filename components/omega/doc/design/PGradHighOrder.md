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

A guiding choice that follows from principles 2 and 3 is that the metric of success *for
simulations* is **low absolute PGF error at affordable (coarse) resolution**, not high
asymptotic order of convergence. Omega will be strongly resource-limited in how fine a
mesh it can run, so a scheme whose error is small at coarse resolution is preferred over
one with a steeper convergence slope but larger coarse-resolution error. Order of
accuracy is configurable and a fourth-order variant is targeted, but it is a means to the
end of coarse-resolution accuracy, not the objective itself.

The measured convergence rate plays a different but equally essential role: it is the
**verification** metric. A scheme that does not converge at its designed order is
mis-implemented, and no amount of favorable absolute error excuses that. The two metrics
therefore answer two different questions and both must be satisfied — *is the code correct?*
(convergence rate at the designed order) and *is it useful at the resolutions we can run?*
(absolute error).

### 1.1 Why this work is delivered in two phases

Realistic global configurations run with `PressureGradCentered` show spurious flow in the
bottom layer over sloping bathymetry, large enough to drive numerical instability. What is
missing there is not high-order accuracy — it is *consistency*: the centered scheme does not
account for the pressure force on a layer's sloping top and bottom interfaces at all, and the
resulting error accumulates downward through the column (§3.7.3). That error is present at
second order and does not require a fourth-order scheme to remove.

The work is therefore split so the consistency fix can be delivered and tested on its own:

- **Phase 1 — a consistent second-order scheme.** The finite-volume control-volume form,
  the sloping-interface integral, the equation-of-state expansion about a state shared across
  each edge, and a mean-preserving *linear* reconstruction of $\Theta$ and $S$ in pressure.
  The horizontal operator is the same two-cell stencil the centered scheme already uses, so
  the scheme remains second order in the horizontal. What Phase 1 buys is the robustness
  property of Requirement 2.3: the pressure gradient is zero to machine precision for any
  resting ocean whose profile varies linearly with pressure, no matter how the layers tilt.
- **Phase 2 — fourth order.** A parabolic vertical reconstruction and a wider horizontal
  stencil, plus the option of a second-order equation-of-state expansion. This raises the
  order of accuracy; it does not change the robustness property, which Phase 1 already
  establishes and Phase 2 must preserve.

The two phases share one implementation, one set of configuration options, and one test
suite; Phase 2 widens settings that Phase 1 puts in place. Section 3 marks which parts of
the formulation belong to which phase, and §4.5 records what each phase delivers and what it
depends on.

## 2 Requirements

### 2.1 Requirement: Higher accuracy than the centered scheme at affordable resolution

The high-order PGF must produce substantially lower absolute error than
`PressureGradCentered` at the coarse-to-moderate resolutions Omega can afford to run,
for representative stratified columns with horizontal gradients of temperature, salinity,
surface pressure, and coordinate slope. Absolute error at affordable resolution is what
determines whether the scheme is *useful*; it is complemented by the separate verification
requirement in §2.6, which determines whether it is *correct*.

### 2.2 Requirement: Bounded TEOS-10 cost

The number of TEOS-10 specific-volume evaluations performed by the PGF must be
independent of the quadrature order and bounded at approximately one evaluation (plus its
first derivatives) per cell per layer per time step — comparable to the cost Omega already
pays to compute the `Eos::SpecVol` field. The scheme must not require evaluating the full
TEOS-10 polynomial at each quadrature point, nor integrating the TEOS-10 polynomial
itself semi-analytically.

### 2.3 Requirement: Robustness for thin, steeply sloped layers via discrete hydrostatic consistency

The scheme must suppress spurious pressure-gradient accelerations when layers are thin and steeply
sloped relative to surfaces of constant pressure. Four things are required.

**2.3.1 No spurious PGF for water columns the scheme resolves exactly.** For any resting ocean
whose temperature and salinity profile the scheme's vertical reconstruction reproduces exactly, the
computed horizontal PGF must be **zero to machine precision** at every edge and layer — for any
tilt of the coordinate surfaces, any layer thickness, and any bathymetry. This is the core
robustness property: it says the scheme adds nothing of its own when there is nothing to add.

**2.3.2 The design must say which water columns those are.** "Water columns the scheme resolves
exactly" is not a fixed set; it depends on the reconstruction, and it widens as the reconstruction
order increases. Each phase of this design must state its own set explicitly (§3.7.3) so that the
guarantee in §2.3.1 can be tested rather than assumed.

**2.3.3 That set must be large enough to matter.** A scheme that is exact only for an ocean of
uniform temperature and salinity would satisfy §2.3.1 trivially and be of no practical use, since
no such ocean exists. The set must at minimum include profiles that vary **linearly with
pressure**, which captures both compressibility and smooth stratification to leading order.

**2.3.4 Everywhere else, substantially smaller error than the centered scheme.** For realistic
profiles, which no reconstruction reproduces exactly, some residual PGF is unavoidable and is
*correct* — the two neighboring columns genuinely do hold slightly different water (§3.7.1). The
requirement is that this residual shrink at least one order faster in layer thickness than
`PressureGradCentered` does, and that it not be swamped by errors the scheme itself introduces
through inexact integration or an inconsistent geopotential.

These properties must hold without reference to any background or reference profile, so that they
hold equally in the warm pool, at the poles, and over a seamount.

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
sub-options (edge stencil width, vertical-reconstruction mode, quadrature) must be configurable,
and the centered scheme must be recoverable as a configuration of the new one. The two phases must
share one set of configuration keys, so that a configuration written for Phase 1 continues to work
unchanged when Phase 2 lands.

### 2.6 Requirement: Verified order of accuracy

The implemented scheme must converge, under refinement of a smooth manufactured or
quasi-analytic reference solution, at the order of accuracy it is configured for
(nominally fourth order for the targeted variant, second order when reduced to the
centered limit). This is a verification requirement: a measured slope that falls short of
the designed order indicates a defect in the implementation — a mis-weighted quadrature
point, an inconsistent reconstruction stencil, a dropped correction term — and must be
diagnosed rather than accepted, even when the absolute error at coarse resolution already
satisfies Requirement 2.1. The two requirements are independent and both are gating:
Requirement 2.1 establishes that the scheme is useful at the resolutions Omega can afford,
Requirement 2.6 establishes that it is the scheme the design describes.

### 2.7 Desired: Extensibility to sixth order and to tidal/geoid geopotential terms

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
interfaces.

```{note}
The signs on the two metric terms are inherited from
{ref}`omega-design-governing-eqns-omega1` and have not been independently re-derived here. Applying
Leibniz' rule to lines 2–4 so that they telescope to $-\int\rho_0\nabla(\alpha p)\,d\tilde z$
appears to require the opposite signs on those two lines; with the signs as printed, the discrete
resting-state cancellation of §3.7 leaves a factor of two rather than zero. This may be a difference
in how $\nabla\tilde z^{\text{top}}$ or the traction normal is defined rather than an error.
`OmegaV1GoverningEqns.md` is not maintained against the code and cannot settle it. The convention
must be re-derived from $\tilde z \equiv -p/(\rho_0 g)$ against the conventions the code actually
fixes — `PressureGradCentered` for the tendency sign, `CellsOnEdge` ordering and `EdgeMask`, and
`VertCoord` for accumulation direction and interface indexing — before implementation, using the
§3.9 reduction to the centered scheme as the check with a known answer.
``` This is the Adcroft et al. (2008) finite-volume route — the net pressure force
on each control volume is obtained by integrating in-situ pressure over the faces, rather
than by forming the pointwise product $\alpha\nabla p$. We adopt this form because it
matches Omega's non-Boussinesq layer-integral momentum equation exactly and is the form
MOM6 uses successfully for the same problem.

Adcroft et al. (2008) were targeting an isopycnal model (constant density within a layer)
and the analytically integrable Wright equation of state. Omega differs on both counts:
density varies within a layer in a general ALE coordinate, and TEOS-10 is far more
expensive and not practical to integrate in closed form. The remainder of this section
adapts the form [](#ho-target) to those two realities.

#### 3.1.1 The four terms are one control volume and must be built together

The four terms of [](#ho-target) are not independent pieces to be discretized separately and added
up. They are the pressure forces on the four faces of a single control volume — the two vertical
faces where the layer meets its neighboring columns, and the sloping top and bottom interfaces —
together with the weight of the water inside. In a genuinely resting ocean these forces balance
exactly. Whether the *discrete* forces balance depends entirely on whether they are built from a
common description of the water column, and that is what §3.7 is about.

Two observations make this practical to enforce.

First, use pressure itself as the vertical integration variable within a column. It is equivalent
to $\tilde z$ by [](#p-linear), but unlike $\tilde z$ it means the same thing in every column
regardless of surface pressure. In terms of $p$, the hydrostatic relation is
$\partial z / \partial p = -\alpha/g$, so the side-wall pressure force and the geometric height
entering the geopotential are integrals of **the same specific volume over pressure**:

$$
\Pi_{i,k} = \frac{1}{g}\int_{p_{i,k}^{\text{top}}}^{p_{i,k}^{\text{bot}}} \hat\alpha_{i,k}(p)\, p \; dp,
\qquad
z_{i}(p) = z_{i}^{\text{anchor}} + \frac{1}{g}\int_{p}^{p_i^{\text{anchor}}} \hat\alpha_{i}(p')\, dp',
$$ (layer-integrals)

where $\hat\alpha_{i,k}(p)$ is the reconstructed specific volume in column $i$, layer $k$
(§3.3–§3.4). Writing both terms this way makes the shared ingredient explicit: it is
$\hat\alpha(p)$, and if the pressure term and the geopotential term use different versions of it,
they cannot balance.

Second, the layer's top and bottom interfaces are *sloping*, so the pressure force on them is an
integral along the slope, not a value at a point. Between two cell centers the interface sweeps
through a range of pressures, and the force per unit slope is the **average of $\hat\alpha p$ over
that pressure range** — not the average of the two cells' interface values. The two agree only when
the tilt is small. Supplying this integral is the piece `PressureGradCentered` omits altogether,
and it is the central content of Phase 1; the explicit form is given in §3.8.

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
TEOS-10 polynomial: the derivatives are analytic derivatives of that same polynomial, taken
at the same normalized state, so no second call to the equation of state occurs (they are the
same quantities computed by `gsw.specvol_first_derivatives`, which the Polaris reference
solution already uses). The arithmetic is not free — $\alpha_p$ reuses the pressure
coefficients assembled for $\alpha_0$, but $\alpha_\Theta$ and $\alpha_S$ need coefficient
sets of their own — though it shares the square root and the normalization, and it is the
*call count*, not the polynomial arithmetic, that Requirement 2.2 bounds. Indeed $\alpha_0$ is
exactly the `Eos::SpecVol` field Omega already computes, so only the three derivative fields
are additional work.

The first-order expansion [](#alpha-taylor) is the default. A second-order expansion (adding
$\alpha_{\Theta\Theta}$, $\alpha_{SS}$, $\alpha_{\Theta S}$, $\alpha_{pp}$, …) is an option
where stratification is strong; it adds derivative evaluations but still no in-integral
TEOS-10 calls.

#### 3.3.1 The expansion point must be shared across the edge

Expanding about the *cell-local* reference state $(\Theta_{i,k}, S_{i,k}, p_{i,k}^{\text{mid}})$
would quietly defeat the robustness property of §3.7. Two neighboring columns would then be using
two slightly different approximations to the *same* equation of state, taken at two different
points. Their specific-volume profiles $\hat\alpha_{L,k}(p)$ and $\hat\alpha_{R,k}(p)$ would then
disagree — by $O(\alpha_{pp}\,\Delta p^{\text{mid}}_e\, \tilde h)$ — even for a resting ocean the
reconstruction reproduces perfectly. The forces on the two vertical faces would no longer balance,
and the scheme would generate spurious flow from nothing but its own EOS approximation.

The design therefore **shares the expansion point across each edge**. The coefficients
[](#alpha-derivs) are computed once per cell per layer, exactly as described above (one TEOS-10
evaluation, Requirement 2.2), and are then averaged to the edge together with the reference state
itself:

$$
\bar\alpha_{0}^{e,k}, \; \bar\alpha_{\Theta}^{e,k}, \; \bar\alpha_{S}^{e,k}, \; \bar\alpha_{p}^{e,k},
\qquad
\left(\bar\Theta^{e,k}, \bar S^{e,k}, \bar p^{e,k}\right) = \tfrac12 \sum_{i \in CE(e)} \left(\Theta_{i,k}, S_{i,k}, p_{i,k}^{\text{mid}}\right),
$$ (edge-ref)

and **both** columns' contributions to edge $e$ are evaluated with this single set, so that the two
sides of the edge see one and the same equation of state.

This has a cost consequence worth recording. The column integral $\Pi_{i,k}$ is no longer a per-cell
quantity that can be computed once and then differenced across each of the cell's edges; it depends
on the edge through [](#edge-ref) and must be evaluated per edge, twice. On a hexagonal TRiSK mesh
this is roughly three times as many integral evaluations as a cell-based formulation. The TEOS-10
call count is unaffected, which is the cost Requirement 2.2 binds; the extra work is polynomial
arithmetic on coefficients already in cache.

Because $\Theta$, $S$ are reconstructed as low-order polynomials in $\tilde z$ (and across
edges) and $p$ is linear in $\tilde z$ by [](#p-linear), the expansion [](#alpha-taylor)
makes $\alpha$ a **low-order polynomial** whose product with $p$ integrates **exactly** by a
low-order Gauss rule. The semi-analytic integration is performed on the cheap Taylor model,
not on TEOS-10 — directly resolving the feasibility and cost concerns with applying
Adcroft-style analytic integration to TEOS-10.

### 3.4 Mean-preserving vertical reconstruction

Within layer $k$ of column $i$ we reconstruct the conservative temperature and absolute
salinity as deviations from the prognostic layer means,

$$
\Theta(p) = \Theta_{i,k} + \Theta'_{i,k}(p), \qquad
S(p) = S_{i,k} + S'_{i,k}(p),
$$ (vert-recon)

where the deviations are polynomials in pressure built from the neighboring layer means and
constrained to **integrate to zero over the layer**
($\int_{p_{i,k}^{\text{top}}}^{p_{i,k}^{\text{bot}}} \Theta'_{i,k}\, dp = 0$). Because pressure
and $\tilde z$ differ within a column only by a fixed linear map [](#p-linear), a polynomial in
one is a polynomial of the same degree in the other; pressure is used here because it is the
variable in which the profiles of two neighboring columns can be compared (§3.7).

The mean-preserving property (Requirement 2.4) guarantees the PGF uses the same layer-mean state
as the rest of the model; the reconstruction only supplies the sub-layer shape needed to integrate
$\alpha p$.

**Phase 1 uses linear deviations; Phase 2 uses parabolic (PPM-style) ones.** The degree sets the
scheme's exact set in §3.7.3 directly: linear deviations make the scheme exact for profiles that
vary linearly with pressure, parabolic ones for profiles that vary quadratically.

That correspondence only holds if the reconstruction actually reproduces those profiles, which
places two requirements on how the deviations are built:

- **The estimator must be exact on a non-uniform vertical grid.** Omega's layers are not of equal
  thickness, and a slope or curvature formula derived assuming equal thickness will not recover a
  linear profile exactly when they are not. The estimator must be built from the layer means and
  the actual interface pressures, and must return the exact slope (Phase 1) or slope and curvature
  (Phase 2) when the underlying profile is of that degree, for any distribution of layer
  thicknesses. This is testable directly and is part of the unit test in §5.2.
- **Any limiter must be inactive on such profiles.** Monotonicity limiters are not required by the
  PGF and are not applied by default here, since the reconstruction feeds an integral rather than
  an advective flux. If one is later added for robustness, it must leave smooth monotone data
  untouched, or it will break the cancellation of §3.7.2 precisely where the profile is well
  resolved.

The isopycnal-model assumption of constant properties within a layer is the degenerate choice
$\Theta'_{i,k} = S'_{i,k} = 0$. This is retained as a **verification-only** configuration — it is
what makes the reduction to the centered scheme in §3.9 possible, and it isolates compressibility
on its own in §5.2 — but it is not a supported production setting in either phase, because it
gives up the exactness for linearly varying profiles that Requirement 2.3.3 asks for.

### 3.5 Analytic layer integral of the side-wall term

Combining [](#p-linear), [](#alpha-taylor), and [](#vert-recon), the side-wall integrand $\alpha p$
within a single column and layer is a polynomial in pressure. With the first-order equation-of-state
expansion, $\alpha$ inherits the degree of the reconstruction: linear in Phase 1, parabolic in
Phase 2. Multiplying by $p$ gives a quadratic (Phase 1) or cubic (Phase 2) integrand, and in both
cases the layer integral

$$
\Pi_{i,k} \equiv \frac{1}{g}\int_{p_{i,k}^{\text{top}}}^{p_{i,k}^{\text{bot}}} \alpha\, p \; dp
\;=\; \frac{\Delta p_{i,k}}{g} \sum_{q} w_q \, \alpha\!\left(p_q\right) p\!\left(p_q\right),
$$ (sidewall-int)

is evaluated **exactly** by a two-point Gauss–Legendre rule, which is exact through cubic. The nodes
$p_q$ and weights $w_q$ are on the layer's pressure interval, and $\alpha(p_q)$ comes from
[](#alpha-taylor). No TEOS-10 calls occur inside the sum.

The quadrature is configurable, but it is not a free accuracy knob: it must be **at least** exact for
the reconstructed integrand, or condition 2 of §3.7.2 fails and the cancellation is lost. Raising it
beyond that changes nothing, since the integrand is a polynomial of known degree. Lowering it — to
midpoint, say — is meaningful only as the verification configuration of §3.9.

### 3.6 The edge operator

The discrete tendency lives at edge $e$ and is the edge-normal projection of the PGF. The side-wall
line integral $\int_{\partial A}(\cdots)\,dl$ in [](#ho-target) becomes, on the TRiSK C-grid, the
difference of the two adjacent columns' contributions across the edge, divided by the cell-center
distance $d_e$:

$$
\left[\nabla_n \Pi\right]_{e,k} = \frac{1}{d_e}\sum_{i \in CE(e)} -n_{e,i}\, \Pi_{i,k}^{(e)},
$$ (edge-grad)

where $\Pi_{i,k}^{(e)}$ is column $i$'s layer integral [](#sidewall-int) evaluated with the
equation-of-state coefficients shared across edge $e$ [](#edge-ref). **Phase 1 uses this operator as
written** — the same two-cell stencil `PressureGradCentered` uses, so the horizontal accuracy is
unchanged from today. Everything Phase 1 gains is in the vertical and in the interface terms.

Nothing is assumed constant across the edge: the horizontal contrast in $\Theta$ and $S$ between the
two columns enters at full strength through the two column integrals. What Phase 1 does *not* do is
raise the order of the horizontal difference, which is Phase 2's job.

#### 3.6.1 Raising the horizontal order without losing the cancellation (Phase 2)

The natural way to reach fourth order — reconstruct $\Theta$ and $S$ from cell centers to the edge
with a cubic stencil, then form $\alpha p$ there — **would forfeit the robustness property of §3.7**.
It is worth being explicit about this, because it is the approach a reader would reasonably assume.
The cancellation is a statement about a *pair* of columns bounding one control volume; a quantity
interpolated to the edge from four or more cells does not belong to any control volume, and there is
no reason for the pressure and geopotential terms built from it to balance. Accuracy at the cell
centers does not carry over to a cancellation at the edge between them (assumption A1, §3.7.6).

The constraint Phase 2 must satisfy instead is:

> The higher-order edge operator must be expressible as a **weighted sum of two-column pair
> contributions**, each pair built exactly as in Phase 1 — with its own shared expansion point and
> its own control volume.

Because each pair contribution is individually zero for a resting ocean in the exact set, any
weighted sum of them is zero as well, and the machine-precision property is inherited rather than
re-derived. A fourth-order edge-normal difference built from the nearest and next-nearest cells
along the normal, with each cell pair contributing its own Phase 1 integral, satisfies this; a
scheme that reconstructs to the edge first does not.

Constructing such a stencil on Omega's unstructured TRiSK mesh is the **principal open design
question of Phase 2** and is deliberately not settled here. The relevant machinery — the wider edge
neighborhoods used for high-order tracer reconstruction, cf.
{ref}`omega-design-governing-eqns-omega1` §10 and White & Adcroft (2008) — exists, but its accuracy
on variable-resolution meshes and the cost of the extra pair evaluations both need assessment before
the form is fixed. Phase 1 does not depend on the answer.

### 3.7 Where the cancellation is exact, and where it is not

This section makes Requirement 2.3 concrete: it explains why the obvious version of the robustness
property cannot be delivered, states the three conditions the implementation must meet, lists the
water columns that cancel exactly under each phase, and records what the design is still assuming
and must therefore test. It does *not* rely on subtracting a background profile.

#### 3.7.1 Why "uniform $\Theta$ and $S$" is the wrong condition

The natural way to state the robustness property is to require the scheme to return zero whenever
the layer means $\Theta_{i,k}$, $S_{i,k}$ are the same in both cells of an edge, for any tilt and
any vertical profile. **That requirement cannot be met, and asking for it would send the
implementation chasing an impossible target.**

When coordinate surfaces tilt relative to surfaces of constant pressure, layer $k$ in one column and
layer $k$ in its neighbor cover *different pressure ranges*. Each column builds its own
specific-volume profile from its own layer means on its own set of interfaces. Unless the vertical
reconstruction happens to reproduce the true profile exactly, the two columns therefore describe
slightly *different water*. Real water columns that differ do exert a real pressure gradient on one
another. A scheme that returned zero there would be hiding an error, not avoiding one.

The same argument runs the other way, and this is the more useful half: take a perfectly uniform
ocean whose temperature and salinity depend only on pressure, and sample it onto tilted layers. The
resulting layer means are **different in the two columns**, because the layers average over
different pressure ranges. So "uniform layer means" is neither necessary nor sufficient for the
cancellation we want. The condition that matters is a property of the reconstructed profiles, not of
the layer means, and it is stated next.

#### 3.7.2 Three conditions the implementation must meet

For each edge and layer:

1. **Both sides of the edge use the same specific-volume profile.** The two columns must evaluate
   $\hat\alpha$ as one and the same function of pressure over the pressure range the layer spans.
   This is what the edge-shared expansion point of §3.3.1 delivers, and it is why that choice is
   load-bearing rather than cosmetic.
2. **Every face integral is exact for that profile.** The side-wall integrals must be exact for
   $\hat\alpha(p)\,p$ (§3.5 chooses the quadrature to guarantee this), and the sloping top and
   bottom interfaces must use the pressure-averaged $\hat\alpha p$ of §3.8, not the average of the
   two cells' interface values.
3. **The geopotential is built from that same profile.** The geometric height entering
   $g\nabla z$ must be the integral of the same $\hat\alpha$, layer by layer, down to a common
   anchor. Satisfying this requires a per-edge correction to the geometric height difference;
   see §3.7.4.

Conditions 2 and 3 are about the discretization alone and are under the implementation's control.
Meeting them means the scheme returns the **exact pressure gradient of the water column it has
reconstructed** — it contributes no error of its own beyond the reconstruction. Condition 1 is
about the state: it holds when the reconstruction reproduces the true profile, and fails, by
however much the reconstruction misses, when it does not. Meeting all three gives a PGF that is
**zero to machine precision**, for any tilt, any layer thickness, and any bathymetry.

The point worth carrying away, because it is what makes this tractable, is that **the reconstructed
profile does not have to be accurate for the cancellation to be exact — it only has to be shared.**
The scheme needs the two columns to describe *one* water column, not the *right* one. Accuracy is a
separate concern, addressed by reconstruction order.

#### 3.7.3 Which water columns cancel exactly

Condition 1 holds whenever the true vertical profile is one the reconstruction reproduces exactly:
two mean-preserving reconstructions of such a profile, built on *different* sets of interfaces from
*different* layer means, recover the identical profile anyway. Together with the edge-shared
expansion point [](#edge-ref) this gives, for each phase:

| Vertical profile of $\Theta$, $S$ | `Centered` | Phase 1 (linear reconstruction) | Phase 2 (parabolic reconstruction) |
|---|---|---|---|
| Uniform ($\alpha$ constant) | exact | exact | exact |
| Uniform $\Theta$, $S$; compressibility only | $O(\tilde h)$ | **exact** | **exact** |
| Linear in pressure | $O(\tilde h)$ | **exact** | **exact** |
| Quadratic in pressure | $O(\tilde h)$ | $O(\tilde h^2)$ | **exact** |
| Realistic (general smooth) | $O(\tilde h)$ | $O(\tilde h^2)$ | $O(\tilde h^3)$ |
| Horizontal structure in $\Theta$, $S$ | $O(d_e^2)$ | $O(d_e^2)$ | $O(d_e^4)$ |

"Exact" means zero to machine precision. Orders are in layer thickness $\tilde h$ at fixed tilt, and
in cell spacing $d_e$ for the last row. The $\tilde h$ entries carry a factor set by how strongly
the profile curves and by how far the interfaces are displaced between the two columns, so they
shrink as either shrinks.

Two consequences deserve emphasis:

- Phase 1 is exact for profiles that are **linear in pressure**, not merely for vertically uniform
  ones. That distinction is the whole reason Phase 1 is worth doing on its own: a vertically uniform
  ocean does not exist, but linear-in-pressure captures compressibility exactly and smooth
  stratification to leading order, which is most of what a real column looks like over a single
  layer. It satisfies Requirement 2.3.3.
- `PressureGradCentered` is only first order in layer thickness for *every* stratified profile,
  including a horizontally uniform one.

That last point is worth showing, since it identifies the error Phase 1 targets. Take specific
volume uniform in the horizontal within each layer, and let the coordinate merely redistribute
thickness between the two columns ($\sum_j \Delta \tilde h_j = 0$). The centered scheme's
edge-normal Montgomery-potential difference is then

$$
\Delta \bar M_k = g \rho_0 \sum_{j>k} \left( \alpha_j - \alpha_k \right) \Delta \tilde h_j,
$$ (centered-error)

which is nonzero whenever $\alpha$ varies with depth and interfaces are displaced, and which
**accumulates downward**: the deeper the layer, the more layers contribute to the sum. The vertical
contrast in *in-situ* $\alpha$ is dominated by compressibility (a few percent over the full ocean
depth, against a few tenths of a percent from $\Theta$ and $S$), and that is exactly the part
Phase 1 cancels.

This first-order behaviour has since been **measured**, and it is no longer an assertion. On the
Polaris `horiz_press_grad` resting-state configurations (uniform $\Theta$, $S$, flat floor, tilted
coordinate), `PressureGradCentered` gives a fitted exponent of $1.0000$ in the coordinate tilt at
three vertical resolutions, and the `ztilde_gradient` variant converges at $\approx 1.1$ in
horizontal resolution once the bottom layer is included in the comparison. Two independent
measurements of the same first-order behaviour, in agreement. The corresponding absolute errors reach
$2\times10^{-5}\ \mathrm{m\,s^{-2}}$ at a coordinate tilt of 50 m/km with 256 m layers — the order of
the bottom-layer error seen in realistic global configurations. Whether this downward accumulation is in fact what produces the bottom-layer noise
seen in realistic global configurations is a plausible diagnosis, not an established one; it is
carried as A4 in §3.7.6 and tested in §5.3.

#### 3.7.4 Condition 3 and the geopotential

Condition 3 constrains how the geopotential term is built. Two questions arise, and they have
different answers; an earlier revision of this design conflated them and drew the wrong conclusion
from the first.

**The quadrature question, which is already settled.** `VertCoord::computeGeomZHeight` builds $z$ by
accumulating $\Delta z = \rho_0\,\alpha_{i,k}\,\tilde h_{i,k}$ — apparently a midpoint rule, where
condition 3 asks for the integral of the reconstructed $\hat\alpha$. For a Phase 1 reconstruction
these are **the same quantity**. Integrating [](#alpha-taylor) with the linear deviations of §3.4
over the layer,

$$
\frac{1}{g}\int_{p^{\text{top}}_{i,k}}^{p^{\text{bot}}_{i,k}} \hat\alpha_{i,k}(p)\,dp
= \frac{\alpha_{0}\,\Delta p_{i,k}}{g} = \rho_0\,\alpha_{i,k}\,\tilde h_{i,k},
$$ (z-increment-exact)

because $\int \Theta'\,dp = \int S'\,dp = 0$ by the mean-preserving constraint, and
$\int (p - p^{\text{mid}})\,dp = 0$ because $p^{\text{mid}}$ is the exact arithmetic midpoint of the
two interface pressures. The midpoint rule *is* the exact layer integral of a Phase 1
$\hat\alpha$. No change to `VertCoord` is required, there is no answer-changing baseline step, and
Requirement 2.4 is satisfied by construction rather than by negotiation — the PGF and the rest of
the model share one $z$ because they compute the same thing.

This is a Phase-1-only result and Phase 2 must re-examine it. Parabolic deviations are fine provided
they remain mean-preserving, but the second-order equation-of-state expansion of §3.3 contributes
$\tfrac12\alpha_{pp}(p-p^{\text{mid}})^2$ and cross terms such as
$\tfrac12\alpha_{\Theta\Theta}\Theta'^2$, none of which integrate to zero over the layer. If that
option is adopted, [](#z-increment-exact) no longer holds and the question returns.

**The sharing question, which is the real constraint.** What condition 3 actually requires is that
the geopotential be built from the *same* $\hat\alpha$ the pressure terms use — and §3.3.1 makes
that an **edge** quantity, since its coefficients and expansion point are averages over the two
cells of the edge. A cell-based $z$ cannot carry an edge-dependent $\hat\alpha$ however it is
integrated. Two mismatches follow, both $O(\alpha_{pp}(\Delta p^{\text{mid}}_e)^2)$ per layer and
both accumulating down the column:

1. `VertCoord`'s $z$ uses the exact per-cell TEOS-10 $\alpha$; the PGF uses the edge-shared
   linearized $\hat\alpha$. The residual is the
   $\tfrac12\alpha_{pp}(p^{\text{mid}}_i - \bar p^e)^2$ term.
2. The geopotential body force is the gradient of $\Phi$ **at fixed pseudo-height**, i.e. at fixed
   pressure. Differencing each column's layer-mean $\Phi$ across the edge compares layer means taken
   over *different* pressure ranges whenever the interfaces tilt, and that difference does not
   vanish for a resting ocean.

**Resolution: the PGF computes a per-edge correction to the `VertCoord` $z$ difference**, built from
the edge-shared $\hat\alpha$ and accumulated down the column. Stated as a correction rather than a
replacement, this does not create the second, slightly different $z$ that Requirement 2.4 exists to
prevent: `GeomZInterface`/`GeomZMid` remain the model's one geometric height and are unchanged, and
the PGF adds a term to a *difference* that is identically zero when the two columns' interfaces are
level. Accumulating only the correction, rather than re-deriving $z$, also preserves most of the
precision headroom §3.7.5 warns is otherwise consumed, since the correction is small by
construction.

The implementation consequence is recorded in §4.1.3: the correction is a column prefix sum and
therefore cannot live inside a per-edge, per-vertical-chunk kernel call.

One related question is left to implementation time: `VertCoord` builds $z$ upward from the
bathymetry while pressure is built downward from the surface, so the two accumulate round-off from
opposite ends of the column. Which end the correction accumulates from is a round-off question
(§3.7.5), not a consistency one — condition 3 is satisfied either way.

#### 3.7.5 Round-off in the deep ocean

Separate from everything above, and not fixed by it, is how much precision the cancellation itself
consumes. At 4000 m the side-wall integrand $\hat\alpha p$ is about
$4\times10^{4}\ \mathrm{m^2\,s^{-2}}$, while the baroclinic signal that survives in the edge
tendency is of order $10^{-6}\ \mathrm{m\,s^{-2}}$ times $d_e$. Roughly ten significant digits are
consumed before any physics appears. Double precision leaves adequate margin; a
`OMEGA_SINGLE_PRECISION` build does not, and the machine-precision cancellation of §3.7.2 would
simply be invisible beneath round-off.

If that turns out to matter, the fix is to compute in **perturbation form**: subtract a local
reference profile from $\hat\alpha$ before integrating and add its contribution back analytically.
This is *not* a Shchepetkin–McWilliams (2003)-style reference-profile subtraction in the usual
sense. The reference here is local to the edge and layer, its contribution cancels identically
rather than approximately, and the accuracy of the scheme does not depend at all on how well it
matches the local column — so it does not degrade in the strong-gradient, steep-layer regions where
a global reference profile would. It changes nothing in exact arithmetic. Whether it is needed is a
measurement (§5.2, run in both precisions), not an assumption.

#### 3.7.6 Assumptions that still need testing

Everything above concerns properties the implementation can be built to have. The following are
assumptions this design is making that only testing can confirm:

- **A1 — Accuracy at cell centers does not imply accuracy at edges.** The cancellation above is a
  statement about a *pair* of columns at an edge. Reconstructing $\Theta$ and $S$ accurately at each
  cell center does not by itself make the combination $\hat\alpha(\Theta,S,p)\,p$ accurate at the
  edge between them, because that combination is nonlinear. The horizontal order claimed in §3.6 is
  therefore a design target to be confirmed by measured convergence (Requirement 2.6), not something
  that follows from the cell-centered reconstruction order.
- **A2 — The EOS expansion is good enough across the edge.** [](#alpha-taylor) is expanded about the
  shared edge state [](#edge-ref), and its error grows with the horizontal contrast in $\Theta$ and
  $S$ across the edge and with the layer's pressure range. Whether the first-order expansion suffices
  in frontal regions, or whether the second-order option of §3.3 is needed, is an open question.
- **A3 — The residual outside the exact set is small enough in practice.** The $O(\tilde h^2)$ and
  $O(\tilde h^3)$ entries in §3.7.3 describe how the error *scales*; how large it actually is at the
  vertical resolutions Omega can afford is unknown.
- **A4 — PGF error is what is driving the observed instability.** Identifying [](#centered-error)
  with bottom-layer noise in realistic global runs is a plausible diagnosis, not a demonstrated one.
  A configuration whose profile falls in Phase 1's exact set gives a direct test: if spurious
  bottom-layer flow persists there, the cause lies elsewhere — most likely in the layer-mean
  treatment in the tracer and remapping operators (§4.3) — and Phase 1 will not cure it. See §5.3.

### 3.8 Geopotential and metric terms

The geopotential body-force term (first line of [](#ho-target)) is
$g\,\nabla z + \nabla(\phi_{TP}+\phi_{SAL})$. The geometric height $z$ and geopotential are
provided by `VertCoord` (`GeomZMid`/`GeomZInterface`, `GeopotentialMid`; computed via
{ref}`omega-design-governing-eqns-omega1` Eqs. discrete-z and the geopotential relation).
Its layer average and edge-normal gradient are evaluated with the same high-order edge
reconstruction [](#edge-grad). The tidal-potential and self-attraction-and-loading
contributions enter through `VertCoord` and are differenced identically (Requirement 2.7);
they are zero in early Omega versions.

The two metric terms (third and fourth lines of [](#ho-target)) are the pressure force on the
layer's sloping top and bottom interfaces. They use the **same** edge-normal operator to evaluate
the interface slopes $\nabla \tilde z_k^{\text{top}}$, $\nabla \tilde z_k^{\text{bot}}$, but the
quantity multiplying each slope needs care. Because the interface is sloping, it sweeps through a
range of pressures between the two cell centers, and the force is an integral along it. The correct
factor is therefore the **average of $\hat\alpha p$ over that pressure range**,

$$
\left\langle \hat\alpha p \right\rangle_{e,k}^{\text{top}}
= \frac{1}{p_{R,k}^{\text{top}} - p_{L,k}^{\text{top}}}
  \int_{p_{L,k}^{\text{top}}}^{p_{R,k}^{\text{top}}} \hat\alpha(p)\, p \; dp,
$$ (metric-divdiff)

and likewise at the bottom interface, reducing to
$\hat\alpha\!\left(p^{\text{top}}\right) p^{\text{top}}$ when the interface is level in the two
columns. Because $\hat\alpha(p)$ is a low-order polynomial (§3.3), this integral is evaluated in
closed form at no meaningful cost.

**Implement this as an integral, not as an average multiplied by a slope.** Since
$\nabla_n\tilde z^{\text{top}}_{e,k} = -(p^{\text{top}}_{R} - p^{\text{top}}_{L})/(\rho_0 g\,d_e)$,
the divisor in [](#metric-divdiff) and the slope it multiplies cancel identically, and the whole
term is

$$
\frac{1}{g\,d_e}\int_{p^{\text{top}}_{L,k}}^{p^{\text{top}}_{R,k}} \hat\alpha(p)\, p \; dp .
$$ (metric-integral)

Forming the average literally divides by a quantity that is zero wherever the interfaces are level,
which is most of the domain; the reduction noted above is a $0/0$ in that form and is well defined
only in [](#metric-integral).

The tempting alternative — averaging the two cells' interface values,
$\tfrac12[(\hat\alpha p)_L + (\hat\alpha p)_R]$ — agrees with [](#metric-divdiff) only when the
tilt is small, and leaves behind a residual that grows with the square of the tilt and does *not*
vanish for the water columns §3.7.3 says should cancel exactly. It would therefore break
condition 2 of §3.7.2 and forfeit the robustness property in exactly the steeply sloped regions
this design is meant to fix.

This integral is what `PressureGradCentered` omits altogether, and supplying it is the central
content of Phase 1. These terms vanish for level interfaces and supply the leading correction for
sloping layers; together with the side-wall integral [](#sidewall-int) and the geopotential term
they complete the force balance on the layer's control volume (§3.1.1).

### 3.9 Reduction to the centered scheme

The new scheme collapses to the implemented `PressureGradCentered` form
($T^p_{e,k} = -\nabla M + \tfrac12(p_0+p_1)\nabla\alpha - \nabla\Phi$, with $M = \alpha p + g z$)
in one configuration: **specific volume constant within each layer**, meaning both
$\Theta' = S' = 0$ and the equation-of-state expansion [](#alpha-taylor) truncated to $\alpha_0$,
combined with the two-cell edge operator [](#edge-grad). This is the isopycnal-model assumption,
and it is the verification-only mode of §3.4.

Two things fall out of that configuration that are worth noting, because they show the reduction is
structural rather than a special case bolted on:

- With $\alpha$ constant in the layer, $\alpha p$ is linear in pressure, so the sloping-interface
  average [](#metric-divdiff) reduces *exactly* to the average of the two cells' interface values.
  The interface term does not need to be switched off to recover the centered scheme; it degenerates
  on its own.
- For the same reason, the layer integral [](#sidewall-int) is exact under any symmetric quadrature
  rule, midpoint included, so the quadrature setting does not enter the reduction either.

What remains is the Montgomery-potential algebra, and the two forms agree. This gives a direct
regression path (§5.5) and confirms the new code reproduces the existing scheme where it should.
The agreement is algebraic, not bit-for-bit — the operations are performed in a different order
(§4.4).

### 3.10 Per-step algorithm summary

Steps 1–3 are per cell and layer; step 4 is per edge and layer. Phase differences are marked.

1. From `VertCoord`: read `PressureInterface`, `PressureMid`, `GeomZInterface`/`GeomZMid`,
   geopotential, and interface pseudo-heights (already computed diagnostically each step). The
   geometric height must satisfy condition 3 of §3.7.2, which requires the per-edge correction of
   §3.7.4; `VertCoord` itself needs no change.
2. Obtain $\alpha_0$ (= the existing `Eos::SpecVol` field) and the derivatives
   $\alpha_\Theta, \alpha_S, \alpha_p$ from one TEOS-10 evaluation ([](#alpha-derivs)). Both phases
   need all four; Phase 2 optionally adds second derivatives (§3.3).
3. Build the mean-preserving deviations $\Theta', S'$ ([](#vert-recon)) — **linear in Phase 1,
   parabolic in Phase 2** — using the actual non-uniform interface pressures (§3.4).
4. For each edge: form the shared expansion point [](#edge-ref) from the two adjacent cells;
   evaluate each column's layer integral [](#sidewall-int) with those shared coefficients; add the
   sloping-interface terms [](#metric-integral) and the corrected geopotential difference (§3.8,
   §3.7.4); assemble
   $T^p_{e,k}$ and accumulate into the tendency with `EdgeMask`. **Phase 1** uses the two-cell
   operator [](#edge-grad); **Phase 2** uses the wider stencil of §3.6.1, built as a weighted sum of
   such two-cell pair contributions.

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

New sub-options for the high-order scheme. Both phases use the same keys; Phase 2 adds values
rather than keys, so no configuration written for Phase 1 needs to change when Phase 2 lands:

```yaml
    PressureGrad:
       PressureGradType: 'FiniteVolume'   # Centered | FiniteVolume
       HorzOrder: 2                       # 2 = two-cell stencil (Phase 1); 4 = wide stencil (Phase 2)
       VerticalReconstruction: 'linear'   # 'linear' (Phase 1) | 'ppm' (Phase 2)
                                          #  | 'constant' (verification only, see below)
       QuadraturePoints: 2                # per-layer Gauss points for the side-wall integral
```

- **Phase 1 default and target:** `HorzOrder: 2`, `VerticalReconstruction: 'linear'`.
- **Phase 2 target:** `HorzOrder: 4`, `VerticalReconstruction: 'ppm'`.
- `QuadraturePoints: 2` is exact for the integrand in both phases (§3.5) and should not normally be
  changed. It is a knob because lowering it is needed for the verification configuration below, not
  because raising it buys accuracy.

`VerticalReconstruction: 'constant'` sets the specific volume constant within each layer — both
$\Theta' = S' = 0$ and the equation-of-state expansion truncated to $\alpha_0$. It is **verification
only**: combined with `HorzOrder: 2` it recovers `PressureGradCentered` (§3.9) and supports the
permanent regression test of §5.5. It is not a supported production setting, because it gives up the
exactness for linearly varying profiles that Requirement 2.3.3 asks for. The implementation should
log a warning if it is selected outside a test.

`'constant'` is **not** the compressibility-only configuration, and the two are easy to confuse.
Truncating the expansion to $\alpha_0$ discards $\alpha_p$, which is exactly the term that makes the
scheme exact for the "uniform $\Theta$, $S$; compressibility only" row of §3.7.3. The configuration
that isolates compressibility is `'linear'` applied to a vertically uniform $\Theta$, $S$ profile:
the reconstruction slopes are then zero on their own, $\alpha_p$ is retained, and no separate
setting is needed. Tests that mean to isolate compressibility — including §5.2 — must select
`'linear'`, not `'constant'`.

`HorzOrder` selects the width of the edge *stencil* — how many cell pairs contribute — not the order
of an interpolation of $\Theta$ and $S$ onto the edge. The distinction is not cosmetic; see §3.6.1.

#### 4.1.2 New EOS support

`PressureGradHighOrder` needs $\alpha$ together with its first derivatives. The `Eos` class
(Eos.h) currently exposes `computeSpecVol`, `computeSpecVolDisp`, and
`computeBruntVaisalaFreqSq`, but no specific-volume derivatives. This design adds one method
and three device arrays:

```c++
// New on the Eos class
Array2DReal SpecVolDCt; ///< d(alpha)/d(ConservTemp) at cell centers
Array2DReal SpecVolDSa; ///< d(alpha)/d(AbsSalinity) at cell centers
Array2DReal SpecVolDP;  ///< d(alpha)/d(Pressure)    at cell centers

// Pressure is relative pressure in Pa; the derivatives are returned per degC,
// per (g/kg), and per Pa respectively, into the members above and SpecVol.
void computeSpecVolAndDerivs(const Array2DReal &ConservTemp,
                             const Array2DReal &AbsSalinity,
                             const Array2DReal &Pressure);
```

`Eos` owns the derivative arrays exactly as it owns `SpecVol`, `SpecVolDisplaced` and
`BruntVaisalaFreqSq`: they are allocated in the constructor and registered as `Field`s in the
`Eos` group, so they can be written to a stream and any number of consumers can read them
without each allocating its own copy. Their valid range must permit negative values —
$\alpha_S < 0$ everywhere, and $\alpha_\Theta < 0$ in cold fresh water — unlike `SpecVol`.

Because `computeSpecVolAndDerivs` fills `SpecVol` as well, it replaces rather than accompanies
a call to `computeSpecVol`; the two are kept separate so that the derivative arithmetic is
paid only where it is needed.

The TEOS-10 derivatives are obtained by differentiating the 75-term polynomial analytically.
The pressure derivative reuses the coefficients `calcPCoeffs` already assembles for
$\alpha$ itself and so is free; the $\Theta$ and $S_A$ derivatives need their own coefficient
sets, built from the same $s$ and $t$ (no additional square root). The marginal cost over
`computeSpecVol` is therefore roughly two extra coefficient assemblies and no additional
TEOS-10 evaluation, which is what Requirement 2.2 bounds. The linear and constant EOS
options supply trivial analytic derivatives.

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
                                const Array2DReal &ConservTemp,  // new
                                const Array2DReal &AbsSalinity,  // new
                                const Array2DReal &SpecVolDCt,   // new
                                const Array2DReal &SpecVolDSa,   // new
                                const Array2DReal &SpecVolDP)    // new
                                const;
```

(The existing centered signature is unchanged.) Additional cached members hold the quadrature nodes
and weights and, in Phase 2, the wide-stencil cell lists and weights. The functor implements
§3.3–§3.8 per edge and vertical chunk, accumulating into `Tend` with `EdgeMask`, exactly as the
centered functor does.

Two aspects of the loop structure follow from §3.3.1 and are worth stating here, because they differ
from the obvious implementation:

- **The layer integral is computed inside the edge loop, not cached per cell.** Both adjacent
  columns' integrals [](#sidewall-int) use equation-of-state coefficients formed at the edge
  [](#edge-ref), so the same cell yields a different integral at each of its edges. The per-cell
  quantities that *can* be cached are the four EOS coefficients and the reconstruction slopes, which
  is where the TEOS-10 cost lives; the per-edge work is polynomial arithmetic on those cached
  values. This is what keeps Requirement 2.2 satisfied despite roughly three times as many integral
  evaluations on a hexagonal mesh.
- **Phase 2's wide stencil is a loop over cell pairs, not a wider interpolation.** Each pair
  contributes a complete Phase 1 evaluation with its own shared expansion point, and the pair
  results are combined with the stencil weights (§3.6.1). Implementing Phase 2 as a wider
  reconstruction feeding a single evaluation would be simpler and would break the property of §3.7.
- **The geopotential correction of §3.7.4 is a column scan and cannot live in this functor.** It is a
  prefix sum down each column with edge-dependent coefficients, so it is not expressible as an
  independent per-vertical-chunk operation. It is computed in a separate kernel that fills a per-edge
  array,

  ```c++
  Array2DReal GeopotCorrection;  ///< (NEdgesAll, NVertLayers), owned by PressureGrad
  ```

  with a `parallelForOuter` over edges and a `parallelScanInner` down the column, in the same shape
  as `VertCoord::computeGeomZHeight`. The functor then reads it chunk-wise like any other input. The
  cost is one edge-sized 2-D array and one column scan per edge per step. This is the one structural
  addition Phase 1 makes beyond the per-edge, per-chunk pattern the centered scheme uses.

The Phase 1 and Phase 2 code paths differ only in the reconstruction degree (§3.4) and in whether
the pair loop has one entry or several. There is one functor, not two.

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

### 4.4 Retention of the centered implementation

`PressureGradCentered` (PGrad.h:25) is deliberately **kept as a separate functor** rather
than reimplemented as the lowest-order configuration of `PressureGradHighOrder`, even though
§3.9 shows the latter reduces to it. The redundancy is small — the centered functor is a
header-only, ~40-line loop body with no supporting machinery of its own — and it buys two
things that a single implementation cannot provide:

1. **An independent cross-check.** The two functors read the mesh, `VertCoord`, and EOS state
   through separately written code. Their agreement to round-off (§5.5) therefore tests the
   shared upstream state — edge masks, interface indexing, `VertCoord` conventions — and not
   just the PGF arithmetic. Collapsing them would make that comparison self-referential: a
   defect upstream of the order switch would appear identically in both limits and cancel.
2. **A stable default.** The algebraic reduction of §3.9 does not imply bit-for-bit agreement,
   because the reduced high-order path performs the same operations in a different order.
   Replacing the default PGF would be an answer-changing change for every existing
   configuration, which this design does not require.

Removing `PressureGradCentered` becomes reasonable once `FiniteVolume` is promoted to the
default and has served out its period as the reference implementation. That is deliberately
left as **follow-up work**, to be taken up as a separate, answer-changing change with its own
baseline step — not as part of this design.

### 4.5 What each phase delivers and depends on

#### 4.5.1 Phase 1 — a consistent second-order scheme

**Delivers.** The finite-volume control-volume form (§3.1.1); the sloping-interface integral
(§3.8), which the centered scheme omits entirely; the equation-of-state expansion about a state
shared across each edge (§3.3.1); mean-preserving linear reconstruction of $\Theta$ and $S$ in
pressure (§3.4). The result is a pressure gradient that is zero to machine precision for any
resting ocean whose profile varies linearly with pressure, at any tilt, thickness, or bathymetry
(§3.7.3) — Requirement 2.3 in full.

**Does not deliver.** Fourth-order accuracy. The horizontal operator is the same two-cell stencil
in use today (§3.6), so horizontal truncation error is unchanged from `PressureGradCentered`.
Requirements 2.1 and 2.6 are met at second order only.

**Depends on.** Nothing unresolved. The `VertCoord` geopotential question of §3.7.4 is settled
there: no change to `VertCoord` is needed and no baseline step is required, because the midpoint rule
is already the exact layer integral of a Phase 1 $\hat\alpha$. What condition 3 does require — the
per-edge geopotential correction — is part of this phase's own implementation (§4.1.3), not a
prerequisite in another module.

**Code and cost.** Three new `Eos` derivative fields and one new method (§4.1.2); the
`PressureGradHighOrder` functor; the per-edge geopotential correction array and its column scan
(§3.7.4, §4.1.3); no new TEOS-10 evaluations per cell and layer (Requirement 2.2), with roughly
three times as many polynomial layer integrals on a hexagonal mesh (§4.1.3).

#### 4.5.2 Phase 2 — fourth order

**Delivers.** Parabolic vertical reconstruction (§3.4), widening the exact set to profiles that vary
quadratically with pressure; a wide horizontal stencil (§3.6.1); optionally a second-order
equation-of-state expansion (§3.3). Requirements 2.1 and 2.6 at fourth order.

**Must preserve.** Everything Phase 1 establishes. In particular the machine-precision property must
survive the wider stencil, which is why §3.6.1 constrains that stencil to be a weighted sum of
two-column pair contributions rather than a reconstruction to the edge. The §5.2 gate is rerun
unchanged for Phase 2.

**Open question.** The form of the wide stencil on Omega's unstructured, variable-resolution TRiSK
mesh is not settled by this design (§3.6.1). Resolving it — including the cost of the additional
pair evaluations — is the first task of Phase 2 and does not block Phase 1.

#### 4.5.3 Suggested order of work

1. Write the discrete form out in full — including the interface-metric sign convention, which
   $\hat\alpha$ the interface integral uses when the two columns' reconstructions differ, and the
   exact form of the per-edge geopotential correction (§3.7.4) — and confirm numerically that it
   returns zero for a profile linear in pressure at large tilt. Doing this in the Polaris
   two-column harness makes it executable and doubles as the reference implementation §5.1 needs.
2. Run the assumption-A4 diagnostic of §5.3, which uses the *existing* centered scheme and so can be
   done immediately and in parallel with step 1. If spurious bottom-layer flow survives a profile
   that Phase 1 would resolve exactly, the cause is elsewhere in the model and the priority of this
   work should be reconsidered before it is built.
3. Implement and verify Phase 1 against §5.1, §5.2, §5.3, and §5.5.
4. Take up Phase 2, starting from the stencil question in §3.6.1. Re-examine [](#z-increment-exact)
   before adopting the second-order equation-of-state expansion.

## 5 Verification and Testing

Testing reuses and extends the Polaris `horiz_press_grad` task family
(`polaris/tasks/ocean/horiz_press_grad`), which is already Omega-only and built around a
quasi-analytic, layer-mean TEOS-10 reference solution (`reference.py`, surface-anchored,
4-point Gauss). That analytic reference remains valid as "truth" for the high-order scheme,
so the principal changes are scheme selection and revised pass criteria, plus two full-model
acceptance tests.

The testing plan applies the two independent gates set out in §1: an **absolute error
tolerance at a representative coarse resolution** (Requirement 2.1), which establishes that
the scheme helps at the resolutions Omega can afford, and a **measured order of convergence**
(Requirement 2.6), which verifies that the implementation is the scheme this design
describes. Neither substitutes for the other, and both must pass.

Which tests gate which phase:

| Test | Phase 1 | Phase 2 |
|---|---|---|
| §5.1 Two-column convergence | gates, at second order | gates, at fourth order |
| §5.2 Machine-precision cancellation | gates | rerun unchanged; must still pass |
| §5.3 Seamount resting state | gates | rerun; error should drop further |
| §5.4 Overflow | gates | rerun |
| §5.5 Reduction to centered | gates | rerun unchanged |
| §5.6 Cost check | gates | gates; stencil width must not change the EOS count |

The A4 diagnostic within §5.3 is run *before* Phase 1 implementation and uses the existing centered
scheme, so it gates nothing but informs whether the work should proceed as prioritized (§4.5.3).

### 5.1 Test: Two-column HPGA convergence (extend existing)

Extend the four existing variants — `temperature_gradient`, `salinity_gradient`,
`surface_pressure_gradient`, `ztilde_gradient` — to run the centered scheme alongside the new one:

- **Scheme selection.** Add `PressureGrad: { PressureGradType: FiniteVolume, … }` to
  `forward.yaml` and parametrize each task over three configurations:
  - `centered` — the legacy `PressureGradCentered` functor, unchanged;
  - `finite_volume_phase1` — `HorzOrder: 2`, `VerticalReconstruction: linear`,
    `QuadraturePoints: 2`;
  - `finite_volume_phase2` — `HorzOrder: 4`, `VerticalReconstruction: ppm` (added when Phase 2
    lands).

  (These variant names are provisional; the final spelling follows Polaris' naming conventions and
  is settled on the Polaris side. What matters here is that the last two are the *same*
  implementation at two orders, distinct from the legacy functor.)

  `finite_volume_phase1` is not a stand-in for the legacy scheme: it is second order in the
  horizontal like `centered`, but it is *consistent*, so its absolute error should be markedly
  lower even though its convergence slope is the same. Both must be run, and the comparison between
  them is the clearest single measure of what Phase 1 buys. Running only the legacy functor under a
  "second order" label would measure the convergence of code this design does not change and leave
  the new code's Phase 1 path unverified.

  The forward step still runs a single time step with only `PressureGradTendencyEnable: true`,
  reading the PGF acceleration from `NormalVelocityTend`.
- **Reference.** `reference.py`/`analysis.py` compare `NormalVelocityTend` against the
  layer-mean analytic HPGA (unchanged). For the high-order scheme the layer-mean comparison
  remains the correct target, since the scheme is a finite-volume, layer-averaged
  discretization.
- **Accuracy gate (new, Requirement 2.1):** at a representative coarse resolution (e.g. the
  coarsest in `horiz_resolutions`), the absolute RMS HPGA error vs. the reference must be below
  a tolerance, **and** the new scheme's RMS error must be below the centered RMS error at that
  same resolution (the scheme must demonstrably help where it matters). This gate applies to
  Phase 1, where it is the primary measure of value, since Phase 1 does not change the
  convergence slope.
- **Verification gate (Requirement 2.6):** the measured slope of RMS error vs. resolution,
  `omega_vs_reference_convergence_rate_*`, must fall within a band around the configured order
  of accuracy — nominally ~2 for `finite_volume_phase1` and ~4 for `finite_volume_phase2`. This band
  is retuned from its present values rather than loosened; a slope outside it fails the test and is
  treated as an implementation defect to be diagnosed, not as a tolerance to be widened.

  **The `centered` bands are *not* uniform across the four variants and must be set per variant from
  measurement.** Measured on chrysalis with the bottom layer included, `PressureGradCentered` gives
  $\approx 1.6$ (`temperature_gradient`), $\approx 1.8$ (`salinity_gradient`), $\approx 2.0$
  (`surface_pressure_gradient`) and $\approx 1.1$ (`ztilde_gradient`). The last is the first-order
  resting-state behaviour of §3.7.3 showing up directly, and it is the reason a single "~2 for
  centered" band would fail three variants out of four. This corrects an earlier statement here that
  the `centered` band was unchanged from its present values.
- **Asymptotic range (Phase 2 implementation-time task):** it is not yet established that the
  existing `horiz_resolutions` sweep spans a range where a fourth-order slope is cleanly
  measurable — the sweep may be too coarse to have entered the asymptotic regime at its fine end,
  or fine enough that the reference solution's own quadrature error and roundoff contaminate the
  slope. Determining the usable window, and extending or tightening the sweep (and, if needed, the
  order of the Gauss quadrature in `reference.py`) so the designed order can be resolved, is part
  of implementing Phase 2. Phase 1's second-order slope is measurable on the existing sweep.
- **Consistency check (retained):** `omega_vs_polaris_rms_threshold` (~1e-10 m/s²) — Omega's
  forward output must still match the Python-computed HPGA, confirming the implementation
  matches the intended discretization.
- **Horizontal-contrast sweep (new, assumption A2):** the `temperature_gradient` and
  `salinity_gradient` variants are run at several amplitudes of the horizontal contrast, up to and
  beyond values typical of ocean fronts. The equation-of-state expansion [](#alpha-taylor) is taken
  about a state shared across the edge, so its error grows with that contrast, and A2 (§3.7.6) is
  the assumption that the first-order expansion remains adequate. **Pass:** the error at fixed
  resolution grows no faster than linearly with the contrast amplitude, and repeating the largest
  amplitude with the second-order expansion of §3.3 changes the answer by less than the accuracy
  gate above. If it does not, the second-order expansion becomes the default rather than an option.
- **Cfg keys.** New keys mirror the existing ones (`horiz_press_grad.cfg`): a coarse-resolution
  absolute tolerance, a `finite_volume_vs_centered` ratio gate, per-scheme expected-rate bands, and
  the contrast amplitudes for the sweep.

**Covers:** Requirements 2.1, 2.2 (the bounded-EOS path is exercised), 2.4, 2.5, 2.6; the last row
of the §3.7.3 table (horizontal structure in $\Theta$, $S$); assumptions A1, A2, A3.

### 5.2 Unit test: Discrete hydrostatic consistency (exact resting state)

This is the primary test of the steep-layer robustness property (Requirement 2.3, §3.7). It follows
the table in §3.7.3 row by row rather than being a single pass/fail.

Set up a two-column (or seamount) configuration with the coordinate interfaces deliberately tilted
between the columns — including steep slopes and thin layers — and initialize the layer means as the
**exact layer averages** of a prescribed continuous profile $\Theta(p)$, $S(p)$. Under tilt those
averages come out *different in the two columns*, and that is the point: a configuration built by
copying identical layer means into both columns would not exercise the property at all (§3.7.1).
Three groups of profiles are run:

- **Profiles the scheme resolves exactly:** $\Theta$, $S$ linear in pressure, including the constant
  case, which isolates compressibility on its own. Run these with
  `VerticalReconstruction: 'linear'`, **not** `'constant'` — see §4.1.1; `'constant'` discards
  $\alpha_p$ and so cannot be exact for compressibility. **Pass:** the PGF is zero at every edge and
  layer to **machine precision** (double-precision builds; the threshold tracks `Real`'s epsilon and
  the size of the hydrostatic terms, not a physical tolerance). Phase 1 and Phase 2 must both pass.
- **Profiles it does not:** $\Theta$, $S$ quadratic in pressure, then a realistic profile.
  **Pass:** the residual shrinks like $\tilde h^2$ (Phase 1) and $\tilde h^3$ (Phase 2) as the
  vertical grid is refined at fixed tilt, matching §3.7.3. A residual that does not shrink at the
  tabulated rate means one of the three conditions in §3.7.2 has been broken somewhere in the
  implementation; it is a bug to find, not a tolerance to widen.
- **Guard tests:** rerun an exactly resolved profile with (a) the endpoint-average interface term in
  place of [](#metric-integral), (b) a cell-local expansion point in place of [](#edge-ref), and
  (c) the per-edge geopotential correction of §3.7.4 disabled. All three must *fail* the
  machine-precision check. Without these, a passing result could just as easily come from a symmetry
  of the test setup as from the scheme being right — and guard (c) in particular is the direct
  evidence that the geopotential correction is load-bearing rather than decorative.

  Guards must be checked against a deliberately broken configuration, not only against a passing
  one. A guard that cannot fire is worse than no guard, because it looks like protection.

A separate and much smaller unit test covers the reconstruction estimator on its own (§3.4): given
layer means sampled from a profile of the reconstruction's own degree on a **deliberately
non-uniform** set of layer thicknesses, the recovered slope (Phase 1) or slope and curvature
(Phase 2) must match the exact values to round-off. This is worth testing separately because it is
the most likely place for the machine-precision gate above to fail, and it localizes the failure
immediately.

The test is also run in a single-precision build, to measure the round-off floor of §3.7.5 and
settle whether the perturbation form is needed.

Condition 3 of §3.7.2 requires the per-edge geopotential correction of §3.7.4, so the
machine-precision check cannot pass until that is in place; this test, together with guard (c)
above, is its acceptance criterion. It is implemented as a fast C++ unit test and also run as a
configuration of the Polaris two-column task.

**Covers:** Requirement 2.3 in full; rows 1–5 of the §3.7.3 table; the round-off floor of §3.7.5.

### 5.3 Test: Seamount resting state (steep-layer robustness)

Use the Polaris seamount task with layers tilted over a seamount and a stratification that is a
function of pressure alone — note that this gives layer means that *differ* between neighboring
columns wherever the layers tilt (§3.7.1), which is the situation the test is meant to create.
Integrate for a fixed period. **Pass:** the maximum spurious velocity stays below a threshold and
is substantially smaller for the new scheme than for the centered scheme. This is the dynamical
counterpart of §5.2 — the same resting state, now run forward in the full solver — and the direct
test of Requirement 2.3.

**Testing assumption A4 (§3.7.6).** Run the seamount case twice: once with a stratification Phase 1
resolves exactly ($\Theta$, $S$ linear in pressure) and once with a realistic profile. Spurious
velocity that survives the first run cannot be PGF truncation error, since the PGF is zero there to
machine precision, and must come from somewhere else in the model — most likely the layer-mean
treatment in the tracer and remapping operators (§4.3). This is the cheapest available check on
whether accelerating this work will actually cure the bottom-layer instability that motivated the
phasing, and it should be run **before** Phase 1 implementation is finished, with the centered
scheme on the realistic profile as the control. A null result would not invalidate the design, but
it would change its priority.

**Covers:** Requirement 2.3 under dynamics; assumptions A3, A4.

### 5.4 Test: Overflow (full non-Boussinesq dynamics)

Use the Polaris overflow task to exercise the PGF within the full non-Boussinesq equations with
strongly sloping layers and active dynamics. Unlike the seamount task, the state here is neither
resting nor close to a profile the scheme resolves exactly, so this is the test of how the scheme
behaves where none of the exact-cancellation results apply and the residual of §2.3.4 is all there
is. **Pass:** the solution remains stable and the down-slope evolution agrees with the reference
behavior; spurious mixing and velocity attributable to PGF error are reduced relative to the
centered scheme.

This is also the most demanding test of assumption A2: the overflow front carries a large horizontal
temperature contrast across edges with steeply sloping layers, which is exactly where the shared
edge expansion point [](#edge-ref) is worked hardest. If the §5.1 contrast sweep indicates the
second-order equation-of-state expansion is needed, this task is where that judgment is confirmed
under dynamics.

**Covers:** Requirements 2.1 and 2.3 under realistic, coupled conditions; the last two rows of the
§3.7.3 table; assumptions A2, A3.

### 5.5 Test: Reduction to the centered scheme (permanent regression)

Configure the new scheme in the verification-only mode of §3.9 — `HorzOrder: 2`,
`VerticalReconstruction: constant` — and confirm it reproduces `PressureGradCentered` to round-off
on the two-column test. This guards Requirement 2.5 and protects the existing default during
refactoring. The test is established with Phase 1 and rerun unchanged for Phase 2, where it also
confirms the wider stencil collapses correctly to the two-cell one.

Note that this configuration is not a supported production setting (§4.1.1); it exists so that this
comparison is possible. Keeping it costs one branch in the reconstruction and one in the
equation-of-state expansion.

**Covers:** Requirement 2.5; §3.9.

### 5.6 Cost check

Requirement 2.2 bounds the number of TEOS-10 evaluations, and nothing above tests it — the scheme
could satisfy every accuracy gate while quietly calling the equation of state inside the quadrature
loop. Two cheap checks close that:

- **Evaluation count.** With an instrumented `Eos`, confirm the number of specific-volume
  evaluations per time step is one per cell per layer and is **unchanged** when `QuadraturePoints`
  and `HorzOrder` are varied. This is the property Requirement 2.2 actually states, and it is a
  counter comparison, not a timing measurement, so it is deterministic and suitable for CI.
- **Wall time.** Record PGF kernel time relative to `PressureGradCentered` on a representative
  configuration, as a performance regression guard. The expected cost is dominated by the per-edge
  layer integrals (§4.1.3), roughly three times as many as a cell-based formulation on a hexagonal
  mesh; a result far above that suggests the per-cell coefficients are being recomputed per edge
  rather than cached.

**Covers:** Requirement 2.2.

### 5.7 Coverage summary

| Requirement / assumption | Verified by |
|---|---|
| 2.1 Accuracy at affordable resolution | §5.1 accuracy gate; §5.3; §5.4 |
| 2.2 Bounded TEOS-10 cost | §5.6 |
| 2.3 Robustness for thin, sloped layers | §5.2 (machine precision); §5.3 (under dynamics) |
| 2.4 Consistency with model state | §5.1 (layer-mean reference); §5.2 (reconstruction unit test) |
| 2.5 Runtime-selectable, backward compatible | §5.1 (three configurations); §5.5 |
| 2.6 Verified order of accuracy | §5.1 verification gate |
| A1 Edge accuracy ≠ cell accuracy | §5.1 verification gate (Phase 2) |
| A2 EOS expansion adequate across an edge | §5.1 contrast sweep; §5.4 |
| A3 Residual small enough in practice | §5.1 accuracy gate; §5.3; §5.4 |
| A4 PGF error causes the instability | §5.3 diagnostic, run before Phase 1 completes |
| §3.7.4 per-edge geopotential correction | §5.2 (cannot pass without it; guard (c) tests it directly) |
| §3.7.5 Round-off floor | §5.2, run in both precisions |

Requirement 2.7 (extensibility) is not testable directly; it is addressed by the configuration
design of §4.1.1 and by the phase structure of §4.5, which is itself the evidence that the framework
extends without restructuring.

This test is retained permanently rather than treated as a one-time transition check, and
that choice is the reason `PressureGradCentered` is kept as a separate implementation rather
than reimplemented as the lowest-order configuration of `PressureGradHighOrder` (§4.4). The
two functors read the mesh, `VertCoord`, and EOS state through independently written code, so
their agreement to round-off is evidence about *shared upstream* state as well as about the
PGF arithmetic: a wrong edge mask, a mis-indexed interface array, or a misinterpreted
`VertCoord` convention shows up as a disagreement. Were the centered scheme replaced by an
order-2 configuration of the new code, this comparison would reduce to comparing an
implementation against itself, and any defect upstream of the order switch would cancel out
of it.
