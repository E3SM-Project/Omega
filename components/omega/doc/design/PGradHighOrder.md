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
sub-options (horizontal reconstruction order, vertical-reconstruction mode, quadrature) must
be configurable, and the centered scheme must be recoverable as the lowest-order limit.

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
TEOS-10 polynomial, since the first derivatives reuse the same polynomial coefficients (the
same quantities computed by `gsw.specvol_first_derivatives`, which the Polaris reference
solution already uses). Indeed $\alpha_0$ is exactly the `Eos::SpecVol` field Omega already
computes, so only the three derivative fields are additional work.

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

Two consequences must be recorded:

- **Cost.** The column integral $\Pi_{i,k}$ is no longer a per-cell quantity that can be computed
  once and then differenced across each of the cell's edges; it depends on the edge through
  [](#edge-ref) and must be evaluated per edge, twice. On a hexagonal TRiSK mesh this is roughly
  three times as many integral evaluations as a cell-based formulation. The TEOS-10 call count is
  unaffected, which is the cost Requirement 2.2 binds; the extra work is polynomial arithmetic on
  coefficients already in cache.
- **The PGF is no longer the gradient of a single scalar.** With the expansion point depending on
  the edge, the discrete PGF is not exactly a discrete gradient, so it is not automatically
  curl-free and could in principle inject spurious vorticity. It is zero at every edge for a resting
  ocean the scheme resolves exactly, so nothing is injected in that case; elsewhere the effect
  enters at truncation order. This trade is deliberate — a robustness property that holds exactly,
  in exchange for a potential-form property that held exactly — but the magnitude of any spurious
  vorticity is something §5 must measure rather than assume away.

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

### 3.7 Where the cancellation is exact, and where it is not

This section makes Requirement 2.3 concrete: it explains why the obvious version of the robustness
property cannot be delivered, states the three conditions the implementation must meet, lists the
water columns that cancel exactly under each phase, and records what the design is still assuming
and must therefore test. It does *not* rely on subtracting a background profile.

#### 3.7.1 Why "uniform $\Theta$ and $S$" is the wrong condition

An earlier draft of this design required the scheme to return zero whenever the layer means
$\Theta_{i,k}$, $S_{i,k}$ are the same in both cells of an edge, for any tilt and any vertical
profile. **That requirement cannot be met, and asking for it would send the implementation chasing
an impossible target.**

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
   anchor. This one is currently *not* satisfied; see §3.7.4.

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
Phase 1 cancels. Whether this downward accumulation is in fact what produces the bottom-layer noise
seen in realistic global configurations is a plausible diagnosis, not an established one; it is
carried as A4 in §3.7.6 and tested in §5.3.

#### 3.7.4 A prerequisite in `VertCoord`

Condition 3 is a constraint on `VertCoord`, not on `PGrad`, and **it is not currently satisfied**.
`VertCoord::computeGeomZHeight` builds $z$ by accumulating
$\Delta z = \rho_0\,\alpha_{i,k}\,\tilde h_{i,k}$ — a midpoint rule using the layer-mean specific
volume, working upward from the bottom. If the pressure term integrates a reconstructed
$\hat\alpha$ while the geopotential term differences a $z$ built this way, the two are working from
**different specific volumes**, and the cancellation fails at $O(\alpha_{pp}\tilde h^2)$ — precisely
the order Phase 1 exists to recover. Worse, because $z$ is a running sum up the column, the mismatch
**accumulates over every layer between the anchor and the layer in question** rather than staying
local.

This is a prerequisite for Phase 1, not an implementation detail. There are three ways to resolve
it:

1. **Align `VertCoord`.** Integrate $z$ using the same reconstruction the PGF uses. Satisfies
   condition 3 outright, but changes answers for every existing configuration and modifies a module
   the PGF does not own.
2. **Build $z$ increments inside the PGF.** Leave `VertCoord` alone and have the PGF form its own
   layer $z$ increments from its own $\hat\alpha$. No answer change elsewhere, but it creates a
   second, slightly different $z$ in the model, which Requirement 2.4 exists to prevent.
3. **Accept the mismatch and document the error floor.** The guarantee of §3.7.2 weakens from "zero
   to machine precision" to "zero to $O(\alpha_{pp}\tilde h^2)$ accumulated over the column", and the
   test in §5.2 cannot pass.

**Recommendation: option 1.** Options 2 and 3 both give up the property this phase exists to
deliver. The change is answer-changing but small and well contained, and it is worth having on its
own merits: it makes the model's $z$ a more accurate integral of the hydrostatic relation. The
decision must be made before Phase 1 implementation begins, and it needs a baseline step (§5).

One related question is deliberately left to implementation time: `VertCoord` builds $z$ upward from
the bathymetry while pressure is built downward from the surface, so the two accumulate round-off
from opposite ends of the column. Whether the PGF should re-anchor its own accumulation is a
round-off question (§3.7.5), not a consistency one — condition 3 is satisfied either way.

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
- **A5 — Spurious vorticity from the edge-shared expansion point is negligible** (§3.3.1).

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

### 5.1 Test: Two-column HPGA convergence (extend existing)

Extend the four existing variants — `temperature_gradient`, `salinity_gradient`,
`surface_pressure_gradient`, `ztilde_gradient` — to run the centered scheme and the new
scheme at two orders:

- **Scheme selection.** Add `PressureGrad: { PressureGradType: FiniteVolume, … }` to
  `forward.yaml` and parametrize each task over three configurations:
  - `centered` — the legacy `PressureGradCentered` functor, unchanged;
  - `finite_volume_order4` — `ReconstructionOrder: 4`, `VerticalReconstruction: ppm` (the
    target);
  - `finite_volume_order2` — the same code in its centered limit (`ReconstructionOrder: 2`,
    `VerticalReconstruction: constant`, midpoint quadrature).

  (These variant names are provisional; the final spelling follows Polaris' naming conventions
  and is settled on the Polaris side. What matters here is that the last two are the *same*
  implementation at two orders, distinct from the legacy functor.)

  The third configuration exists because the order-2 verification gate below must exercise the
  *new* implementation. Running the legacy functor under a "second order" label would measure
  the convergence of code this design does not change, and would leave the new code's
  lowest-order path unverified. It also supplies the round-off comparison of §5.5 at every
  resolution rather than at a single configuration. The forward step still runs a single time
  step with only `PressureGradTendencyEnable: true`, reading the PGF acceleration from
  `NormalVelocityTend`.
- **Reference.** `reference.py`/`analysis.py` compare `NormalVelocityTend` against the
  layer-mean analytic HPGA (unchanged). For the high-order scheme the layer-mean comparison
  remains the correct target, since the scheme is a finite-volume, layer-averaged
  discretization.
- **Accuracy gate (new, Requirement 2.1):** at a representative coarse resolution (e.g. the
  coarsest in `horiz_resolutions`), the absolute RMS HPGA error vs. the reference must be below
  a tolerance, **and** the high-order RMS error must be below the centered RMS error at that
  same resolution (the scheme must demonstrably help where it matters).
- **Verification gate (Requirement 2.6):** the measured slope of RMS error vs. resolution,
  `omega_vs_reference_convergence_rate_*`, must fall within a band around the configured order
  of accuracy — nominally ~4 for `finite_volume_order4` and ~2 for `finite_volume_order2`
  (and for `centered`, whose band is unchanged from today). This band is retuned from its
  present values rather than loosened; a slope outside it fails the test and is treated as an
  implementation defect to be diagnosed, not as a tolerance to be widened.
- **Asymptotic range (implementation-time task):** it is not yet established that the existing
  `horiz_resolutions` sweep spans a range where a fourth-order slope is cleanly measurable —
  the sweep may be too coarse to have entered the asymptotic regime at its fine end, or fine
  enough that the reference solution's own quadrature error and roundoff contaminate the slope.
  Determining the usable window, and extending or tightening the sweep (and, if needed, the
  order of the Gauss quadrature in `reference.py`) so the designed order can be resolved, is
  part of implementing this test.
- **Consistency check (retained):** `omega_vs_polaris_rms_threshold` (~1e-10 m/s²) — Omega's
  forward output must still match the Python-computed HPGA, confirming the implementation
  matches the intended discretization.
- **Cfg keys.** New keys mirror the existing ones (`horiz_press_grad.cfg`): a coarse-resolution
  absolute tolerance, a `finite_volume_vs_centered` ratio gate, and per-scheme expected-rate
  bands.

Tests Requirements 2.1, 2.2 (via the bounded-EOS implementation exercised), 2.4, 2.5, 2.6.

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
  case, which isolates compressibility on its own. **Pass:** the PGF is zero at every edge and layer
  to **machine precision** (double-precision builds; the threshold tracks `Real`'s epsilon and the
  size of the hydrostatic terms, not a physical tolerance). Phase 1 and Phase 2 must both pass.
- **Profiles it does not:** $\Theta$, $S$ quadratic in pressure, then a realistic profile.
  **Pass:** the residual shrinks like $\tilde h^2$ (Phase 1) and $\tilde h^3$ (Phase 2) as the
  vertical grid is refined at fixed tilt, matching §3.7.3. A residual that does not shrink at the
  tabulated rate means one of the three conditions in §3.7.2 has been broken somewhere in the
  implementation; it is a bug to find, not a tolerance to widen.
- **Guard tests:** rerun an exactly resolved profile with (a) the endpoint-average interface term in
  place of [](#metric-divdiff), and (b) a cell-local expansion point in place of [](#edge-ref).
  Both must *fail* the machine-precision check. Without these, a passing result could just as easily
  come from a symmetry of the test setup as from the scheme being right.

The test is also run in a single-precision build, to measure the round-off floor of §3.7.5 and
settle whether the perturbation form is needed.

Condition 3 of §3.7.2 depends on `VertCoord` (§3.7.4), so the machine-precision check cannot pass
until whichever resolution is chosen there is in place; this test therefore doubles as the
acceptance criterion for that prerequisite. It is implemented as a fast C++ unit test and also run
as a configuration of the Polaris two-column task.

### 5.3 Test: Seamount resting state (steep-layer robustness)

Use the Polaris seamount task with tilted layers over a seamount and a horizontally uniform
$T$/$S$ stratification, integrated for a fixed period. **Pass:** the maximum spurious velocity
stays below a threshold and is substantially smaller for the high-order scheme than for the
centered scheme. This is the dynamical counterpart of §5.2 (the same resting state, now run
forward in the full solver) and the direct test of Requirement 2.3.

**Testing assumption A4 (§3.7.6).** Run the seamount case twice: once with a stratification Phase 1
resolves exactly ($\Theta$, $S$ linear in pressure) and once with a realistic profile. Spurious
velocity that survives the first run cannot be PGF truncation error, since the PGF is zero there to
machine precision, and must come from somewhere else in the model — most likely the layer-mean
treatment in the tracer and remapping operators (§4.3). This is the cheapest available check on
whether accelerating this work will actually cure the bottom-layer instability that motivated the
phasing, and it should be run **before** Phase 1 implementation is finished, with the centered
scheme on the realistic profile as the control. A null result would not invalidate the design, but
it would change its priority.

### 5.4 Test: Overflow (full non-Boussinesq dynamics)

Use the Polaris overflow task to exercise the PGF within the full non-Boussinesq equations
with strongly sloping layers and active dynamics. **Pass:** the solution remains stable and
the down-slope evolution agrees with the reference behavior; spurious mixing/velocity
attributable to PGF error is reduced relative to the centered scheme. This tests Requirements
2.1 and 2.3 under realistic, coupled conditions.

### 5.5 Test: Reduction to the centered scheme (permanent regression)

Configure the high-order option in its lowest-order limit (§3.9:
`ReconstructionOrder: 2`, `VerticalReconstruction: constant`, midpoint quadrature) and confirm
it reproduces `PressureGradCentered` to round-off on the two-column test. This guards
Requirement 2.5 and protects the existing default during refactoring.

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
