(omega-design-pressure-grad-high-order)=
# Higher-Order Horizontal Pressure Gradient

**Table of Contents**
1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Algorithmic Formulation](#3-algorithmic-formulation)
4. [Design](#4-design)
5. [Verification and Testing](#5-verification-and-testing)
6. [Supporting evidence](#6-supporting-evidence)

## 1 Overview

This document specifies a more accurate discretization of the horizontal pressure gradient force
(PGF) in Omega's momentum equation, delivered in two phases: a consistent second-order scheme, then
a fourth-order one. It is the "future design document" promised in
{ref}`Pressure Gradient <omega-design-pressure-grad>` (§2.3, §3.2), which commits Omega to a
high-order option "similar to
[Adcroft et al. 2008](https://doi.org/10.1016/j.ocemod.2008.02.001)" but defers the details. The
existing centered, second-order Montgomery-potential scheme (`PressureGradCentered`) remains the
default and the reference implementation; this design fills in the second `PressureGradType` entry,
which the `PressureGrad` class already anticipates as a placeholder.

That entry, and the functor behind it, are named **`FiniteVolume`** and
**`PressureGradFiniteVolume`** — not anything containing "high order" (§4.1.1). The choice is
deliberate. Phase 1 uses the same two-cell horizontal stencil as the centered scheme and is second
order, not higher order, so naming the class for an order it does not have would misdescribe what
ships first; and what Phase 1 does deliver is *consistency* (§1.1), which is a separate property from
order of accuracy. `FiniteVolume` refers to the layer-integrated control-volume form of §3.1, which
both phases share.

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

One property of Omega's vertical coordinate simplifies the *pressure* part of this problem,
and §3.1.1–§3.2 turn on it. Pseudo-height is defined as $\tilde z \equiv -p/(\rho_0 g)$ with
**no offset** ({ref}`omega-design-governing-eqns-omega1` §5.1, which considers adding one and
rejects it so that "$\tilde z$ varies identically to pressure"). A surface of constant
$\tilde z$ is therefore a surface of constant $p$ in every column and at every time, so
$\nabla_{\tilde z}\,p \equiv 0$: the pressure-gradient term vanishes when the gradient is taken
along a coordinate surface, and the entire horizontal PGF is the geopotential compared at fixed
pressure.

**This is a statement about coordinate surfaces, not about layer interfaces, and the distinction is
the whole subject of this document.** Layer interfaces sit at $\tilde z^{\text{top}}_{i,k}(x,y,t)$ —
wherever the ALE coordinate puts them. They are tilted, they move in time, and layer $k$ occupies a
*different pressure range* in every column. Nothing here assumes otherwise, and a scheme that did
would be worthless: it is exactly that difference in pressure range between neighbouring columns that
produces the pressure-gradient error this design exists to remove.

The two facts coexist because $\nabla_{\tilde z}$ differentiates along a coordinate surface, while the
layer interfaces enter the layer-integrated form of §3.1 as the *limits of integration*. Their
horizontal variation is what generates that form's two metric terms, through Leibniz' rule, and those
terms are individually large. What §3.1.1 shows is that four individually nonzero terms **sum** to
zero — not that any one of them vanishes, and not that the interfaces are level. Every quantity the
resulting scheme actually computes is a tilt term: if interfaces were isobaric, the cross-edge
interface pressure differences would vanish, the scheme would collapse to a plain difference of
`GeomZInterface`, and there would be nothing left to design.

The reduction accordingly holds for any ALE layering *of pseudo-height* — p-star, sigma-p,
z-tilde-star, or any hybrid — because all of them relocate interfaces within $\tilde z$ and none of
them redefines $\tilde z$. It would fail only if Omega replaced pseudo-height with a coordinate that
is not a function of pressure alone: geometric $z$, a true $z$-star, or a genuine isopycnal
coordinate. §3.8 states the scope precisely and records what would have to come back in that event.

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
missing there is not high-order accuracy — it is *consistency*. §3.9 shows that
`PressureGradCentered` is, exactly and identically, the **first-order** conversion of a
height difference taken at fixed layer index into one taken at fixed pressure. Its error is
therefore the truncation of that conversion, and it accumulates downward through the column
(§3.7.3). That error is present at second order and does not require a fourth-order scheme
to remove.

The work is therefore split so the consistency fix can be delivered and tested on its own:

- **Phase 1 — a consistent second-order scheme.** The exact fixed-pressure geopotential
  difference: the equation-of-state expansion about a state shared across each edge, a
  mean-preserving *linear* reconstruction of $\Theta$ and $S$ in pressure, and the remainder
  by which the exact fixed-pressure height difference exceeds the centered scheme's
  first-order estimate of it (§3.5). The horizontal operator is the same two-cell stencil the
  centered scheme already uses, so the scheme remains second order in the horizontal. What
  Phase 1 buys is the robustness property of Requirement 2.3: the pressure gradient is zero to
  machine precision for any resting ocean whose profile varies linearly with pressure, no
  matter how the layers tilt.
- **Phase 2 — fourth order.** A parabolic vertical reconstruction and a wider horizontal
  stencil, plus the option of a second-order equation-of-state expansion. This raises the
  order of accuracy; it does not change the robustness property, which Phase 1 establishes and
  Phase 2 inherits by construction rather than having to re-establish (§3.6.1).

The two phases share one implementation, one set of configuration options, and one test
suite; Phase 2 widens settings that Phase 1 puts in place. Section 3 marks which parts of
the formulation belong to which phase, and §4.5 records what each phase delivers and what it
depends on.

## 2 Requirements

### 2.1 Requirement: Higher accuracy than the centered scheme at affordable resolution

The new PGF must produce substantially lower absolute error than
`PressureGradCentered` at the coarse-to-moderate resolutions Omega can afford to run,
for representative stratified columns with horizontal gradients of temperature, salinity,
surface pressure, and coordinate slope. Absolute error at affordable resolution is what
determines whether the scheme is *useful*; it is complemented by the separate verification
requirement in §2.6, which determines whether it is *correct*.

### 2.2 Requirement: Bounded TEOS-10 cost

The number of TEOS-10 specific-volume evaluations performed by the PGF must be
independent of the reconstruction order and the horizontal stencil width, and bounded at
approximately one evaluation (plus its first derivatives) per cell per layer per time step —
comparable to the cost Omega already pays to compute the `Eos::SpecVol` field. The scheme
must not require evaluating the full TEOS-10 polynomial at sub-layer sample points, nor
integrating the TEOS-10 polynomial itself semi-analytically.

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
through an inconsistent geopotential or a truncated fixed-pressure comparison.

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

The new scheme must be selectable at runtime through the existing `PressureGrad`
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

The scheme starts from the layer-integrated, finite-volume form of the pressure
and geopotential terms derived in
{ref}`the governing equations <omega-design-governing-eqns-omega1>`
(the {ref}`full volume-integral form <omega-v1-vh-momentum-reynolds2>`). Neglecting the
turbulent correlations and dropping the resolved-component notation, the layer tendency is
({ref}`Pressure Gradient <omega-design-pressure-grad>` §3.2):

$$
T^p &= - \int_A \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \rho_0 \, \left( \nabla \Phi \right) \, d\tilde{z} \, dA \\
& - \int_{\partial A} \left( \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \rho_0 \left(\alpha p \right) \, d\tilde{z} \right) dl \\
& + \int_A \rho_0 \left[ \alpha \, p \, \nabla \tilde{z}_k^{\text{top}} \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} \, dA \\
& - \int_A \rho_0 \left[ \alpha \, p \, \nabla \tilde{z}_k^{\text{bot}} \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \, dA \\
& + \int_A \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \rho_0 \, p \, \nabla \alpha \, d\tilde{z} \, dA.
$$ (ho-target)

The five terms are: the geopotential (gravity) body force integrated over the layer volume;
the side-wall integral of $\alpha p$ (the pressure traction on the cell faces); two metric
terms accounting for the pressure traction on the sloping top and bottom layer interfaces;
and a $p\nabla\alpha$ term.

**Two of those terms differ from [](#ho-target) as inherited from
{ref}`omega-design-governing-eqns-omega1`, and the differences are corrections.** The fifth term is
**absent** there; omitting it leaves a scheme that does not converge at all. The signs on lines 3 and 4
are **opposite** to those printed there; with the inherited signs the cancellation of §3.1.1 doubles
rather than vanishing. Both follow from applying Leibniz' rule directly, immediately below, and neither
is a convention difference in how the traction normal or $\nabla\tilde z^{\text{top}}$ is defined.
`OmegaV1GoverningEqns.md` is not maintained against the code and should be corrected there rather than
overridden here; §6.2 states what is owed to it and §6.1 the supporting measurements.

The derivation is short enough to give in full, because everything in §3 rests on it.
Differentiating $\tilde z \equiv -p/(\rho_0 g)$ and applying the hydrostatic relation gives the mass
element $\rho\,dz = \rho_0\,d\tilde z$ ({ref}`omega-design-governing-eqns-omega1` §5). The layer
integral of the pointwise acceleration $-(\alpha\nabla p + \nabla\Phi)$ at fixed $\tilde z$ is
therefore
$-\int\rho_0(\alpha\nabla p + \nabla\Phi)\,d\tilde z$. Writing
$\alpha\nabla p = \nabla(\alpha p) - p\nabla\alpha$ and applying Leibniz' rule to the first
piece,

$$
\int_{\tilde z_k^{\text{bot}}}^{\tilde z_k^{\text{top}}} \nabla(\alpha p)\,d\tilde z
= \nabla \int_{\tilde z_k^{\text{bot}}}^{\tilde z_k^{\text{top}}} \alpha p \,d\tilde z
- (\alpha p)_{\text{top}}\,\nabla \tilde z_k^{\text{top}}
+ (\alpha p)_{\text{bot}}\,\nabla \tilde z_k^{\text{bot}},
$$ (leibniz)

gives exactly lines 2–5 of [](#ho-target), with the divergence theorem turning the first
term on the right into the side-wall line integral. This is the Adcroft et al. (2008)
finite-volume route — the net pressure force on each control volume obtained by integrating
in-situ pressure over its faces, rather than by forming the pointwise product
$\alpha\nabla p$ — and it is the form MOM6 uses successfully for the same problem.

#### 3.1.1 In pseudo-height the four pressure terms cancel identically

Lines 2–5 of [](#ho-target) are, by construction, nothing but
$-\int_A\int\rho_0\,\alpha\nabla_{\tilde z} p\,d\tilde z\,dA$ rewritten. In Omega's
coordinate that quantity is **identically zero**.

Pseudo-height is $\tilde z \equiv -p/(\rho_0 g)$ with no offset, per
{ref}`omega-design-governing-eqns-omega1` §5.1, which considers adding one and rejects it so
that "$\tilde z$ varies identically to pressure". The map $p \leftrightarrow \tilde z$ is
therefore *universal*: it is the same in every column, at every time, irrespective of surface
pressure. A surface of constant $\tilde z$ is a surface of constant $p$, so

$$
\nabla_{\tilde z}\, p \;\equiv\; 0,
$$ (grad-p-zero)

and lines 2–5 of [](#ho-target) sum to zero — not approximately, and not only for a resting
ocean, but as an algebraic identity for any state whatsoever. This is the standard result that
in an isobaric coordinate the horizontal pressure gradient force is $-\nabla_p\Phi$.

**What [](#grad-p-zero) does and does not say.** It says the coordinate *surfaces* are isobars.
It says nothing about the layer *interfaces*, which are at $\tilde z^{\text{top}}_{i,k}(x,y,t)$
and are tilted, moving, and at different pressures in different columns — as they must be in an
ALE coordinate. The horizontal variation of those limits of integration is precisely what
produces the two metric terms through [](#leibniz), and each of the four pressure terms is
individually large: measured on a two-column state at 50 m/km tilt, the side-wall term alone is
$0.55$ m s$^{-2}$. The cancellation is a property of their *sum*. Read the other way, it is a
statement that the four terms carry no information the geopotential term does not already carry,
so discretizing them can only add error — not a statement that the tilt they describe is absent.

Three things follow, and they set the shape of the rest of §3.

1. **The entire horizontal PGF is the geopotential compared at fixed pressure.** §3.2 states
   the resulting target. Comparing "at fixed pressure" is hard exactly because the interfaces
   are not at fixed pressure; the difficulty of the original problem is preserved, relocated
   into a single term.
2. **Adcroft's problem is not Omega's problem.** Adcroft et al. (2008) needed the face-integral
   form because their coordinate was *not* isobaric, so $\alpha\nabla_r p$ was a real term
   requiring careful treatment. Omega's coordinate is isobaric by definition, and remains so
   under ALE: p-star, sigma-p and any other layering Omega adopts move interfaces *within*
   $\tilde z$ without redefining $\tilde z$, so [](#grad-p-zero) is untouched. The remaining
   problem is not a force-balance problem at all — it is the problem of integrating the
   hydrostatic relation accurately enough that two columns' geopotentials can be differenced at
   a common pressure.
3. **Discretizing the four terms separately is worse than not discretizing them.** Because they
   must cancel, any discretization of them contributes only error — and the cancellation is
   delicate in a way that is not apparent from [](#ho-target). It requires the $p\nabla\alpha$
   term to be integrated over the mean of the *two columns' own* layer pressure ranges rather
   than over the edge-mean range, which is the choice an edge-centred implementation invites;
   mispaired, the residual is $O(\tilde h^0)$ and does not converge under vertical refinement at
   all, while remaining invisible to every gate in §5 — on the exact set the mispaired term is
   identically zero. The reduced target of §3.2 has no such requirement, because it has no terms
   to pair.

[](#grad-p-zero) has been confirmed to hold discretely as well as continuously: the four pressure
terms were measured to sum to $3\times10^{-15}$ of the size of the side-wall term alone, at 50 m/km
tilt, for every profile tried. §6.1 gives the numbers and the harness.

### 3.2 The exact target

Pressure is an exactly linear function of $\tilde z$,

$$
p(\tilde z) = p_{k}^{\text{top}} - \rho_0 g \left( \tilde z - \tilde z_k^{\text{top}} \right),
\qquad \tilde z_k^{\text{bot}} \le \tilde z \le \tilde z_k^{\text{top}},
$$ (p-linear)

anchored at the interface pressures $p_k^{\text{top}}$, $p_k^{\text{bot}}$ that `VertCoord`
already provides (`PressureInterface`). No reconstruction of $p$ within a layer is needed or
permitted — using [](#p-linear) is exactly the `VertCoord` anchoring required by
Requirement 2.4. Within a column, pressure is used from here on as the vertical integration
variable in place of $\tilde z$; the two are interchangeable by [](#p-linear), and pressure is
the variable in which two columns' profiles can be compared.

Combining [](#grad-p-zero) with [](#ho-target) and mass-weighting the layer average with
$\rho_0\,d\tilde z = -dp/g$, the exact edge-normal acceleration at edge $e$, layer $k$ is

$$
T^p_{e,k} = -\,\frac{g}{\Delta p_{e,k}} \int_{p_{e,k}^{\text{top}}}^{p_{e,k}^{\text{bot}}}
\left[\nabla_n z\right]_{p} \; dp
\;-\; \left[\nabla_n \left(\phi_{TP} + \phi_{SAL}\right)\right]_{e},
$$ (ho-exact)

in which $[\nabla_n z]_p$ is the edge-normal gradient of geometric height taken **at fixed
pressure**, $\Delta p_{e,k} = p^{\text{bot}}_{e,k} - p^{\text{top}}_{e,k}$ is the edge control
volume's pressure range, and the geometric height of column $i$ as a function of pressure is
the hydrostatic integral

$$
z_{i}(p) = z_{i}^{\text{anchor}} + \frac{1}{g}\int_{p}^{p_i^{\text{anchor}}}
\hat\alpha_{i}(p')\, dp',
$$ (z-of-p)

with $\hat\alpha_{i}(p)$ the reconstructed specific volume in column $i$ (§3.3–§3.4).
Equation [](#ho-exact) is the whole scheme. Everything that follows is either an ingredient of
$\hat\alpha$ (§3.3, §3.4), a statement about how [](#ho-exact) is evaluated without losing
precision (§3.5), or a statement about when it returns zero (§3.7).

Two properties of [](#ho-exact) are worth noticing immediately, because between them they do the work
that a four-term force balance would need several delicate arguments to secure.

- **The robustness property is structural, not a cancellation.** If the two columns' reconstructed
  profiles agree as functions of pressure, $\hat\alpha_L(p) \equiv \hat\alpha_R(p)$, then by
  [](#z-of-p) the two $z_i(p)$ differ by a constant fixed by the anchors, and $[\nabla_n z]_p$ is
  *pointwise* zero at every pressure. The integral in [](#ho-exact) is then zero regardless of
  the tilt, the layer thicknesses, the bathymetry, or the quadrature used to evaluate it. There
  are no terms to pair and no signs to get right. For a resting ocean the anchors agree by the
  inverse-barometer relation — $z^{\text{surf}}_i$ tilts exactly so as to compensate a horizontal
  gradient in $p^{\text{surf}}_i$ — which is the state the `surface_pressure_gradient` Polaris
  variant initializes.
- **The vertical quadrature is not an accuracy knob.** $[\nabla_n z]_p$ varies across the layer
  only through $\frac{d}{dp}[\nabla_n z]_p = -\frac{1}{g}[\nabla_n\hat\alpha]_p$, the horizontal
  contrast in specific volume at fixed pressure — a quantity that is *small* and that vanishes
  identically on the exact set. The layer average in [](#ho-exact) is therefore a nearly-constant
  integrand, handled in closed form for a polynomial $\hat\alpha$. There is no quadrature setting
  (§4.1.1).

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
reconstruction reproduces perfectly. By [](#z-of-p) the two columns' $z_i(p)$ would then differ by
more than a constant, $[\nabla_n z]_p$ would be nonzero, and the scheme would generate spurious
flow from nothing but its own EOS approximation. This choice remains load-bearing under the reduced
target of §3.2 — it is the *only* remaining way for the exactness of [](#ho-exact) to fail on a
resolved profile — and guard test (b) of §5.2 confirms it fires.

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

This has a cost consequence worth recording. The column quantity built from $\hat\alpha$ — under §3.2
the hydrostatic integral [](#z-of-p) — is not a per-cell quantity that can be computed once and
differenced across each of the cell's edges. It depends on the edge through [](#edge-ref) and must be
evaluated per edge, twice. On a hexagonal TRiSK mesh this is roughly three times as many column evaluations as a
cell-based formulation. The TEOS-10 call count is unaffected, which is the cost Requirement 2.2
binds; the extra work is polynomial arithmetic on coefficients already in cache. §3.5 shows how to
arrange it so that only a *small* correction is accumulated per edge, which bounds both the cost
and the round-off.

Because $\Theta$, $S$ are reconstructed as low-order polynomials in pressure and the expansion
[](#alpha-taylor) is linear in $\Theta$, $S$ and $p$, $\hat\alpha$ is a **low-order polynomial in
pressure** and every integral of it required by §3.2 and §3.5 is available in closed form. No
TEOS-10 evaluation occurs inside any integral — which is what resolves the feasibility and cost
concerns with applying Adcroft-style analytic integration to TEOS-10.

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

### 3.5 The remainder relative to the centered scheme

Equation [](#ho-exact) is not evaluated directly. It is evaluated as the centered scheme plus a
remainder, for one reason: the fixed-pressure height difference is a small residual of two large
quantities, and grouping the large parts so that they cancel *analytically* rather than in floating
point is what preserves the precision §3.7.5 is about.

**`PressureGradCentered` is the first-order fixed-pressure shift, exactly.** This is an algebraic
identity about the code in `PGrad.h`, not an approximation, and it is the pivot of this section.
Write $\Delta_e f \equiv \sum_{i\in CE(e)} -n_{e,i}\,f_i$ so that $[\nabla_n f]_e = \Delta_e f/d_e$,
and let $q_{i,k}$ = `PressureInterface`, $Z_{i,k}$ = `GeomZInterface`, $\alpha_{i,k}$ = `SpecVol`,
with layer $k$ bounded by interfaces $k$ and $k+1$. Define

$$
\mathcal{S}_{e,k} \;\equiv\; \tfrac{1}{2}\left(\Delta_e Z_{k} + \Delta_e Z_{k+1}\right)
\;+\; \frac{\bar\alpha_{e,k}}{2g}\left(\Delta_e q_{k} + \Delta_e q_{k+1}\right),
\qquad
\bar\alpha_{e,k} = \tfrac{1}{2}\!\!\sum_{i \in CE(e)}\!\! \alpha_{i,k}.
$$ (centered-shift)

Then $T^{p,\text{ctr}}_{e,k} = -\frac{g}{d_e}\,\mathcal{S}_{e,k} - \left[\nabla_n(\phi_{TP}+\phi_{SAL})\right]_e$. The proof is two lines: with
$M_{i,k\pm} = \alpha_{i,k}\,q_{i,\cdot} + g Z_{i,\cdot}$, the identity
$\Delta_e(\alpha q) = \bar\alpha\,\Delta_e q + \bar q\,\Delta_e \alpha$ splits the Montgomery
difference, and the $\bar q\,\Delta_e\alpha$ halves cancel the scheme's $\bar p\,\nabla\alpha$ term
**term for term**, because $\tfrac12(\bar q_{k} + \bar q_{k+1}) = \bar p^{\text{mid}}_{e,k}$ exactly
— `PressureMid` is the exact arithmetic midpoint of the two interface pressures
(`VertCoord.cpp:1067`).

What [](#centered-shift) *is*, physically: Taylor-shifting each column's height from its own
interface pressure to a pressure common to both gives
$z_i(p) = Z_{i,k} - \frac{\alpha}{g}\left(p - q_{i,k}\right) + O\!\left((p-q_{i,k})^2\right)$, hence

$$
g\left[\nabla_n z\right]_{p} \;=\; g\left[\nabla_n Z\right]_{k}
\;+\; \alpha \left[\nabla_n q\right]_{k} \;+\; O\!\left(\frac{\alpha_p}{g}(\Delta_e q)^2\right).
$$ (first-order-shift)

So $\mathcal{S}_{e,k}$ is the trapezoidal average over the layer's two interfaces of the
**first-order** conversion from fixed layer index to fixed pressure, with the layer-mean $\alpha$
standing in for $\hat\alpha$ at the interface. The centered scheme's entire Montgomery-potential
apparatus is that conversion and nothing else, and its error is the conversion's truncation. This is
why §3.7.3's first-order behaviour is what it is, and it is the sharpest available statement of what
Phase 1 fixes: **Phase 1 computes the centered scheme's error and removes it.**

**The remainder.** Define

$$
\mathcal{R}_{e,k} \;\equiv\; \frac{1}{\Delta p_{e,k}}
\int_{p^{\text{top}}_{e,k}}^{p^{\text{bot}}_{e,k}} \Delta_e z(p)\; dp \;-\; \mathcal{S}_{e,k},
\qquad
T^p_{e,k} = -\frac{g}{d_e}\left(\mathcal{S}_{e,k} + \mathcal{R}_{e,k}\right)
- \left[\nabla_n(\phi_{TP}+\phi_{SAL})\right]_e ,
$$ (ho-discrete)

with $z_i(p)$ from [](#z-of-p) built on the edge-shared $\hat\alpha^{(e)}$. Equation
[](#ho-discrete) is exactly [](#ho-exact) rewritten; `VertCoord`'s $Z$ enters through
$\mathcal{S}$ and leaves again through $\mathcal{R}$, and cancels in exact arithmetic. That is what
keeps Requirement 2.4 satisfied — the scheme does not construct a second geometric height, it
constructs a correction to a difference — and it is why the split is a numerical device rather than
a change of target.

$\mathcal{R}$ has four sources. Three are local to the layer; the fourth is not.

1. **Interface value versus layer mean.** [](#centered-shift) weights $\Delta_e q$ by the layer-mean
   $\bar\alpha_{e,k}$ where [](#first-order-shift) calls for $\hat\alpha^{(e)}$ evaluated *at the
   interface*. Contributes $\propto \hat\alpha'\,\Delta p_{k}\,\Delta_e q_{k}$ per interface, with
   opposite sign at the layer's top and bottom.
2. **Second and higher order in the shift.** The $O((\Delta_e q)^2)$ term of
   [](#first-order-shift). Its column-symmetric part cancels in the edge difference; what survives
   is $-\frac{(\Delta_e q_k)^2}{8g}\,\Delta_e \hat\alpha'$ and higher, i.e. it is driven by the
   *horizontal contrast* in compressibility, not by compressibility itself.
3. **Layer average versus interface trapezoid.** [](#ho-discrete) averages $\Delta_e z(p)$ over the
   edge control volume's pressure range; [](#centered-shift) takes the mean of two interface values.
   The gap is set by the curvature $\frac{d^2}{dp^2}\Delta_e z = -\frac{1}{g}\Delta_e\hat\alpha'$.
4. **The `VertCoord` height is built from a different $\hat\alpha$.** `VertCoord` accumulates
   $\Delta Z = \rho_0\,\alpha_{i,k}\,\tilde h_{i,k}$, which by §3.7.4's [](#z-increment-exact) is the
   exact layer integral of a *cell-local* Phase 1 $\hat\alpha$. The scheme's water column is the
   *edge-shared* one, whose layer mean is

   $$
   \left\langle\hat\alpha^{(e)}_{i,k}\right\rangle = \bar\alpha_0^{e,k} + \bar\alpha_\Theta^{e,k}\left(\Theta_{i,k}-\bar\Theta^{e,k}\right) + \bar\alpha_S^{e,k}\left(S_{i,k}-\bar S^{e,k}\right) + \bar\alpha_p^{e,k}\left(p^{\text{mid}}_{i,k}-\bar p^{e,k}\right),
   $$ (alpha-edge-layer-mean)

   the first-order Taylor estimate of $\alpha_{i,k}$ from the edge state. The per-layer mismatch

   $$
   \varepsilon_{i,k} = \rho_0\,\tilde h_{i,k}\left(\alpha_{i,k} - \left\langle\hat\alpha^{(e)}_{i,k}\right\rangle\right)
   $$ (eos-remainder)

   is therefore precisely the **second-order Taylor remainder of the equation of state across the
   edge**, and its edge difference **accumulates down the column**. This is the one non-local source,
   and it is why Phase 1 needs a column scan (§4.1.3). It also identifies source 4 with assumption A2
   (§3.7.6): the same second-order EOS expansion that would answer A2 shrinks this source directly.

Three requirements on how $\mathcal{R}$ is built, all of them testable:

- **R1 — built from the edge-shared $\hat\alpha^{(e)}$** of [](#edge-ref), so that condition 1 of
  §3.7.2 is met. Guard test (b) of §5.2 checks this.
- **R2 — evaluated analytically, never as a numerical difference.** Each of the four sources above
  is a small quantity with a closed form. Forming $\mathcal{R}$ by computing
  $\frac{1}{\Delta p}\int \Delta_e z\,dp$ and $\mathcal{S}$ separately and subtracting them would
  throw away exactly the precision the split exists to protect.
- **R3 — $\mathcal{R}_{e,k} = -\mathcal{S}_{e,k}$ to round-off whenever
  $\hat\alpha^{(e)}_L(p) \equiv \hat\alpha^{(e)}_R(p)$.** This is the machine-precision property of
  Requirement 2.3.1 restated in the form the implementation must satisfy, and it is what §5.2 tests.
  Note that it makes the test *sharper* than a check that the tendency is small: $\mathcal{S}$ is
  the centered scheme's answer on the same state, which §3.7.3 measures at $2\times10^{-5}$
  m s$^{-2}$, so R3 asserts a cancellation of ten orders of magnitude against a known nonzero
  reference rather than against zero.

#### 3.5.1 Accumulating without cancellation

$\mathcal{S}_{e,k}$ as written in [](#centered-shift) is the sum of two terms of order 200 m, at a
coordinate tilt of 50 m/km over a 4 km edge, that cancel to a few millimetres. That is
`PressureGradCentered`'s existing round-off exposure, not something the new scheme introduces
(§3.7.5), but the same algebra that produced [](#centered-shift) removes it, and the design
recommends doing so.

Let $\Gamma_{e,k} = \Delta_e Z_{k} + \frac{\bar\alpha_{e,k}}{g}\Delta_e q_{k}$ and
$\Gamma^{+}_{e,k} = \Delta_e Z_{k+1} + \frac{\bar\alpha_{e,k}}{g}\Delta_e q_{k+1}$, so that
$\mathcal{S}_{e,k} = \tfrac12(\Gamma_{e,k} + \Gamma^{+}_{e,k})$. Using
$Z_{i,k+1} - Z_{i,k} = -\rho_0\alpha_{i,k}\tilde h_{i,k}$ and
$q_{i,k+1} - q_{i,k} = \rho_0 g \tilde h_{i,k}$ together with the same product rule as before,

$$
\Gamma^{+}_{e,k} - \Gamma_{e,k} = -\,\rho_0 \, \overline{\tilde h}_{e,k}\, \Delta_e \alpha_{k},
\qquad
\Gamma_{e,k+1} - \Gamma^{+}_{e,k}
= \frac{\Delta_e q_{k+1}}{g}\left(\bar\alpha_{e,k+1} - \bar\alpha_{e,k}\right).
$$ (gamma-increments)

Both increments are products of a small factor with a bounded one — the *horizontal* contrast in
specific volume within a layer, and the *vertical* contrast between adjacent layers — so
$\mathcal{S}$ can be walked down the column with no cancellation of large numbers anywhere, starting
from $\Gamma_{e,1}$, which is the inverse-barometer residual at the sea surface and is itself small
for a state at rest.

[](#gamma-increments) should be reconciled against §3.7.3's [](#centered-error), which describes the
same quantity in a different grouping, as a check on both. In the limit [](#centered-error) is written
for — specific volume horizontally uniform within each layer, thickness merely redistributed — the
first increment vanishes identically and the second one carries the whole signal, reproducing
[](#centered-error)'s downward accumulation through the horizontal displacement of the interfaces.
The two expressions are expected to agree up to a $k$-independent constant fixed by the anchor, and
confirming that (including the constant, and without assuming
$\sum_j \Delta_e\tilde h_j = 0$) is one of the algebraic checks assigned to the reference derivation
in §4.5.3.

Consequences worth stating, because they change several things downstream:

- Both $\mathcal{S}$ and source 4 of $\mathcal{R}$ are column prefix sums. They should be
  accumulated in the same scan, from the same end of the column (§3.7.4 closing note), and cost one
  scan rather than two.
- The perturbation form §3.7.5 holds in reserve becomes unnecessary in double precision and probably
  in single as well. That is now a measurement (§5.2 in both precisions) with a plausible expected
  answer rather than an open question.
- The reduction to `PressureGradCentered` (§3.9) becomes *algebraic* rather than bit-for-bit,
  because [](#gamma-increments) evaluates $\mathcal{S}$ in a different order from `PGrad.h`. This is
  the trade the design accepts: §4.4 already expects agreement to round-off rather than bit-for-bit,
  and it is a better cross-check for being arithmetically independent.

### 3.6 The edge operator

The discrete tendency lives at edge $e$ and is the edge-normal projection of the PGF. Every
horizontal difference in §3.5 is the two-cell TRiSK operator

$$
\left[\nabla_n f\right]_{e} = \frac{1}{d_e}\sum_{i \in CE(e)} -n_{e,i}\, f_{i}^{(e)},
$$ (edge-grad)

applied to column quantities $f_i^{(e)}$ evaluated with the equation-of-state coefficients shared
across edge $e$ [](#edge-ref). **Phase 1 uses this operator as written** — the same two-cell stencil
`PressureGradCentered` uses, so the horizontal accuracy is unchanged from today. Everything Phase 1
gains is in the vertical and in the fixed-pressure comparison.

Nothing is assumed constant across the edge: the horizontal contrast in $\Theta$ and $S$ between the
two columns enters at full strength through the two columns' $\hat\alpha^{(e)}_i(p)$ and hence
through $z_i(p)$. What Phase 1 does *not* do is raise the order of the horizontal difference, which
is Phase 2's job.

#### 3.6.1 Raising the horizontal order (Phase 2)

**The robustness property of §3.7 places no constraint on the wide stencil.** This is worth stating
explicitly, because the opposite would be true of a four-term force-balance formulation, where a
quantity interpolated to the edge from four or more cells belongs to no control volume and so has no
reason to balance. That would force the stencil to be a weighted sum of two-column pair
contributions, and settling its form on an unstructured, variable-resolution mesh would be a
prerequisite for Phase 2. Under [](#ho-exact) it is not.

The reason is that the quantity being differenced is $z_i(p)$ — a *per-column function of pressure*,
not a nonlinear product formed at the edge. For a resting ocean in the exact set, [](#z-of-p) makes
$z_i(p)$ the same function of $p$ in every cell of the neighbourhood, up to the constant set by each
column's anchor, and for a resting ocean those constants agree as well. **Any** consistent
edge-normal gradient operator — any operator that annihilates a constant, which every consistent
gradient operator does — therefore returns zero. The machine-precision property is inherited by
construction, and there is no pair-decomposition requirement.

What remains for Phase 2 is an ordinary accuracy question rather than a structural one:

- The **stencil width and weights** are chosen for order of accuracy on Omega's unstructured,
  variable-resolution TRiSK mesh, drawing on the wider edge neighbourhoods used for high-order
  tracer reconstruction (cf. {ref}`omega-design-governing-eqns-omega1` §10 and White & Adcroft
  (2008)). Assumption A1 (§3.7.6) still applies: the achieved order must be measured, not assumed.
- The **expansion point** of [](#alpha-taylor) must be shared across whatever set of cells the
  stencil couples, for the same reason §3.3.1 gives at a two-cell edge. Sharing one point across the
  whole neighbourhood is the simplest option and is sufficient for exactness; building the stencil as
  a weighted sum of pairs, each with its own shared point, keeps the expansion error smaller and is
  therefore preferable on assumption-A2 grounds. Both work; the choice is now a trade-off to be
  measured rather than a constraint imposed by the robustness property.

Phase 1 does not depend on either answer.

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
   This is what the edge-shared expansion point of §3.3.1 delivers, and under the reduced target of
   §3.2 it is the *only* remaining mechanism by which exactness can fail on a resolved profile.
2. **The two columns are compared at a common pressure, to all orders retained.** Layer $k$ spans
   different pressure ranges in the two columns, so $[\nabla_n z]_p$ requires shifting each
   column's height from its own interface pressures to a pressure shared with its neighbour.
   Truncating that shift at first order is what `PressureGradCentered` does, and its truncation is
   the whole of the error Phase 1 removes; the shift must therefore be carried to the order of the
   reconstruction, which is what §3.5's sources 1–3 supply.
3. **The height being shifted is the hydrostatic integral of that same profile.** `VertCoord`'s
   $Z$ is built from each cell's own exact TEOS-10 $\alpha$, not from the edge-shared
   $\hat\alpha^{(e)}$, so the difference between the two must be accounted for. This is §3.5's
   source 4, and it is a column prefix sum; see §3.7.4.

Conditions 2 and 3 are about the discretization alone and are under the implementation's control.
Meeting them means the scheme returns the **exact pressure gradient of the water column it has
reconstructed** — it contributes no error of its own beyond the reconstruction. Condition 1 is
about the state: it holds when the reconstruction reproduces the true profile, and fails, by
however much the reconstruction misses, when it does not. Meeting all three gives a PGF that is
**zero to machine precision**, for any tilt, any layer thickness, and any bathymetry.

Note that conditions 2 and 3 are, between them, exactly the requirement R3 of §3.5 — that
$\mathcal{R}$ cancel $\mathcal{S}$ to round-off on a resolved profile. Stating them separately is
useful because they fail in distinguishable ways: a defect in condition 2 leaves a residual local
to the layer, while a defect in condition 3 leaves one that grows with depth. §5.2's guard tests
are arranged to tell them apart.

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

#### 3.7.4 Condition 3 and the column scan

Condition 3 constrains how the height being shifted is built. Two questions arise, and they have
different answers; they are easy to conflate, and doing so leads to the conclusion that `VertCoord`
must be changed, which it need not be.

**The quadrature question, which is settled.** `VertCoord::computeGeomZHeight` builds $z$ by
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
the height be the hydrostatic integral of the *same* $\hat\alpha$ the fixed-pressure shift uses —
and §3.3.1 makes that an **edge** quantity, since its coefficients and expansion point are averages
over the two cells of the edge. A cell-based $Z$ cannot carry an edge-dependent $\hat\alpha$ however
it is integrated. This is §3.5's source 4, and it is $O(\alpha_{pp}(\Delta p^{\text{mid}}_e)^2)$ per
layer, accumulating down the column. §3.5 gives its closed form: the edge-shared profile's layer mean
is [](#alpha-edge-layer-mean), the per-layer mismatch against `VertCoord`'s increment is the
second-order equation-of-state remainder [](#eos-remainder), and what the scheme needs is $\Delta_e$
of its running sum down the column.

**Resolution: the PGF accumulates $\Delta_e$ of $\varepsilon_{i,k}$ down each edge's column**, built
from the edge-shared $\hat\alpha$, and adds it to the local terms of $\mathcal{R}$. This does not
create the second, slightly different $z$ that Requirement 2.4 exists to prevent:
`GeomZInterface`/`GeomZMid` remain the model's one geometric height and are unchanged, and what the
PGF accumulates is a per-edge *difference* of second-order EOS remainders, identically zero when the
two columns' states coincide.

Two things about it are easy to get wrong.

- **This piece is genuinely small; the scheme as a whole is not a small correction.** Source 4 is
  second order in the cross-edge pressure contrast, so accumulating it consumes little precision. But
  sources 1–3, and $\mathcal{S}$ itself, are the same size as the terms they correct — at 50 m/km both
  are $\approx 0.5$ m s$^{-2}$ before cancellation. Writing something as a "correction" therefore buys
  no precision on its own; what buys it is [](#gamma-increments) and §3.5.1, which remove the
  cancellation algebraically.
- **The scan is not optional and not a detail.** §4.1.3 records the implementation consequence: a
  column prefix sum cannot live inside a per-edge, per-vertical-chunk kernel call. Both this scan and
  $\mathcal{S}$'s accumulation via [](#gamma-increments) are prefix sums and should share one pass.

One related question is left to implementation time: `VertCoord` builds $z$ upward from the
bathymetry while pressure is built downward from the surface, so the two accumulate round-off from
opposite ends of the column. Which end the scan of this section and the accumulation of
[](#gamma-increments) run from is a round-off question (§3.7.5), not a consistency one — condition 3
is satisfied either way. It is now a single decision rather than two, since both should share one
pass.

#### 3.7.5 Round-off in the deep ocean

Separate from everything above, and not fixed by it, is how much precision the cancellation
consumes. At 4000 m the hydrostatic height difference across a steeply tilted edge is of order
$10^{2}$ m, while the part of it that survives in the edge tendency is of order
$10^{-3}$ m — the millimetre-scale residual of [](#centered-shift). Roughly five significant digits
are consumed there, and roughly ten by the time the tendency is formed against a signal of order
$10^{-6}\ \mathrm{m\,s^{-2}}$ times $d_e$. Double precision leaves adequate margin; a
`OMEGA_SINGLE_PRECISION` build does not, and the machine-precision cancellation of §3.7.2 would
simply be invisible beneath round-off.

Three points.

- **This exposure belongs to `PressureGradCentered`, not to the new scheme.** It is entirely in
  $\mathcal{S}_{e,k}$ as [](#centered-shift) writes it — two terms of order 200 m cancelling to
  millimetres — which is what the default PGF has always computed. $\mathcal{R}$ adds none of its
  own: requirement R2 of §3.5 exists precisely so that every one of its four sources is formed from
  quantities that are individually small. Any conclusion reached here therefore applies to the
  centered scheme as well, and it would be wrong to treat a single-precision round-off floor as an
  objection to `FiniteVolume` specifically.
- **§3.5.1 removes it, in both schemes.** Accumulating $\mathcal{S}$ through
  [](#gamma-increments) never forms the two large terms at all. If that form is adopted — and the
  design recommends it — the ten-digit estimate above does not apply to the implemented scheme, and
  the expected answer to the single-precision question changes from "probably fails" to "probably
  passes". That is a prediction to be measured (§5.2, run in both precisions), not an assumption.
- **The perturbation form remains in reserve and is probably not needed.** It would subtract a local
  reference profile from $\hat\alpha$ before integrating and add its contribution back analytically.
  This is *not* a Shchepetkin–McWilliams (2003)-style reference-profile subtraction in the usual
  sense: the reference is local to the edge and layer, its contribution cancels identically rather
  than approximately, and the accuracy of the scheme does not depend at all on how well it matches
  the local column — so it does not degrade in the strong-gradient, steep-layer regions where a
  global reference profile would. It changes nothing in exact arithmetic. With
  [](#gamma-increments) in place it addresses a problem that has already been removed
  algebraically, so it should be adopted only against a measurement.

#### 3.7.6 Assumptions that still need testing

Everything above concerns properties the implementation can be built to have. The following are
assumptions this design is making that only testing can confirm:

- **A1 — Accuracy at cell centers does not imply accuracy at edges.** Reconstructing $\Theta$ and
  $S$ accurately at each cell center does not by itself make $z_i(p)$, and hence
  $[\nabla_n z]_p$, accurate at the edge between them, because $\hat\alpha(\Theta,S,p)$ is a
  nonlinear function of the reconstructed quantities. The horizontal order claimed in §3.6 is
  therefore a design target to be confirmed by measured convergence (Requirement 2.6), not something
  that follows from the cell-centered reconstruction order. Note that A1 is now purely an *accuracy*
  assumption: the reformulation of §3.2 removes its consequences for the machine-precision property,
  which is why §3.6.1's stencil constraint could be dropped.
- **A2 — The EOS expansion is good enough across the edge.** [](#alpha-taylor) is expanded about the
  shared edge state [](#edge-ref), and its error grows with the horizontal contrast in $\Theta$ and
  $S$ across the edge and with the layer's pressure range. Whether the first-order expansion suffices
  in frontal regions, or whether the second-order option of §3.3 is needed, is an open question.
  A2 is now measurable more directly than by an accuracy sweep alone: §3.5's source 4 *is* the
  second-order remainder of that expansion, so the magnitude of the accumulated part of
  $\mathcal{R}$ relative to its local part is a running diagnostic of how hard the expansion is being
  worked. Reporting that ratio is cheap and should be done.
- **A3 — The residual outside the exact set is small enough in practice.** The $O(\tilde h^2)$ and
  $O(\tilde h^3)$ entries in §3.7.3 describe how the error *scales*; how large it actually is at the
  vertical resolutions Omega can afford is unknown.
- **A4 — PGF error is what is driving the observed instability.** Identifying [](#centered-error)
  with bottom-layer noise in realistic global runs is a plausible diagnosis, not a demonstrated one.
  A configuration whose profile falls in Phase 1's exact set gives a direct test: if spurious
  bottom-layer flow persists there, the cause lies elsewhere — most likely in the layer-mean
  treatment in the tracer and remapping operators (§4.3) — and Phase 1 will not cure it. See §5.3.
- **A5 — The treatment at the top and bottom of the column is adequate.** [](#ho-exact) averages over
  the edge control volume's pressure range $[p^{\text{top}}_{e,k}, p^{\text{bot}}_{e,k}]$, which in
  general is not either column's own range, so each column's reconstruction must be evaluated
  slightly outside the layer it was built for. In the interior this is a small extrapolation across a
  smooth reconstruction and is harmless. At the surface, where the two columns have different
  $p^{\text{surf}}$, and at the sea floor, where they may have different `maxLevelCell` and where a
  partial cell carries the whole signal for the Polaris `bathymetry_step` configuration, it is not
  obviously harmless and **has not been tested** — the reference derivation was run with equal layer
  counts in both columns. The extrapolation must be defined explicitly, and whichever choice is made
  must reduce to a single function when the two reconstructions coincide, or condition 1 of §3.7.2
  fails at the column ends. §5.1's `bathymetry_step` and `surface_pressure_gradient` variants are
  where it is measured.

### 3.8 Scope: the reduction and the vertical coordinate

The reduction of §3.1.1 rests on one property of Omega — that $\tilde z$ is a function of pressure
alone — so it is worth recording precisely how much freedom that leaves and where the boundary is.

**It holds for any ALE layering of pseudo-height.** p-star, sigma-p, z-tilde-star and hybrid schemes
all relocate layer interfaces within $\tilde z$, in space and in time, and none of them redefines
$\tilde z$. Interfaces may be arbitrarily tilted, arbitrarily thin, and at wholly different pressures
in neighbouring columns; [](#grad-p-zero) is untouched, because it is a statement about coordinate
surfaces (§3.1.1). Every term the scheme computes is a tilt term, so there is no configuration in
which the reduction quietly turns the scheme into a no-op.

**It fails only if Omega abandons pseudo-height.** A coordinate $r$ that is not a function of pressure
alone — geometric $z$, a true $z$-star, or a genuine isopycnal coordinate — gives $\nabla_r p \neq 0$,
and then all four pressure terms of [](#ho-target) become real work rather than a discretization of
zero. Should such a change ever be contemplated, [](#ho-target) as corrected here and the derivation in
§3.1 are the starting point, and §3.1.1's third consequence is the pitfall to carry forward with them.

**The geopotential body force** — the first line of [](#ho-target) — is the whole target. It is
$g\,\nabla z + \nabla(\phi_{TP}+\phi_{SAL})$; $z$ enters as `VertCoord`'s `GeomZInterface` through
[](#centered-shift), and the tidal-potential and self-attraction-and-loading contributions are
supplied by `VertCoord` and differenced with [](#edge-grad) (Requirement 2.7). They are zero in early
Omega versions and are unaffected by anything in §3.5.

### 3.9 Reduction to the centered scheme

Setting $\mathcal{R}_{e,k} \equiv 0$ in [](#ho-discrete) recovers `PressureGradCentered` exactly.
That is the content of [](#centered-shift), which is an algebraic identity about the implemented
code rather than a limit reached by tuning options, so the reduction is **structural in the strict
sense**: there is one term, the reduction is switching it off, and it holds at any tilt, any layer
thickness and any stratification.

The reduction is therefore exercised by a single configuration switch, `RemainderEnable: false`
(§4.1.1), which doubles as guard test (c) of §5.2. The agreement is algebraic rather than bit-for-bit
if §3.5.1's accumulation is used, since that evaluates $\mathcal{S}$ in a different order from
`PGrad.h` (§4.4); it would be bit-for-bit if [](#centered-shift) were assembled literally, and that is
a deliberate trade against round-off.

Because the reduction holds at any tilt, §5.5 should be run *at* nonzero tilt and over a sweep of
tilts. That is not merely permitted, it is the stronger test: the whole of the difference between the
two schemes is $\mathcal{R}$, so a tilt sweep confirms that $\mathcal{S}$ alone reproduces
`PressureGradCentered` and that both schemes are in fact responding to tilt. A reduction demonstrated
only on level interfaces would establish much less.

### 3.10 Per-step algorithm summary

Steps 1–3 are per cell and layer; step 4 is a per-edge column scan; step 5 is per edge and layer.
Phase differences are marked.

1. From `VertCoord`: read `PressureInterface`, `PressureMid`, `GeomZInterface`, `PseudoThickness`,
   and the tidal/SAL potentials (already computed diagnostically each step). `VertCoord` itself
   needs no change and no baseline step (§3.7.4).
2. Obtain $\alpha_0$ (= the existing `Eos::SpecVol` field) and the derivatives
   $\alpha_\Theta, \alpha_S, \alpha_p$ from one TEOS-10 evaluation ([](#alpha-derivs)). Both phases
   need all four; Phase 2 optionally adds second derivatives (§3.3).
3. Build the mean-preserving deviations $\Theta', S'$ ([](#vert-recon)) — **linear in Phase 1,
   parabolic in Phase 2** — using the actual non-uniform interface pressures (§3.4).
4. For each edge, in one column scan: form the shared expansion point [](#edge-ref) layer by layer;
   accumulate $\mathcal{S}_{e,k}$ through [](#gamma-increments) and the source-4 part of
   $\mathcal{R}_{e,k}$ through §3.7.4, from the same end of the column (§3.7.4 closing note).
5. For each edge and layer: add the local sources 1–3 of $\mathcal{R}_{e,k}$ (§3.5), assemble
   [](#ho-discrete) with the tidal/SAL difference, and accumulate into `Tend` with `EdgeMask`.
   **Phase 1** uses the two-cell operator [](#edge-grad); **Phase 2** uses the wider stencil of
   §3.6.1.

## 4 Design

The `FiniteVolume` scheme extends the existing `PressureGrad` class in
`components/omega/src/ocn/PGrad.{h,cpp}`. The class manager, creation/retrieval/removal
methods, named-instance map, and dispatch logic are unchanged; the work is to flesh out the
`PressureGradFiniteVolume` functor (currently a no-op placeholder) and to widen the inputs the
manager hands to it.

### 4.1 Data types and parameters

#### 4.1.1 Parameters

Configuration lives under the existing `PressureGrad` YAML group. The `PressureGradType`
enum (PGrad.h:22) is updated to give the new scheme a descriptive name and to drop
the unimplemented second placeholder, leaving a commented stub for a future variant:

```c++
enum class PressureGradType {
   Centered,      // existing 2nd-order Montgomery scheme
   FiniteVolume   // layer-integrated finite-volume scheme (this design)
   // , <FutureVariant>   // e.g. a 6th-order option, added when implemented
};
```

New sub-options for `FiniteVolume`. Both phases use the same keys; Phase 2 adds values
rather than keys, so no configuration written for Phase 1 needs to change when Phase 2 lands:

```yaml
    PressureGrad:
       PressureGradType: 'FiniteVolume'   # Centered | FiniteVolume
       HorzOrder: 2                       # 2 = two-cell stencil (Phase 1); 4 = wide stencil (Phase 2)
       VerticalReconstruction: 'linear'   # 'linear' (Phase 1) | 'ppm' (Phase 2)
       RemainderEnable: true              # false = verification only, reduces to Centered (§3.9)
```

- **Phase 1 default and target:** `HorzOrder: 2`, `VerticalReconstruction: 'linear'`,
  `RemainderEnable: true`.
- **Phase 2 target:** `HorzOrder: 4`, `VerticalReconstruction: 'ppm'`.

`RemainderEnable: false` switches off $\mathcal{R}_{e,k}$, which by [](#centered-shift) leaves
exactly `PressureGradCentered` (§3.9). It is **verification only**: it supports the permanent
regression test of §5.5 and serves as guard test (c) of §5.2. It is not a supported production
setting — it *is* the default scheme, reached by a more expensive path — and the implementation
should log a warning if it is selected outside a test.

`HorzOrder` selects the width of the edge stencil used by [](#edge-grad). It is purely an accuracy
setting: it carries no constraint from the robustness property (§3.6.1).

There is no vertical quadrature setting: the one integrand in [](#ho-discrete) is nearly constant
across the layer and is evaluated in closed form (§3.2). `FiniteVolume` names the layer-integrated
control-volume form of §3.1 rather than a face-by-face assembly of the pressure traction, which is
the narrower sense Adcroft et al. (2008) use.

The configuration that isolates compressibility on its own is `'linear'` applied to a vertically uniform
$\Theta$, $S$ profile — the reconstruction slopes are then zero and $\alpha_p$ is retained — and it
needs no setting of its own. Tests that mean to isolate compressibility, including §5.2, select
`'linear'`.

#### 4.1.2 New EOS support

`PressureGradFiniteVolume` needs $\alpha$ together with its first derivatives. The `Eos` class
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

#### 4.1.3 `PressureGradFiniteVolume` functor

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

(The existing centered signature is unchanged.) Additional cached members hold, in Phase 2, the
wide-stencil cell lists and weights. The functor implements step 5 of §3.10 per edge and vertical
chunk, accumulating into `Tend` with `EdgeMask`, exactly as the centered functor does.

Three aspects of the loop structure differ from the obvious implementation:

- **Column quantities built from $\hat\alpha$ are computed per edge, not cached per cell.** Both
  adjacent columns' contributions use equation-of-state coefficients formed at the edge
  [](#edge-ref), so the same cell yields a different result at each of its edges. The per-cell
  quantities that *can* be cached are the four EOS coefficients and the reconstruction slopes, which
  is where the TEOS-10 cost lives; the per-edge work is polynomial arithmetic on those cached
  values. This is what keeps Requirement 2.2 satisfied despite roughly three times as many column
  evaluations on a hexagonal mesh.
- **The column scan of §3.10 step 4 cannot live in this functor.** $\mathcal{S}_{e,k}$ via
  [](#gamma-increments) and source 4 of $\mathcal{R}_{e,k}$ are both prefix sums down each column
  with edge-dependent coefficients, so neither is expressible as an independent per-vertical-chunk
  operation. They are computed together in a separate kernel filling two per-edge arrays,

  ```c++
  Array2DReal ShiftFirstOrder;  ///< (NEdgesAll, NVertLayers) — S_{e,k}, eq. (centered-shift)
  Array2DReal RemainderAccum;   ///< (NEdgesAll, NVertLayers) — source 4 of R_{e,k}, §3.7.4
  ```

  with a `parallelForOuter` over edges and a `parallelScanInner` down the column, in the same shape
  as `VertCoord::computeGeomZHeight`. The functor then reads them chunk-wise like any other input.
  The cost is two edge-sized 2-D arrays and one column scan per edge per step. This is the one
  structural addition Phase 1 makes beyond the per-edge, per-chunk pattern the centered scheme uses,
  and it is not optional: §3.5's requirement R2 rules out forming either quantity by subtracting
  large numbers inside the chunk loop.
- **$\mathcal{R}$'s local sources are the only part that is a plain per-chunk computation.** Sources
  1–3 of §3.5 depend on the layer's own interface pressures, thicknesses and reconstruction slopes,
  and belong in the functor.

The Phase 1 and Phase 2 code paths differ only in the reconstruction degree (§3.4) and in the width
of the edge stencil (§3.6.1). There is one functor, not two.

### 4.2 Methods

`PressureGrad::computePressureGrad` keeps its role of selecting the configured option, but
its input list grows so the `FiniteVolume` branch can reach $\Theta$, $S$, and the EOS
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
`FiniteVolume` branch dispatches to `PressureGradFiniteVolume` through the same
`parallelForOuter`/`parallelForInner` team pattern used today (PGrad.cpp). The call site in
`Tendencies.cpp` is updated to pass $\Theta$, $S$ (from the tracer state) and the `Eos`
instance, which it already references for `SpecVol`.

Creation, retrieval, and removal (`init`, `create`, `get`, `getDefault`, `clear`, `erase`)
are unchanged. The `FiniteVolume` constructor additionally caches the reconstruction stencil and
quadrature weights from config.

### 4.3 Consistency and follow-up work

This PGF is deliberately more accurate in the vertical than the rest of Omega's
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
than reimplemented as the lowest-order configuration of `PressureGradFiniteVolume`, even though
§3.9 shows the latter reduces to it. The redundancy is small — the centered functor is a
header-only, ~40-line loop body with no supporting machinery of its own — and it buys two
things that a single implementation cannot provide:

1. **An independent cross-check.** The two functors read the mesh, `VertCoord`, and EOS state
   through separately written code. Their agreement to round-off (§5.5) therefore tests the
   shared upstream state — edge masks, interface indexing, `VertCoord` conventions — and not
   just the PGF arithmetic. Collapsing them would make that comparison self-referential: a
   defect upstream of the order switch would appear identically in both limits and cancel.
2. **A stable default.** The algebraic reduction of §3.9 does not imply bit-for-bit agreement,
   because the reduced `FiniteVolume` path performs the same operations in a different order.
   Replacing the default PGF would be an answer-changing change for every existing
   configuration, which this design does not require.

Removing `PressureGradCentered` becomes reasonable once `FiniteVolume` is promoted to the
default and has served out its period as the reference implementation. That is deliberately
left as **follow-up work**, to be taken up as a separate, answer-changing change with its own
baseline step — not as part of this design.

### 4.5 What each phase delivers and depends on

#### 4.5.1 Phase 1 — a consistent second-order scheme

**Delivers.** The exact fixed-pressure geopotential difference [](#ho-exact), evaluated as
`PressureGradCentered` plus the remainder $\mathcal{R}_{e,k}$ (§3.5); the equation-of-state expansion
about a state shared across each edge (§3.3.1); mean-preserving linear reconstruction of $\Theta$ and
$S$ in pressure (§3.4); the per-edge column scan that supplies $\mathcal{S}$ and $\mathcal{R}$'s
accumulated part (§3.7.4, §4.1.3). The result is a pressure gradient that is zero to machine
precision for any resting ocean whose profile varies linearly with pressure, at any tilt, thickness,
or bathymetry (§3.7.3) — Requirement 2.3 in full.

**Does not deliver.** Fourth-order accuracy. The horizontal operator is the same two-cell stencil
in use today (§3.6), so horizontal truncation error is unchanged from `PressureGradCentered`.
Requirements 2.1 and 2.6 are met at second order only.

**Depends on.** Nothing unresolved. The `VertCoord` question of §3.7.4 is settled there: no change to
`VertCoord` is needed and no baseline step is required, because the midpoint rule is already the exact
layer integral of a Phase 1 $\hat\alpha$. What condition 3 does require — the column scan — is part
of this phase's own implementation (§4.1.3), not a prerequisite in another module.

**Open at implementation time.** Two items, neither of them blocking. Which end of the column the
scan accumulates from (§3.7.4 closing note) is a round-off question, settled by measurement. And the
treatment at the top and bottom of the column, assumption A5 of §3.7.6, has to be defined explicitly
rather than inherited — it is the one part of the formulation the reference derivation did not cover,
since it ran with equal layer counts in both columns.

**Code and cost.** Three new `Eos` derivative fields and one new method (§4.1.2); the
`PressureGradFiniteVolume` functor; two per-edge arrays and one column scan (§3.7.4, §4.1.3); no new
TEOS-10 evaluations per cell and layer (Requirement 2.2), with roughly three times as many column
evaluations on a hexagonal mesh (§4.1.3). Most of the arithmetic is in the column scan.

#### 4.5.2 Phase 2 — fourth order

**Delivers.** Parabolic vertical reconstruction (§3.4), widening the exact set to profiles that vary
quadratically with pressure; a wide horizontal stencil (§3.6.1); optionally a second-order
equation-of-state expansion (§3.3). Requirements 2.1 and 2.6 at fourth order.

**Must preserve.** Everything Phase 1 establishes. The machine-precision property is inherited by any
consistent edge-normal stencil (§3.6.1) rather than having to be re-established, so what Phase 2 must
not break is narrower than it appears: it must keep the expansion point shared across whatever cells
the stencil couples, and it must re-examine [](#z-increment-exact) before adopting the second-order
equation-of-state expansion, which would otherwise reopen the §3.7.4 question. The §5.2 gate is rerun
unchanged.

**No blocking open question.** The form of the wide stencil on Omega's unstructured,
variable-resolution TRiSK mesh is an ordinary accuracy-and-cost choice, settled by measurement against
Requirement 2.6, because exactness does not constrain it (§3.6.1).

#### 4.5.3 Suggested order of work

1. **Write the discrete form out in full, in the Polaris two-column harness, before writing C++.**
   Doing it there makes it executable — exactness becomes a pytest that runs in seconds on a login
   node — and doubles as the independent reference implementation §5.1 needs, which it would not be
   if written by reading the finished C++. *Substantially done*: the derivation of §3.1, the
   corrected [](#ho-target), the cancellation [](#grad-p-zero) and the numerical zero for a profile
   linear in pressure at large tilt have all been established this way. What remains is to re-express
   the confirmed result in the reduced form of §3.5, and to settle four things this design leaves to
   it:

   - the closed form of each of $\mathcal{R}$'s four sources, with requirement R2 of §3.5 satisfied —
     no source formed by subtracting large numbers;
   - the reconciliation of [](#gamma-increments) against [](#centered-error), including the anchor
     constant (§3.5.1);
   - which end of the column the scan accumulates from (§3.7.4);
   - the top- and bottom-of-column treatment, assumption A5 (§3.7.6), which the existing derivation
     did not cover.
2. Run the assumption-A4 diagnostic of §5.3, which uses the *existing* centered scheme and so can be
   done immediately and in parallel with step 1. If spurious bottom-layer flow survives a profile
   that Phase 1 would resolve exactly, the cause is elsewhere in the model and the priority of this
   work should be reconsidered before it is built.
3. Implement and verify Phase 1 against §5.1, §5.2, §5.3, and §5.5.
4. Take up Phase 2. Re-examine [](#z-increment-exact) before adopting the second-order
   equation-of-state expansion; that, and not the stencil, is the item that could reopen a settled
   question.

## 5 Verification and Testing

Testing reuses and extends the Polaris `horiz_press_grad` task family
(`polaris/tasks/ocean/horiz_press_grad`), which is already Omega-only and built around a
quasi-analytic, layer-mean TEOS-10 reference solution (`reference.py`, surface-anchored,
4-point Gauss). That analytic reference remains valid as "truth" for `FiniteVolume`,
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
  - `finite_volume_phase1` — `HorzOrder: 2`, `VerticalReconstruction: linear`;
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
  layer-mean analytic HPGA (unchanged). For `FiniteVolume` the layer-mean comparison
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
  `VerticalReconstruction: 'linear'` (§4.1.1). **Pass:** the PGF is zero at every edge and layer to
  **machine precision** (double-precision builds; the threshold tracks `Real`'s epsilon and the size
  of the hydrostatic terms, not a physical tolerance). Phase 1 and Phase 2 must both pass.

  Assert this in the sharper form of §3.5's requirement R3 — that
  $\mathcal{R}_{e,k} = -\mathcal{S}_{e,k}$ to round-off — and report both quantities, not only their
  sum. $\mathcal{S}$ is the centered scheme's answer on the same state and is nonzero and
  tilt-dependent; checking against it turns a comparison with zero into a cancellation of ten orders
  of magnitude against a known reference, and it is what distinguishes a working scheme from one that
  has stopped responding to tilt at all.
- **Profiles it does not:** $\Theta$, $S$ quadratic in pressure, then a realistic profile.
  **Pass:** the residual shrinks like $\tilde h^2$ (Phase 1) and $\tilde h^3$ (Phase 2) as the
  vertical grid is refined at fixed tilt, matching §3.7.3. A residual that does not shrink at the
  tabulated rate means one of the three conditions in §3.7.2 has been broken somewhere in the
  implementation; it is a bug to find, not a tolerance to widen.
- **Guard tests.** Rerun an exactly resolved profile with each of the following, and require each to
  *fail* the machine-precision check. Without them, a passing result could as easily come from a
  symmetry of the test setup, or from the scheme having become insensitive to tilt, as from the scheme
  being right.

  | guard | breaks | expected to fire because |
  |---|---|---|
  | (a) **tilt sensitivity** — assert $\mathcal{S}_{e,k}$ and $\mathcal{R}_{e,k}$ are *individually* nonzero and grow with tilt | nothing; this one runs on the passing configuration | a bug that zeroed the tilt response would satisfy every other check in this section perfectly. This is the guard against mistaking "no tilt terms" for "correct tilt terms" |
  | (b) cell-local expansion point in place of [](#edge-ref) | condition 1, §3.7.2 | the two columns then use different approximations to one equation of state, so $z_i(p)$ differ by more than a constant (§3.3.1) |
  | (c) `RemainderEnable: false` | conditions 2 and 3 | this is exactly `PressureGradCentered` (§3.9), whose error on this state is the thing Phase 1 removes |
  | (d) source 4 of $\mathcal{R}$ omitted, local sources retained | condition 3 only | isolates the column scan: a residual that grows with depth rather than being flat down the column |
  | (e) $\mathcal{R}$'s local sources truncated at first order | condition 2 only | isolates the fixed-pressure shift, leaving a residual local to each layer |

  Guards (d) and (e) are the pair that makes a failure diagnosable: condition 2 and condition 3 fail
  in distinguishable ways (§3.7.2), and running both tells the implementer which one is broken without
  another round trip. Guard (a) is the one to write first, because it is the only one that can fire on
  a configuration where every other check passes.

  Guards must be checked against a deliberately broken configuration, not only against a passing
  one. A guard that cannot fire is worse than no guard, because it looks like protection.

A separate and much smaller unit test covers the reconstruction estimator on its own (§3.4): given
layer means sampled from a profile of the reconstruction's own degree on a **deliberately
non-uniform** set of layer thicknesses, the recovered slope (Phase 1) or slope and curvature
(Phase 2) must match the exact values to round-off. This is worth testing separately because it is
the most likely place for the machine-precision gate above to fail, and it localizes the failure
immediately.

The test is also run in a single-precision build, to measure the round-off floor of §3.7.5. Note that
§3.7.5 now predicts this should *pass* if [](#gamma-increments) is used for $\mathcal{S}$, and that
the same measurement applies to `PressureGradCentered`; running the centered scheme alongside is what
makes the result interpretable.

Condition 3 of §3.7.2 requires the column scan of §3.7.4, so the machine-precision check cannot pass
until that is in place; this test, together with guards (c) and (d) above, is its acceptance
criterion. It is implemented as a fast C++ unit test and also run as a configuration of the Polaris
two-column task.

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

Configure the new scheme with `HorzOrder: 2`, `RemainderEnable: false` and confirm it reproduces
`PressureGradCentered` to round-off on the two-column test. This guards Requirement 2.5 and protects
the existing default during refactoring. The test is established with Phase 1 and rerun unchanged for
Phase 2, where it also confirms the wider stencil collapses correctly to the two-cell one.

**Run it at nonzero tilt, and at several tilts.** Per §3.9 the reduction is exact at any tilt, so a
level-interface configuration would be a weaker test than the design supports — and a tilt sweep here
is what confirms that $\mathcal{S}$ alone is reproducing the centered scheme rather than the two
schemes agreeing because neither is responding to tilt.

Note that `RemainderEnable: false` is not a supported production setting (§4.1.1); it exists so that
this comparison is possible. Keeping it costs one branch.

This test is retained permanently rather than treated as a one-time transition check, and that choice
is the reason `PressureGradCentered` is kept as a separate implementation rather than reimplemented as
the lowest-order configuration of `PressureGradFiniteVolume` (§4.4). The two functors read the mesh,
`VertCoord`, and EOS state through independently written code, so their agreement to round-off is
evidence about *shared upstream* state as well as about the PGF arithmetic: a wrong edge mask, a
mis-indexed interface array, or a misinterpreted `VertCoord` convention shows up as a disagreement.
Were the centered scheme replaced by an order-2 configuration of the new code, this comparison would
reduce to comparing an implementation against itself, and any defect upstream of the order switch
would cancel out of it.

**Covers:** Requirement 2.5; §3.9.

### 5.6 Cost check

Requirement 2.2 bounds the number of TEOS-10 evaluations, and nothing above tests it — the scheme
could satisfy every accuracy gate while quietly calling the equation of state inside the column loop.
Two cheap checks close that:

- **Evaluation count.** With an instrumented `Eos`, confirm the number of specific-volume
  evaluations per time step is one per cell per layer and is **unchanged** when `HorzOrder` and
  `VerticalReconstruction` are varied. This is the property Requirement 2.2 actually states, and it is
  a counter comparison, not a timing measurement, so it is deterministic and suitable for CI.
- **Wall time.** Record PGF kernel time relative to `PressureGradCentered` on a representative
  configuration, as a performance regression guard. The expected cost is dominated by the per-edge
  column scan (§4.1.3), roughly three times as many column evaluations as a cell-based formulation on
  a hexagonal mesh; a result far above that suggests the per-cell coefficients are being recomputed
  per edge rather than cached.

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
| A2 EOS expansion adequate across an edge | §5.1 contrast sweep; §5.4; the source-4 diagnostic of §3.7.6 |
| A3 Residual small enough in practice | §5.1 accuracy gate; §5.3; §5.4 |
| A4 PGF error causes the instability | §5.3 diagnostic, run before Phase 1 completes |
| A5 Top- and bottom-of-column treatment | §5.1 `bathymetry_step` and `surface_pressure_gradient` |
| §3.5 remainder is load-bearing | §5.2 guards (c), (d), (e); §5.5 at nonzero tilt |
| §3.5 the scheme responds to tilt at all | §5.2 guard (a); §5.5 tilt sweep |
| §3.7.4 column scan | §5.2 (cannot pass without it; guard (d) isolates it) |
| §3.7.5 Round-off floor | §5.2, run in both precisions, alongside `Centered` |

Requirement 2.7 (extensibility) is not testable directly; it is addressed by the configuration
design of §4.1.1 and by the phase structure of §4.5, which is itself the evidence that the framework
extends without restructuring.

## 6 Supporting evidence

Two things §3 asserts are not self-evident from the algebra alone: that the cancellation of §3.1.1
survives discretization, and that the form of [](#ho-target) used here is the correct one. This
section supplies the evidence for both.

### 6.1 Numerical confirmation of the reduction and the exact set

The claims below were established in a purpose-built two-column Python harness before any C++ was
written, and they are the reason §3.2 takes [](#ho-exact) as its target rather than assembling the
four pressure terms. The harness builds its state on the Polaris `horiz_press_grad` initialization:
interface pressures prescribed directly, $\Theta$ and $S$ sampled as `Init._interpolate_t_s` does,
equation-of-state coefficients from one `gsw.specvol_first_derivatives` call per cell per layer, $z$
accumulated as `VertCoord` does, geometry from `hydrostatic_consistency.cfg`
($d_e = 4$ km, $\tilde z_{\text{bot}} = -4096$ m). Where a truth value is needed it is
$-g\,\partial z/\partial x|_p$ by direct segment-wise Gauss integration of the true TEOS-10
$\alpha(p)$, sharing no code with the scheme under test.

| claim | §ref | measurement |
|---|---|---|
| the four pressure terms of [](#ho-target) sum to zero | §3.1.1 | $1.4$–$5.0\times10^{-15}$ m s$^{-2}$ against a side-wall term of $0.55$; relative $3\times10^{-15}$, at 50 m/km tilt, for profiles linear in $p$, quadratic in $p$, and with a horizontal temperature contrast |
| the $p\nabla\alpha$ term of [](#ho-target) is required | §3.1 | with it, error $1.45\times10^{-9} \to 2.26\times10^{-11}$ at rate $2.00$ per refinement; without it, error pinned at $4.55\times10^{-6}$ at rate $0.00$ — the scheme returns $9.0\times10^{-7}$ against a true $5.42\times10^{-6}$ |
| exactness on the exact set, independent of tilt | §3.7.3 | residual at the round-off floor, $9\times10^{-16}$ to $9\times10^{-15}$, **not growing** across 0.5–200 m/km at 16 and 64 layers. `Centered` on the same states is exactly first order in tilt |
| $O(\tilde h^2)$ off the exact set | §3.7.3 | quadratic-in-$p$ profile at 50 m/km: rates 2.21, 2.14, 2.09, 2.05 under vertical refinement |
| the edge-shared expansion point is load-bearing | §3.3.1 | replacing it with a cell-local one costs $2.2\times10^{-9}$ to $3.0\times10^{-5}$ across the tilt range — §5.2 guard (b) |
| the fixed-pressure comparison is the whole scheme, not a correction | §3.5, §3.7.4 | omitting it costs $2.8\times10^{-3}$ to $1.1$ m s$^{-2}$ — three to five orders larger than any other guard, and the reason §3.5 presents $\mathcal{R}$ as a component rather than a refinement |

Three limits on what this establishes, all of them live.

- **Differing `maxLevelCell` is untested.** The harness gave both columns the same layer count, so the
  `bathymetry_step` case — where a bottom partial cell carries the whole signal — is not covered. This
  is assumption A5 of §3.7.6.
- **Nothing has been run in Omega**, or through the real Polaris task machinery, so the interaction
  with p-star initialization and partial cells is untested.
- **The rows above were measured on an assembly of all four pressure terms**, not on [](#ho-discrete).
  Given that those terms cancel to $3\times10^{-15}$ relative, what the exactness and $O(\tilde h^2)$
  rows measured *was* the geopotential term alone, so they carry over — but by that argument rather
  than by direct measurement, and confirming it against [](#ho-discrete) is the first task of §4.5.3
  step 1.

### 6.2 Corrections owed to `OmegaV1GoverningEqns.md`

[](#ho-target) is inherited from {ref}`omega-design-governing-eqns-omega1`, which is not maintained
against the code. Three corrections are owed to it, all following from [](#leibniz), and they should
be made there rather than left as local overrides here.

1. **The $+\int_A\int\rho_0\,p\,\nabla\alpha\,d\tilde z\,dA$ term is missing.** Row 2 of §6.1 is what
   it costs: a scheme built on the four-term form without it does not converge at all. The omission is
   invisible in every resting-state configuration, because $\partial_x\alpha|_{\tilde z} = 0$ when
   $\Theta$ and $S$ are horizontally uniform, which is why it survived review.
2. **The signs on the two metric terms are reversed.** [](#leibniz) gives
   $+(\alpha p)_{\text{top}}\nabla\tilde z^{\text{top}}$ and
   $-(\alpha p)_{\text{bot}}\nabla\tilde z^{\text{bot}}$. With the signs as printed there, the
   cancellation of §3.1.1 doubles rather than vanishing.
3. **The small-slope interface normals are not needed.** Its §6 reaches [](#ho-target) through
   $\mathbf n^{\text{top}} \approx (-\nabla\tilde z^{\text{top}}, 1)$, unnormalized. The derivation in
   §3.1 here uses Leibniz' rule directly, so neither the corrected [](#ho-target) nor [](#ho-exact)
   depends on that approximation.
