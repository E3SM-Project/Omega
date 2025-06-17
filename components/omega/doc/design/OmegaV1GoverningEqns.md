# Omega V1: Governing Equations

<!--
Add later, if it seems necessary. There is a toc on the left bar.
**Table of Contents**
1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Algorithmic Formulation](#3-algorithmic-formulation)
4. [Design](#4-design)
5. [Verification and Testing](#5-verification-and-testing)
-->


## 1. Overview

This design document describes the governing equations for Omega, the Ocean Model for E3SM Global Applications. Overall, Omega is an unstructured-mesh ocean model based on TRiSK numerical methods ([Thuburn et al. 2009](https://www.sciencedirect.com/science/article/pii/S0021999109004434)) that is specifically designed for modern exascale computing architectures. The algorithms in Omega will be mostly identical to those in MPAS-Ocean, but it will be written in c++ rather than Fortran in order to take advantage of the Kokkos performance portability library to run on GPUs ([Trott et al. 2022](https://ieeexplore.ieee.org/document/9485033)). Significant differences between MPAS-Ocean and Omega are:

1. Omega is non-Boussinesq. This means that the full 3D density is used everywhere, and results in a mass-conserving model. MPAS-Ocean and POP were Boussinesq, so that a reference density $\rho_0$ is used in the pressure gradient term, and were therefore volume-conserving models. In Omega the layered mass-conservation equation is in terms of pressure-thickness ($h=\rho g \Delta z$). In MPAS-Ocean the simple thickness ($h=\Delta z$) is the prognostic volume variable (normalized by horizontal cell area).
1. Omega will use the updated equation of state TEOS10, while MPAS-Ocean used the Jackett-McDougall equation of state.

The planned versions of Omega are:

- **Omega-0: Shallow water equations with identical vertical layers and inactive tracers.** In his first version, there is no vertical transport or advection. The tracer equation is horizontal advection-diffusion, but tracers do not feed back to dynamics. Pressure gradient is simply gradient of sea surface height. Capability is similar to [Ringler et al. 2010](https://www.sciencedirect.com/science/article/pii/S0021999109006780)
- **Omega-1.0: Layered ocean, idealized, no surface fluxes.** This adds active temperature, salinity, and density as a function of pressure in the vertical. Vertical advection and diffusion terms are added to the momentum and tracer equations. An equation of state and simple vertical mixing, such as constant coefficient, are needed. Capability and testing are similar to [Petersen et al. 2015](http://www.sciencedirect.com/science/article/pii/S1463500314001796). Tests include overflow, internal gravity wave, baroclinic channel, seamount, and vertical merry-go-round.
- **Omega-1.1: Layered ocean, idealized, with surface fluxes.** Addition of simple vertical mixing scheme such as Pacanowski & Philander; nonlinear equation of state (TEOS10); tracer surface restoring to a constant field; constant-in-time wind forcing; and flux-corrected transport for horizontal advection. Testing will be with the baroclinic gyre and single column tests of surface fluxes and vertical mixing.
- **Omega-2.0: Coupled within E3SM, ocean only.** Ability to run C cases (active ocean only) within E3SM. Requires addition of E3SM coupling infrastructure; simple analysis (time-averaged output of mean, min, max); split baroclinic-barotropic time; global bounds checking on state. Testing and analysis similar to [Ringler et al. 2013](https://www.sciencedirect.com/science/article/pii/S1463500313000760), except Omega uses a non-Boussinesq formulation.
- **Omega-2.1: E3SM fully coupled** Ability to run G cases (active ocean and sea ice) and B cases (all components active) within E3SM. This will include: a full vertical mixing scheme, such as KPP; frazil ice formation in the ocean; and a submosescale parameterization. Omega 2.1 will mostly be run at eddy-resolving resolutions, and will not include parameterizations for lower resolutions such as Gent-McWilliams and Redi mixing.  Simulations will be compared to previous E3SM simulations, including [Caldwell et al. 2019](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019MS001870) and [Petersen et al. 2019](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2018MS001373).

This document describes the governing equations for the layered ocean model, which are applicable for Omega-1.0 onwards. Specific terms, such as the pressure gradient, vertical mixing, and parameterizations, are left in their general form here but are described in more detail in other design documents.

## 2. Requirements

The requirements in the [Omega-0 design document](OmegaV0ShallowWater) still apply. Additional design requirements for Omega-1 are:

### Omega will be a hydrostatic, non-Boussinesq ocean model.
Omega will adopt a non-Boussinesq formulation, meaning it retains the full, spatially and temporally varying fluid density $\rho({\bf x}, t)$ in all governing equations. This ensures exact mass conservation, as opposed to volume conservation in Boussinesq models that approximate density as a constant reference value $\rho_0$ outside the pressure gradient. The equations are derived in terms of layer-integrated mass per unit area (see [layered equations](#layered-equations)), and no approximation is made that filters out compressibility or density variations.

The model also assumes hydrostatic balance in the vertical momentum equation, which is a standard and well-justified simplification for large-scale geophysical flows. In such regimes, vertical accelerations are typically small compared to the vertical pressure gradient and gravitational forces. This assumption simplifies the dynamics and removes most sound waves while retaining fidelity for the mesoscale to planetary-scale ocean circulation that Omega is designed to simulate.

### Omega will use TEOS10 for the equation of state.
See additional [EOS design document](EOS)

### Omega-1.0 will add new terms for the pressure gradient, vertical mixing, and vertical advection.
See forthcoming design documents on the pressure gradient, vertical mixing, and vertical advection.

## 3. Conservation Equations

### Control Volume Formulation

We begin with the continuous, control-volume form of the conservation equations. Consider an arbitrary control volume $V(t)$ with a bounding control surface $\partial V(t)$. The fluid has density $\rho({\bf x},t)$ and velocity ${\bf v}({\bf x},t)$. The control surface is moving at a velocity ${\bf v}_r({\bf x},t)$.  The vector ${\bf n}$ denotes the outward-facing unit normal vector on the control surface $\partial V(t)$, and is used in all surface integrals to represent the direction of fluxes across the boundary. The conservation equations are

mass:

$$
\frac{d}{dt} \int_{V(t)} \rho({\bf x},t) \, dV
+ \int_{\partial V(t)}\rho({\bf x},t)\left({\bf v}({\bf x},t) - {\bf v}_r \right) \cdot {\bf n} \, dA
= 0
$$ (continuous-mass)

tracers:

$$
\frac{d}{dt} \int_{V(t)} \rho({\bf x},t) \, \varphi({\bf x},t) \, dV
+ \int_{\partial V(t)}\rho({\bf x},t)\, \varphi({\bf x},t) \left({\bf v}({\bf x},t) - {\bf v}_r \right) \cdot {\bf n} \, dA
= 0
$$ (continuous-tracer)

momentum:

$$
&\frac{d}{dt} \int_{V(t)} \rho({\bf x},t)\,  {\bf v}({\bf x},t) \, dV
+ \int_{\partial V(t)}\rho({\bf x},t)\, {\bf v}({\bf x},t) \left({\bf v}({\bf x},t) - {\bf v}_r \right) \cdot {\bf n} \, dA
\\ & \; \; \; =
\int_{V(t)} \rho({\bf x},t) \, {\bf b}({\bf x},t)\, dV
+ \int_{\partial V(t)} {\bf f}({\bf n},{\bf x},t)  \, dA
$$ (continuous-momentum)

The operator $\frac{d}{dt}$ used here denotes the rate of change within a moving control volume, sometimes referred to as a Reynolds transport derivative. It differs from the partial derivative $\frac{\partial}{\partial t}$, which represents the local rate of change at a fixed point in space (Eulerian frame), and from the material derivative $\frac{D}{Dt}$, which follows an individual fluid parcel (Lagrangian frame). The use of $\frac{d}{dt}$ allows for conservation laws to be expressed in a general framework that includes both stationary and moving control volumes, consistent with the Reynolds transport theorem.

These equations are taken from [Kundu et al. 2016](https://doi.org/10.1016/C2012-0-00611-4) (4.5) for mass conservation and (4.17) for momentum conservation.
All notation is identical to Kundu, except that we use ${\bf v}$ for the three-dimensional velocity (${\bf u}$ will be used below for horizontal velocities) and ${\bf v}_r$ and $\partial V$ match the notation in [Ringler et al. 2013](https://www.sciencedirect.com/science/article/pii/S1463500313000760) Appendix A.2.


The tracer equation is simply mass conservation, where the conserved quantity is the tracer mass $\rho \varphi$, as $\varphi$ is the tracer concentration per unit mass.
In all three equations, the first term is the change of the quantity within the control volume; the second term is the flux through the moving boundary.
If the control surface moves with the fluid in a Lagrangian fashion, then ${\bf v}_r={\bf v}$ and the second term (flux through the boundary) vanishes---there is no net mass, momentum, or tracer transport across the moving surface.


The momentum equation is an expression of Newton's second law and includes two types of external forces.
The first is the body force, represented here as ${\bf b}({\bf x}, t)$, which encompasses any volumetric force acting throughout the fluid, such as gravitational acceleration or the Coriolis force. In some contexts, body forces may be expressible as the gradient of a potential, ${\bf b} = -\nabla_{3D} \Phi$, but this is not assumed in general.
The second is the surface force ${\bf f}$, which acts on the boundary of the control volume and includes pressure and viscous stresses. These forces appear as surface integrals over the boundary and drive momentum exchange between adjacent fluid parcels or between the fluid and its environment.
The derivation of the momentum equation may also be found in [Leishman 2025](https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/conservation-of-momentum-momentum-equation/#chapter-260-section-2), Chapter 21, equation 10.

### Horizontal \& Vertical Separation

In geophysical flows, the vertical and horizontal directions are treated differently due to rotation and stratification, which leads to different characteristic scales of motion.
To make this distinction explicit, we reformulate the conservation equations in a geometry that separates horizontal and vertical fluxes.
We partition the surface integral into contributions from the fixed side walls $\partial V^{\text{side}}$ and the time-varying upper and lower surfaces, $\partial V^{\text{top}}(t)$ and $\partial V^{\text{bot}}(t)$.
The top and bottom surfaces are not necessarily flat, so we retain their full geometry.

As an example, we consider the tracer equation and drop the explicit notation for spatial and temporal dependence for clarity. We write the control-volume form as:

$$
\frac{d}{dt} \int_{V(t)} \rho \, \varphi \, dV
&+
\int_{\partial V^{\text{side}}} \rho \varphi \left({\bf v} - {\bf v}_r \right) \cdot {\bf n} \, dA \\
&+
\int_{\partial V^{\text{top}}(t)} \rho \varphi \left({\bf v} - {\bf v}_r \right) \cdot {\bf n} \, dA
+
\int_{\partial V^{\text{bot}}(t)} \rho \varphi \left({\bf v} - {\bf v}_r \right) \cdot {\bf n} \, dA
= 0
$$ (tr-v-h-split)

The unit normals on the top and bottom surfaces are given by:

$$
{\bf n}^{\text{top}} = \frac{(-\nabla z^{\text{top}}, 1)}{\sqrt{1 + |\nabla z^{\text{top}}|^2}}, \quad
{\bf n}^{\text{bot}} = \frac{(-\nabla z^{\text{bot}}, 1)}{\sqrt{1 + |\nabla z^{\text{bot}}|^2}}
$$ (top-bot-normal)

In typical Omega configurations, the slope of the top and bottom surfaces will be small, i.e., $|\nabla z^{\text{top}}| \ll 1$ and $|\nabla z^{\text{bot}}| \ll 1$. Under this small-slope approximation, we neglect the square root in the denominator of the unit normal, and write:

$$
{\bf n}^{\text{top}} \approx (-\nabla z^{\text{top}}, 1), \quad
{\bf n}^{\text{bot}} \approx (-\nabla z^{\text{bot}}, 1)
$$

This allows sloping-surface contributions such as $\nabla z^{\text{top}}$ to be retained while avoiding more complex metric factors. The approximation is accurate to leading order in slope and consistent with hydrostatic and layered modeling frameworks.

Thus, the flux integrals across sloping boundaries retain contributions from both the vertical and horizontal components of ${\bf v} - {\bf v}_r$, and include the slope terms $\nabla z^{\text{top}}$ and $\nabla z^{\text{bot}}$.

This formulation keeps the geometry fully general and allows for future manipulation or approximations. Subsequent approximations---such as retaining only the vertical component---can be applied explicitly where appropriate in later sections.

To facilitate integration over a fixed horizontal domain, we introduce the following notation:

- $A$ is the horizontal footprint (in the $x$–$y$ plane) of the control volume $V(t)$.
- $dA$ is the horizontal area element.
- $\partial A$ is the boundary of $A$, and $dl$ is the line element along this boundary.
- ${\bf n}_\perp$ is the outward-pointing unit normal vector in the horizontal plane, defined on $\partial A$. It lies in the $x$–$y$ plane and is orthogonal to $dl$.

We now project the tracer equation [](#tr-v-h-split) onto a horizontal domain $A$ (the footprint of the control volume), over which the top and bottom boundaries vary in height. The side walls remain fixed in time and space. Using this projection, we obtain:

$$
\frac{d}{dt} \int_{A} \int_{z^{\text{bot}}(x,y,t)}^{z^{\text{top}}(x,y,t)} \rho \, \varphi \, dz \, dA
&+
\int_{\partial A} \left( \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \varphi \, {\bf u} \, dz \right) \cdot {\bf n}_\perp \, dl \\
&+
\int_A \rho \varphi \left[
  (w - w_r) - {\bf u} \cdot \nabla z^{\text{top}}
\right]_{z = z^{\text{top}}} dA \\
&-
\int_A \rho \varphi \left[
  (w - w_r) - {\bf u} \cdot \nabla z^{\text{bot}}
\right]_{z = z^{\text{bot}}} dA
= 0
$$ (tr-v-h-separation)


Here ${\bf v} = ({\bf u},w)$ separates the three-dimensional velocity into horizontal velocity ${\bf u}$ and vertical velocity $w$.
In the final equation, only velocity components aligned with the boundary normals contribute to each integral, so perpendicular components drop out.
Since the side boundary $\partial V^{\text{side}}$ is fixed in space, ${\bf v}_r = 0$ there and drops out of the corresponding term.
The domain $A$ is the fixed horizontal footprint of the control volume in the $x$–$y$ plane, over which the top and bottom surfaces $z^{\text{top}}(x,y,t)$ and $z^{\text{bot}}(x,y,t)$ may vary in space and time.

Using this procedure of separating the horizontal from the vertical, the governing equations are

mass:

$$
\frac{d}{dt} \int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \, dz \, dA
&+
\int_{\partial A}\left( \int_{z^{\text{bot}}}^{z^{\text{top}}}\rho\, {\bf u} \, dz \right) \cdot {\bf n}_\perp \, dl \\
&+
\int_{A} \rho \left[ (w - w_r) - {\bf u} \cdot \nabla z^{\text{top}} \right]_{z=z^{\text{top}}} \, dA
-
\int_{A} \rho \left[ (w - w_r) - {\bf u} \cdot \nabla z^{\text{bot}} \right]_{z=z^{\text{bot}}} \, dA
= 0
$$ (vh-mass)

tracers:

$$
\frac{d}{dt} \int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \, \varphi \, dz \, dA
&+
\int_{\partial A}\left( \int_{z^{\text{bot}}}^{z^{\text{top}}}\rho\, \varphi \, {\bf u} \, dz \right) \cdot {\bf n}_\perp \, dl \\
&+
\int_{A} \rho \varphi \left[ (w - w_r) - {\bf u} \cdot \nabla z^{\text{top}} \right]_{z=z^{\text{top}}} \, dA
-
\int_{A} \rho \varphi \left[ (w - w_r) - {\bf u} \cdot \nabla z^{\text{bot}} \right]_{z=z^{\text{bot}}} \, dA
= 0
$$ (vh-tracer)

momentum:

$$
\frac{d}{dt} \int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \, {\bf v} \, dz \, dA
&+
\int_{\partial A} \left( \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \, {\bf v} \otimes {\bf u} \, dz \right) \cdot {\bf n}_\perp \, dl \\
&+
\int_A \rho \, {\bf v} \left[ (w - w_r) - {\bf u} \cdot \nabla z^{\text{top}} \right]_{z=z^{\text{top}}} \, dA
-
\int_A \rho \, {\bf v} \left[ (w - w_r) - {\bf u} \cdot \nabla z^{\text{bot}} \right]_{z=z^{\text{bot}}} \, dA \\
&=
\int_A \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \, {\bf b} \, dz \, dA
+
\int_{\partial A} \left( \int_{z^{\text{bot}}}^{z^{\text{top}}} {\bf f} \, dz \right) dl \\
&\quad
+ \int_A \left[ {\bf f} \right]_{z = z^{\text{top}}} \, dA
- \int_A \left[ {\bf f} \right]_{z = z^{\text{bot}}} \, dA
$$ (vh-momentum)

The momentum advection term contains ${\bf v} \otimes {\bf u}$, the outer (or tensor) product of the full velocity ${\bf v} = ({\bf u}, w)$ with the horizontal velocity ${\bf u}$. This object is a 3×2 tensor: the three rows correspond to the components of momentum being advected (in $x$, $y$, and $z$), and the two columns correspond to the directions of horizontal transport. This structure naturally arises in the surface integral over $\partial A$, where the tensor is contracted with the horizontal unit normal vector ${\bf n}_\perp$ to yield a vector flux through the vertical sides of the control volume.

### Hydrostatic Approximation

The momentum equation [](#vh-momentum) consists of three components for the $(x,y,z)$ directions. Thus we may rewrite it as

horizontal momentum:

$$
& \frac{d}{dt} \int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \, {\bf u} \, dz \, dA
+
\int_{\partial A} \left( \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \, {\bf u} \otimes {\bf u} \, dz \right) \cdot {\bf n}_\perp \, dl \\
& +
\int_{A} \rho \, {\bf u} \left[ (w - w_r) - {\bf u} \cdot \nabla z^{\text{top}} \right]_{z = z^{\text{top}}} \, dA
-
\int_{A} \rho \, {\bf u} \left[ (w - w_r) - {\bf u} \cdot \nabla z^{\text{bot}} \right]_{z = z^{\text{bot}}} \, dA \\
& =
\int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \, {\bf b}_\perp \, dz \, dA
+
\int_{\partial A} \left( \int_{z^{\text{bot}}}^{z^{\text{top}}} {\bf f}_\perp \, dz \right) \, dl \\
& +
\int_{A} \left[ {\bf f}_\perp \right]_{z = z^{\text{top}}} \, dA
-
\int_{A} \left[ {\bf f}_\perp \right]_{z = z^{\text{bot}}} \, dA
$$ (h-momentum)

vertical momentum:

$$
& \frac{d}{dt} \int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \, w \, dz \, dA
+
\int_{\partial A} \left( \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \, w \, {\bf u} \, dz \right) \cdot {\bf n}_\perp \, dl \\
& +
\int_{A} \rho \, w \left[ (w - w_r) - {\bf u} \cdot \nabla z^{\text{top}} \right]_{z = z^{\text{top}}} \, dA
-
\int_{A} \rho \, w \left[ (w - w_r) - {\bf u} \cdot \nabla z^{\text{bot}} \right]_{z = z^{\text{bot}}} \, dA \\
& =
\int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \, b_z \, dz \, dA
+
\int_{\partial A} \left( \int_{z^{\text{bot}}}^{z^{\text{top}}} f_z \, dz \right) \, dl \\
& +
\int_{A} \left[ f_z \right]_{z = z^{\text{top}}} \, dA
-
\int_{A} \left[ f_z \right]_{z = z^{\text{bot}}} \, dA
$$ (v-momentum)

where the potential gradient vector is ${\bf b} = ({\bf b}_\perp, b_z)$ and the surface forces are ${\bf f} = ({\bf f}_\perp, f_z)$.

The hydrostatic approximation applies to the vertical component of the momentum equation and assumes that the leading-order balance is between the vertical pressure gradient and the gravitational body force. All other terms---such as vertical acceleration, advection, and viscous or turbulent stresses---are assumed to be negligible in comparison.
Applying this assumption to [](#v-momentum),

$$
\int_A \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \, b_z \, dz \, dA
+ \int_A \left[ f_z \right]_{z = z^{\text{top}}} \, dA
- \int_A \left[ f_z \right]_{z = z^{\text{bot}}} \, dA = 0
$$ (v-hydrostatic1)

Assuming the vertical body force is gravity, $b_z = -g$, and the vertical surface stress is from pressure, $f_z = -p$, the equation becomes:

$$
- \int_A \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho g \, dz \, dA
- \int_A \left[ p \right]_{z = z^{\text{top}}} \, dA
+ \int_A \left[ p \right]_{z = z^{\text{bot}}} \, dA = 0
$$ (v-hydrostatic2)

Although the top and bottom surfaces may be sloping, the vertical component of the pressure force simplifies to $\pm p \, dA$ to leading order. This is because the projection of the pressure force onto the vertical direction introduces a factor of $\hat{\bf n} \cdot \hat{\bf z} \approx 1 - \tfrac{1}{2}|\nabla z|^2$, while the sloping surface area element adds a compensating factor of $\sqrt{1 + |\nabla z|^2} \approx 1 + \tfrac{1}{2}|\nabla z|^2$. These cancel to second order, and the net vertical pressure force is simply the pressure value multiplied by the horizontal area element $dA$. Thus, no explicit slope terms appear in the hydrostatic balance.

Rewriting:

$$
\int_A \left( \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho g \, dz
+ \left[ p \right]_{z = z^{\text{top}}}
- \left[ p \right]_{z = z^{\text{bot}}} \right) dA = 0
$$ (v-hydrostatic)

Because this equation holds for any horizontal region $A$, the integrand must vanish, yielding the hydrostatic pressure relation:

$$
 \left[  p \right]_{z=z^{\text{bot}}}
 =
 \left[  p \right]_{z=z^{\text{top}}}
 + \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \, g \, dz.
$$ (integral-hydrostatic)

Taking the limit as $dz \rightarrow 0$, we arrive at the typical form of the hydrostatic approximation,

$$
 \frac{\partial p}{\partial z}  = - \rho g
$$ (hydrostatic)


## 4. Favre Averaging

The most common approach to determine the structure of the small scale stresses in ocean modeling is through Reynolds' averaging.  In this approach, a generic field $\phi$ is broken into a mean and deviatoric component, i.e.

$$
\phi = \overline{\phi} + \phi^\prime
$$

When deriving the Reynolds' averaged equations, averages of terms with a single prime are discarded by construction.  This is an attractive approach for Boussinesq ocean models since the density is assumed constant in all equations.  When the ocean model is non Boussinesq, this leads to difficulties.  For example, consider the first term in equation [](#h-momentum), if a Reynolds' decomposition and averaging is performed,

$$
\frac{d}{dt} \int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \,  {\bf u}  \, dz \, dA = \frac{d}{dt} \int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}} (\overline{\rho {\bf u}} +  \,  \overline{\rho^\prime {\bf u}^\prime})  \, dz \, dA
$$

In this equation, the products of prime and average drop out by construction.  The $\overline{\rho^\prime {\mathbf u}^\prime}$ term is an unnecessary complication and difficult to parameterize.  To circumvent this complication, Omega will adopt Favre averaging [(Pope 2000)](https://elmoukrie.com/wp-content/uploads/2022/04/pope-s.b.-turbulent-flows-cambridge-university-press-2000.pdf), which for the generic variable $\phi$ is

$$
\phi = \hat{\phi} + \phi^"
$$

Where $\hat{\phi} \equiv \frac{\overline{\rho \phi}}{\overline{\rho}}$, and the double prime indicates deviations from this density weighted mean.  In the definition, the overbar is an averaging operator with identical properties to a Reynolds' average.  Using this relation, we can relate Reynolds' average to Favre average by considering the average of $\rho \phi$.  The standard Reynolds' approach gives

$$
\overline{\rho \phi} = \overline{\rho}\overline{\phi} + \overline{\rho^\prime \phi^\prime}
$$

Isolating $\overline{\phi}$,

$$
\overline{\phi} = \frac{\overline{\rho \phi}}{\overline{\rho}} + \frac{\overline{\rho^\prime \phi^\prime}}{\overline{\rho}}
$$

The first term on the right side of the equation is the definition of a Favre average, which yields

$$
\overline{\phi} = \hat{\phi} + \frac{\overline{\rho^\prime \phi^\prime}}{\overline{\rho}}
$$

Throughout much of the ocean, we expect the second term to be $O(10^{-3})$ smaller than the first, but could be large in highly turbulent regions. With this adoption, all prognostic and diagnostic variables in Omega are interpreted as Favre averages.
This choice ensures that the governing equations are closed in terms of density-weighted means, avoiding the need to model second-order density fluctuations like $\overline{\rho' \phi'}$ that would otherwise arise in a Reynolds framework.

## 5. Layered Equations

### Pseudo-Height

In our non-Boussinesq hydrostatic framework, we adopt a vertical coordinate based on pseudo-height,

$$
\tilde{z}(\hat{p}) = -\frac{1}{\rho_0 g} \, \hat{p}
$$ (def-pseudo-height)

The pseudo-height is simply the pressure normalized by two constants, a reference density $\rho_0$ and gravitational acceleration $g$.
A pseudo-height coordinate is effectively a pressure coordinate, but comes with the intuition and units of distance that people are familiar with.
It is convenient to relate these variables using differentials as

$$
d\tilde{z} = -\frac{1}{\rho_0 g}\, d\hat{p} = \frac{\hat{\rho}}{\rho_0} \, dz
$$ (def-dtildez)

where the second equality uses the hydrostatic balance $d\hat{p} = -\hat{\rho} g dz$.  We note that hydrostatic balance applies equivalently to Favre averaged variables.

This shows that the pseudo-height is nearly the same as the physical height.

In order to convert from $z$ to ${\tilde z}$ the pressure must be computed,

$$
\tilde{z}(z) = -\frac{1}{\rho_0 g} \, p(z) = -\frac{1}{\rho_0 g}\left( p^\text{surf} + \int_{z}^{z^\text{surf}} \rho(z') g dz'\right).
$$ (formula-pseudo-height)

Here, $z'$ is a dummy variable of integration.

The pseudo-velocity in the vertical is

$$
\tilde{w} = \frac{\rho}{\rho_0} \, w.
$$ (def-pseudo-velocity)

This simply falls out of the definitions above, as

$$
\tilde{w} = \frac{d{\tilde z}}{dt} = \frac{\rho}{\rho_0}\frac{dz}{dt} = \frac{\rho}{\rho_0} \, w.
$$ (def-pseudo-velocity)

As above, $\tilde{w}$ has identical units and very similar values to $w$. But in a non-Boussinesq model, it is the vertical *mass* transport that is the physically relevant quantity, not the volume transport. To this end, $\rho w$ is the mass transport per unit area in kg/m$^2$/s. The pseudo-velocity *is* the Eulerian mass transport, but with a convenient normalization of $\rho_0$.

[Griffies 2018](https://doi.org/10.2307/j.ctv301gzg) p. 37  argues for the use of pseudo-velocities, which he calls the density-weighted velocity, for non-Boussinesq models.
Griffies recommends a value of $\rho_0=1035$ kg/m$^3$, following p. 47 of [Gill (1982)](https://doi.org/10.1016/S0074-6142(08)60028-5), because ocean density varies less than 2% from that value.

The use of a constant $\rho_0$ in defining pseudo-height does not imply the Boussinesq approximation. In Boussinesq models, $\rho$ is set to $\rho_0$ everywhere except in the buoyancy term (i.e., the vertical pressure gradient or gravitational forcing). Here, by contrast, we retain the full $\rho$ in all terms, and use $\rho_0$ only as a normalization constant—for example, so that $d\tilde{z} \approx dz$ when $\rho \approx \rho_0$. This preserves full mass conservation while making vertical units more intuitive.

Here we explain the reasoning for the choice of defining $\tilde{z}$ as directly proportional to pressure in [](#def-pseudo-height).
The differential form of the hydrostatic balance $dp = -\rho g dz$ implies that we could choose an arbitrary offset. One could set the offset such that $\tilde{z}=0$ at $z=0$, so that the equilibrium sea surface height matches. Or one could set $\tilde{z}^{\text{floor}} = z^{\text{floor}}$ at some reference depth. Our definition [](#def-pseudo-height) was made so that $\tilde{z}$ varies in space and time in the same way as pressure, and the additional normalization by $\rho_0 g$ was included so that units and values for $\tilde{z}$,  $\tilde{h}$,  and $\tilde{w}$ are intuitive and easy to work with. A major advantage of [](#def-pseudo-height) is that $\tilde{z}^{\text{floor}}$ is proportional to the bottom pressure, and can be used directly for the barotropic pressure gradient in time-split methods.

### Vertical Discretization

The previous equation set [](#vh-mass) to [](#vh-momentum) is for an arbitrary layer bounded by $z^{\text{top}}$ above and $z^{\text{bot}}$ below.
We now provide the details of the vertical discretization.
The ocean is divided vertically into $K_{max}$ layers, with $k=0$ at the top and increasing downwards (opposite from $z$).
Layer $k$ is bounded between $z_k^{\text{top}}$ above and $z_{k+1}^{\text{top}}$ below (i.e. $z_k^{\text{bot}} = z_{k+1}^{\text{top}}$).

The layer thickness of layer k, used in MPAS-Ocean, is

$$
h_k = \int_{z_{k+1}^{\text{top}}}^{z_k^{\text{top}}} dz.
$$ (def-thickness)

In Omega we will use the pseudo-thickness,

$$
{\tilde h}_k(x,y,t)
&= \int_{{\tilde z}_{k+1}^{\text{top}}}^{{\tilde z}_k^{\text{top}}} d\tilde{z} \\
&= \frac{1}{\rho_0} \int_{z_{k+1}^{\text{top}}}^{z_k^{\text{top}}} \rho \, dz \\
&= \frac{1}{\rho_0 g} \left( \hat{p}_{k+1}^{\text{top}} - \hat{p}_k^{\text{top}} \right)
$$ (def-pseudo-thickness)

which is the mass per unit area in the layer, normalized by $\rho_0$. This pseudo-thickness and layer-averaging will be used to express conservation laws in a mass-weighted coordinate system.  Pseudo-thickness, rather than geometric
thickness will be the prognostic variable in Omega.

The density-weighted average of any variable $\phi({\bf x},t)$ in layer $k$ is

$$
{\overline \phi}^{\tilde{z}}_k(x,y,t) =
\frac{\int_{z_{k+1}^{\text{top}}}^{z_k^{\text{top}}} \rho \phi dz}
     {\int_{z_{k+1}^{\text{top}}}^{z_k^{\text{top}}} \rho dz}
     =
\frac{\frac{1}{\rho_0}\int_{z_{k+1}^{\text{top}}}^{z_k^{\text{top}}} \rho \phi dz}
     {\frac{1}{\rho_0}\int_{z_{k+1}^{\text{top}}}^{z_k^{\text{top}}} \rho dz}
     =
\frac{\int_{{\tilde z}_{k+1}^{\text{top}}}^{{\tilde z}_k^{\text{top}}} \phi d{\tilde z}}
     {\int_{{\tilde z}_{k+1}^{\text{top}}}^{{\tilde z}_k^{\text{top}}} d{\tilde z}}
=
\frac{1}{{\tilde h}_k}\int_{{\tilde z}_{k+1}^{\text{top}}}^{{\tilde z}_k^{\text{top}}} \phi d{\tilde z}.
$$ (def-layer-average)

Rearranging, is it useful to note that

$$
\frac{1}{\rho_0}\int_{z_{k+1}^{\text{top}}}^{z_k^{\text{top}}} \rho \phi dz
 =
{\tilde h}_k {\overline \phi}^{\tilde{z}}_k(x,y,t).
$$ (h-phi)

This relation is frequently used in discretized fluxes and conservation equations to replace integrals with layer-mean quantities.

### Layered Tracer & Mass

Substituting [](#def-pseudo-velocity) into the tracer equation [](#vh-tracer) for layer $k$, we have

$$
\frac{d}{dt} \int_{A} \int_{z_{k+1}^{\text{top}}}^{z_k^{\text{top}}} \rho \, \varphi \, dz \, dA
+
\int_{\partial A} \left( \int_{z_{k+1}^{\text{top}}}^{z_k^{\text{top}}} \rho \, \varphi \, {\bf u} \, dz \right) \cdot {\bf n}_\perp \, dl & \\
+
\int_{A} \left[ \rho_0 \varphi (\tilde{w} - \tilde{w}_r) - \rho \varphi {\bf u} \cdot \nabla z^{\text{top}} \right]_{z = z^{\text{top}}} \, dA & \\
-
\int_{A} \left[ \rho_0 \varphi (\tilde{w} - \tilde{w}_r) - \rho \varphi {\bf u} \cdot \nabla z^{\text{bot}} \right]_{z = z^{\text{bot}}} \, dA
& = 0
$$ (Aintegral-tracer)

where we converted to pseudo-height for the last two terms using [](#formula-pseudo-height).

While working in pseudo-height coordinates, vertical transport across a sloping layer interface must still be expressed in terms of the slope of the geometric height, $\nabla z^{\text{top}}$, of that interface. Although it may be tempting to rewrite this slope as $\nabla \tilde{z}^{\text{top}}$, the two are not exactly equivalent due to the nonlinear transformation between pressure and height. Differentiating the pseudo-height definition under hydrostatic balance yields:

$$
\nabla z^{\text{top}} = \frac{\rho_0}{\rho(z^{\text{top}})} \nabla \tilde{z}^{\text{top}} + \frac{1}{\rho(z^{\text{top}}) g} \nabla p^{\text{surf}}
$$

This shows that the geometric slope differs from the pseudo-height slope by a correction term involving the horizontal gradient of surface pressure. Under typical oceanographic conditions—where density is close to the reference $\rho_0$ and surface pressure gradients are modest—this correction may be small. However, we retain $\nabla z^{\text{top}}$ explicitly in the vertical transport terms to avoid introducing assumptions that may not hold in all regimes.

Substituting [](#h-phi) into the first two terms in [](#Aintegral-tracer) and dividing by $\rho_0$,

$$
\frac{d}{dt} \int_A \tilde{h}_k \, \overline{\varphi}^{\tilde{z}}_k \, dA
+
\int_{\partial A} \left( \tilde{h}_k \, \overline{\varphi {\bf u}}^{\tilde{z}}_k \right) \cdot {\bf n}_\perp \, dl & \\
+
\int_A \left[ \varphi (\tilde{w} - \tilde{w}_r) - \frac{\rho}{\rho_0} \varphi {\bf u} \cdot \nabla z^{\text{top}} \right]_{z = z^{\text{top}}} \, dA & \\
-
\int_A \left[ \varphi (\tilde{w} - \tilde{w}_r) - \frac{\rho}{\rho_0} \varphi {\bf u} \cdot \nabla z^{\text{bot}} \right]_{z = z^{\text{bot}}} \, dA
& = 0
$$ (Aintegral-tracer2)

The horizontal tracer flux term is expanded using $\varphi \equiv \overline{\varphi}^{\tilde{z}}_k + \delta \varphi$ and ${\bf u} \equiv \overline{\bf u}^{\tilde{z}}_k + \delta {\bf u}$ to yield

$$
\frac{d}{dt} \int_{A} {\tilde h}_k {\overline \varphi}^{\tilde{z}}_k   \, dA
+
   \int_{\partial A}\left( {\tilde h}_k \, \left(\left(\overline{\varphi}^{\tilde{z}}_k + \delta \varphi\right)\left(\overline{\bf u}^{\tilde{z}}_k + \delta {\bf u}\right)\right) \, \right) \cdot {\bf n} \, dl
 + \int_{A}\left[ \varphi \left({\tilde w} - {\tilde w}_r \right) \right]_{{\tilde z}={\tilde z}_k^{\text{top}}} \, dA
 - \int_{A}\left[ \varphi \left({\tilde w} - {\tilde w}_r \right) \right]_{{\tilde z}={\tilde z}_{k+1}^{\text{bot}}} \, dA
= 0.
$$

$$
\frac{1}{\rho} \frac{d}{dt} \int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}} \overline{\rho} \, \hat{\varphi} \, dz \, dA
+
   \frac{1}{\rho} \int_{\partial A}\left( \int_{z^{\text{bot}}}^{z^{\text{top}}} \overline{\rho}\, (\hat{\varphi}\hat{\bf u} + \overline{\varphi^{\prime \prime} {\bf u}^{\prime \prime}}) \, dz \right) \cdot {\bf n} \, dl
 + \int_{A}\left[ \hat{\varphi} \hat{\tilde w}_{tr} + \overline{\varphi^{\prime \prime} {\tilde w}_{tr}^{\prime \prime}} \right]_{{\tilde z}={\tilde z}_k^{\text{top}}} \, dA
 - \int_{A}\left[ \hat{\varphi} \hat{\tilde w}_{tr} + \overline{\varphi^{\prime \prime} {\tilde w}_{tr}^{\prime \prime}} \right]_{{\tilde z}={\tilde z}_{k+1}^{\text{bot}}} \, dA
= 0
$$ (Aintegral-tracer2)

Applying [](#formula-pseudo-height) to the first two terms and taking the limit as $A \rightarrow 0$ and using Gauss's theorem for the horizontal advection,

$$
\frac{d{\tilde h}_k {\overline \varphi}^{\tilde{z}}_k  }{dt}
+
   \nabla \cdot \left({\tilde h}_k \overline{\varphi {\bf u}}^{\tilde{z}}_k  \right)
 + \left[ \varphi {\tilde w}_{tr} \right]_{{\tilde z}={\tilde z}_k^{\text{top}}}
 - \left[ \varphi {\tilde w}_{tr} \right]_{{\tilde z}={\tilde z}_{k+1}^{\text{bot}}}
= -\left(\nabla \cdot \left({\tilde h}_k \overline{\varphi {\bf u}}^{\tilde{z}}_k  \right) + \left[\overline{\varphi^{\prime \prime} {\tilde w}_{tr}^{\prime \prime}}\right]_{{\tilde z}={\tilde z}_k^{\text{top}}} - \left[\overline{\varphi^{\prime \prime} {\tilde w}_{tr}^{\prime \prime}} \right]_{{\tilde z}={\tilde z}_{k+1}^{\text{bot}}}\right).
$$ (layer-tracer-pre)

Here ${\tilde w}_{tr}={\tilde w} - {\tilde w}_r$ is the normalized vertical mass transport through the layer interface.
It is a pseudo-transport velocity in units of m/s.
Equation [](#layer-tracer-pre) is identical to [Ringler et al. 2013](https://www.sciencedirect.com/science/article/pii/S1463500313000760) (A.24), and further explanation can be found in that appendix.

As a final step, we substitute the approximation $\overline{\varphi {\bf u}}^{\tilde{z}}_k \approx \overline{\varphi }^{\tilde{z}}_k \overline{{\bf u}}^{\tilde{z}}_k$.
The mass equation is identical to the tracer equation with $\varphi=1$. Dropping the overlines for simpler notation, the layered version of mass and tracer conservation are

mass:

$$
\frac{d{\tilde h}_k }{dt}
+
   \nabla \cdot \left({\tilde h}_k{\bf u}_k \right)
 + \left[ {\tilde w}_{tr} \right]_{{\tilde z}={\tilde z}_k^{\text{top}}}
 - \left[ {\tilde w}_{tr} \right]_{{\tilde z}={\tilde z}_{k+1}^{\text{bot}}}
= 0
$$ (layer-mass)

tracer:

$$
\frac{d{\tilde h}_k \varphi_k  }{dt}
+
   \nabla \cdot \left({\tilde h}_k \varphi_k {\bf u}_k \right)
 + \left[ \varphi {\tilde w}_{tr} \right]_{{\tilde z}={\tilde z}_k^{\text{top}}}
 - \left[ \varphi {\tilde w}_{tr} \right]_{{\tilde z}={\tilde z}_{k+1}^{\text{bot}}}
= 0
$$ (layer-tracer)

### Layered Momentum

We now derive the horizontal momentum equation in our non-Boussinesq, hydrostatic framework, following the same finite-volume approach used for mass and tracer conservation. We work with a pseudo-height vertical coordinate $\tilde{z}$ as defined in [Pseudo-Height Coordinate Section](pseudo-height).

We begin by specifying the forces in the full three-dimensional momentum equation [](#continuous-momentum),

$$
\frac{d}{dt} \int_{V(t)} \rho\,  {\bf v} \, dV
+ \int_{\partial V(t)}\rho\, {\bf v} \left({\bf v} - {\bf v}_r \right) \cdot {\bf n} \, dA
= \mathbf{F}_\text{total}[V(t)],
$$ (continuous-momentum2)

with total forces given by:

$$
\mathbf{F}_\text{total}[V(t)] =
- \int_{V(t)} \rho\, \mathbf{f} \times \mathbf{u} \, dV
- \int_{V(t)} \rho\, \nabla_{3D} \Phi \, dV
- \int_{\partial V(t)} p \, \mathbf{n} \, dA
+ \int_{\partial V(t)} \boldsymbol{\tau} \cdot \mathbf{n} \, dA
$$ (momentum-Ftotal)

Each term on the right-hand side corresponds to a physically distinct force acting on the fluid within the control volume:

- The first term is the **Coriolis force**, where $ \mathbf{f} $ is the vector Coriolis parameter (e.g., $ f \hat{\mathbf{z}} $ on the sphere).
- The second term represents the **gravitational force**, expressed in terms of the gradient of the gravitational potential $ \Phi(x, y, z, t) $, which may include effects such as tides and self-attraction and loading.
- The third term is the **pressure force**, which acts on the boundary surfaces and is naturally expressed as a surface integral. It gives rise to both horizontal pressure gradients and contributions from sloping surfaces.
- The fourth term represents the surface stresses. They include **horizontal stress forces** across the vertical sides of the control volume, such as those due to lateral friction or subgrid momentum transfer.
They also include **vertical stress forces** (e.g., wind stress at the ocean surface and bottom drag at the seafloor), projected into horizontal momentum and acting across the top and bottom surfaces.

#### Pressure Term

The pressure force term may be converted from the boundary to the interior with Gauss' divergence theorem (see [Kundu et al. 2016](https://doi.org/10.1016/C2012-0-00611-4) p. 119),

$$
- \int_{\partial V(t)} p \, \mathbf{n} \, dA
= - \int_{V(t)} \nabla_{3D} p \, dV .
$$ (gradp-Gauss)

Considering only the horizontal components, we have

$$
 - \int_{V(t)} \nabla_{\perp} p \, dV
= - \int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}}
  \nabla_{\perp} p \, dz \, dA
= - \int_{A}   \overline{\nabla_{\perp} p}^{z}  \, dA.
$$ (gradp-h)

#### Stress Term

Likewise, the stress tensor integrated over the surface may be converted to a volume integral with Gauss' theorem ([Kundu et al. 2016](https://doi.org/10.1016/C2012-0-00611-4) p. 125 eqn 4.20b),

$$
\int_{\partial V(t)} \boldsymbol{\tau} \cdot \mathbf{n} \, dA
= \int_{V(t)} \nabla_{3D} \cdot \boldsymbol{\tau} \, dV
$$ (stress-Gauss)

for clarity, this can be written in index notation as

$$
\int_{\partial V(t)} n_i \tau_{ij} \, dA
= \int_{V(t)} \frac{\partial}{\partial x_i} \left( \tau_{ij} \right) \, dV
$$ (stress-Gauss-index)

Taking only the horizontal ($j$=1,2),

$$
 \int_{V(t)} \frac{\partial}{\partial x_i} \left( \tau_{ij} \right) \, dV
= \int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}}
 \frac{\partial}{\partial x_i} \left( \tau_{ij} \right)
\, dz \, dA
=  \int_{A}
  \overline{
 \frac{\partial}{\partial x_i} \left( \tau_{ij} \right) }^z dA
$$ (gradp-h)

#### Horizontal momentum

Putting the pressure and stress term into [](h-momentum), and using Gauss' Theorem on the advection, the horizontal momentum equation is

$$
& \frac{d}{dt} \int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \, {\bf u} \, dz \, dA
+
   \int_{A}\nabla_\perp \cdot \left( \int_{z^{\text{bot}}}^{z^{\text{top}}}\rho\, {\bf u} \otimes {\bf u} \, dz \right) dA
 + \int_{A}\left[ \rho\, {\bf u} \left(w - w_r \right) \right]_{z=z^{\text{top}}} \, dA
 - \int_{A}\left[ \rho\, {\bf u} \left(w - w_r \right) \right]_{z=z^{\text{bot}}} \, dA
\\ & \; =
- \int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \,  \mathbf{f} \times \mathbf{u} \, dz \, dA
- \int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}} \rho \, \nabla_\perp \Phi \, dz \, dA
 - \int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}}
  \nabla_{\perp} p \, dz \, dA
+ \int_{A} \int_{z^{\text{bot}}}^{z^{\text{top}}}
 \frac{\partial}{\partial x_i} \left( \tau_{ij} \right)
\, dz \, dA.
$$ (h-momentum-p-tau)

Dividing by $\rho_0$ and taking vertical averages using [](#def-thickness) and [](#def-pseudo-thickness),

$$
\frac{d}{dt} \int_{A} \tilde{h}\, \overline{  {\bf u} }^{\tilde z} \, dA
+
   \int_{A}\nabla_\perp \cdot \left( \tilde{h}\, \overline{{\bf u} \otimes {\bf u} }^{\tilde z} \right) dA
& + \int_{A}\left[ {\bf u}\, \tilde{w}_{tr} \right]_{z=z^{\text{top}}} \, dA
 - \int_{A}\left[ {\bf u}\, \tilde{w}_{tr} \right]_{z=z^{\text{bot}}} \, dA
\\ & \; =
- \int_{A} \tilde{h} \,\overline{  \mathbf{f} \times \mathbf{u} }^{\tilde z} \, dA
- \int_{A} \tilde{h} \,\overline{  \nabla_\perp \Phi }^{\tilde z} \, dA
- \int_{A} \frac{1}{\rho_0} \overline{\nabla_{\perp} p}^{z}  \, dA
+ \int_{A} \frac{1}{\rho_0} \overline{ \frac{\partial}{\partial x_i} \left( \tau_{ij} \right) }^z dA
$$ (h-momentum-p-tau)

Taking the limit as $A \rightarrow 0$, we arrive at the local, horizontally continuous form:

$$
\frac{\partial \tilde{h}\, \overline{  {\bf u} }^{\tilde z} }{\partial t}
& +
   \nabla_\perp \cdot \left( \tilde{h}\, \overline{{\bf u} \otimes {\bf u} }^{\tilde z} \right)
 + \left[ {\bf u}\, \tilde{w}_{tr} \right]_{z=z^{\text{top}}}
 - \left[ {\bf u}\, \tilde{w}_{tr} \right]_{z=z^{\text{bot}}}
\\ & =
-  \tilde{h} \,\overline{  \mathbf{f} \times \mathbf{u} }^{\tilde z}
-  \tilde{h} \,\overline{  \nabla_\perp \Phi }^{\tilde z}
-  \frac{1}{\rho_0} \overline{\nabla_{\perp} p}^{z}
+  \frac{1}{\rho_0} \overline{ \frac{\partial}{\partial x_i} \left( \tau_{ij} \right) }^z
$$ (momentum-layered-differential-1)

The first two terms are the material derivative, confined within the horizontal layer. Using the product rule, and layered mass conservation [](#layer-mass),

$$
\frac{\partial \tilde{h}\, \overline{  {\bf u} }^{\tilde z} }{\partial t}
+ \nabla_\perp \cdot \left( \tilde{h}\, \overline{{\bf u} \otimes {\bf u} }^{\tilde z} \right)
&=
\frac{D_\perp \tilde{h}\, \overline{  {\bf u} }^{\tilde z} }{D t} \\
&=
\tilde{h}\frac{D_\perp  \overline{  {\bf u} }^{\tilde z} }{D t}
+
\overline{ {\bf u} }^{\tilde z} \frac{D_\perp \tilde{h}}{D t} \\
&=
\tilde{h}\left(
  \frac{\partial  \overline{  {\bf u} }^{\tilde z} }{\partial t}
+ \nabla_\perp \cdot \left( \overline{{\bf u} \otimes {\bf u} }^{\tilde z} \right)  \right)
- \overline{ {\bf u} }^{\tilde z} \left(
   \left[ {\tilde w}_{tr} \right]_{{\tilde z}={\tilde z}_k^{\text{top}}}
 - \left[ {\tilde w}_{tr} \right]_{{\tilde z}={\tilde z}_{k+1}^{\text{bot}}}
\right)
$$ (2D-material-der-product)

We now substitute [](#2D-material-der-product) into [](#momentum-layered-differential-1) and divide by $\tilde h$ to get

$$
 \frac{\partial \overline{  {\bf u} }^{\tilde z} }{\partial t}
+ \nabla_\perp \cdot \left(  \overline{{\bf u} \otimes {\bf u} }^{\tilde z} \right)
& + \frac{\left[ {\bf u}\, \tilde{w}_{tr} \right]_{z=z^{\text{top}}}
 - \left[ {\bf u}\, \tilde{w}_{tr} \right]_{z=z^{\text{bot}}} }{\tilde h}
-   \frac{\overline{ {\bf u} }^{\tilde z}}{\tilde{h}}  \left(
   \left[ {\tilde w}_{tr} \right]_{{\tilde z}={\tilde z}_k^{\text{top}}}
 - \left[ {\tilde w}_{tr} \right]_{{\tilde z}={\tilde z}_{k+1}^{\text{bot}}}\right)
\\ &  =
-  \overline{  \mathbf{f} \times \mathbf{u} }^{\tilde z}
-  \overline{  \nabla_\perp \Phi }^{\tilde z}
-  \frac{1}{\rho_0 \, \tilde{h}} \overline{\nabla_{\perp} p}^{z}
+  \frac{1}{\rho_0 \, \tilde{h}} \overline{ \frac{\partial}{\partial x_i} \left( \tau_{ij} \right) }^z
$$ (momentum-layered-differential-2)

Note that the coefficient of the pressure gradient and stress tensor can be rewritten as

$$
\frac{1}{\rho_0 \, \tilde{h}} =
\frac{1}{\rho_0 \,\frac{1}{\rho_0} \int_{z_{k+1}^{\text{top}}}^{z_k^{\text{top}}} \rho dz}
= \frac{1} {\overline{\rho}^z}
= \overline{\alpha}^z,
$$ (grad-p-coeff)

where $\alpha$ is the specific volume.
Gathering the vertical advection terms, [](#momentum-layered-differential-2) becomes

$$
 \frac{\partial \overline{  {\bf u} }^{\tilde z} }{\partial t}
 + \nabla_\perp \cdot \left(  \overline{{\bf u} \otimes {\bf u} }^{\tilde z} \right)
&
+ \left( \nabla_\perp \cdot  {\bf u}  \right) {\bf u}
+ \frac{1}{\tilde h}\left(
 \left[ {\bf u}\right]_{z=z^{\text{top}}}-\overline{ {\bf u} }^{\tilde z}\right) \, \left[\tilde{w}_{tr} \right]_{z=z^{\text{top}}}
 -
  \frac{1}{\tilde h}\left(
 \left[ {\bf u}\right]_{z=z^{\text{bot}}}-\overline{ {\bf u} }^{\tilde z}\right) \, \left[\tilde{w}_{tr} \right]_{z=z^{\text{bot}}}
\\ &  =
-  \overline{  \mathbf{f} \times \mathbf{u} }^{\tilde z}
-  \overline{  \nabla_\perp \Phi }^{\tilde z}
-  \overline{\alpha}^z \overline{\nabla_{\perp} p}^{z}
+  \overline{\alpha}^z \overline{ \frac{\partial}{\partial x_i} \left( \tau_{ij} \right) }^z
$$ (momentum-layered-differential-3)

In addition, the horizontal advection may be rewritten as

$$
\nabla_\perp \cdot \left( {\bf u} \otimes {\bf u}  \right)
&=
\nabla_\perp \cdot \left( {\bf u} {\bf u}^T  \right) \\
&= \left( \nabla_\perp \cdot  {\bf u}  \right) {\bf u}
+  {\bf u} \cdot \nabla_\perp {\bf u}
$$ (adv2d-prod)

The term ${\bf u} \cdot \nabla_\perp {\bf u}$ may be replaced with the vector identity

$$
\begin{aligned}
{\bf u} \cdot \nabla_\perp{\bf u}
&= (\nabla_\perp\times {\bf u}) \times {\bf u} + \nabla_\perp\frac{|{\bf u}|^2}{2} \\
&= \left( \boldsymbol{k} \cdot (\nabla_\perp\times {\bf u})\right)
\left( \boldsymbol{k} \times {\bf u} \right) + \nabla_\perp\frac{|{\bf u}|^2}{2} \\
&= \zeta {\bf u}^{\perp} + \nabla_\perp K,
\end{aligned}
$$ (advection-identity)

where $\zeta$ is relative vorticity and $K$ is kinetic energy.
This step separates the horizontal advection into non-divergent and non-rotational components, which is useful in the final TRiSK formulation.

Now [](#momentum-layered-differential-3) becomes

$$
 \frac{\partial \overline{  {\bf u} }^{\tilde z} }{\partial t}
 +
   \left(  \zeta + f  \right) {\overline{  {\bf u} }^{\tilde z}}^\perp
&
+ \left( \nabla_\perp \cdot  {\bf u}  \right) {\bf u}
+ \frac{1}{\tilde h}\left(
 \left[ {\bf u}\right]_{z=z^{\text{top}}}-\overline{ {\bf u} }^{\tilde z}\right) \, \left[\tilde{w}_{tr} \right]_{z=z^{\text{top}}}
 -
  \frac{1}{\tilde h}\left(
 \left[ {\bf u}\right]_{z=z^{\text{bot}}}-\overline{ {\bf u} }^{\tilde z}\right) \, \left[\tilde{w}_{tr} \right]_{z=z^{\text{bot}}}
\\ &  =
- \nabla_\perp K
-  \overline{  \nabla_\perp \Phi }^{\tilde z}
-  \overline{\alpha}^z \overline{\nabla_{\perp} p}^{z}
+  \overline{\alpha}^z \overline{ \frac{\partial}{\partial x_i} \left( \tau_{ij} \right) }^z
$$ (momentum-layered-differential-4)

Notes from Mark:

1. I don't know what to do with the z-average in the advection term of [](#momentum-layered-differential-2) and here. In the mean time, I'll just skip to the layer-averaged terms.
2. I don't know what to make of the $\left( \nabla_\perp \cdot  {\bf u}  \right) {\bf u}$ term, and whether it cancels with some of those odd vertical transport terms.



### Xylar: Horizontal Velocity Equation

We now derive a prognostic equation for the layer-averaged horizontal velocity, starting from the horizontal momentum and mass conservation equations. The goal is to express the equation in terms of the velocity field $\overline{\mathbf{u}}^{\tilde{z}}$ and to isolate the nonlinear advection term in a standard form that can be further decomposed using vector identities.

Let us define the layer-averaged velocity as

$$
\mathbf{U}(x, y, t) \equiv \overline{\mathbf{u}}^{\tilde{z}}.
$$ (layer-avg-velocity)

We begin with the horizontally continuous, layer-integrated horizontal momentum equation [](momentum-layered):

$$
\begin{aligned}
\frac{\partial}{\partial t} \left(  \tilde{h}\, \mathbf{U} \right)
&+ \nabla_z \cdot \left(  \tilde{h}\, \overline{\mathbf{u} \otimes \mathbf{u}}^{\tilde{z}} \right) \\
&+ \left[ \mathbf{u} \tilde{w}_{tr} \right]_{\tilde{z}^{\text{top}}}
- \left[ \mathbf{u} \tilde{w}_{tr} \right]_{\tilde{z}^{\text{bot}}} \\
&=
- \int_{V(t)} \rho\, \mathbf{f} \times \mathbf{u} \, dV
- \tilde{h}\, \overline{ \mathbf{f} \times \mathbf{u} + \nabla_z \Phi }^{\tilde{z}} \\
&\quad - \frac{1}{\rho_0}\nabla_z \left( \int_{{z}^{\text{bot}}}^{{z}^{\text{top}}} p \, d{z} \right)
- \frac{1}{\rho_0}\left[ p \nabla_z \tilde{z} \right]_{\tilde{z}^{\text{top}}}
+ \frac{1}{\rho_0}\left[ p \nabla_z \tilde{z} \right]_{\tilde{z}^{\text{bot}}} \\
&\quad + \frac{1}{\rho_0}\nabla_z \cdot \left( \int_{{z}^{\text{bot}}}^{{z}^{\text{top}}} \boldsymbol{\tau}_h \, d{z} \right)
+ \frac{1}{\rho_0}\left[ \boldsymbol{\tau}_h^z \right]_{\tilde{z}^{\text{top}}}
- \frac{1}{\rho_0}\left[ \boldsymbol{\tau}_h^z \right]_{\tilde{z}^{\text{bot}}}
\end{aligned}
$$ (momentum-velocity-start)

We also have the mass conservation equation [](layer-mass):

$$
\frac{\partial  \tilde{h}}{\partial t}
+ \nabla_z \cdot ({\tilde h} \mathbf{U})
+ \left[ \tilde{w}_{tr} \right]_{\tilde{z}^{\text{top}}}
- \left[ \tilde{w}_{tr} \right]_{\tilde{z}^{\text{bot}}}
= 0.
$$ (mass-conservation-repeated)

Multiplying the mass conservation equation by $\mathbf{U}$ and subtracting it from the momentum equation eliminates time derivatives of $h$ and isolates the evolution of velocity. After simplification, we obtain:

$$
\begin{aligned}
\tilde{h} \left(
\frac{\partial \mathbf{U}}{\partial t}
+ \mathbf{U} \cdot \nabla_z \mathbf{U}
\right)
&+ \left[ \mathbf{u} \tilde{w}_{tr} \right]_{\tilde{z}^{\text{top}}}
- \mathbf{U} \left[ \tilde{w}_{tr} \right]_{\tilde{z}^{\text{top}}} \\
&- \left[ \mathbf{u} \tilde{w}_{tr} \right]_{\tilde{z}^{\text{bot}}}
+ \mathbf{U} \left[ \tilde{w}_{tr} \right]_{\tilde{z}^{\text{bot}}}
= \text{(forces)}
\end{aligned}
$$ (velocity-eq-general)

The term $ \mathbf{U} \cdot \nabla_z \mathbf{U} $ may be rewritten using a standard vector identity:

$$
\mathbf{U} \cdot \nabla_z \mathbf{U}
= \nabla_z \left( \frac{1}{2} |\mathbf{U}|^2 \right)
+ (\nabla_z \times \mathbf{U}) \times \mathbf{U} = \nabla_z K + \zeta \times \mathbf{U},
$$ (velocity-vector-identity)
where $K$ is the kinetic energy per unit mass and $\zeta$ is the relative voriticity.

This decomposition splits the nonlinear advection into a **potential** component (gradient of kinetic energy) and a **rotational** component (advection by relative vorticity).

We now combine the relative vorticity term $\zeta \times \mathbf{U}$ with the vertically averaged Coriolis force term $\overline{\mathbf{f} \times \mathbf{u}}^{\tilde{z}}$. If we assume that the Coriolis vector $\mathbf{f}$ is uniform within the layer or varies slowly enough to justify replacing $\mathbf{u}$ with its vertical average $\mathbf{U}$, then the two terms may be added directly. This gives rise to a term of the form $\boldsymbol{\eta} \times \mathbf{U}$, where the **absolute vorticity** is defined as

$$
\boldsymbol{\eta} = \zeta + \mathbf{f}.
$$ (absolute-vorticity)

This substitution makes the rotational part of the velocity equation appear in its familiar form from geophysical fluid dynamics, emphasizing the role of absolute vorticity in governing the curvature and rotation of horizontal flow. Substituting into [](velocity-eq-general) and returning to the $\overline{\mathbf{u}}^{\tilde{z}}$ notation yields:

$$
\begin{aligned}
\frac{\partial \overline{\mathbf{u}}^{\tilde{z}}}{\partial t}
+ \nabla_z K
+ \boldsymbol{\eta} \times \overline{\mathbf{u}}^{\tilde{z}}
&= \overline{ \nabla_z \Phi }^{\tilde{z}} \\
&\quad - \frac{1}{\rho_0{\tilde h}} \nabla_z \left( \int_{{z}^{\text{bot}}}^{{z}^{\text{top}}} p\, d{z} \right)
- \frac{1}{\rho_0{\tilde h}} \left[ p \nabla_z \tilde{z} \right]_{\tilde{z}^{\text{top}}}
+ \frac{1}{\rho_0{\tilde h}} \left[ p \nabla_z \tilde{z} \right]_{\tilde{z}^{\text{bot}}} \\
&\quad + \frac{1}{\rho_0{\tilde h}} \nabla_z \cdot \left( \int_{{z}^{\text{bot}}}^{{z}^{\text{top}}} \boldsymbol{\tau}_h \, d{z} \right)
+ \frac{1}{\rho_0{\tilde h}} \left[ \boldsymbol{\tau}_h^z \right]_{\tilde{z}^{\text{top}}}
- \frac{1}{\rho_0{\tilde h}} \left[ \boldsymbol{\tau}_h^z \right]_{\tilde{z}^{\text{bot}}} \\
&\quad + \frac{1}{\tilde h} \sum_{i = \text{top, bot}} \left[ (\mathbf{u} - \overline{\mathbf{u}}^{\tilde{z}}) \tilde{w}_{tr} \right]_i
\end{aligned}
$$ (velocity-layered)

This is a prognostic equation for the layer-averaged horizontal velocity, written in a form familiar from geophysical fluid dynamics. The right-hand side includes the horizontal pressure gradient, gravitational body force, divergence of horizontal and vertical stresses, and corrections due to unresolved vertical momentum exchange.

### General Vertical Coordinate

The vertical layer interfaces $(z_k^{bot}, z_k^{top})$ (or equivalently $(p_k^{bot}, p_k^{top})$) can vary as a function of $(x,y,t)$. Thus, these equations describe a general vertical coordinate, and these interface surfaces may be chosen arbitrarily by the user. See [Adcroft and Hallberg 2006](https://www.sciencedirect.com/science/article/pii/S1463500305000090) Section 2 and [Griffies et al (2000)](http://sciencedirect.com/science/article/pii/S1463500300000147) Section 2.  It is convenient to introduce a new variable, the coordinate $r(x,y,z,t)$, where $r$ is constant at the middle of each layer. To be specific, we could design $r$ to be the layer index $k$ at the mid-depth of the layer. That is, define the mid-depth as

$$
z_k^{mid}(x,y,t) = \frac{z_k^{bot}(x,y,t) +  z_k^{top}(x,y,t)}{2}
$$ (def-z)
and generate the function $r$ such that

$$
r(x,y,z_k^{mid}(x,y,t),t) = k
$$ (def-r)
and interpolate linearly in the vertical between mid-layers. When we write the layered form of the equations, we must take into account of the tilted layers using the chain rule.

We now rewrite derivatives in order to convert from horizontal coordinates to tilted coordinates. Let $(x,y)$ be the original horizontal coordinates, which are perpendicular to $z$, and the horizontal gradient be written as $\nabla_z=(\partial/\partial x, \partial/\partial y)$, as above. Now define a tilted coordinate system using the layers defined by $r$, where the within-layer horizontal coordinates are $(x',y')$ and the along-layer gradient is written as $\nabla_r=(\partial/\partial x', \partial/\partial y')$. We construct $r$ to be monotonic in $z$, so we can invert it as $z(x',y',r,t)$. Now horizontal derivatives along the tilted direction $x'$ for any field $\varphi(x,y,z,t)$ can be expanded using the chain rule as

$$
\frac{\partial }{\partial x'} \left[ \varphi(x(x'),y(y'),z(x',y',r,t),t) \right]
= \frac{\partial \varphi}{\partial x}\frac{\partial x}{\partial x'} +  \frac{\partial \varphi}{\partial z} \frac{\partial z}{\partial x'}
$$ (dvarphidx)

We may define the tilted horizontal variable $x'$ as we please. The simplest definition is $x'(x)\equiv x$. Then $\partial x / \partial x'=1$. Rearranging [](#dvarphidx) and repeating for $y$,

$$
\begin{aligned}
\frac{\partial \varphi}{\partial x} &=
\frac{\partial \varphi}{\partial x'}
- \frac{\partial \varphi}{\partial z} \frac{\partial z}{\partial x'}\\
\frac{\partial \varphi}{\partial y} &=
\frac{\partial \varphi}{\partial y'}
- \frac{\partial \varphi}{\partial z} \frac{\partial z}{\partial y'}.
\end{aligned}
$$ (dvarphidxy)

This may be written in vector form as

$$
\nabla_z \varphi = \nabla_r \varphi - \frac{\partial \varphi}{\partial z} \nabla_r z.
$$ (dvarphidnabla)

### Pressure Gradient

For most terms, we can safely assume $\nabla_z \approx \nabla_r$ because the vertical to horizontal aspect ratio even at very high horizontal resolution is on the order of $\epsilon ~ 10^{-3}$ (thought may need to re-assess this assumption if we decide to use strongly sloped layers). This applies, for example, to the curl operator uset to compute the relative vorticity and the gradient applied to the kinetic energy in [](z-integration-momentum).

The exception is the pressure gradient term.  This is becasue strong vertical and week horizontal pressure gradients mean that both terms in the chain rule [](dvarphidnabla) are of the same order and must be retained.  Substituting pressure for $\varphi$ in [](#dvarphidnabla),

$$
\begin{aligned}
-\frac{1}{\rho} \nabla_z p
&=-\frac{1}{\rho}  \nabla_r p +\frac{1}{\rho}  \frac{\partial p}{\partial z} \nabla_r z \\
&=-\frac{1}{\rho} \nabla_r p - g \nabla_r z \\
&=-v \nabla_r p - \nabla_r \Phi_g\\
\end{aligned}
$$ (gradp)

where we have substituted hydrostatic balance [](hydrostatic-balance), specific volume $\alpha\equiv 1/\rho$, and $\Phi_g=gz$.

The general form of the geopotential may include the Earth's gravity, tidal forces, and self attraction and loading (SAL), and may be written as

$$
\Phi = \Phi_g + \Phi_{tides} + \Phi_{SAL} + c
$$ (def-geopotential)

where $c$ is an arbitrary constant. Therefore, the pressure gradient and geopotential gradient may be written together as

$$
\begin{aligned}
-v \nabla_z p - \nabla_z \Phi
&= -v \nabla_z p - \nabla_z \Phi_{tides} - \nabla_z \Phi_{SAL}\\
&=-v  \nabla_r p - \nabla_r \Phi_g - \nabla_r \Phi_{tides} - \nabla_r \Phi_{SAL}\\
&=-v  \nabla_r p - \nabla_r \Phi. \\
\end{aligned}
$$ (gradp-gradphi)

On the first line, note that $\nabla_z \Phi_g=\nabla_z gz=0$. For tides and SAL we assume that these forces do not vary in the vertical due to the small aspect ratio of the ocean, so that the vertical derivative in the expansion [](dvarphidnabla) is zero. This means that $\nabla_z \Phi_{tides}=\nabla_r \Phi_{tides}$ and $\nabla_z \Phi_{SAL}=\nabla_r \Phi_{SAL}$.
For versions 1.0 and 2.0 of Omega we only consider a constant gravitational force, and will not include tides and SAL.  Further details will be provided in the forthcoming pressure gradient design document.

See [Adcroft and Hallberg 2006](https://www.sciencedirect.com/science/article/pii/S1463500305000090) eqn. 1 and [Griffies et al](http://sciencedirect.com/science/article/pii/S1463500300000147) eqn 2 for additional examples of the pressure gradient in tilted coordinates. The additional terms due to the expansion of $\nabla_z$ to $\nabla_r$ in the rest of the equations are small and are ignored.

Some publications state that the transition from Boussinesq to non-Boussinesq equations is accompanied by a change from z-coordinate to pressure-coordinates. However, we use a general vertical coordinate, so the vertical may be referenced to $z$ or $p$. In a purely z-coordinate model like POP, only the $\nabla p$ term is used in [](gradp). In a purely p-coordinate model, only $\nabla z$ remains, as described in  [de Szoeke and Samelson 2002](https://journals.ametsoc.org/view/journals/phoc/32/7/1520-0485_2002_032_2194_tdbtba_2.0.co_2.xml). In a general vertical coordinate model the layer interface placement is up to the user's specification, and so both terms are kept.

## 5. Discrete Equations

The horizontally discretized layered equations are as follows. We have dropped the $r$ in $\nabla_r$ for conciseness, and the operator $\nabla$ from here on means within-layer.

$$
\frac{\partial u_{e,k}}{\partial t} + \left[ \frac{{\bf k} \cdot \nabla \times u_{e,k} +f_v}{[h_{i,k}]_v}\right]_e\left([h_{i,k}]_e u_{e,k}^{\perp}\right)
+ \frac{ \left[ \omega_{i,k}^{top}   \right]_e u_{e,k}^{top}
-        \left[ \omega_{i,k+1}^{top} \right]_e u_{e,k+1}^{top}}{ \left[h_{i,k}\right]_e }
=
- \left[ \alpha_{i,k} \right]_e \nabla p_{i,k} - \nabla \Phi_{i,k}
- \nabla K_{i,k} +  { \bf D}^u_{e,k} + {\bf F}^u_{e,k}
$$ (discrete-momentum)

$$
\frac{\partial h_{i,k}}{\partial t} + \nabla \cdot \left([h_{i,k}]_e u_{e,k}\right)
+ \omega_{i,k}^{top} - \omega_{i,k+1}^{top}
= Q^h_{i,k}
$$ (discrete-thickness)

$$
\frac{\partial h_{i,k} \varphi_{i,k}}{\partial t} + \nabla \cdot \left(u_{e,k} [h_{i,k} \varphi_{i,k}]_e \right)
+ \varphi_{i,k}^{top} \omega_{i,k}^{top} - \varphi_{i,k+1}^{top}\omega_{i,k+1}^{top}
= D^\varphi_{i,k} + Q^\varphi_{i,k}
$$ (discrete-tracer)

$$
p_{i,k} = p_{i}^{surf} + \sum_{k'=1}^{k-1} g h_{i,k'} + \frac{1}{2} g h_{i,k}
$$ (discrete-pressure)

$$
\alpha_{i,k} = f_{eos}(p_{i,k},\Theta_{i,k},S_{i,k})
$$ (discrete-eos)

$$
z_{i,k}^{top} = z_{i}^{floor} + \sum_{k'=k}^{K_{max}} \alpha_{i,k'}h_{i,k'}
$$ (discrete-z)

The subscripts $i$, $e$, and $v$ indicate cell, edge, and vertex locations and subscript $k$ is the layer.  Square brackets $[\cdot]_e$ and $[\cdot]_v$ are quantities that are interpolated to edge and vertex locations. For vector quantities, $u_{e,k}$ denotes the normal component at the center of the edge, while $u_{e,k}^\perp$ denotes the tangential component. We have switched from $\varphi_{i,k}^{bot}$ to the identical $\varphi_{i,k+1}^{top}$ for all variables in order for the notation to match the array names in the code.  The superscripts $surf$ and $floor$ are the surface and floor of the full ocean column. All variables without these superscripts indicate that they are layer-averaged, as defined in [](def-mass-thickness-average), and can be considered to represent a mid-layer value in the vertical. The mid-layer location is equivalently the average in $z$, $p$, or $h$ (mass), since density $\rho_{i,k}$ is considered constant in the cell.

We refer to these as the discrete equations, but time derivatives remain continuous. The time discretization is described in the [time stepping design document](TimeStepping.md). The velocity, mass-thickness, and tracers are solved prognostically using [](discrete-momentum), [](discrete-thickness), [](discrete-tracer). At the new time, these variables are used to compute pressure [](discrete-pressure), specific volume [](discrete-eos), and z-locations [](discrete-z). Additional variables are computed diagnostically at the new time: $u^{\perp}$, $K$, $\omega$, $z^{mid}$, $\Phi$, etc. The initial geopotential is simply $\Phi=gz$, but additional gravitational terms may be added later.

The horizontal operators $\nabla$, $\nabla\cdot$, and $\nabla \times$ are now in their discrete form. In the TRiSK design, gradients ($\nabla$) map cell centers to edges; divergence ($\nabla \cdot$) maps edge quantities to cells; and curl ($\nabla \times$) maps edges to vertices. The exact form of operators and interpolation stencils remain the same as those given in [Omega-0 design document](OmegaV0ShallowWater.md#operator-formulation). The discrete version of terms common with Omega-0, such as advection, potential vorticity, and $\nabla K$, can be found in [Omega-0 Momentum Terms](OmegaV0ShallowWater.md#momentum-terms) and [Omega-0 Thickness and Tracer Terms](OmegaV0ShallowWater.md#thickness-and-tracer-terms).


### Momentum Dissipation

The discretized momentum dissipation ${ \bf D}^u_{e,k}$ may include these terms, which are detailed in the subsections below.

$$
{ \bf D}^u_{e,k} =  \nu_2 \nabla^2 u_{e,k} - \nu_4 \nabla^4 u_{e,k} +
\frac{\partial }{\partial z} \left( \nu_v \frac{\partial u_{e,k}}{\partial z} \right)
$$ (discrete-mom-del2)

#### Laplacian dissipation (del2)

$$
 \nu_2 \nabla^2 u_{e,k} = \nu_2 \left( \nabla D_{i,k} - \nabla^{\perp} \zeta_{v,k} \right)
$$ (discrete-mom-del2)

where $D$ is divergence and $\zeta$ is relative vorticity. See [Omega V0 Section 3.3.4](OmegaV0ShallowWater.md#del2-momentum-dissipation)

#### Biharmonic dissipation (del4)
As in [Omega V0 Section 3.3.5](OmegaV0ShallowWater.md#del4-momentum-dissipation), biharmonic momentum dissipation is computed with two applications of the Del2 operator above.

$$
 - \nu_4 \nabla^4 u_{e,k}
= - \nu_4 \nabla^2 \left( \nabla^2 u_{e,k} \right)
$$ (discrete-mom-del4)

#### Vertical momentum diffusion
Vertical derivatives may be computed with either $z$ or $p$ as the independent variable,

$$
\frac{\partial }{\partial z} \left( \nu_v \frac{\partial u}{\partial z} \right)
= \frac{\partial }{\partial p}\frac{\partial p}{\partial z} \left( \nu_v \frac{\partial u}{\partial p} \frac{\partial p}{\partial z}\right)
= \rho g^2\frac{\partial }{\partial p} \left( \nu_v \rho \frac{\partial u}{\partial p} \right).
$$ (mom-vert-diff-z-p)

We choose to use $z$ values for simplicity. A single vertical derivative of an arbitrary variable $\varphi$ at mid-layer is

$$
\frac{\partial \varphi_k}{\partial z}
= \frac{\varphi_k^{top} - \varphi_k^{bot} }{z_k^{top} - z_k^{bot}}
$$ (vertderiv1)

and a second derivative is

$$
\frac{\partial }{\partial z} \left(
\frac{\partial \varphi_k}{\partial z} \right)
=
\frac{1}{z_{k}^{top} - z_{k+1}^{top}} \left(
\frac{\varphi_{k-1} - \varphi_k }{z_{k-1}^{mid} - z_k^{mid}}
 -
\frac{\varphi_{k} - \varphi_{k+1} }{z_{k}^{mid} - z_{k+1}^{mid}}
\right)
$$ (vertderiv2)

Thus, the vertical momentum diffusion is

$$
\frac{\partial }{\partial z} \left( \nu_v \frac{\partial u_{e,k}}{\partial z} \right)
=
\frac{1}{z_{e,k}^{top} - z_{e,k+1}^{top}} \left(
\nu_{e,k}^{top}
\frac{u_{e,k-1} - u_k }{z_{e,k-1}^{mid} - z_k^{mid}}
 -
\nu_{e,k+1}^{top}
\frac{u_{e,k} - u_{e,k+1} }{z_{e,k}^{mid} - z_{e,k+1}^{mid}}
\right)
$$ (discrete-mom-vert-diff)

This stencil is applied as an implicit tri-diagonal solve at the end of the time step. See details in the [tridiagonal solver design document](TridiagonalSolver) and forthcoming vertical mixing design document.

### Momentum Forcing
The discretized momentum forcing ${ \bf F}^u_{e,k}$ may include:

#### Wind Forcing

The wind forcing is applied as a top boundary condition during implicit vertical mixing as

$$
\frac{\tau_{e}}{[ h_{i,k}]_e}
$$

where $\tau$ is the wind stress in Pa. Since the mass-thickness $h$ is in kg/s/m$^2$, this results in the desired units of m/s$^2$ for a momentum tendency term.

#### Bottom Drag

Bottom Drag is applied as a bottom boundary condition during implicit vertical mixing as

$$
- C_D \frac{u_{e,k}\left|u_{e,k}\right|}{[\alpha_{i,k}h_{i,k}]_e} .
$$ (discrete-mom-bottom)

The units of specific volume times mass-thickness $\alpha h$ are length (m), so that the full term has units of m/s$^2$.

#### Rayleigh Drag

Rayleigh drag is a simple drag applied to every cell.  It is used to ensure stability during spin-up.

$$
- Ra \, u_{e,k}
$$ (discrete-mom-Ra)

### Tracer Diffusion

The discretized tracer diffusion $ D^\varphi_{i,k}$ may include these terms, which are detailed below. Here $\kappa_2$ and $\kappa_4$ are written in front of the operator for simplicity.

$$
D^\varphi_{i,k} =  \kappa_2 \nabla^2 \varphi_{i,k} - \kappa_4 \nabla^4 \varphi_{i,k} +
\frac{\partial }{\partial z} \left( \kappa_v \frac{\partial \varphi_{i,k}}{\partial z} \right)
$$ (discrete-tracer-diff)

#### Laplacian diffusion (del2)
The Laplacian may be written as the divergence of the gradient,

$$
 h_{i,k} \nabla \cdot \left( \kappa_{2,e,k} \nabla \varphi_{i,k} \right).
$$ (discrete-tracer-del2)

See [Omega V0 Section 3.3.2](OmegaV0ShallowWater.md#del2-tracer-diffusion) for details of this calculation.

#### Biharmonic diffusion (del4)
The biharmonic is a Laplacian operator applied twice,

$$
 - h_{i,k} \nabla \cdot \left( \kappa_{4,e,k} \nabla
\right[
\nabla \cdot \left(  \nabla \varphi_{i,k} \right)
\left]
 \right).
$$ (discrete-tracer-del4)

Each of these operators are written as horizontal stencils in the [Omega V0 Operator Formulation Section](OmegaV0ShallowWater.md#operator-formulation)

#### Vertical tracer diffusion
As discussed above in the [momentum section](#vertical-momentum-diffusion), vertical derivatives may be written in terms of $z$ or $p$,

$$
\frac{\partial }{\partial z} \left( \kappa_v \frac{\partial {\bf \varphi}}{\partial z} \right)
= \rho g^2 \frac{\partial }{\partial p} \left( \kappa_v \rho \frac{\partial {\bf \varphi}}{\partial p} \right)
$$ (discrete-tracer-vertdiff)
and $z$ is chosen. The second derivative stencil is

$$
h_{i,k} \frac{\partial }{\partial z} \left( \kappa_v \frac{\partial \varphi_{i,k}}{\partial z} \right)
=
\frac{h_{i,k}}{z_{i,k}^{top} - z_{i,k+1}^{top}} \left(
\kappa_{i,k}^{top}
\frac{\varphi_{i,k-1} - \varphi_k }{z_{i,k-1}^{mid} - z_k^{mid}}
 -
\kappa_{i,k+1}^{top}
\frac{\varphi_{i,k} - \varphi_{i,k+1} }{z_{i,k}^{mid} - z_{i,k+1}^{mid}}
\right).
$$ (discrete-tracer-vert-diff)

Like the momentum term, this is applied using a tridiagonal solver in the
[tridiagonal solver](TridiagonalSolver) in the implicit vertical mixing step.

### MPAS-Ocean Equations of Motion

The MPAS-Ocean layered formulation are provided here for reference. MPAS-Ocean solves for momentum, thickness, and tracers at layer $k$. These are continuous in the horizontal and discrete in the vertical.

$$
\frac{\partial {\bf u}_k}{\partial t}
+ \frac{1}{2}\nabla \left| {\bf u}_k \right|^2
+ ( {\bf k} \cdot \nabla \times {\bf u}_k) {\bf u}^\perp_k
+ f{\bf u}^{\perp}_k
+ \frac{w_k^{bot}{\bf u}_k^{bot} - w_k^{top}{\bf u}_k^{top}}{h_k}
=
- \frac{1}{\rho_0}\nabla p_k
- \frac{\rho g}{\rho_0}\nabla z^{mid}_k
+ \nu_h\nabla^2{\bf u}_k
+ \frac{\partial }{\partial z} \left( \nu_v \frac{\partial {\bf u}_k}{\partial z} \right),
$$ (mpaso-continuous-momentum)

$$
\frac{\partial h_k}{\partial t} + \nabla \cdot \left( h_k^e {\bf u}_k \right) + w_k^{bot} - w_k^{top} = 0,
$$ (mpaso-continuous-thickness)

$$
\frac{\partial h_k\varphi_k}{\partial t} + \nabla \cdot \left( h_k^e\varphi_k^e {\bf u}_k \right)
+ \varphi_k^{bot} w_k^{bot} - \varphi_k^{top} w_k^{top}
= \nabla\cdot\left(h_k^e \kappa_h \nabla\varphi_k \right)
+ h_k \frac{\partial }{\partial z} \left( \kappa_v \frac{\partial \varphi_k}{\partial z} \right).
$$ (mpaso-continuous-tracer)

The layer thickness $h$, vertical velocity $w$, pressure $p$, and tracer $\varphi$, are cell-centered quantities, while the horizontal velocity ${\bf u}$ and $e$ superscript are variables interpolated to the cell edges.


## 6. Variable Definitions

Table 1. Definition of variables. Geometric variables may be found in the [Omega V0 design document, Table 1](OmegaV0ShallowWater.md#variable-definitions)

| symbol  | name   | units    | location | name in code | notes  |
|---------------------|-----------------------------|----------|-|---------|-------------------------------------------------------|
|$D_{i,k}$   | divergence | 1/s      | cell | Divergence  |$D=\nabla\cdot\bf u$ |
|${\bf D}^u_{k} $, $ D^u_{e,k} $ | momentum dissipation terms | m/s$^2$ | edge | |see [Momentum Dissipation Section](#momentum-dissipation) |
|$ D_{e,k}^\varphi$ | tracer diffusion terms | | cell | |see [Tracer Diffusion Section](#tracer-diffusion) |
|$f_v$       | Coriolis parameter| 1/s      | vertex   | FVertex  |  $f = 2\Omega sin(\phi)$, $\Omega$ rotation rate, $\phi$ latitude|
|${\bf F}^u_{k} $, $ F^u_{e,k} $      | momentum forcing | m/s$^2$    | edge     |   | see [Momentum Forcing Section](#momentum-forcing) |
|$f_{eos}$ | equation of state | -  | any | function call | |
|$g$ | gravitational acceleration | m/s$^2$ | constant  | Gravity |
|$h_{i,k}$ | layer mass-thickness | kg/m$^2$  | cell | LayerThickness | see [](def-h) |
|$k$ | vertical index |  |
|${\bf k}$ | vertical unit vector |  |
|$K_{min}$ | shallowest active layer |  |
|$K_{max}$ | deepest active layer |  |
|$K_{i,k}$  | kinetic energy    | m$^2$/s$^2$  | cell     | KineticEnergyCell  |$K = \left\| {\bf u} \right\|^2 / 2$ |
|$p_{i,k}$ | pressure | Pa | cell | Pressure | see [](discrete-pressure) |
|$p^{floor}_i$ | bottom pressure | Pa | cell | PFloor | pressure at ocean floor
|$p^{surf}_i$ | surface pressure | Pa | cell | PSurface | due to atm. pressure, sea ice, ice shelves
|$q_{v,k}$ | potential vorticity         | 1/m/s    | vertex   | PotentialVorticity  |$q = \left(\zeta+f\right)/h$ |
|$Q^h_{i,k}$ | mass source and sink terms| kg/s/m$^2$ | cell |   |
|$Q^\varphi_{i,k}$ | tracer source and sink terms|kg/s/m$^2$ or similar| cell |   |
|$Ra$      | Rayleigh drag coefficient   | 1/s      | constant |   |  |
|$S_{i,k}$ | salinity | PSU | cell | Salinity | a tracer $\varphi$  |
|$t$       | time    | s        | none     |   |  |
|${\bf u}_k$   | velocity, vector form       | m/s      | - |   |  |
|$u_{e,k}$   | velocity, normal to edge      | m/s      | edge     | NormalVelocity  | |
|$u^\perp_{e,k}$   | velocity, tangential to edge      | m/s      | edge     | TangentialVelocity  |${\bf u}^\perp = {\bf k} \times {\bf u}$|
|$\alpha_{i,k}$ | specific volume | m$^3$/kg | cell  | SpecificVolume | $v = 1/\rho$ |
|$w_{i,k}$ | vertical velocity | m/s | cell  | VerticalVelocity | volume transport per m$^2$ |
|$z$ | vertical coordinate | m | - | | positive upward |
|$z^{top}_{i,k}$ | layer top z-location | m | cell | ZTop | see [](discrete-z) |
|$z^{mid}_{i,k}$ | layer mid-depth z-location | m | cell | ZMid |
|$z^{surf}_{i}$ | ocean surface, i.e. sea surface height  | m | cell | ZSurface | same as SSH in MPAS-Ocean |
|$z^{floor}_{i}$ | ocean floor z-location | m | cell | ZFloor | -bottomDepth from MPAS-Ocean |
|$\zeta_{v,k}$   | relative vorticity| 1/s      | vertex   |  RelativeVorticity |$\zeta={\bf k} \cdot \left( \nabla \times {\bf u}\right)$ |
|$\Theta_{i,k}$ | conservative temperature | C | cell  | Temperature  | a tracer $\varphi$ |
|$\kappa_2$| tracer diffusion  | m$^2$/s    | cell     |   |  |
|$\kappa_4$| biharmonic tracer diffusion | m$^4$/s    | cell     |   |  |
|$\kappa_v$| vertical tracer diffusion | m$^2$/s    | cell     |   |  |
|$\nu_2$   | horizontal del2 viscosity         | m$^2$/s    | edge     |   | |
|$\nu_4$   | horizontal biharmonic (del4) viscosity        | m$^4$/s    | edge     |   |  |
|$\nu_v$| vertical momentum diffusion | m$^2$/s    | edge       |   |  |
|$\varphi_{i,k}$ | tracer | kg/m$^3$ or similar | cell | | e.g. $\Theta$, $S$ |
|$\rho_{i,k}$ | density | kg/m$^3$ | cell  | Density |
|$\rho_0$ | Boussinesq reference density | kg/m$^3$ | |  constant |
|$\tau_i$ | wind stress | Pa=N/m$^2$ | edge |  SurfaceStress |
|$\Phi_{i,k}$ | geopotential| | cell | Geopotential |$\partial \Phi / \partial z = g$ for gravity |
|$\omega$   | mass transport | kg/s/m^2      | cell | VerticalTransport |$\omega=\rho w$|


## 7. Verification and Testing

Capability and testing are similar to [Petersen et al. 2015](http://www.sciencedirect.com/science/article/pii/S1463500314001796). The following tests are in idealized domains and do not require surface fluxes or surface restoring. For the following tests to show results comparable to those published with other models, the full dynamic sequence of density, pressure, momentum, and advection must work correctly. The successful completion of the following tests is a validation of the primitive equation functions in Omega 1.0. All of the following tests may exercise a linear equation of state or the nonlinear TEOS10. The first four tests quantify the anomalous mixing caused by the numerical schemes. The first five are on cartesian planes with regular hexagon meshes.

### Lock Exchange (Optional)
The Lock Exchange is the simplest possible test of a primitive equation model. There is an analytic formulation for the wave propagation speed. It is listed as optional because the Overflow tests the same dynamics.
See [Petersen et al. 2015](http://www.sciencedirect.com/science/article/pii/S1463500314001796) and the compass `lock_exchange` case.

### Overflow
The Overflow test case adds bathymetry to the Lock Exchange. It is a particularly effective test of vertical mass and tracer advection, and vertical mixing. It is useful to compare different vertical coordinates, like level (z- or p-level) versus terrain-following (sigma).
See [Petersen et al. 2015](http://www.sciencedirect.com/science/article/pii/S1463500314001796) and the compass `overflow` case.


### Internal Gravity Wave
The internal gravity wave tests horizontal and vertical advection.
See [Petersen et al. 2015](http://www.sciencedirect.com/science/article/pii/S1463500314001796) and the `internal_wave` case in both compass and polaris.

### Baroclinic Channel
This is the first test to add the Coriolis force and uses a three-dimensional domain. It is designed to result in an eddying simulation at sufficiently high resolution. This tests the combination of Coriolis and pressure gradient forces that produce geostrophic balance, as well as horizontal advection and dissipation for numerical stability.
See [Petersen et al. 2015](http://www.sciencedirect.com/science/article/pii/S1463500314001796) and the `baroclinic_channel` case in both compass and polaris.

### Seamount with zero velocity.
This is a 3D domain with a seamount in the center, where temperature and salinity are stratified in the vertical and constant in the horizontal. The test is simply that an initial velocity field of zero remains zero. For z-level layers the velocity trivially remains zero because the horizontal pressure gradient is zero. For tilted layers, this is a test of the pressure gradient error and the velocity is never exactly zero. This is a common test for sigma-coordinate models like ROMS because the bottom layers are extremely tilted along the seamount, but it is a good test for any model with tilted layers. Omega will use slightly tilted layers in p-star mode (pressure layers oscillating with SSH) and severely tilted layers below ice shelves, just like MPAS-Ocean. See [Ezer et al. 2002](https://www.sciencedirect.com/science/article/pii/S1463500302000033), [Haidvogel et al. 1993](https://journals.ametsoc.org/view/journals/phoc/23/11/1520-0485_1993_023_2373_nsofaa_2_0_co_2.xml), [Shchepetkin and McWilliams 2003](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2001JC001047), and previous MPAS-Ocean [confluence page](https://acme-climate.atlassian.net/wiki/spaces/OCNICE/blog/2015/11/19/40501447/MPAS-O+Sigma+coordinate+test+sea+mount).

### Cosine Bell on the Sphere
This test uses a fixed horizontal velocity field to test horizontal tracer advection. It is repeated from [Omega-0 design document](OmegaV0ShallowWater) and is important to conduct again as we convert Omega to a layered primitive-equation model. See `cosine_bell` case in both compass and polaris.

### Merry-Go-Round
This is an exact test for horizontal and vertical tracer advection. A fixed velocity field is provided, and a tracer distribution is advected around a vertical plane. See the `merry_go_round` test in compass, and the results on the [merry-go-round pull request](https://github.com/MPAS-Dev/compass/pull/108) and [compass port pull request](https://github.com/MPAS-Dev/compass/pull/452).


## References
This section is for references without webpage links. These are mostly textbooks.

- Cushman‐Roisin, B., & Beckers, J.M. (2011). Introduction to Geophysical Fluid Dynamics: Physical and Numerical Aspects. Academic Press.
- Gill, A. E. (2016). Atmosphere—Ocean dynamics. Elsevier.
- Kundu, P.K., Cohen, I.M., Dowling D.R. (2016) Fluid Mechanics 6th Edition, Academic Press.
- Pedlosky, J. (1987). Geophysical Fluid Dynamics (Vol. 710). Springer.
- Vallis, G. K. (2017). Atmospheric and oceanic fluid dynamics. Cambridge University Press.

## OLD: Continuous Equations

The continuous form of the conservation equations are as follows. See [Kundu et al. 2016](https://doi.org/10.1016/C2012-0-00611-4), chapter 4, eqns 4.7 and 4.22 or the [MOM5 manual](https://mom-ocean.github.io/assets/pdfs/MOM5_manual.pdf) eqn 7.7. This is before any assumptions are made, so this is a compressible, non-hydrostatic, non-Boussinesq fluid. Here all variables are a function of $(x,y,z)$, ${\bf u}_{3D}$ denotes the three-dimensional velocity vector, ${\bf u}_{3D} \otimes {\bf u}_{3D} = {\bf u}_{3D}{\bf u}_{3D}^T$ is the tensor product, $\nabla_{3D}$ is the three-dimensional gradient, $D/Dt$ is the material derivative, and other variables defined in the [Variable Definition Section](#variable-definitions) below.

momentum:

$$
\frac{D \rho {\bf u}_{3D} }{D t} \equiv
\frac{\partial \rho {\bf u}_{3D}}{\partial t}
 + \nabla_{3D} \cdot \left( \rho {\bf u}_{3D} \otimes {\bf u}_{3D}  \right)
  = - \nabla_{3D} p
   - \rho \nabla_{3D} \Phi
+ \rho {\bf D}^u_{3D} + \rho {\bf F}^u_{3D}
$$ (continuous-momentum)

mass:

$$
\frac{D \rho}{D t} \equiv
\frac{\partial \rho }{\partial t}
 + \nabla_{3D} \cdot \left( \rho  {\bf u}_{3D}  \right)
= 0
$$ (continuous-mass)

tracers:

$$
\frac{D \rho \varphi }{D t} \equiv
\frac{\partial \rho \varphi}{\partial t}
 + \nabla_{3D} \cdot \left( \rho \varphi {\bf u}_{3D}  \right)
= D^\varphi + Q^\varphi
$$ (continuous-tracer)

Here we have express the following terms as a general operators, with examples of specific forms provided below: the dissipation ${\bf D}^u$, momentum forcing ${\bf F}^u$, tracer diffusion $D^\varphi$, and tracer sources and sinks $Q^\varphi$. The graviational potential, $\Phi$, is written in a general form, and may include Earth's gravity, tidal forces, and self attraction and loading.

## OLD: Momentum Equations

Geophysical fluids such as the ocean and atmosphere are rotating and stratified, and horizontal velocities are orders of magnitude larger than vertical velocities. It is therefore convenient to separate the horizontal and vertical as ${\bf u}_{3D} = \left( {\bf u}, w \right)$ and $\nabla_{3D} = \left( \nabla_z, d/dz \right)$ where $z$ is the vertical direction in a local Cartesian coordinate system aligned with gravity (approximately normal to Earth's surface), and $w$ is the vertical velocity. The $z$ subscript on $\nabla_z$ is to remind us that this is the true horizontal gradient (perpendicular to $z$), as opposed to gradients within tilted layers used in the following section. The Earth's gravitational force is included as $\Phi_g = gz $ so that $ \nabla_{3D} \Phi_g =  g{\bf k}$. The rotating frame of reference results in the Coriolis force $f {\bf k} \times {\bf u} \equiv f {\bf u}^\perp$, where $f$ is the Coriolis parameter and ${\bf u}^\perp$ is the horizontal velocity vector rotated $90^\circ$ counterclockwise from $\bf u$ in the horizontal plane. See any textbook in the [References](#references) for a full derivation.


#### Assumptions

For a primitive equation ocean model, we assume the fluid is hydrostatic. For Omega we are not making the Boussinesq assumption, so all density-dependent terms use the full density. In particular, the density coefficient of the pressure gradient in [](#continuous-momentum) is not a constant, as it is in primitive equation models like POP and MPAS-Ocean.

**Hydrostatic:** Beginning with the vertical momentum equation,

$$
\frac{D \rho w }{D t}
  =
  -  \frac{\partial p}{\partial z} - \rho g + \rho {\bf k} \cdot {\bf D}^u_{3D} + \rho {\bf k} \cdot {\bf F}^u_{3D}
$$ (continuous-vert-mom)

assume that advection of vertical momentum $Dw/Dt$, dissipation, and forcing are small, and that the first order balance is between pressure gradient and buoyancy,

$$
\frac{\partial p}{\partial z}
  = - \rho g.
$$ (hydrostatic-balance)
We then integrate from $z$ to the surface $z^{surf}$ to obtain the hydrostatic pressure equation,

$$
p(x,y,z) = p^{surf}(x,y) + \int_{z}^{z^{surf}} \rho g dz'.
$$ (continuous-hydrostatic-pressure)

The constitutive equation is the equation of state,

$$
\rho = f_{eos}(p,\Theta,S).
$$ (continuous-eos)

where conservative temperature, $\Theta$, and absolute salinity, $S$, are examples of tracers $\varphi$.

The Boussinesq primitive equations also make an incompressibility assumption, which is identical to an assumption of constant density. Non-Boussinesq models do not make that assumption and are not explicitly incompressible. However, the mass conservation equation [](continuous-mass), along with an equation of state for sea water where density only varies slightly, results in a fluid that is nearly incompressible.

A concern when using the full, compressible continuity equation is that this might support acoustic waves with wave speeds on the order of 1500 m/s, requiring an extremely small time step. According to [Griffies and Adcroft (2008)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/177GM18) and [de Szoeke and Samelson (2002)](https://doi.org/10.1175/1520-0485(2002)032%3C2194:TDBTBA%3E2.0.CO;2), the hydrostatic approximation removes vertical sound waves, leaving only barotropic acoustic modes called Lamb waves.  Fortunately, the Lamb waves can be "subsumed" into the external gravity mode because the scale height of the ocean is much larger (200 km) than its depth (~5 km).  This suggests that Lamb waves should not produce any additional constraints on our barotropic time step (though we should keep an eye on this).  For more details on Lamb waves, see [Dukowicz (2013)](https://doi.org/10.1175/MWR-D-13-00148.1)

#### Momentum Advection

Here we expand the momentum advection in continuous form. This will be a useful reference when we expand the layered version in a later section. The momentum advection,

$$
\nabla_{3D} \cdot \left( \rho {\bf u}_{3D} \otimes {\bf u}_{3D}  \right)
 = \nabla_{3D} \cdot \left( \rho {\bf u}_{3D}  {\bf u}_{3D}^T  \right),
$$ (advection)

may be written out fully as the three $(x,y,z)$ cartesian components coordinates as

$$
&\partial_x \left( \rho u u\right) + \partial_y \left( \rho v u\right) + \partial_z \left( \rho w u\right) \\
&\partial_x \left( \rho u v\right) + \partial_y \left( \rho v v\right) + \partial_z \left( \rho w v\right) \\
&\partial_x \left( \rho u w\right) + \partial_y \left( \rho v w\right) + \partial_z \left( \rho w w\right) .
$$ (advection-written-out)

The third line is the vertical component, and was assumed to be small in the previous section. The first two lines are the horizontal components may be written as

$$
\nabla_z \cdot \left( \rho {\bf u} \otimes {\bf u}  \right) + \partial_z \left( w \rho {\bf u}\right)
$$ (adv2d)

where ${\bf u} = (u,v)$ is the horizontal velocity vector and $\nabla_z=(\partial_x,\partial_y)$ is the horizontal gradient, and the tensor product ${\bf u} \otimes {\bf u}={\bf u}  {\bf u}^T$. Using the product rule, this can be expanded as

$$
\nabla_z \cdot \left( \rho {\bf u} \otimes {\bf u}  \right) + \partial_z \left( w \rho {\bf u}\right)
&= \left( \nabla_z \cdot  {\bf u}  \right) \rho {\bf u}
+  {\bf u} \cdot \nabla_z \left( \rho {\bf u} \right)
+ \partial_z \left( w \rho {\bf u}\right) \\
&= \left( \nabla_z \cdot   {\bf u}  \right) \rho {\bf u}
+ \left( {\bf u} \cdot \nabla_z  \rho  \right) {\bf u}
+ \left( {\bf u} \cdot \nabla_z  {\bf u} \right)  \rho
+ \partial_z \left( w \rho {\bf u}\right)
$$ (adv2d-prod)

The term ${\bf u} \cdot \nabla_z {\bf u}$ may be replaced with the vector identity

$$
\begin{aligned}
{\bf u} \cdot \nabla_z {\bf u}
&= (\nabla_z \times {\bf u}) \times {\bf u} + \nabla_z \frac{|{\bf u}|^2}{2} \\
&= \left( \boldsymbol{k} \cdot (\nabla_z \times {\bf u})\right)
\left( \boldsymbol{k} \times {\bf u} \right) + \nabla_z \frac{|{\bf u}|^2}{2} \\
&= \zeta {\bf u}^{\perp} + \nabla_z K,
\end{aligned}
$$ (advection-identity)

where $\zeta$ is relative vorticity and $K$ is kinetic energy. This step separates the horizontal advection into non-divergent and non-rotational components, which is useful in the final TRiSK formulation.

#### Final Continuous Equations

The final form of the continuous conservation equations for a non-Boussinesq, hydrostatic ocean are

momentum:

$$
\frac{\partial (\rho \mathbf{u})}{\partial t}
+ \left( \nabla_z \cdot  {\bf u}  \right) \rho {\bf u}
+  {\bf u} \cdot \nabla_z \left( \rho {\bf u} \right)
+ \partial_z \left( w \rho \mathbf{u} \right)
+ f \rho {\bf u}^\perp
  = -  \nabla_z p
  - \rho\nabla_z \Phi
+ \rho{\bf D}^u + \rho{\bf F}^u
$$ (continuous-momentum-final)

mass:

$$
\frac{\partial \rho }{\partial t}
 + \nabla_z \cdot \left( \rho  {\bf u}  \right)
 + \frac{\partial}{\partial z} \left(  \rho  w \right)
= 0
$$ (continuous-mass-final)

tracers:

$$
\frac{\partial \rho \varphi}{\partial t}
 + \nabla_z \cdot \left( \rho \varphi {\bf u}  \right)
 + \frac{\partial}{\partial z} \left(  \rho \varphi w \right)
= D^\varphi + Q^\varphi
$$ (continuous-tracer-final)

equation of state:

$$
\rho = f_{eos}(p,\Theta,S).
$$ (continuous-eos-final)

hydrostatic pressure:

$$
p(x,y,z) = p^{surf}(x,y) + \int_{z}^{z^{surf}} \rho g dz'.
$$ (continuous-hydrostatic-pressure-final)

Here the $\nabla_z$ operators are exactly horizontal; we expand terms for tilted layers in the following section. The momentum diffusion terms include Laplacian (del2), biharmonic (del4), and vertical viscosity,

$$
{\bf D}^u
= \nu_2 \nabla_z^2 {\bf u} - \nu_4 \nabla_z^4 {\bf u}
 + \frac{\partial }{\partial z} \left( \nu_v \frac{\partial {\bf u}}{\partial z} \right)
$$ (continuous-h_mom_diff)
and may also include a Rayleigh drag and eventually parameterizations. Momentum forcing is due to wind stress and bottom drag. Similarly, the tracer diffusion terms include Laplacian (del2), and vertical viscosity,

$$
D^\varphi =
 \nabla_z\cdot\left(\rho \kappa_2 \nabla_z\varphi \right)
+ \rho \frac{\partial }{\partial z}
  \left( \kappa_v \frac{\partial \varphi}{\partial z} \right),
$$ (continuous-v_tr_diff)

and may also include a biharmonic (del4) term and parameterizations such as Redi mixing. Sources and sinks include surface fluxes from the atmosphere and land, and bio-geo-chemical reactions.
All of the diffusion and forcing terms are written in more detail with the [Discrete Equations](#discrete-equations) below.

## OLD: Layered Equations

Here we derive the layered equations by discretizing in the vertical, while the horizontal remains continuous. We discretize by integrating in the vertical from the lower surface $z=z_k^{bot}(x,y)$ to $z=z_k^{top}(x,y)$ for the layer with index $k$, as described in [Ringler et al. 2013](https://www.sciencedirect.com/science/article/pii/S1463500313000760) Appendix A.2. Equivalently, we can vertically integrate from a deeper pressure surface $p=p_k^{bot}(x,y)$ (higher pressure) to $p=p_k^{top}(x,y)$ where $p$ and $z$ are related by the hydrostatic pressure equation [](#continuous-hydrostatic-pressure).

### Layer Integration

For non-Boussinesq layered equations we begin by defining the pressure-thickness of layer $k$ as

$$
h_k(x,y,t)
 \equiv \int_{z_k^{bot}}^{z_k^{top}} \rho dz = \frac{1}{g} \int_{p_k^{top}}^{p_k^{bot}} dp.
$$ (def-h)

The letter $h$ is used because this is the familiar variable for thickness in Boussinesq primitive equations. In the Boussinesq case the thickness equation describes conservation of volume because $h_k^{Bouss}$ is volume normalized by horizontal area, resulting in a height,

$$
h_k^{Bouss}(x,y,t)
 \equiv \int_{z_k^{bot}}^{z_k^{top}} dz.
$$ (def-h-bouss)

In this document we remain with the more general non-Boussinesq case, where $h_k$ is mass per unit area (kg/m$^2$). Since horizontal cell area remains constant in time, the thickness equation is a statement of conservation of mass.

Throughout this derivation we can write all equations equivalently in $z$-coordinates (depth), or in $p$-coordinates (pressure). From the hydrostatic equation, any quantity $\varphi$ may be integrated in $z$ or $p$ as

$$
\int_{z_k^{bot}}^{z_k^{top}} \varphi \rho dz
= \frac{1}{g} \int_{p_k^{top}}^{p_k^{bot}} \varphi dp.
$$(depth-pressure-integral-conversion)

We can convert from an interfacial depth surface $z^{top}$ to a pressure surface $p^{top}$ with the hydrostatic equation [](continuous-hydrostatic-pressure-final):

$$
p^{top}(x,y) =  p^{surf}(x,y) + \int_{z^{top}(x,y)}^{z^{surf}(x,y)} \rho g dz
$$ (def-p-surf)

where $z^{surf}$ is the sea surface height and $p^{surf}$ is the surface pressure at $z^{surf}$ imposed by the atmosphere or floating ice. Note that pressure increases with depth.  This means that positive $p$ points downward, so that the $top$ and $bot$ extents of the integration limits are flipped in [](#depth-pressure-integral-conversion).

For any three-dimensional quantity $\varphi(x,y,z,t)$, the mass-thickness-averaged quantity in layer $k$ is defined as

$$
\varphi_k(x,y,t)
\equiv \frac{\int_{z_k^{bot}}^{z_k^{top}} \rho \varphi dz}{\int_{z_k^{bot}}^{z_k^{top}} \rho dz}
= \frac{\int_{z_k^{bot}}^{z_k^{top}} \rho \varphi dz}{h_k}
$$(def-mass-thickness-average)

At this point our derivation has not made any assumptions about density, and may be used for both Boussinesq and non-Boussinesq fluids. A Boussinesq derivation would now assume small variations in density and replace $\rho(x,y,z,t)$ with a constant $\rho_0$ everywhere but the pressure gradient coefficient. In that case $\rho$ divides out in [](#def-mass-thickness-average) the Boussinesq layer quantities would simply be thickness-weighted averages.

We can now derive the layered equations. Integrate the continuous equations [](continuous-momentum-final), [](continuous-mass-final), [](continuous-tracer-final) in $z$ from $z_k^{bot}$ to $z_k^{top}$,

$$
\int_{z_k^{bot}}^{z_k^{top}} \frac{\partial \rho {\bf u}}{\partial t}  dz
+  \int_{z_k^{bot}}^{z_k^{top}}  \left( \nabla_z \cdot  {\bf u}  \right) \rho {\bf u} dz
+  \int_{z_k^{bot}}^{z_k^{top}}  {\bf u} \cdot \nabla_z \left( \rho {\bf u} \right) dz
+  \int_{z_k^{bot}}^{z_k^{top}} \frac{\partial }{\partial z} \left(  \rho {\bf u} w \right) dz
+  \int_{z_k^{bot}}^{z_k^{top}}
 f \rho {\bf u}^\perp dz
  =  \int_{z_k^{bot}}^{z_k^{top}}  \left[- \nabla_z p
  - \rho \nabla_z \Phi
+ \rho {\bf D}^u + \rho {\bf F}^u \right] dz
$$ (z-integration-momentum)

$$
\frac{\partial }{\partial t}  \int_{z_k^{bot}}^{z_k^{top}} \rho dz
 +\nabla_z \cdot \left( \int_{z_k^{bot}}^{z_k^{top}} \rho  {\bf u} dz \right)
 + \int_{z_k^{bot}}^{z_k^{top}}\frac{\partial}{\partial z} \left(  \rho  w \right) dz
= 0
$$ (z-integration-mass)

$$
\frac{\partial }{\partial t} \int_{z_k^{bot}}^{z_k^{top}} \rho \varphi  dz
 + \nabla_z \cdot \left( \int_{z_k^{bot}}^{z_k^{top}} \rho \varphi {\bf u} dz  \right)
 + \int_{z_k^{bot}}^{z_k^{top}} \frac{\partial}{\partial z} \left(  \rho \varphi w \right) dz
= \int_{z_k^{bot}}^{z_k^{top}} \left( D^\varphi + Q^\varphi \right) dz
$$ (z-integration-tracers)

This results in conservation equations that are valid over the layer. The momentum variables are simply vertically averaged. Tracer variables are vertically mass-averaged, as defined in [](def-mass-thickness-average).

In order to deal with nonlinear terms where we take the integrals of products, we may assume that the variables are piecewise constant in the vertical within each layer, i.e.

$$
a(x,y,z,t) = a_k(x,y,t) \in [z_k^{bot}, z_k^{top}).
$$ (discrete-a)

Then variables may come out of the integral as needed. For example, we can handle the advection terms by making this assumptoin for ${\bf u}$.

$$
 \int_{z_k^{bot}}^{z_k^{top}} \nabla_z \cdot \left( \rho \varphi {\bf u}  \right) dz
&= \nabla_z \cdot \left( \int_{z_k^{bot}}^{z_k^{top}} \rho \varphi {\bf u} dz  \right)\\
&= \nabla_z \cdot \left( {\bf u}_k \int_{z_k^{bot}}^{z_k^{top}} \rho \varphi  dz  \right)\\
&= \nabla_z \cdot \left( {\bf u}_k  h_k \varphi_k \right).\\
$$ (tracer-adv)

Likewise, for the horizontal momentum advection,

$$
\int_{z_k^{bot}}^{z_k^{top}}  \left( \nabla_z \cdot  {\bf u}  \right) \rho {\bf u} dz
+  \int_{z_k^{bot}}^{z_k^{top}}  {\bf u} \cdot \nabla_z \left( \rho {\bf u} \right) dz
&=   \left( \nabla_z \cdot  {\bf u}_k  \right) \int_{z_k^{bot}}^{z_k^{top}}  \rho {\bf u} dz
+    {\bf u}_k \cdot \nabla_z \left( \int_{z_k^{bot}}^{z_k^{top}}\rho {\bf u} dz\right)  \\
&=   \left( \nabla_z \cdot  {\bf u}_k  \right) h_k {\bf u}_k
+    {\bf u}_k \cdot \nabla_z \left( h_k {\bf u}_k \right)  \\
$$ (mom-adv)

Our governing equations are now discrete in the vertical, but remain continuous in the horizontal and in time.


This results in the layered conservation equations. This system is now discrete in the vertical, but remains continuous in the horizontal and in time.

$$
\frac{\partial h_k {\bf u}_k}{\partial t}
+   \left( \nabla_z \cdot  {\bf u}_k  \right) h_k {\bf u}_k
+    {\bf u}_k \cdot \nabla_z \left( h_k {\bf u}_k \right)
+ \left[ \rho_k w_k {\bf u}_k \right]^{top} - \left[ \rho_k w_k {\bf u}_k \right]^{bot}
+ f h_k {\bf u}_k^{\perp}
=
- \nabla_z p_k
- h_k \nabla_z \Phi_k
+ h_k {\bf D}_k^u + h_k {\bf F}_k^u
$$ (layered-momentum-1)

$$
\frac{\partial h_k}{\partial t} + \nabla_z \cdot \left(h_k {\bf u}_k\right) + \left[ \rho_k w_k \right]^{top} - \left[ \rho_k w_k \right]^{bot}= Q^h_k
$$ (layered-mass-1)

$$
\frac{\partial h_k \varphi_k}{\partial t} + \nabla_z \cdot \left(h_k {\bf u}_k \varphi_k\right)
+ \left[ \varphi_k \rho_k w_k \right]^{top} - \left[ \varphi_k \rho_k w_k \right]^{bot}
=  D^\varphi_k + Q^\varphi_k.
$$ (layered-tracer-1)

The term $\left( \nabla_z \cdot  {\bf u}_k  \right) h_k {\bf u}_k$ is assumed to be small for seawater, which is nearly incompressible. In the Boussinesq approximation this term is formally zero due to the incompressibility assumption. The remaining terms from the momentum material derivative  may be rewritten as

$$
\frac{\partial h_k {\bf u}_k}{\partial t}
+    {\bf u}_k \cdot \nabla_z \left( h_k {\bf u}_k \right)
&= h_k \frac{\partial {\bf u}_k}{\partial t}
+ {\bf u}_k\frac{\partial h_k} {\partial t}
+  h_k  {\bf u}_k \cdot \nabla_z \left(  {\bf u}_k \right)
+    {\bf u}_k \cdot \nabla_z \left( h_k  \right) {\bf u}_k\\
&= h_k \frac{\partial {\bf u}_k}{\partial t}
+ {\bf u}_k\frac{\partial h_k} {\partial t}
+  h_k  {\bf u}_k \cdot \nabla_z {\bf u}_k
+    {\bf u}_k \cdot \nabla_z \left( h_k  \right) {\bf u}_k \\
&= h_k \left[ \frac{\partial {\bf u}_k}{\partial t}
+   {\bf u}_k \cdot \nabla_z {\bf u}_k \right]
+ {\bf u}_k \left[\frac{\partial h_k} {\partial t}
+    {\bf u}_k \cdot \nabla_z  h_k   \right] \\
&= h_k \left[ \frac{\partial {\bf u}_k}{\partial t}
+ \zeta {\bf u}^{\perp} + \nabla_z K
   \right]
+ {\bf u}_k \left[\frac{\partial h_k} {\partial t}
+    {\bf u}_k \cdot \nabla_z  h_k   \right] \\
$$ (mom-h-adv)

The terms in the second square brackets drop out due to conservation of mass (note: I need to check on vertical thickness flux and source terms). We now divide the momentum equation by $h_k$ to obtain

$$
\frac{\partial {\bf u}_k}{\partial t}
+ \zeta {\bf u}^{\perp} + \nabla_z K
+ \frac{\left[ \rho_k w_k {\bf u}_k \right]^{top} - \left[ \rho_k w_k {\bf u}_k \right]^{bot}}{h_k}
+ f {\bf u}_k^{\perp}
=
- \nabla_z p_k
- \nabla_z \Phi_k
+ {\bf D}_k^u + {\bf F}_k^u
$$ (layered-momentum-2)

One could derive the layered equations by integrating in pressure rather than $z$. In that case, one would have multiplied  [](z-integration-momentum), [](z-integration-mass), [](z-integration-tracers) by $g$ and expressed the integrals in terms of pressure, with integration limits from $p_k^{top}$ to $p_k^{bot}$.
Derivations of ocean model equations with pressure as the vertical variable may be found in [de Szoeke and Samelson 2002](https://journals.ametsoc.org/view/journals/phoc/32/7/1520-0485_2002_032_2194_tdbtba_2.0.co_2.xml) and [Losch et al. 2003](https://journals.ametsoc.org/view/journals/phoc/34/1/1520-0485_2004_034_0306_hsacgc_2.0.co_2.xml).

The vertical advection terms use the fundamental theorem of calculus to integrate a derivative in $z$, resulting in boundary conditions at the layer interfaces. This makes intuitive sense as a mass balance. For example, in the absence of horizontal advection and sources, the mass equation is simply

$$
\frac{\partial h_k}{\partial t}
= \left[ \rho_k w_k \right]^{bot} - \left[ \rho_k w_k \right]^{top}.
$$ (layered-mass-vert-adv)

This states that the change in mass in the layer is the incoming mass from below minus the outgoing mass above, since the vertical velocity $w$ is positive upwards. The vertical advection of momentum and tracers have a similar interpretation.

### Vertical Transport
The integration in [](#z-integration-momentum) to [](#z-integration-tracers) changes the vertical velocity $w$ in m/s to a vertical mass-thickness transport $\omega=\rho w$ in kg/s/m$^2$. Here $w$ is the Latin letter and $\omega$ is the Greek letter omega.  One can think of fluid velocity $w$ as a volume transport, normalized by area, in units of length per time. Analogously, $\omega$ is mass-thickness transport, which is mass transport per unit area (kg/s/m$^2$). The variables $w$ and $\omega$ have the same sign convention of upward (positive $z$) for positive transport.

### Final Layered Equations

mrp temp:

The momentum equation can be rewritten using the product rule on $\rho {\bf u}$, mass conservation, and dividing by $\rho$, as:

$$
\frac{D {\bf u}_{3D} }{D t} \equiv
\frac{\partial {\bf u}_{3D}}{\partial t}
+ {\bf u}_{3D}\cdot \nabla_{3D} {\bf u}_{3D}
  = - \frac{1}{\rho} \nabla_{3D} p
   -  \nabla_{3D} \Phi
+ {\bf D}^u + {\bf F}^u
$$ (continuous-momentum-rho)

mrp temp:

The advection term may be separated into horizontal and vertical parts as

$$
{\bf u}_{3D}\cdot \nabla_{3D} {\bf u}_{3D}
=
{\bf u} \cdot \nabla_z {\bf u} + w \frac{\partial {\bf u}}{\partial z}.
$$ (advection-3d2d)

The horizontal component may be replaced with the vector identity

$$
\begin{aligned}
{\bf u} \cdot \nabla_z {\bf u}
&= (\nabla_z \times {\bf u}) \times {\bf u} + \nabla_z \frac{|{\bf u}|^2}{2} \\
&= \left( \boldsymbol{k} \cdot (\nabla_z \times {\bf u})\right)
\left( \boldsymbol{k} \times {\bf u} \right) + \nabla_z \frac{|{\bf u}|^2}{2} \\
&= \zeta {\bf u}^{\perp} + \nabla_z K,
\end{aligned}
$$ (advection-identity)

where $\zeta$ is relative vorticity and $K$ is kinetic energy. This step separates the horizontal advection into non-divergent and non-rotational components, which is useful in the final TRiSK formulation.

momentum:

$$
\frac{\partial {\bf u}_k}{\partial t}
+ q_k h_k {\bf u}_k^{\perp}
+ \frac{\left[ \alpha_k \omega_k {\bf u}_k \right]^{top} - \left[ \alpha_k \omega_k {\bf u}_k \right]^{bot}}{h_k}
=
- \alpha_k \nabla_r p_k - \nabla_r \Phi_k
- \nabla_r K_k
+ {\bf D}_k^u + {\bf F}_k^u
$$ (layered-momentum)

mass:

$$
\frac{\partial h_k}{\partial t} + \nabla_r \cdot \left(h_k {\bf u}_k\right) + \omega_k^{top} - \omega_k^{bot}= Q^h_k
$$ (layered-mass)

tracers:

$$
\frac{\partial h_k \varphi_k}{\partial t} + \nabla_r \cdot \left(h_k {\bf u}_k \varphi_k\right)
+ \varphi_k^{top} \omega_k^{top} - \varphi_k^{bot}\omega_k^{bot}
=  D^\varphi_k + Q^\varphi_k.
$$ (layered-tracer)

The superscripts $top$ and $bot$ mean that the layered variable is interpolated to the top and bottom layer interface, respectively. In [](layered-momentum) the non-divergent momentum advection and Coriolis term were combined and expressed in terms of the potential vorticity $q_k$,

$$
 \zeta_k {\bf u}_k^\perp + f {\bf u}_k^\perp = \frac{\zeta_k + f }{h_k}h_k{\bf u}_k^\perp \equiv q_k h_k {\bf u}^\perp_k.
$$ (potential-vort-adv)

Surface fluxes $Q^h_k$ have been added to the mass equation for precipitation, evaporation, and river runoff.  These fluxes, like mass transport $\omega$, are in units of kg/s/m$^2$.


<!--
Comments for later consideration:
## Outstanding questions
Treatment of outcropping isobars and boundary conditions. Treatment of topography/partial bottom cells.
Is there an equivalent p* coordinate? This would imply stretching/squeezing for both the top (external forcing from atmosphere and sea ice) and the bottom ($p_{surf}$ variations from forcing and barotropic thickness variations) pressure layers. Other examples of coordinate choices include:  a normalized pressure-$\sigma$ coordinates ([Huang et al, 2001](https://link.springer.com/article/10.1007/s00376-001-0001-9) is the only ocean model example I could find), where $\sigma = \frac{p - p_{surf}}{p_{bt}}$, and $p_{bt} = p_b - p_{surf}$. Another possibility is to use the coordinate $\sigma = p/p_{surf}$.

## Main References
This document is based on [de Szoeke and Samelson, 2002](https://journals.ametsoc.org/view/journals/phoc/32/7/1520-0485_2002_032_2194_tdbtba_2.0.co_2.xml), [Losch et al. 2004](https://journals.ametsoc.org/downloadpdf/view/journals/phoc/34/1/1520-0485_2004_034_0306_hsacgc_2.0.co_2.pdf), and [Adcroft's 2004 lecture notes](https://ocw.mit.edu/courses/12-950-atmospheric-and-oceanic-modeling-spring-2004/resources/lec25). Other references of interest include: [Hallberg et al., 2024](https://www.cesm.ucar.edu/sites/default/files/2024-02/2024-OMWG-B-Hallberg.pdf), [Losch et al poster](http://mitgcm.org/~mlosch/NBposter.pdf), [MOM6 code](https://github.com/mom-ocean/MOM6/blob/main/src/core/MOM_PressureForce_Montgomery.F90)... MOM6 has 2 different Pressure Gradient implementations: Finite Volume (including an analytical integration based on [Adcroft et al., 2007](https://www.sciencedirect.com/science/article/pii/S1463500308000243)) and Montgomery-potential form (based on [Hallberg, 2005](https://www.sciencedirect.com/science/article/pii/S1463500304000046)), both of which can be used in non-Boussinesq mode.

--->
