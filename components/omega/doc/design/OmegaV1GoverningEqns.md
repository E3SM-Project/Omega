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

1. Omega is non-Boussinesq. This means that the full 3D density is used everywhere, and results in a mass-conserving model. MPAS-Ocean and POP were Boussinesq, so that a reference density $\rho_0$ is used in the pressure gradient term, and were therefore volume-conserving models. In Omega the layered mass-conservation equation is in terms of mass-thickness ($h=\rho \Delta z$). In MPAS-Ocean the simple thickness ($h=\Delta z$) is the prognostic volume variable (normalized by horizontal cell area).
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
See discussion in introduction. The non-Boussinesq formulation uses the full density throughout resulting in governing equations that conserve mass rather than volume. The substantial change is that the "thickness" variable, $h$, is now a mass-thickness, $h = \rho \Delta z$. This is explained in the derivation of the [layered equations below](#layered-equations).

### Omega will use TEOS10 for the equation of state.
See additional [EOS design document](EOS)

### Omega-1.0 will add new terms for the pressure gradient, vertical mixing, and vertical advection.
See forthcoming design documents on the pressure gradient, vertical mixing, and vertical advection.

## 3. Continuous Equations

The continuous form of the conservation equations are as follows. See [Kundu et al. 2024](https://www.amazon.com/Fluid-Mechanics-Pijush-K-Kundu/dp/012405935X), chapter 4, eqns 4.7 and 4.22 or the [MOM5 manual](https://mom-ocean.github.io/assets/pdfs/MOM5_manual.pdf) eqn 7.7. This is before any assumptions are made, so this is a compressible, non-hydrostatic, non-Boussinesq fluid. Here all variables are a function of $(x,y,z)$, ${\bf u}_{3D}$ denotes the three-dimensional velocity vector, ${\bf u}_{3D} \otimes {\bf u}_{3D} = {\bf u}_{3D}{\bf u}_{3D}^T$ is the tensor product, $\nabla_{3D}$ is the three-dimensional gradient, $D/Dt$ is the material derivative, and other variables defined in the [Variable Definition Section](#variable-definitions) below.

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

Geophysical fluids such as the ocean and atmosphere are rotating and stratified, and horizontal velocities are orders of magnitude larger than vertical velocities. It is therefore convenient to separate the horizontal and vertical as ${\bf u}_{3D} = \left( {\bf u}, w \right)$ and $\nabla_{3D} = \left( \nabla_z, d/dz \right)$ where $z$ is the vertical direction in a local Cartesian coordinate system aligned with gravity (approximately normal to Earth's surface), and $w$ is the vertical velocity. The $z$ subscript on $\nabla_z$ is to remind us that this is the true horizontal gradient (perpendicular to $z$), as opposed to gradients within tilted layers used in the following section. The Earth's gravitational force is included as $\Phi_{gravity} = gz $ so that $ \nabla_{3D} \Phi_{gravity} =  g{\bf k}$. The rotating frame of reference results in the Coriolis force $f {\bf k} \times {\bf u} \equiv f {\bf u}^\perp$, where $f$ is the Coriolis parameter and ${\bf u}^\perp$ is the horizontal velocity vector rotated $90^\circ$ counterclockwise from $\bf u$ in the horizontal plane. See any textbook in the [References](#references) for a full derivation.


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
\nabla_z \cdot \left( \rho {\bf u} \otimes {\bf u}  \right) + \partial_z \left( \rho w {\bf u}\right)
$$ (adv2d)

where ${\bf u} = (u,v)$ is the horizontal velocity vector and $\nabla_z=(\partial_x,\partial_y)$ is the horizontal gradient, and the tensor product ${\bf u} \otimes {\bf u}={\bf u}  {\bf u}^T$. Using the product rule, this can be expanded as

$$
\nabla_z \cdot \left( \rho {\bf u} \otimes {\bf u}  \right) + \partial_z \left( \rho w {\bf u}\right)
&= \left( \nabla_z \cdot  {\bf u}  \right) \rho {\bf u} 
+  {\bf u} \cdot \nabla_z \left( \rho {\bf u} \right) 
+ \partial_z \left( \rho w {\bf u}\right) \\
&= \left( \nabla_z \cdot   {\bf u}  \right) \rho {\bf u} 
+ \left( {\bf u} \cdot \nabla_z  \rho  \right) {\bf u} 
+ \left( {\bf u} \cdot \nabla_z  {\bf u} \right)  \rho
+ \partial_z \left( \rho w {\bf u}\right) 
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
+ \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u})
+ \partial_z (\rho \mathbf{u} w)
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

## 4. Layered Equations

Here we derive the layered equations by discretizing in the vertical, while the horizontal remains continuous. We discretize by integrating in the vertical from the lower surface $z=z_k^{bot}(x,y)$ to $z=z_k^{top}(x,y)$ for the layer with index $k$, as described in [Ringler et al. 2013](https://www.sciencedirect.com/science/article/pii/S1463500313000760) Appendix A.2. Equivalently, we can vertically integrate from a deeper pressure surface $p=p_k^{bot}(x,y)$ (higher pressure) to $p=p_k^{top}(x,y)$ where $p$ and $z$ are related by the hydrostatic pressure equation [](#continuous-hydrostatic-pressure).

### Layer Integration

For non-Boussinesq layered equations we begin by defining the mass-thickness of layer $k$ as

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
\frac{\partial }{\partial t} \int_{z_k^{bot}}^{z_k^{top}} \rho {\bf u} dz
+  \int_{z_k^{bot}}^{z_k^{top}} \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) dz
+  \int_{z_k^{bot}}^{z_k^{top}} \frac{\partial }{\partial z} \left(  \rho w{\bf u} \right) dz
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

Then variables may come out of the integral as needed. For example,

$$
\int_{z_k^{bot}}^{z_k^{top}}   \zeta {\bf u}^\perp dz
&= \int_{z_k^{bot}}^{z_k^{top}}   \zeta_k {\bf u}^\perp_k dz \\
&= \zeta_k {\bf u}^\perp_k \int_{z_k^{bot}}^{z_k^{top}}  dz \\
&= \zeta_k {\bf u}^\perp_k \Delta z_k,
$$ (nonlinear-z)

where $\Delta z_k\equiv z_k^{top}-z_k^{bot}$. In the momentum equation this results in a $\Delta z_k$ factor in every term except vertical advection, which is divided out.

Our governing equations are now discrete in the vertical, but remain continuous in the horizontal and in time.


This results in conservation equations that are valid over the layer. The momentum variables are simply vertically averaged. Tracer variables are vertically mass-averaged, as defined in [](def-mass-thickness-average). We assume that variables such as ${\bf u}$, $\zeta$, and $K$ are constant in the vertical so that they may be taken out of the integral where products occur.  The subscript $k$ denotes a variable average over layer $k$. This system is now discrete in the vertical, but remains continuous in the horizontal and in time.

$$
\frac{\partial h_k {\bf u}_k}{\partial t}
+ \nabla_z \cdot (h_k \mathbf{u}_k \otimes \mathbf{u}_k)
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


Equivalently, one could have multiplied  [](z-integration-momentum), [](z-integration-mass), [](z-integration-tracers) by $g$ and expressed the integrals in terms of pressure, with integration limits from $p_k^{top}$ to $p_k^{bot}$.
Derivations of ocean model equations with pressure as the vertical variable may be found in [de Szoeke and Samelson 2002](https://journals.ametsoc.org/view/journals/phoc/32/7/1520-0485_2002_032_2194_tdbtba_2.0.co_2.xml) and [Losch et al. 2003](https://journals.ametsoc.org/view/journals/phoc/34/1/1520-0485_2004_034_0306_hsacgc_2.0.co_2.xml).

The vertical advection terms use the fundamental theorem of calculus to integrate a derivative in $z$, resulting in boundary conditions at the layer interfaces. This makes intuitive sense as a mass balance. For example, in the absence of horizontal advection and sources, the mass equation is simply

$$
\frac{\partial h_k}{\partial t}
= \left[ \rho_k w_k \right]^{bot} - \left[ \rho_k w_k \right]^{top}.
$$ (layered-mass-vert-adv)

This states that the change in mass in the layer is the incoming mass from below minus the outgoing mass above, since the vertical velocity $w$ is positive upwards. The vertical advection of momentum and tracers have a similar interpretation.

### Questions

1. How do integrals in  [](#z-integration-momentum) to [](#z-integration-tracers) turn into the layered form [](layered-momentum-1) to [](layered-tracer-1) when there are products in the terms? For example, $\zeta_k {\bf u}^\perp_k$, $\nabla K$ and $\rho \varphi {\bf u}$.
1. Using the definition of the layer average, the momentum equation [](layered-momentum-1) should also be thickness-weighted. Perhaps that goes away with conservation of mass?
1. How do I justify changing from $\nabla_z$ to $\nabla_r$ for terms other than the pressure gradient? Particularly $\nabla K$?

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
&=-v \nabla_r p - \nabla_r \Phi_{gravity}\\
\end{aligned}
$$ (gradp)

where we have substituted hydrostatic balance [](hydrostatic-balance), specific volume $v\equiv 1/\rho$, and $\Phi_{gravity}=gz$.

The general form of the geopotential may include the Earth's gravity, tidal forces, and self attraction and loading (SAL), and may be written as

$$
\Phi = \Phi_{gravity} + \Phi_{tides} + \Phi_{SAL} + c
$$ (def-geopotential)

where $c$ is an arbitrary constant. Therefore, the pressure gradient and geopotential gradient may be written together as

$$
\begin{aligned}
-v \nabla_z p - \nabla_z \Phi
&= -v \nabla_z p - \nabla_z \Phi_{tides} - \nabla_z \Phi_{SAL}\\
&=-v  \nabla_r p - \nabla_r \Phi_{gravity} - \nabla_r \Phi_{tides} - \nabla_r \Phi_{SAL}\\
&=-v  \nabla_r p - \nabla_r \Phi. \\
\end{aligned}
$$ (gradp-gradphi)

On the first line, note that $\nabla_z \Phi_{gravity}=\nabla_z gz=0$. For tides and SAL we assume that these forces do not vary in the vertical due to the small aspect ratio of the ocean, so that the vertical derivative in the expansion [](dvarphidnabla) is zero. This means that $\nabla_z \Phi_{tides}=\nabla_r \Phi_{tides}$ and $\nabla_z \Phi_{SAL}=\nabla_r \Phi_{SAL}$.
For versions 1.0 and 2.0 of Omega we only consider a constant gravitational force, and will not include tides and SAL.  Further details will be provided in the forthcoming pressure gradient design document.

See [Adcroft and Hallberg 2006](https://www.sciencedirect.com/science/article/pii/S1463500305000090) eqn. 1 and [Griffies et al](http://sciencedirect.com/science/article/pii/S1463500300000147) eqn 2 for additional examples of the pressure gradient in tilted coordinates. The additional terms due to the expansion of $\nabla_z$ to $\nabla_r$ in the rest of the equations are small and are ignored.

Some publications state that the transition from Boussinesq to non-Boussinesq equations is accompanied by a change from z-coordinate to pressure-coordinates. However, we use a general vertical coordinate, so the vertical may be referenced to $z$ or $p$. In a purely z-coordinate model like POP, only the $\nabla p$ term is used in [](gradp). In a purely p-coordinate model, only $\nabla z$ remains, as described in  [de Szoeke and Samelson 2002](https://journals.ametsoc.org/view/journals/phoc/32/7/1520-0485_2002_032_2194_tdbtba_2.0.co_2.xml). In a general vertical coordinate model the layer interface placement is up to the user's specification, and so both terms are kept.

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
+ v_k^{top}\omega_k^{top} {\bf u}_k^{top} - v_k^{bot}\omega_k^{bot}{\bf u}_k^{bot}
=
- v_k \nabla_r p_k - \nabla_r \Phi_k
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


## 5. Discrete Equations

The horizontally discretized layered equations are as follows. We have dropped the $r$ in $\nabla_r$ for conciseness, and the operator $\nabla$ from here on means within-layer.

$$
\frac{\partial u_{e,k}}{\partial t} + \left[ \frac{{\bf k} \cdot \nabla \times u_{e,k} +f_v}{[h_{i,k}]_v}\right]_e\left([h_{i,k}]_e u_{e,k}^{\perp}\right)
+ \left[ v_{i,k}^{top}\omega_{i,k}^{top} \right]_e u_{e,k}^{top} - \left[ v_{i,k+1}^{top}\omega_{i,k+1}^{top} \right]_e u_{e,k+1}^{top}
=
- \left[ v_{i,k} \right]_e \nabla p_{i,k} - \nabla \Phi_{i,k}
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
v_{i,k} = f_{eos}(p_{i,k},\Theta_{i,k},S_{i,k})
$$ (discrete-eos)

$$
z_{i,k}^{top} = z_{i}^{floor} + \sum_{k'=k}^{K_{max}} v_{i,k'}h_{i,k'}
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
- C_D \frac{u_{e,k}\left|u_{e,k}\right|}{[v_{i,k}h_{i,k}]_e} .
$$ (discrete-mom-bottom)

The units of specific volume times mass-thickness $vh$ are length (m), so that the full term has units of m/s$^2$.

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
+ w_k^{bot}{\bf u}_k^{bot}
- w_k^{top}{\bf u}_k^{top} =
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
|$f_{eos}$ | equation of state | -  | any | function call | may produce density or specific volume |
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
|$v_{i,k}$ | specific volume | m$^3$/kg | cell  | SpecificVolume | $v = 1/\rho$ |
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


<!--
Comments for later consideration:
## Outstanding questions
Treatment of outcropping isobars and boundary conditions. Treatment of topography/partial bottom cells.
Is there an equivalent p* coordinate? This would imply stretching/squeezing for both the top (external forcing from atmosphere and sea ice) and the bottom ($p_{surf}$ variations from forcing and barotropic thickness variations) pressure layers. Other examples of coordinate choices include:  a normalized pressure-$\sigma$ coordinates ([Huang et al, 2001](https://link.springer.com/article/10.1007/s00376-001-0001-9) is the only ocean model example I could find), where $\sigma = \frac{p - p_{surf}}{p_{bt}}$, and $p_{bt} = p_b - p_{surf}$. Another possibility is to use the coordinate $\sigma = p/p_{surf}$.

## Main References
This document is based on [de Szoeke and Samelson, 2002](https://journals.ametsoc.org/view/journals/phoc/32/7/1520-0485_2002_032_2194_tdbtba_2.0.co_2.xml), [Losch et al. 2004](https://journals.ametsoc.org/downloadpdf/view/journals/phoc/34/1/1520-0485_2004_034_0306_hsacgc_2.0.co_2.pdf), and [Adcroft's 2004 lecture notes](https://ocw.mit.edu/courses/12-950-atmospheric-and-oceanic-modeling-spring-2004/resources/lec25). Other references of interest include: [Hallberg et al., 2024](https://www.cesm.ucar.edu/sites/default/files/2024-02/2024-OMWG-B-Hallberg.pdf), [Losch et al poster](http://mitgcm.org/~mlosch/NBposter.pdf), [MOM6 code](https://github.com/mom-ocean/MOM6/blob/main/src/core/MOM_PressureForce_Montgomery.F90)... MOM6 has 2 different Pressure Gradient implementations: Finite Volume (including an analytical integration based on [Adcroft et al., 2007](https://www.sciencedirect.com/science/article/pii/S1463500308000243)) and Montgomery-potential form (based on [Hallberg, 2005](https://www.sciencedirect.com/science/article/pii/S1463500304000046)), both of which can be used in non-Boussinesq mode.

--->
