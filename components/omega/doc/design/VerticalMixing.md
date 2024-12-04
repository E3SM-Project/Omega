## 1. Overview

 

Representation of unresolved vertical fluxes of momentum, heat, salt, and biogeochemical tracers in ocean models is essential to simulations fidelity.  Models of turbulent fluxes spans a wide range of complexity, but the models generally fall into a few categories: simply polynomial relationships, equilibrium turbulence models, and prognostic turbulence models.  

 

## 2. Requirements

### 2.1 Requirement: Vertical mixing interface should be modular

To prepare for additional mixing closures and later improvements, the interface to vertical mixing must allow for easy connection of closures.  It is assumed that the strength of vertical mixing from each closure is additive, such that the final vertical diffusion coefficient is the sum of the coefficients and potential non local terms from each closure.

### 2.2 Requirement: Vertical mixing models should be have local and/or gradient free terms

For optimal performance, initial models of vertical turbulent fluxes for Omega should be cast as an implicit, down gradient mixing and an explicit, gradient free, component, e.g.,

$$
\overline{w' \phi'}\approx \kappa(\overline{\rho},\overline{u},\overline{v}) \left(\frac{\partial \overline{\phi}}{\partial z} + \gamma(\overline{\rho},\overline{u},\overline{v}) \right)
$$

Here, $\phi$ is a generic tracer, $\overline{w'\phi'}$ is the vertical turbulent flux of that tracer, and $\gamma$ is the gradient free portion of the flux (often referred to as 'non-local').  The vertical diffusivity can be as simple as a constant value, to complex functions of quantities like shear and stratification.  A similar equation can be written for turbulent momentum fluxes.

The use of a flux term proportional to the local gradient allows this problem to be cast implicitly, allowing for larger time steps, reducing the cost of the mixing model.  When more complex schemes are considered, other techniques such as subcycling can be considered to improve performance.

### 2.3 Requirement: Vertical mixing models cannot ingest single columns at a time

Many standard vertical mixing libraries (e.g. CVMix) receive and operate on one column at a time.  This results in poor performance, especially on accelerators.  Any utilized model of vertical mixing must ingest and operate on multiple columns concurrently.

## 3. Algorithmic Formulation

For Omega-0 a few simple vertical mixing algorithms will be used.  

### 3.1 Richardson number dependent mixing

One of the simplest class of models of vertical turbulence are polynomial functions of the gradient Richardson number.  As in MPAS-Ocean, we choose to define the Richardson number as

$$
Ri = \frac{N^2}{\left|\frac{\partial \mathbf{U}}{\partial z}\right|^2}
$$

where $N^2$ is the Brunt Vaisala Frequency, which for non Boussinesq flows is defined as 

$$ 
N^2 = \frac{g}{\rho}\frac{\partial \rho}{\partial z}
$$

This term is discretized as

$$
N^2(k) = g \frac{\ln \rho_{DD}(k) - \ln \rho(k)}{z_m(k-1) - z_m(k)}
$$

Here $\rho_{DD}$ is the density of the fluid of layer k-1 displaced to layer k adiabatically.

The shear in the denominator is discretized as 

$$
\left| \frac{\partial \mathbf{U}}{\partial z}\right|^2 = \left(\frac{U_n(k-1)-U_n(k)}{z_m(k-1) - z_m(k)}\right)^2 + \left(\frac{U_t(k-1)-U_t(k)}{z_m(k-1) - z_m(k)}\right)^2
$$

Where the subscripts *n* and *t* are the normal and tangential velocities respectively.

With the definition of the gradient richardson number, the viscosity and diffusivity are defined as

$$ 
\nu = \frac{\nu_o}{(1+\alpha Ri)^n} + \nu_b
$$

$$
\kappa = \frac{\nu}{(1+\alpha Ri)} + \kappa_b
$$

in these formulae, $\alpha$ and *n* are tunable parameters, most often chosen as 5 and 2 respectively.  $\nu_b$ and $\kappa_b$ are the background viscosity and diffusivity respectively.

### 3.2 Convective Instability Mixing

Commonly, mixing due to convective instability is treated as a step function for the diffusivity and viscosity.  This is often represented as

$$
\kappa
= 
\begin{cases}
\kappa_{conv} \quad \text{ if } N^2 \leq 0\\
0 \quad \text{ if } N^2 > 0
\end{cases}
$$

A similar expression is utilized for viscosity.  The effect of this formula is to homogenize T, S, and normal velocity, for unstable stratification, but also in neutral stratification.  The behavior in unstable stratification is likely correct, but for neutral stratification, homogenization of momentum is questionable.  To allow for a different behavior, we slightly modifty the algorithm above to, 

$$
\kappa
= 
\begin{cases}
\kappa_{conv} \quad \text{ if } N^2 < 0\\
0 \quad \text{ if } N^2 \geq 0
\end{cases}
$$

