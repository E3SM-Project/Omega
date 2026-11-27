(omega-tech-vertical-mixing)=

# Vertical Mixing

Omega computes vertical viscosity and diffusivity as the sum of three
contributions: a constant background value, convective mixing in statically
unstable water, and mixing driven by shear instability. The vertical mixing
module also calculates the gradient Richardson number and smooths it with a
1-2-1 filter before using it in the shear-instability-driven mixing. The
coefficients are computed in `src/ocn/VertMix.cpp`, and the options are
described in {ref}`omega-user-model-vertical-mixing`. Boundary-layer mixing
is described in {ref}`omega-tech-kpp`.

## Background mixing

A constant background mixing value that represents small-scale mixing processes not explicitly resolved or modeled. Typically, this is chosen to represent low values of vertical mixing happening in the ocean's stratified interior and is assumed roughly equivalent to the globally averaged interior mixing from all sources (e.g., from internal wave breaking).

## Convective mixing

Enhanced convective adjustment mixing that occurs in statically unstable regions of the water column to parameterize convection and homogenize properties. In Omega this mixing is defaulted to occur when the squared Brunt Vaisala Frequency is less than 0.0 (unstable), and is off when equal to (neutral) or greater than (stable) 0.0. The convective diffusivity is thus added to the background diffusivity ($\kappa_b$) to form the total diffusivity coefficient $\kappa$.

$$
\kappa =
\begin{cases}
\kappa_b + \kappa_{conv} \quad \text{ if } N^2 < N^2_{crit}\\
\kappa_b \quad \text{ if } N^2 \geq N^2_{crit}
\end{cases}
$$

This is different than some current implementations (i.e. in MPAS-Ocean and the CVMix package), where convective adjustment occurs both with unstable and neutral conditions ($N^2 \leq N^2_{crit}$). $\kappa_{conv}$ and $N^2_{crit}$ are constant parameters set in the `VertMix` section of the yaml file (`Diffusivity` and `TriggerBVF` under the `Convective` header).

## Shear-instability-driven mixing

Mixing induced by vertical pseudo-velocity shear, implemented using the scheme of {cite:t}`Large1994` (LMD94), through the gradient Richardson number (ratio of buoyancy to shear).

$$
\nu_{shear} = \kappa_{shear} =
\begin{cases}
\nu_o \quad \text{ if } Ri_g < 0\\
\nu_o \left[1 - \left( \frac{Ri_g}{Ri_{crit}} \right)^2 \right]^p \text{ if } 0 \leq Ri_g < Ri_{crit}\\
0.0 \quad \text{ if } Ri_{crit} \leq Ri_g
\end{cases}
$$

where $\nu_o$, $Ri_{crit}$, and $p$ are constant parameters set in the `VertMix` section of the yaml file (`BaseShearValue`, `RiCrit`, and `Exponent` under the `Shear` header). $Ri_g$ is defined as:

$$
Ri_g = \frac{N^2}{\left|\frac{\partial \mathbf{U}}{\partial z}\right|^2}\,,
$$

where $N^2$ is calculated by the EOS based on the ocean state and $\mathbf{U}$ is the magnitude of the horizontal velocity. $Ri_g$ is calculated by the vertical mixing module and then smoothed with a 1-2-1 (vertical) filter before being used to calculate the shear-instability-driven mixing. $N^2$, $\left| \frac{\partial \mathbf{U}}{\partial z}\right|^2$ and $Ri_g$ of `k` are all defined at the cell center, top interface of layer `k`. $N^2$, $\nu$ and $\kappa$ are set to zero at the surface.
