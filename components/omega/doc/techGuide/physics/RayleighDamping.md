(omega-tech-rayleigh-damping)=

# Rayleigh Damping

Rayleigh damping adds a linear drag on the normal velocity to the momentum equation:

$$
\left. \frac{\partial u}{\partial t} \right|_{\mathrm{Rayleigh}} = -c_R u ,
$$

where $c_R$ (s$^{-1}$) is the damping coefficient and $1/c_R$ is the e-folding time scale of the damping. The damping acts on every active layer of every edge, but not on pseudo-thickness or tracers.

## Discrete form

The damping is treated implicitly, with backward Euler, as part of the
implicit vertical mixing of velocity, so it is stable for any value of
$c_R \Delta t$. The tridiagonal system solved for the new normal velocity $Y_k$ in each active layer is

$$
-G_{k-1} Y_{k-1} + \left( G_{k-1} + G_k + H_k \right) Y_k - G_k Y_{k+1} = X_k .
$$

The backward-Euler Rayleigh term adds $\Delta t \, c_R \, \tilde{h}_k$ to the diagonal-only coefficient $H_k$, where $\tilde{h}_k$ is the pseudo-thickness of layer $k$ interpolated to the edge. With bottom drag, $H_k$ becomes

$$
H_k = \tilde{h}_k \left( 1 + \Delta t \, c_R \right)
+ \delta_{k,k_{\mathrm{bot}}} \, \Delta t \, C_D \, \frac{\rho}{\rho_0} \left| \mathbf{u} \right| .
$$

The off-diagonal coefficients $G_k$ and the right-hand side $X_k = \tilde{h}_k u^n_k$ are unchanged.

With no vertical viscosity or bottom drag, every active layer therefore decays
by exactly $1/(1 + \Delta t \, c_R)$ in a time step. The damping is set up by
the `VelVertMixSetupOnEdge` functor in `src/ocn/VertMix.h`, described in
{ref}`omega-dev-rayleigh-damping`, and its options in
{ref}`omega-user-model-rayleigh-damping`. The formulation was designed in the
{ref}`Rayleigh damping design document <omega-design-rayleigh-damping>`.
