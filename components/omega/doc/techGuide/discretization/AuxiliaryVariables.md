(omega-tech-auxiliary-variables)=

# Auxiliary Variables

Quantities computed from the prognostic state at each time step and used by
the tendency terms. The variables available for output are listed in
{ref}`omega-user-io-fields-auxiliary`.

## Kinetic energy on cells

In {cite:t}`Ringler2010`, the cell-centered kinetic energy ($K_i$) is a geometry-weighted combination of the squared edge-normal velocities surrounding cell ($i$), constructed so that its discrete gradient enters the vector-invariant momentum equation as part of the Bernoulli-gradient term, ($-\nabla(K_i+\Phi_i)$). Although only one velocity component is stored at each C-grid edge, the differently oriented edges collectively represent the two-dimensional velocity magnitude; for uniform flow on an isotropic cell, the construction recovers ($K_i=\tfrac12|\mathbf{u}|^2$). Importantly, this definition is chosen for algebraic compatibility with the edge-based kinetic-energy norm and the discrete continuity equation, enabling the nondissipative momentum terms to conserve total energy to within time-discretization error rather than providing an arbitrary pointwise reconstruction of the velocity magnitude.

It should be noted that this is the irrotational kinetic energy and is not the total kinetic energy. In principle, this implies that that $K_i$ could be an underestimate of the total kinetic energy. Even for a regular hexagon in an irrotational constant flow, this formulation underestimates the total kinetic energy by 13-16\% depeneding on the flow orientation. However, it has been shown with MPAS-Ocean at standard resolution (Icos30) that $K_i$ *exceeds* $K = 0.5 (u_{cell}^2 + v_{cell}^2)$, likely due to velocity noise at the grid scale that is filtered in the course of reconstructing velocities at cell-centers. This offers a justification for employing $K_i$ in other terms of the momentum equation such as bottom drag.
