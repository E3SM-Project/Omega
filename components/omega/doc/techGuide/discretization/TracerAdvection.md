(omega-tech-tracer-advection)=

# Tracer Advection

How Omega advects tracers horizontally, including the monotone flux-corrected
limiter. The options that select these schemes are described in
{ref}`omega-user-model-advection`. Both are implemented by the
`TracerHorzAdvOnCell` term in `src/ocn/TendencyTerms.h` and
`src/ocn/TendencyTerms.cpp`.

## Second order horizontal advection

The horizontal advection is done independently within each ocean layer
so a two dimensional scheme is required to advect mixing ratios.
The second order horizontal advection scheme is described by
{cite:t}`Skamarock2011`, and only a summary of a few equations from that paper
are reproduced here.

The conservative form of the scalar transport equation
within a fluid is given in equation (1) of Skamarock and Gassmann,
$$
 \frac{\partial(\rho\psi)}{\partial t} = -\nabla\cdot\mathbf{V}\rho\psi
$$
where $\rho$ is the fluid density, $\psi$ is the mixing ratio, and $\mathbf{V}$ is
the fluid velocity. This is simplified by assuming that $\rho$ is constant and cancels out.
The advection scheme for a finite volume method involves approximating this equation along an element
edge at a certain time step and multiplying by the time step length to get an approximate scalar transport
from element to element during the time step.
The left hand side of this equation is a time derivative and for Omega this is handled by standard Runge—Kutta methods and
not discussed here but can be found in the Skamarock and Gassmann paper.
The right hand side of this equation is the spatial derivative that needs to be evaluated by the Omega horizontal advection scheme.
Omega ocean uses unstructured Voronoi meshes and a
finite volume scheme with $\psi$ defined on cell centers. A first order approximation
to the spatial derivative of the left hand termacross a single edge of an element would be the difference of $\psi$ at
element centroids divided by the distance from cell centroid to cell centroid.  This is expressed in equation (5) of the paper.
For edge $i$ between elements $i+1/2$ and $i-1/2$,
$$
 \frac{\partial(u\psi_i)}{\partial x} = \frac{1}{\Delta x}[F_{i+1/2}(u\psi) - F_{i-1/2}(u\psi)] + O(\Delta x^2).
$$
Where $u$ is $\mathbf{V}\cdot\mathbf{n}$ where $\mathbf{n}$ is the unit normal along the edge being evaluated and $F$ is the evaluation
of $\psi$ for the elements on each side of the edge being evaluated.
This is used when the TracerHorzAdvOnCell user option is active and HorzTracerFluxOrder is 2.
To make the comparison to higher order methods explicit, this first order method is equivalent to defining $\psi$ as a linear function,
$\psi = c_0 + c_x x + c_y y$, between the two elements sharing the edge $i$ and taking the directed derivative along the edge.


For higher order flux calculations, higher order derivatives of $\psi$ are needed and for the user option of
TracerHorzAdvOnCell along with HorzTracerFluxOrder set to 3, a quadratic approximation of $\psi$ is used,
$$
\psi = c_0 + c_x x + c_y y + c_{xx} x^2 + c_{xy} xy + c_{yy} y^2,
$$
defined for each edge in the mesh. The six coefficients in this quadratic equation are calculated from
a least squares fit of $\psi$ in a neighborhood of elements around an edge. The
least squares fit is defined in equation (12) of the paper,
$$
 \mathbf{f = (P^T P)^{-1}P^T s = Bs}
$$
where the vector $\mathbf{f} = [c_0,c_x,c_y,c_{xx},c_{xy},c_{yy}]$ and $\mathbf{P}$ is an $m\times 6$ matrix with each row,
$(1,x_i,y_i,x_i^2,x_iy_i,y_i^2)$, based on the cell center coordinates where there is one row for each of the $m$ elements in the neighborhood.
This results in a $6\times m$ matrix $\mathbf{B}$. The values of $\mathbf{s}$ are the mixing ratios,
$\mathbf{s}=[\psi_0, \psi_1,...,\psi_m]$, of the elements in the neighborhood at the time step being evaluated.
Note that $\mathbf{B}$ depends on the mesh geometry only, so only needs to
be computed once for each edge during initialization. Since only the second order terms are used to determine the
second order derivatives needed for the higher order flux corrections to the lower order linear flux given above
and the derivatives are taken along a line connecting
cell centers, the amount of data that needs to be stored is minimized and $\mathbf{B}$ reduces to a vector.

The second order derivatives of $\psi$ are then combined as shown in equation (11) of Skamarock and Gassmann to get
a second order horizontal transport algorithm.

The computation of these second order fluxes is complicated by having to do these computations on a sphere.
Figure 2 from the Skamarock and Gassmann paper describs the slight modifications needed to compute the entries
of the $\mathbf{P}$ matrix for a spherical geometry where the distances need to be measured as arc lengths along
the surface.

In practice, the convergence of this higher order algorithm on a sphere is slightly below
second order. For the advection of a $\psi$ in the shape of a cosine bell
being advected around a sphere and compared with the analytic solution, the spatial convergence for grid refinement
is about order 1.7 as shown in  {numref}`tracer-higher-order-convergence`:

```{figure} ../images/higher_order_tracer_convergence_on_sphere.jpeg
:name:  tracer-higher-order-convergence
:align: center
:width: 600 px
Tracer higer order convergence example of a cosine bell advected on a sphere showing an order 1.71 convergence rate
```

## Monotonic flux-corrected transport

It is well known that higher order (order 2 and above) scheme for numerically integrating fluxes suffer
from dispersive "ripples" in the results particularly near steep gradients.  Lower order schemes
produce no ripples but suffer from excessive numerical diffusion. Flux-corrected transport (FCT) is a
technique which embodies the best of both schemes.

The implementation of a monotonic flux-corrected transport (FCT) in Omega follows a standard algorithm as given in
{cite:t}`Zalesak1979`.

The procedure is as follows given the notation above:
1. Compute $F^L_{i+1/2}(u\psi)$, the transportive flux by some low order scheme guarenteed to give
monotonic (ripple-free) results.
2. Compute  $F^H_{i+1/2}(u\psi)$, the transportative flux by some high order scheme.
3. Define the "antififfusive flux":
$$
A_{i+1/2}(u\psi) = F^H_{i+1/2}(u\psi) - F^L_{i+1/2}(u\psi)
$$
4. Compute the updated low order ("transported and diffused") "td" solution:
$$
 w^{td}_i = w^n_i - \frac{1}{\Delta x}[F_{i+1/2}(u\psi) - F_{i-1/2}(u\psi)]
$$
5. Limit he $A_{i+1/2}(u\psi)$ in a manner such that $w^{n+1}$ as computed below is free of extrema
not found in $w^{td}$ = $w^n$;
$$
A^C_{i+1/2}(u\psi) =  C_{i+1/2}(u\psi) A_{i+1/2}(u\psi), \quad  \quad   0 \leq C_{i+1/2} \leq 1
$$
6. Finally apply hte limited antidiffusive fluxes:
$$
 w^{n+1}_i = w^{td}_i - \frac{1}{\Delta x}[A^C_{i+1/2}(u\psi) - A^C_{i-1/2}(u\psi)]
$$


```{figure} ../images/higher_order_tracer_convergence_on_sphere_FCT.jpeg
:name:  tracer-higher-order-convergence_FCT
:align: center
:width: 600 px
Tracer higer order convergence example of a cosine bell advected on a sphere with FCT showing an order 2.23 convergence rate
```


### Example: limiting the ripple effect

The effect of monotonic higher order FCT horizontal transport verses just higher order horizontal transport
can be seen in the Polaris test of the slotted cylinder. The cylinder is advected for one revolution under
a constant velocity field. The result is advecting the solution back to the original configuration as
an exact solution. The difference between the initial condition and the advected solution is the error
in the advection scheme. This difference is shown in the following two plots for non-monotonic and monotonic
advection.

```{figure} ../images/higher_order_slotted_cylinder_non_monotonic_diff.jpg
:name:  tracer-error-in-higher-order-convergence_non-monotonic
:align: center
:width: 600 px
Error in a tracer higher order convergence example of a slotted cylinder advected on a sphere with the standard third order advection
algorithm with no FCT flux correction. Notice the ringing due to the sharp edges of the tracer distribution.
```

```{figure} ../images/higher_order_slotted_cylinder_FCT_diff.jpg
:name:  tracer-error-in-higher-order-convergence_FCT
:align: center
:width: 600 px
Error Tracer higher order convergence example of a slotted cylinder advected on a sphere with FCT monotone
flux correction. Notice the reduction in ringing around the sharp edges of the tracer distribution.
```
