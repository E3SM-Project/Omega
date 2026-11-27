(omega-tech-pressure-gradient)=

# Pressure Gradient

The pressure gradient term in the momentum equation represents the force per
unit mass due to horizontal variations in pressure and geopotential. Omega has
a centered and a finite-volume discretization of it, implemented in
`src/ocn/PGrad.cpp`, `src/ocn/PGradFiniteVolume.h` and `src/ocn/PGradRecon.h`.
How to choose between them is described in
{ref}`omega-user-model-pressure-gradient`.

## Physical background

In the layered non-Boussinesq momentum equation solved in Omega, the pressure
gradient tendency for each edge and layer includes three contributions:

1. **Montgomery potential gradient**: The horizontal gradient of the Montgomery
   potential ($\alpha p + g z$), averaged across the top and bottom interfaces of
   each layer. The Montgomery potential combines the pressure gradient and the
   geopotential, and its gradient along coordinate surfaces accounts for both the
   direct pressure force and the effect of tilted layer interfaces that arise when
   using a general vertical coordinate.

2. **Specific volume correction**: A correction term proportional to the gradient
   of specific volume (inverse density) at each edge. This term ensures that
   horizontal density variations between the two cells sharing an edge are properly
   represented in the pressure gradient force.

3. **External geopotential forcing**: Contributions from the tidal potential and
   the self-attraction and loading (SAL) terms. These represent gravitational
   forcing from astronomical tides and the deformation of the solid Earth and ocean
   surface in response to the ocean mass distribution. These terms are currently set
   to zero and will be provided by a future tidal forcing module.

(omega-tech-pressure-gradient-fv)=

## Finite-volume scheme

The finite-volume scheme compares the geopotential of the two columns sharing
an edge at a **common pressure**, rather than at a common layer index. This
gives a pressure gradient that is zero to machine precision for any resting
ocean whose temperature and salinity vary linearly with pressure, at any
coordinate tilt or layer thickness, provided the two columns' sea floors are
at the same pressure. Over a sloping sea floor a small residual remains,
because the two columns' sea-floor heights come from `VertCoord`'s sums over
different layers. It is smallest where partial bottom cells cut the floor
beneath level interior layers, and with terrain-following layers it shrinks as
the layers are refined. The centered scheme, by contrast, carries an error
that is first order in the coordinate tilt and accumulates downward through
the column.

:::{admonition} Not yet written
:class: warning
The derivation of the finite-volume scheme is still to be written here, from
the {ref}`design document <omega-design-pressure-grad-high-order>`. Progress is
tracked in
[E3SM-Project/Omega#726](https://github.com/E3SM-Project/Omega/issues/726).
:::
