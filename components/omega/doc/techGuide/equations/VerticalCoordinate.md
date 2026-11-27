(omega-tech-vertical-coordinate)=

# Vertical Coordinate

Omega uses pseudo height, $\tilde{z} = \frac{p}{\rho_0 g}$,  as the vertical coordinate,
where $p$ is **relative pressure** (gauge pressure, i.e., absolute pressure minus the standard atmosphere of 101325 Pa).
With this definition, $\tilde{z} = 0$ at the sea surface under standard atmospheric forcing.
In practice, $\tilde{z}$ is not computed directly.
Instead, the model tracks and evolves the pseudo-thickness, $\tilde{h} = \Delta \tilde{z}$, defined as the difference between adjacent layer interfaces.
The pseudo height is essentially a normalized pressure coordinate, with the advantage that it has units of meters.

The governing equations in this coordinate were derived in the
{ref}`V1 governing equation design document <omega-design-governing-eqns-omega1>`.
The vertical coordinate is implemented in `src/ocn/VertCoord.cpp`; the options
and variables it provides are described in
{ref}`omega-user-model-vertical-coordinate`.
