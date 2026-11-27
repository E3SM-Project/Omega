(omega-tech-time-stepping)=

# Time Stepping

:::{admonition} Partly written
:class: warning
This page does not yet describe the forward-backward, RK2, RK4 and
split-explicit time steppers. They are to be written here from the
{doc}`time stepping </design/TimeStepping>` and
{doc}`split time stepping </design/SplitTimeStepping>` design documents. Progress is tracked in
[E3SM-Project/Omega#726](https://github.com/E3SM-Project/Omega/issues/726).
:::

The time steppers are chosen with the options in
{ref}`omega-user-config-time-stepping`, and implemented as described in
{ref}`omega-dev-time-stepping` and {ref}`omega-dev-split-time-stepping`.

## Split-explicit RK2

(omega-tech-time-stepping-split-forcing)=

### Barotropic forcing

The barotropic forcing is the thickness-weighted column mean of the
provisional full-step baroclinic velocity divided by the time step. For each
wet edge column, with $h_k$ from `MeanPseudoThickEdge`, $u^n_{bcl,k}$ from `NormalBclVelCur`,
$R_k$ from `IterVelocityTend`, and split factor $s$,

$$
F = \frac{s}{\Delta t}
\frac{\sum_k h_k (u^n_{bcl,k} + \Delta t R_k)}{\sum_k h_k}.
$$

This includes any column mean remaining in the old baroclinic velocity under
the current thickness weights.
