(omega-design)=
# Design Documents

:::{admonition} Frozen records
:class: note
Each design document records what we were thinking when a part of Omega was
designed. It is not updated as the code changes, so it may no longer match
the code. The {ref}`omega-user-guide`, {ref}`omega-dev-guide` and
{ref}`omega-tech-guide` describe the code as it is.
:::

New designs start from the {doc}`template <Template>`; see
{ref}`omega-dev-docs` for how to write one.

```{toctree}
:maxdepth: 1

OmegaV0ShallowWater
OmegaV1GoverningEqns
Analysis
Broadcast
Config
Coupling
DataTypes
Decomp
Driver
EOS
Error
FillValues
Frazil
Halo
HorzMeshClass
Logging
MachEnv
Metadata
PGrad
PGradHighOrder
IO
IOStreams
RayleighDamping
Reductions
ReorgDocs
State
StateValidation
SubmesoscaleEddies
Tendency
Tendencies
AuxiliaryVariables
AuxiliaryState
TimeMgr
Timers
TimeStepping
SplitTimeStepping
Tracers
TridiagonalSolver
KPPMix
VertAdv
VertCoord
VerticalMixingCoeff

Template
```
