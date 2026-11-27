(omega-dev-tendencies)=

# Tendencies

The `Tendencies` class groups together all of the OMEGA [tendency terms](#omega-dev-tend-terms).
It is responsible for managing their lifetime and providing functions that compute the
tendency terms and the [auxiliary state](#omega-dev-aux-state).
It is possible to have multiple `Tendencies` instances in OMEGA. Every instance has a name
and is tracked in a static C++ map called `AllTendencies`.

## Initialization

An instance of the `Tendencies` class requires a [`HorzMesh`](#omega-dev-horz-mesh), so
the mesh class and all of its dependencies need to be initialized before the `Tendencies` class
can be. The static method:
```c++
OMEGA::Tendencies::init();
```
initializes the default `Tendencies`. A pointer to it can be retrieved at any time using:
```c++
OMEGA::Tendencies* DefTendencies = OMEGA::Tendencies::getDefault();
```
The constructor of the `Tendencies` class initializes member instances of each of the
tendency terms, which each store constant mesh information as private member variables.

## Creation of non-default tendencies

A non-default tendency group can be created with or without custom tendencies.
Without custom tendencies, it is created from a string `Name`, horizontal mesh `Mesh`,
vertical coordinate `VCoord`, vertical advection `VAdv`, pressure gradient `PGrad`,
equation of state `EqState`, vertical mixing `VMix`, number of tracers `NTracers`,
time step `TimeStep`, and a configuration `Options`:
```c++
OMEGA::Tendencies* NewTendencies = OMEGA::Tendencies::create(
    Name, Mesh, VCoord, VAdv, PGrad, EqState, VMix, NTracers, TimeStep, Options);
```
For convenience, this returns a pointer to the newly created instance.
To allow the user to provide custom tendencies, the `create` function can take two additional arguments
`CustomThicknessTend` and `CustomVelocityTend`
```c++
OMEGA::Tendencies* NewTendencies = OMEGA::Tendencies::create(
    Name, Mesh, VCoord, VAdv, PGrad, EqState, VMix, NTracers, TimeStep, Options,
    CustomThicknessTend, CustomVelocityTend);
```
The two custom tendency arguments need to be callable objects that take a Kokkos array `Tend`, ocean state `State`,
auxiliary state `AuxState`, two integers: `ThickTimeLevel` and `VelTimeLevel`, and time instant `Time`.

Given its name, a pointer to a named tendency group
can be obtained at any time by calling the static `get` method:
```c++
OMEGA::Tendencies* NewTendencies = OMEGA::Tendencies::get(Name);
```

## Computation of tendencies
To compute all tendencies for pseudo-thickness, normal velocity, and tracer equations,
given ocean state, `State`, a group of auxiliary variables, `AuxState`, an array of tracers,
`TracerArray`, thickness time level `ThickTimeLevel`, velocity time level `VelTimeLevel`,
and time instant `Time` do:
```c++
Tendencies.computeAllTendencies(State, AuxState, TracerArray, ThickTimeLevel, VelTimeLevel, Time);
```
To call only the pseudo-thickness tendency terms:
```c++
Tendencies.computePseudoThicknessTendencies(State, AuxState, ThickTimeLevel, VelTimeLevel, Time);
```
To call only the normal velocity tendency terms:
```c++
Tendencies.computeVelocityTendencies(State, AuxState, ThickTimeLevel, VelTimeLevel, Time);
```
To call only the tracer tendency terms:
Tendencies.computeTracerTendencies(State, AuxState, TracerArray, ThickTimeLevel, VelTimeLevel);

## Frazil tendencies

Frazil tendencies are calculated by the Frazil class, not the TendencyTerms.
The tendencies are integrated in the tracer tendency compute path and read in from
the frazil computation. When
`Omega.Tendencies.FrazilTendencyEnable` is true, `computeTracerTendenciesOnly`
invokes `FrazilOnCell`, which calls the default `Frazil` object and adds frazil
contributions to:

- `PseudoThicknessTend`
- `TracerTend` for `Temperature`
- `TracerTend` for `Salinity`

The frazil implementation choice and parameters are read from the
`Omega.Frazil` block. Detailed algorithm/physics notes are in
[Frazil](#omega-dev-frazil).

## Removal of tendencies
To erase a specific named tendencies instance use `erase`
```c++
OMEGA::Tendencies::erase(Name);
```
To clear all instances do:
```c++
OMEGA::Tendencies::clear();
```

(omega-dev-tend-terms)=

## Tendency terms

The tendency terms are each implemented as a separate functor. They contain
private member variables to store any constant data needed for computation.
This data is initialized during construction and comes from either an instance
of the `HorzMesh` class or a `Config` object, which are passed as arguments to
the constructor:
```c++
   auto Mesh = OMEGA::HorzMesh::getDefault();
   auto Options = OMEGA::Config::getOmegaConfig();
   OMEGA::PseudoThicknessFluxDivOnCell ThickFluxDivOnC(Mesh, Options));
```

Each functor is called like a function to compute tendency values for a range
of vertical layers at a given mesh element. So, the functors all take as
input the mesh element index and the index for the vertical chunk, along with
arrays computed by [Auxillary Variables](#omega-dev-aux-vars) classes. In each
case, the first argument is a tendency array that gets updated by the
contributon to the tendency for the given term. The functors are designed to be
used inside a Kokkos parallel loop:
```c++
   OMEGA::parallelFor(
       {mesh->NCellsOwned, NVertLayers}, KOKKOS_LAMBDA(int ICell, int KChunk) {
          ThickFluxDixOnC(ThicknessTend, ICell, KChunk, ThickFluxEdge);
       });
```

The following tendency terms for the shallow water equations are currently
implemented:
- `PseudoThicknessFluxDivOnCell`
- `PotentialVortHAdvOnEdge`
- `KEGradOnEdge`
- `SSHGradOnEdge`
- `VelocityDiffusionOnEdge`
- `VelocityHyperDiffOnEdge`
- `SfcStressForcingOnEdge`
- `BottomDragOnEdge`
- `TracerHorzAdvOnCell`
- `TracerDiffOnCell`
- `TracerHyperDiffOnCell`
- `SfcThicknessForcingOnCell`
- `SfcTracerForcingOnCell`
- `SurfaceTracerRestoringOnCell`

### See Also

Additional information on forcing (surface stress, surface mass and tracer flux forcing, and
surface tracer restoring) is detailed in [](omega-dev-forcing).
