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
Without custom tendencies, it is created from a string `Name`, horizontal mesh `Mesh`, number of
vertical levels `NVertLevels`, and a configuration `Options`:
```c++
OMEGA::Tendencies*  NewTendencies = OMEGA::Tendencies::create(Name, Mesh, NVertLevels, Options);
```
For convenience, this returns a pointer to the newly created instance.
To allow the user to provide custom tendencies, the `create` function can take two additional arguments
`CustomThicknessTend` and `CustomVelocityTend`
```c++
OMEGA::Tendencies*  NewTendencies = OMEGA::Tendencies::create(Name, Mesh, NVertLevels, Options, CustomThicknessTend, CustomVelocityTend);
```
The two custom tendency arguments need to be callable objects that take a Kokkos array `Tend`, ocean state `State`,
auxiliary state `AuxState`, two integers: `ThickTimeLevel` and `VelTimeLevel`, and time instant `Time`.

Given its name, a pointer to a named tendency group
can be obtained at any time by calling the static `get` method:
```c++
OMEGA::Tendencies* NewTendencies = OMEGA::Tendencies::get(Name);
```

## Computation of tendencies
To compute all tendencies for both layer thickness and normal velocity equations,
given ocean state, `State`, a group of auxiliary variables, `AuxState`, thickness
time level 'ThickTimeLevel', velocity time level `VelTimeLevel`, and time instant `Time` do:
```c++
Tendencies.computeAllTendencies(State, AuxState, ThickTimeLevel, VelTimeLevel, Time);
```
To call only the layer thickness tendency terms:
```c++
Tendencies.computeThicknessTendencies(State, AuxState, ThickTimeLevel, VelTimeLevel, Time);
```
To call only the normal velocity tendency terms:
```c++
Tendencies.computeVelocityTendencies(State, AuxState, ThickTimeLevel, VelTimeLevel, Time);
```

## Removal of tendencies
To erase a specific named tendencies instance use `erase`
```c++
OMEGA::Tendencies::erase(Name);
```
To clear all instances do:
```c++
OMEGA::Tendencies::clear();
```
