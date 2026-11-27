(omega-dev-aux-state)=

# Auxiliary State

The `AuxiliaryState` class groups together all of Omega [auxiliary
variables](#omega-dev-aux-vars). It is responsible for managing their lifetime,
registering them with IOStreams, and providing functions that compute them,
to be used in the tendency evaluation. It is possible
to have multiple `AuxiliaryState` instances in Omega. Every instance has a name
and is tracked in a static C++ map called `AllAuxStates`.

## Initialization

An instance of the `AuxiliaryState` class requires a [`HorzMesh`](#omega-dev-horz-mesh), so
the mesh class and all of its dependencies need to be initialized before the `AuxiliaryState` class
can be. The static method:
```c++
OMEGA::AuxiliaryState::init();
```
initializes the default `AuxiliaryState`. A pointer to it can be retrieved at any time using:
```c++
OMEGA::AuxiliaryState* DefAuxState = OMEGA::AuxiliaryState::getDefault();
```

## Creation of non-default auxiliary states

A non-default auxiliary state can be created from a string `Name`, horizontal mesh `Mesh`, halo layer `MeshHalo`, vertical coordinate `VCoord`, and a number of Tracers `NTracers`:
```c++
OMEGA::AuxiliaryState*  NewAuxState = OMEGA::AuxiliaryState::create(Name, Mesh, MeshHalo, VCoord, NTracers);
```
For conveniece, this returns a pointer to the newly created state. Given its name, a pointer to a named auxiliary state
can be obtained at any time by calling the static `get` method:
```c++
OMEGA::AuxiliaryState* NewAuxState = OMEGA::AuxiliaryState::get(Name);
```

## Computation of auxiliary variables
To compute all auxiliary variables stored in an auxiliary state `AuxState`,
given ocean state `State`, an array of tracers `TracerArray`, and time level `TimeLevel` do:
```c++
AuxState.computeAll(State, TracerArray, TimeLevel);
```

## Removal of auxiliary states
To erase a specific named auxiliary state use `erase`
```c++
OMEGA::AuxiliaryState::erase(Name);
```
To clear all instances do:
```c++
OMEGA::AuxiliaryState::clear();
```

(omega-dev-aux-vars)=

## Auxiliary variables

Auxiliary variables are variables that need to be computed at every model time
step in order to advance the model state. Similar auxiliary variables are
grouped together in auxiliary variable groups. An auxiliary variable group is
implemented as a C++ class.

### Basic class structure

Each auxiliary group class has the same basic structure. It is a light-weight
class, since it needs to be copied to the GPU to perform computations. Each
class contains:
 - a constructor
 - public arrays for the auxiliary variables
 - public member function(s) that compute the variables
 - member functions for adding metadata and registering with IO
 - static strings that contain the names used for metadata and IO

Additionally, there may be enum members that store user-configurable options.

As an example, the class `KineticAuxVars` that groups kinetic energy and
velocity divergence has the following structure:
```c++
class KineticAuxVars {
 public:
   Array2DReal KineticEnergyCell;
   Array2DReal VelocityDivCell;

   KineticAuxVars(const HorzMesh *mesh, int NVertLayers);
   void addMetaData() const;
   void defineIOFields() const;

   KOKKOS_FUNCTION void
   computeVarsOnCell(int ICell, int KChunk,
                     const Array2DReal &NormalVelEdge) const;

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DR8 EdgeSignOnCell;
   Array1DR8 DcEdge;
   Array1DR8 DvEdge;
   Array1DR8 AreaCell;

   // names used in defining MetaData and IOFields
   static const std::string KineticEnergyCellName;
   static const std::string VelocityDivCellName;
};
```

### Configuration enums

Some auxiliary variables have user-configurable options. In that case, for each
option an enum is defined in the auxiliary group header file. For example, the
flux thickness choice is represented using
```c++
enum FluxThickEdgeOption { Center, Upwind };
```

### Constructor

The constructor takes in a mesh, and does the following:
- allocates the auxiliary variable arrays
- grabs the neccessary data from the mesh
- retrieves configuration options
- calls the `addMetaData()` and `defineIOFields()` helper functions

### Compute member functions

The compute functions operate on a vertical chunk of mesh elements and take in
the element index, the chunk index, and input arrays, which may contain state
variables or other auxiliary variables. There is no output argument, they write
their results directly to the member arrays. There might be more than one
compute function, if variables located on different mesh elements are grouped
together. For example, the `VorticityAuxVars` class provides both
```c++
void computeVarsOnVertex(int IVertex, int KChunk, const Array2DReal &PseudoThickCell, const Array2DReal &NormalVelEdge);
void computeVarsOnEdge(int IEdge, int KChunk);
```
because the vorticity variables are first computed on vertices and then
averaged to edges. Note that the edge compute function doesn't take any array
arguments, because it operates only on the class member arrays.

### Implemented groups
The following auxiliary variable groups are currently implemented:
| Group | Auxiliary Variable | Available options |
| ----- | ------------------- | ------- |
| KineticAuxVars | KineticEnergyCell ||
|| VelocityDiv ||
| PseudoThicknessAuxVars | FluxPseudoThickEdge | Center or Upwind|
|| MeanPseudoThickEdge ||
| VorticityAuxVars |  RelVortVertex ||
||  NormRelVortVertex ||
||  NormPlanetVortVertex ||
||  NormRelVortEdge ||
||  NormPlanetVortEdge ||
| VelocityDel2AuxVars |  Del2Edge ||
||  Del2DivCell ||
||  Del2RelVortVertex ||
| TracerAuxVars | HTracersEdge | Center or Upwind|
|| Del2TracersCell ||
| SurfTracerRestAuxVars | TracersMonthlySurfClimoCell ||
| MixedLayerAuxVars | DenMixLayerDepth ||
|| DenMixLayerIndex ||
