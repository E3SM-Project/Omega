(omega-design-vert-coord)=
# Vertical Coordinate

## 1 Overview
The vertical coordinate module will be responsible for computing and storing information related to the vertical mesh in Omega.
Since Omega will be a non-Boussinesq model, the vertical independent variable will be pressure- as opposed to $z$-based.
Similar to MPAS-Ocean's support for tilted height ($z^\star$) coordinates, Omega V1 will support a general ($p^\star$) vertical coordinate, where mass thickness $h$ varies in proportion to the total mass thickness of the water column. 
The prognostic variable in the layered continuity equation is mass thickness, which can be used to calculate the total pressure at a given vertical level.
The vertical height of a given level is still necessary in the model to compute sea level, the geopotential, and derivatives with respect to $z$ (vertical mixing, etc).
The height of a layer ($\Delta z$) can be found using the mass thickness, density and known bottom elevation relative to the geoid (positive up).
The vertical coordinate module will also serve as a container for information related to the variable number of vertical levels for a given ocean column due to variations in bottom elevation and ice shelf cavities.


## 2 Requirements

### 2.1 Requirement: Compute and store the total pressure at a given level

The total pressure at the bottom of each level will be computed and stored within the vertical coordinate module.
This involves a surface-down summation of the mass thickness variable times the gravitational acceleration, $g$, starting with the surface pressure.
The pressure variable will be needed in the computation of the TEOS-10 equation of state and the pressure gradient.
The variable will be registered with the `IOStreams` framework to enable it to be output during runtime.

### 2.2 Requirement: Compute and store the $z$ location of a given level

The $z$ location at the top of each level will be computed and stored within the vertical coordinate module.
In the model, $z$ is defined as being positive up from the geoid.
The $z$ location is calculated in a bottom-up summation of the mass thickness times the specific volume, starting with the known bottom elevation.
The $z$ location of the interfaces is needed to compute the geopotential gradient and the sea level.
The variable will be registered with the `IOStreams` framework to enable it to be output during runtime.

### 2.3 Requirement: Compute and store the vertical level bounds for each horizontal cell

All tendency terms require a variable vertical loop range to account for variation in bottom elevation and ice shelf cavities.
These bounds will be computed in the vertical coordinate module and used for setting masking arrays that will be used to determine where calculations are valid in the vertical.

### 2.4 Requirement: Compute and store the geopotential

The geopotential will be computed and stored within the vertical coordinate module.
Initially this will be calculated by summing together the $z$ location times $g$.
In the future, tidal potential and self attraction and loading will also have contributions to the geopotential.

### 2.5 Requirement: Compute and store the desired thickness used in the $p^\star$ vertical coordinate

The $p^\star$ vertical coordinate specifies a desired layer thickness that stretches the initial thicknesses by the bottom pressure anomaly.
This involves adding a vertically-weighted bottom pressure perturbation to the reference (initial) thickness.
This desired thickness will be computed and stored within the vertical coordinate module.
The vertical advection algorithm will use this desired thickness to compute the vertical transport.

### Desired: General Arbitrary Lagrangian Eulerian (ALE) coordinate support

In future versions of Omega, the $p^\star$ coordinate will be extended to a more general ALE coordinate

### Desired: Vertical Lagrangian remapping (VLR) support

The vertical coordinate module should retain flexibility to support VLR in future versions of Omega.

## 3 Algorithmic Formulation
The non-Boussinesq momentum equation is

$$ 
\frac{D \mathbf{u}_h}{D t } + f\boldsymbol{k}\times \mathbf{u}_h + \left(\alpha\nabla_A p + \nabla_A \Phi \right) = \boldsymbol{\mathcal{F}},
$$

where $\mathbf{u}_h$ is the horizontal velocity, $f$ is the Coriolis parameter, $\alpha = \frac{1}{\rho}$ is the specific volume, $\rho = \rho(\Theta,S,p)$ is the density, $p$ is the hydrostatic pressure, $\Phi$ is the geopotential, and $\boldsymbol{\mathcal{F}}$ are the dissipative terms.
The operator $\nabla_A$ is the gradient along a constant surface, $A$, and the total derivative is

$$ \frac{D \mathbf{u}_h}{D t } = \left( \frac{\partial}{\partial t} \right)_A  + \mathbf{u}_h\cdot \nabla_A + \omega\frac{\partial}{\partial A}, $$

where $\omega$ is the cross coordinate flow.
In the layered non-Boussinesq equations, the prognostic variable for cell $i$ and layer $k$ is the mass thickness, $h_{i,k}$, so that the geometric thickness (in meters) is a diagnostic variable defined as:

$$ \Delta z_{i,k} = \alpha_{i,k} h_{i,k}. $$

The pressure at vertical cell interfaces is the found by summing the mass thicknesses:

$$  p_{i,k+1/2} = p_{i}^{surf} + g\sum_{k^\prime=1}^k h_{i,k^\prime}. $$

and at the midpoint by

$$ p_{i,k} = p_{i}^{surf} + \sum_{k'=1}^{k-1} g h_{i,k'} + \frac{1}{2} g h_{i,k} $$

The $z$ location of cell interfaces is found by summing the mass thicknesses from the bottom multiplied by the specific volume:

$$ z_{i,k+1/2} = z_i^{floor} + \sum_{k=k^\prime}^{N} \alpha_{i,k^\prime} h_{i,k^\prime}, $$

where $z_i^{floor}$ is the (positive-up) bottom elevation.
The $z$ location of cell centers is given by: 

$$ z_{i,k} = z_i^{floor} + \frac{1}{2} \alpha_{i,k} h_{i,k} + \sum_{k^\prime= k+1}^{N} \alpha_{i,k^\prime} h_{i,k^\prime}. $$

Initially, the geopotential is the sum of the $z$ height times $g$. 
In the future, it will include contributions from the tidal potential ($\Phi_{tp}$) and self attraction and loading ($\Phi_{SAL}$):

$$ \Phi_{i,k} = \left( gz_{i,k} + \Phi_{tp} + \Phi_{SAL} \right). $$

The desired mass thickness used for the $p^\star$ coordinate is:

$$h_k^{p^\star} = h_k^{rest} + h_k^{p}, $$

where the pressure alteration, $h_k^{p}$, is given by:

$$ h_k^{p} =\frac{(p_B-p_{surf})}{g} \frac{W_kh_k^{rest}}{\sum_{k^\prime=1}^K W_{k^\prime}h_{k^\prime}^{rest}}.$$

## 4 Design

### 4.1 Data types and parameters

The `VertCoord` class will be used to compute pressure, $z$ height, and geopotential at each layer.
It will also contain the vertical loop bounds for cells, edges, and vertices.
```c++
class VertCoord {
    public:
        I4 NVertLevels;
        I4 NVertLevelsP1;

        // Variables computed 
        Array2DReal PressureInterface;
        Array2DReal PressureMid;        
        Array2DReal ZInterface;
        Array2DReal ZMid;
        Array2DReal GeopotentialMid;
        Array2DReal LayerThicknessPStar;

        // Vertical loop bounds
        Array1DI4 MinLevelCell;
        Array1DI4 MaxLevelCell;
        Array1DI4 MinLevelEdgeTop;
        Array1DI4 MaxLevelEdgeTop;
        Array1DI4 MinLevelEdgeBot;     
        Array1DI4 MaxLevelEdgeBot;
        Array1DI4 MinLevelVertexTop;
        Array1DI4 MaxLevelVertexTop;     
        Array1DI4 MinLevelVertexBot;
        Array1DI4 MaxLevelVertexBot;

        // p star coordinate variables
        Array2DReal VertCoordMovementWeights;
        Array2DReal RefLayerThickness;

        // Variables from HorzMesh
        Array2DReal BottomDepth;
    private:

        // Variables from HorzMesh
        I4 NCells;
        I4 NEdges;
        I4 NVertices;
        Array2DReal CellsOnEdge;
        Array2DReal CellsOnVertex;

        static VertCoord *DefaultVertCoord;
        static std::map<std::string, std::unique_ptr<VertCoord>> AllVertCoords;
}
```

### 4.2 Methods

```c++
class VertCoord{
    public:
        static VertCoord *create();
        void computePressure();
        void computeZHeight();
        void computeGeopotential();
        void computePStarThickness();
    private:
        void minMaxLevelEdge();
        void minMaxLevelVertex();

}
```
#### 4.2.1 Creation

The constructor will be responsible for storing any static mesh information as private variables and handling the selection any user-specified options.
```c++
VertCoord::VertCoord(const HorzMesh *Mesh,
                     int NVertLevels,
                     Config *Options);
```

The create method will take the same arguments as the constructor plus a name.
It calls the constructor to create a new vertical coordinate instance, and put it in the static map of all vertical coordinate objects.
It will return a pointer to the newly created object.
```c++
VertCoord *VertCoord::create(const std::string &Name,
                             const HorzMesh *Mesh,
                             int NVertLevels,
                             Config *Options);
```

#### 4.2.2 Initialization

The init method will create the default vertical coordinate and return an error code:
```c++
int VertCoord::init();
```

#### 4.2.3 Retrieval

There will be methods for getting the default and non-default vertical coordinate instances:
```c++
VertCoord *VertCoord::getDefault();
VertCoord *VertCoord::get(const std::string &Name);
```

#### 4.2.4 Computation

The public `computePressure` will sum the mass thicknesses times $g$ from the top layer down, starting with the surface pressure.
This will be done with a `parallel_scan` inside a `parallel_for` over the mesh cells using hierarchical parallelism.
```c++
void VertCoord::computePressure(const Array2DReal &PressureInterface,
                                const Array2DReal &PressureMid,
                                const Array2DReal &LayerThickness,
                                const Array2DReal &SurfacePressure) {

}                                
```

The public `computeZHeight` will sum the mass thicknesses times $\alpha$ from the bottom later up, starting with the bottom elevation
This will be done with a `parallel_scan` inside a `parallal_for over the mesh cells using hierarchical parallelism.
```c++
void VertCoord::computeZHeight(const Array2DReal &ZInterface,
                               const Array2DReal &ZMid,
                               const Array2DReal &LayerThickness,
                               const Array2DReal &SpecVol,
                               const Array2DReal &BottomDepth) {

}    
```

The public `computeGeopotential` will sum together the $z$ height times $g$, the tidal potential, and self attraction and loading:
```c++
void VertCoord::computeGeopotential(const Array2DReal &GeopotentialMid,
                                    const Array2DReal &ZMid,
                                    const Array2DReal &TidalPotential
                                    const Array2DReal &SelfAttractionLoading) {
                                  
} 
```
Tidal potential forcing and self attraction and loading will be default-off features.
The will be added to (or excluded from) the geopotential based on config flags.

The public `computePStarThickness` will determine the desired mass thickness used for the $p^\star$ vertical coordinate:
```c++
void VertCoord::computePStarThickness(const Array2DReal &LayerThicknessPStar,
                                      const Array2DReal &VertCoordMovementWeights,
                                      const Array2DReal &RefLayerThickness) {
                                  
} 
```

The private `minMaxLevelEdge` will determine the various vertical loop bounds on edges:
```c++
void VertCoord::minMaxLevelEdge(const Array2DI4 &MinLevelEdgeTop,
                                const Array2DI4 &MaxLevelEdgeTop,
                                const Array2DI4 &MinLevelEdgeBot,
                                const Array2DI4 &MaxLevelEdgeBot,
                                const Array2DI4 &MinLevelCell,
                                const Array2DI4 &MaxLevelCell,
                                const Array2DI4 &CellsOnEdge) {

}
```

The private `minMaxLevelVertex` will determine the various vertical loop bounds on vertices:
```c++
void VertCoord::minMaxLevelVertex(const Array2DI4 &MinLevelVertexTop,
                                  const Array2DI4 &MaxLevelVertexTop,
                                  const Array2DI4 &MinLevelVertexBot,
                                  const Array2DI4 &MaxLevelVertexBot,
                                  const Array2DI4 &MinLevelCell,
                                  const Array2DI4 &MaxLevelCell,
                                  const Array2DI4 &CellsOnVertex) {
                                  
} 
```



#### 4.2.5 Destruction and removal

No operations are needed in the destructor.
The erase method will remove a named vertical coordinate instance, whereas the clear method will remove all of
them. 
Both will call the destructor in the process.
```c++
void VertCoord::erase(const std::string &Name);
void VertCoord::clear();
```

## Verification and Testing

### Test: 
Unit tests will be used to test each of the computations (computePressure, computeZHeight, computeGeopotential, computePStarThickness) for a given mass thickness array. Comparison will be made to a ''truth'' vertical mesh that includes spatially varying `minLevelCell` and `maxLevelCell`.
