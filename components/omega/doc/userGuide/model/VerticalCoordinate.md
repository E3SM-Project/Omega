(omega-user-model-vertical-coordinate)=

# Vertical Coordinate

## Overview

Omega uses pseudo height, a normalized pressure with units of meters, as its
vertical coordinate, and tracks the pseudo-thickness of each layer; see
{ref}`omega-tech-vertical-coordinate` for the definitions.

The `VertCoord` class contains variables and functions relevant to keeping track of:
 - the maximum possible number of layers in a cell, i.e. the extent of the vertical dimension (read in from the mesh file)
 - the bottom depth of each cell (read in from the mesh file)
 - the location of active vertical layers (used to set extents of vertical loop bounds):
   - the max and min indices of active layers in each cell (read in from the mesh file)
   - the max and min indices of active layers for edges and vertices at the bottom and top of the water column (computed from min/max cell layers)
 - pressure (computed from the pseudo-thickness and surface pressure)
 - geometric height (computed from the bottom depth, specific volume, and pseudo-thickness)
 - geopotential (computed from the geometric height and tidal forcing)
 - desired vertical interface locations (computed from pressure, reference layer pseudo-thickness, and user-specified weights)

Multiple instances of the vertical coordinate class can be created and accessed by a unique name.

## Variables

| Variable Name | Description | Units |
| ------------- | ----------- | ----- |
| NVertLayers   | maximum number of vertical layers | - |
| NVertLayersP1 | maximum number of vertical layers plus 1 | - |
| PressureInterface | relative pressure (gauge pressure) at layer interfaces | Pa |
| PressureMid | relative pressure (gauge pressure) at layer midpoints | Pa |
| GeomZInterface | geometric height of layer interfaces | m |
| GeomZMid | geometric height of layer midpoint | m |
| GeopotentialMid | geopotential at layer mid points | m$^2$/s$^2$|
| PseudoThicknessTarget | desired pseudo-thickness based on total perturbation from the reference pseudo-thickness | - |
| MinLayerCell | first active layer for cell | - |
| MaxLayerCell | last active layer for cell | - |
| MinLayerEdgeTop | min of the first active layers for cells on edge | - |
| MaxLayerEdgeTop | min of the last active layer for cells on edge | - |
| MinLayerEdgeBot | max of the first active layer for cells on edge | - |
| MaxLayerEdgeBot | max of the last active layer for cells on edge | - |
| MinLayerVertexTop | min of the first active layer for cells on vertex | - |
| MaxLayerVertexTop | min of the last active layer for cells on vertex | - |
| MinLayerVertexBot | max of the first active layer for cells on vertex | - |
| MaxLayerVertexBot | max of the last active layer for cells on vertex | - |
| VertCoordMovementWeights | weights to specify how total column thickness changes are distributed across layers | - |
| RefPseudoThickness | reference pseudo-thickness used to distribute total column thickness changes | m |
| BottomGeomDepth | positive down distance from the reference geoid to the bottom | m |
| SurfacePressure | relative pressure (gauge pressure) at the top of the ocean column | Pa |

## Surface pressure

`SurfacePressure` is the relative pressure (gauge pressure, i.e., absolute pressure minus the
standard atmosphere of 101325 Pa) at the top of the ocean column. It is the top boundary
condition used when integrating layer pressures downward in `computePressure`, and is
effectively $\tilde{z}$ at the top of the water column (times a scaling factor).
Unlike the other pressure variables, which are recomputed each timestep, `SurfacePressure` is
registered in the `State` and `Restart` field groups, so it is written to restart files and read
from the initial-condition and restart files when present. The read is optional: if
`SurfacePressure` is not present in the file, it defaults to zero rather than causing an error.
For a free-surface ocean run with standard atmospheric forcing, the surface pressure is
zero (0 Pa), so it may simply be omitted from the initial-condition file.
For a coupled ocean–atmosphere run, `SurfacePressure` holds the anomaly of the atmospheric
surface pressure from the standard atmosphere; eventually it will be updated each timestep
via the coupler as a weighted sum of atmosphere, sea-ice, and land-ice pressure.
