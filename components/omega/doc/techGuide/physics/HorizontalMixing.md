(omega-tech-horizontal-mixing)=

# Horizontal Mixing

The Laplacian and biharmonic horizontal mixing of normal velocity and tracers
use coefficients that scale with the local mesh size, so that a single set of
coefficients suits a mesh whose resolution varies. The options are described
in {ref}`omega-user-model-horizontal-mixing`.

## Scaling with mesh size

The width of a cell is the diameter of a circle with the same area, and the
reference width $W_{med}$ (`MedianCellWidth`) is the global median of the
cell widths, computed during mesh initialization.

At each edge $e$, an effective cell width is computed from the mean area of
the two adjacent cells $c_1$ and $c_2$ (`CellsOnEdge`):

$$
w_e = 2 \sqrt{\frac{A_{c_1} + A_{c_2}}{2 \pi}},
$$

where $A_c$ is `AreaCell`. The Laplacian mixing coefficient (`ViscDel2`,
`EddyDiff2`) is multiplied by

$$
s^{(2)}_e = \frac{w_e}{W_{med}}
$$

(`MeshScalingDel2`), and the biharmonic coefficient (`ViscDel4`,
`EddyDiff4`) by

$$
s^{(4)}_e = \left(s^{(2)}_e\right)^3
$$

(`MeshScalingDel4`). With `MeshScalingEnable` set to false, both factors
are one.

The factors are computed in `HorzMesh::computeMeshScaling` and
`HorzMesh::computeMedianCellWidth` in `src/ocn/HorzMesh.cpp`, and applied by
the mixing terms in `src/ocn/TendencyTerms.cpp`.
