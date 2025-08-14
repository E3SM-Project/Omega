(omega-user-decomp)=

## Domain Decomposition (Decomp)

Omega is designed to be run in parallel across multiple nodes and cores
in a clustered architecture. As in many Earth System Models, we decompose
the horizontal domain into subdomains that are distributed across the
system. Communication between the domains is accomplished using the
Message Passing Interface standard. To reduce the communication between
subdomains, a halo of points is filled with information from remote neighbor
cells, edges and vertices so that a time step can largely be completed
without the need to communicate.

The partitioning itself is performed using the
[METIS/ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
libraries that can partition unstructured domains in a way that also optimizes
communication. More details on the mesh, connectivity and partitioning can be
found in the [Developer's Guide](#omega-dev-decomp).

There are two parameters that are set by the user in the input configuration
file. These are:
```yaml
Decomp:
   HaloWidth: 3
   DecompMethod: MetisKWay
```
The HaloWidth is set to be able to compute all of the baroclinic terms in a
timestep without communication and for higher-order tracer advection terms,
this currently must be at least 3.

METIS and ParMETIS support a number of partitioning schemes, but the MetisKWay
and its parallel equivalent ParMetisKWay are the only currently supported
decomposition methods for Omega. The serial form (MetisKWay) is typically
faster for small meshes and decompositions. The ParMetisKWay option is a
parallel form that is faster and requires less memory than the serial option.
Note that the two methods will not result in the same partition.

An input mesh file must be provided that contains at a minimum
  - the total number of cells, edges and vertices (NCells, NEdges, NVertices)
  - the mesh connectivity contained in the arrays CellsOnCell, EdgesOnCell
    VerticesOnCell, CellsOnEdge, EdgesOnEdge, CellsOnVertex, EdgesOnVertex.
Again, a full description of the mesh is given in the
[Developer's Guide](#omega-dev-decomp).
The mesh information is read via parallel IO into an initial linear domain
decomposition and then is partitioned by METIS and rearranged into the
final METIS parallel decomposition.

Once the mesh is decomposed, all of the mesh index arrays are stored in
a Decomp named Default which can be retrieved as described in the
Developer guide. In the future, additional decompositions associated
with processor subsets (as described in MachEnv) should be possible, but
this capability is not yet supported.
