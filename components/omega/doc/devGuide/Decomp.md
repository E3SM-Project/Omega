(omega-dev-decomp)=

## Domain Decomposition (Decomp)

In order to run across nodes in a parallel computer, OMEGA subdivides the
horizontal domain into subdomains that are distributed across the machine
and communicate using message passing via the Message Passing Interface (MPI).
To decompose the domain, we utilize the
[METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) library.
Mesh information is first read using parallel IO into a equal-spaced
linear decomposition, the partitioned by METIS into a more optimal
decomposition. The parallel METIS implementation (ParMETIS) will eventually
be used to reduce memory footpring for high-resolution configurations, but
currently we use the serial METIS library for the partitioning.

METIS requires information about the connectivity in the mesh. In particular,
it needs the total number of cells and edges in the mesh and the connectivity
given by the CellsOnCell array (described below). Once the cells have been
partitioned, we determine multiple layers of halo cells that will be needed
to avoid excessive communcation with neighbors. The CellsOnCell, EdgesOnCell
and VerticesOnCell arrays are then distributed according to this cell
partitioning. Edges and Vertices are then partitioned. Ownership of each
edge/vertex is determined by the first (valid) cell index in the CellsOnEdge
or CellsOnVertex for that edge/vertex. The remaining connectivity arrays
(CellsOnEdge, EdgesOnEdge, VerticesOnEdge, CellsOnVertex, EdgesOnVertex)
are the redistributed to match the edge and vertex distributions. Halos
are filled to ensure all necessary edge and vertex information for the
cell decomposition (and cell halos) are present in the subdomain.

The mesh is fully described by the
[MPAS Mesh Specification](https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf)
which we will reproduce here eventually.
