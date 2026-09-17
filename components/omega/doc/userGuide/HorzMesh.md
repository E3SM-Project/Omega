(omega-user-horz-mesh)=

## Horizontal Mesh

Omega uses the {ref}`MPAS mesh specification <mpas-mesh-specification>`. The
names of the mesh variables have been retained, with the caveat that they now
begin with a capital letter (e.g. `angleEdge` becomes `AngleEdge`).

The Mesh class is meant to be a container for all mesh variables local to a
decomposed sub-domain that can be easily passed among the dycore routines. It
depends on a given [Decomp](#omega-user-decomp) and reproduces the
cell/edge/vertex totals and connectivity information from that class. The Mesh
class also creates parallel I/O decompositions that are used to read in the
additional mesh variables, which are not required for Decomp.

Currently, the Mesh class reads in every field required by the
{ref}`MPAS mesh specification <mpas-mesh-specification>` except those read by
Decomp. The input contents are defined by the
`HorzMeshIn` FieldGroup. Refer to the specification for the
definition, units, and dimensions of each field.

These are all read using the input [IOStream](#omega-user-iostreams) HorzMeshIn
```yaml
  IOStreams:
    HorzMeshIn:
      UsePointerFile: false
      Filename: OmegaMesh.nc
      Mode: read
      Precision: double
      Freq: 1
      FreqUnits: OnStartup
      UseStartEnd: false
      Contents:
        - HorzMeshIn
```
Only the Filename should be changed by the user to point to the relevant input
mesh file. The mesh Filename is sometimes overridden by the driver routine in
the case of unit tests using an optional argument to the
[decomposition](#omega-dev-decomp).

In addition, some mesh metadata are read from the file and stored as mesh
variables. These include ``OnSphere`` (or ``on_a_sphere`` for backcompatibility),
``IsPeriodic`` (``is_periodic``), ``SphereRadius`` (``sphere_radius``),
``XPeriod`` (``x_period``), and ``YPeriod`` (``y_period``). The two flags
OnSphere and IsPeriodic are stored as YES/NO strings in the metadata but as
boolean flags in the code.

In the future, the Mesh class will optionally compute internally those fields
the specification marks as *could be computed internally* — the areas, lengths,
angles, and weights needed for the TRiSK discretization.
