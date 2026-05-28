(omega-user-dynamicinputstreams)=

## Dynamic Input Streams

A regular [IOStream](#omega-user-iostreams) requires every field listed in its
`Contents` to be pre-registered in Omega's field registry before the stream is
read.  Dynamic input streams lift this restriction.  When a stream is marked
with `DynamicFields: true`, Omega inspects the input file at read time,
discovers each field's dimensions and type, allocates storage, and registers the
field automatically.

The primary use case is weight fields used by analysis operators — region masks,
transect edge-sign arrays, and similar arrays whose names and secondary dimensions
are not fixed at compile time.  After being registered, dynamic fields are
indistinguishable from any other Omega field and can be accessed by name through
the standard field registry.

### Configuration

Dynamic streams are placed in the same `IOStreams` section of the Omega input
configuration file as all other streams.  The only additional option is
`DynamicFields: true`:

```yaml
Omega:
  IOStreams:
    MocMasksAndTransects:
      UsePointerFile: false
      Filename: /path/to/oQU240_mocBasinsAndTransects.nc
      Mode: read
      Freq: 1
      FreqUnits: OnStartup
      DynamicFields: true
      Contents:
        - MocCellMasks
        - MocEdgeSigns
```

The names in `Contents` must match the variable names in the netCDF file
exactly.  These names become the Omega-internal names used to retrieve the
fields after reading.

### Constraints on dynamic fields

Each field in a dynamic stream must satisfy the following constraints.
Violations abort initialization with a descriptive error message.

- **Exactly one mesh dimension**: the field must have exactly one dimension that
  is a distributed Omega mesh dimension (`NCells`, `NEdges`, or `NVertices`).
  Fields that lack a mesh dimension (e.g. a 1-D region-only array) are not
  supported.

- **At most one secondary dimension**: beyond the mesh dimension, at most one
  additional non-mesh dimension (e.g. `NMocBasins`) is allowed.

- **Unique field names**: a dynamic field name must not already exist in
  Omega's field registry, whether from another dynamic stream or from model
  state variables.

- **Consistent secondary dimension sizes**: if two streams reference a
  secondary dimension by the same name (e.g. `NMocBasins`), the dimension
  length must be identical in both files.  Use descriptive, unique dimension
  names to avoid unintended conflicts between unrelated streams.

- **Storage type**: all dynamic fields are stored as 64-bit floating-point
  (`R8`) regardless of the native type in the file.  Integer (`I4`, `I8`) and
  single-precision (`R4`) fields are promoted to `R8` when read.

- **Re-read on every initialization**: dynamic fields are not written to restart
  files and must be re-read from the original input file on every model start,
  including restarts.

### Initialization ordering

Dynamic streams are read automatically during `ocnInit` via
`IOStream::readAllDynamic()`, which is called after mesh initialization.
No explicit read call is needed in the configuration or in model code.
Adding a new `DynamicFields: true` stream to `omega.yml` is sufficient for
it to be discovered and read on the next model start.
