(omega-dev-dynamicinputstreams)=

## Dynamic Input Streams

Dynamic input streams extend the [IOStream](#omega-dev-iostreams) mechanism to
support reading fields whose names and secondary dimensions are not known at
compile time.  When a stream is configured with `DynamicFields: true` in the
input YAML, fields in that stream do not need to be pre-registered in the field
registry.  Instead, each field's dimensions and type are discovered from the
open file at read time, the necessary `Dimension` and `Field` objects are created
automatically, storage is allocated, and the data is read.

Any module requiring this functionality must include the `IOStream.h` header
file.  The new `IOStream::readAllDynamic()` call described below replaces
any per-stream `IOStream::read()` calls that would otherwise require the
calling code to know stream names.

### The `DynamicFields` flag

The flag is stored as a private `bool DynamicFields` member on each `IOStream`
instance and defaults to `false`.  It is set in `IOStream::create()` when the
stream's YAML block contains `DynamicFields: true`; if the key is absent the
stream behaves as a standard IOStream.

Dynamic streams skip the normal field-existence check in `IOStream::validate()`
because the fields will not exist until the stream is read.  The stream is
marked validated immediately so that `IOStream::validateAll()` does not abort.

### Reading a dynamic stream

When `IOStream::readStream()` processes a field from a dynamic stream, it calls
the private method `registerAndReadDynamicField()` instead of the standard
`Field::get()` + `readFieldData()` path:

```c++
Error DynErr = registerAndReadDynamicField(InFileID, FieldName);
CHECK_ERROR_ABORT(DynErr,
    "IOStream::readStream: Failed to register/read dynamic field {} "
    "in stream {}", FieldName, Name);
```

Any error from `registerAndReadDynamicField` causes an immediate abort with a
descriptive message.

`registerAndReadDynamicField()` performs the full registration pipeline for one
field:

1. Call `IO::getVarInfo()` to retrieve the field's dimension names, global
   dimension lengths, and native data type from the open file.
2. Classify each dimension.  A dimension is a mesh dimension if
   `Dimension::exists(name) && Dimension::isDistributedDim(name)` is true;
   otherwise it is a secondary dimension.
3. Require exactly one mesh dimension and at most one secondary dimension;
   return an error otherwise.
4. If a secondary dimension is present: look it up in `Dimension::AllDims`.
   If absent, create it with `Dimension::create(name, length)`.  If it already
   exists with the same length, reuse it silently.  If it exists with a
   different length, return an error.
5. Return an error if `Field::exists(fieldName)` is true (name collision).
6. Create the field with `Field::create()` using `TimeDependent=false`.
   Because the first dimension is a distributed mesh dimension,
   `Field::isDistributed()` will automatically return `true`.
7. Allocate a `HostArray1DR8` or `HostArray2DR8` and attach it to the field
   via `Field::attachFieldData()`.
8. Obtain a SCORPIO decomposition for the read.  For 1-D mesh fields a
   temporary decomposition is built from the mesh dimension's global offsets
   and destroyed after use.  For 2-D fields `getOrCreateDynamicDecomp()` is
   called to obtain a cached decomposition.
9. Read the data via `IO::readArray()`.  SCORPIO converts from the native file
   type to `R8` during the read.
10. Copy the resulting `R8` buffer into the field's Kokkos host array.

### Dynamic decomposition cache

2-D decompositions are cached in the private static member
```c++
static std::map<std::tuple<std::string, I4>, int> DynamicDecomps;
```
keyed on `(meshDimName, nSecondary)`.  All dynamic decompositions use
`IOTypeR8`.  The global offset for local mesh index `j` (0-based global index
`globalJ`) at secondary index `r` is `globalJ * nSecondary + r`, which matches
the row-major layout of `HostArray2DR8`.

Cached decompositions survive across restarts within a single process lifetime
and are freed when `IOStream::finalize()` is called.

### `IO::getVarInfo()`

The helper function
```c++
Error IO::getVarInfo(int FileID,
                     const std::string &VarName,
                     int &NVarDims,
                     std::vector<std::string> &DimNames,
                     std::vector<I4> &DimLengths,
                     IO::IODataType &NativeType);
```
wraps `PIOc_inq_varndims`, `PIOc_inq_vartype`, `PIOc_inq_vardimid`,
`PIOc_inq_dimname`, and `PIOc_inq_dimlen` to retrieve all metadata needed to
classify and register a dynamic field.  It is declared in `IO.h` and follows
the same error-return conventions as the other IO functions.

### `IOStream::readAllDynamic()`

All dynamic streams are read through a single call:
```c++
Error DynErr = IOStream::readAllDynamic(ModelClock);
CHECK_ERROR_ABORT(DynErr, "Error reading dynamic input streams");
```
The method iterates `AllStreams` and calls `readStream()` for every stream with
`DynamicFields=true`.  Users can define any number of dynamic streams in
`omega.yml` and they will all be read automatically without modifying source
code.

In `ocnInit`, this call is placed in `initOmegaModulesImpl()` between
`HorzMesh::init()` and `VertCoord::init()`.  Application code that uses
`initOmegaModules()` automatically gets this behavior.  Standalone applications
that perform their own initialization must call `readAllDynamic()` explicitly
after `HorzMesh::init()`.

### Initialization ordering

`registerAndReadDynamicField()` classifies mesh dimensions by checking
`Dimension::isDistributedDim()`, which returns `true` only for dimensions that
were registered during mesh initialization.  Dynamic stream reads must therefore
occur after `HorzMesh::init()` and before any code that depends on the
resulting fields.  The recommended sequence is:

```c++
Decomp::init();
HorzMesh::init();           // registers NCells, NEdges, NVertices
Error DynErr = IOStream::readAllDynamic(ModelClock); // reads all dynamic streams
CHECK_ERROR_ABORT(DynErr, "Error reading dynamic input streams");
// ... construct analysis operators or other consumers ...
```
