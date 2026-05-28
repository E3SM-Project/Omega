(omega-design-DynamicInputStreams)=
# Dynamic Input Streams

**Table of Contents**
1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Algorithmic Formulation](#3-algorithmic-formulation)
4. [Design](#4-design)
5. [Verification and Testing](#5-verification-and-testing)

## 1 Overview

Omega's IOStreams mechanism requires all fields to be pre-registered in the Metadata system
before they can be read. This is appropriate for model state variables whose names, types, and
dimensions are fixed at compile time. However, some analysis capabilities require reading
**weight fields** — region masks, transect edge-sign masks, and similar arrays — whose names
and secondary dimensions are not known until runtime. These fields are user-supplied in
external input files and are not part of the Omega variable set.

DynamicInputStreams extends the existing IOStreams mechanism to support reading such fields.
When a stream is marked with `dynamicFields: true`, Omega discovers each requested field's
dimensions and type directly from the input file, dynamically registers the necessary metadata
and dimensions, allocates storage, and reads the data. The resulting fields are registered in
the same global field registry as all other Omega fields, making them transparently available
to Analysis operators without any special handling.

The primary use case is the Atlantic Meridional Overturning Circulation (MOC) analysis, which
requires:
- A region mask on cells: which cells belong to each ocean basin.
- A transect edge-sign mask on edges: which edges lie on the basin's southern boundary, and
  in which direction the edge normal points relative to the transect.

The design is general, however, and supports any field indexed on cells, edges, or vertices
with an optional secondary non-mesh dimension.

## 2 Requirements

### 2.1 Requirement: Read fields not pre-registered in Metadata

A stream marked with `dynamicFields: true` must be able to read fields that have not been
pre-registered in Omega's Metadata system. Field metadata, dimensions, and types are
discovered from the input file at initialization time. This avoids the need to hard-code
field names or dimensions in Omega source code.

### 2.2 Requirement: Runtime dimension discovery

For each dynamic field, Omega must inspect the input file to determine the field's dimensions
and native data type. Secondary dimensions (e.g., NMocBasins) that are not standard
Omega mesh dimensions are registered in the global MetaDim system so that downstream output
streams can reference them.

### 2.3 Requirement: Exactly one mesh dimension per field

Each dynamic field must have exactly one mesh dimension (NCells, NEdges, or NVertices) and at
most one secondary non-mesh dimension. Fields lacking a mesh dimension (e.g., scalar or 1D
region-only arrays) are outside the scope of this design. This constraint ensures that
existing SCORPIO decompositions can be used for the distributed mesh dimension.

### 2.4 Requirement: Dimension name conflict rules

When a secondary dimension name is encountered during reading:
- If no MetaDim with that name exists: create it with the size from the file.
- If a MetaDim with that name exists and the size matches: reuse it (silent deduplication).
- If a MetaDim with that name exists but the size differs: exit with a hard error.

Input files should use unique, descriptive dimension names (e.g., `NMocBasins`) to avoid
unintended collisions between unrelated streams.

### 2.5 Requirement: Field name collision is a hard error

If a dynamic stream attempts to register a field whose name already exists in
`MetaData::AllFields` (whether from another dynamic stream or from the pre-registered Omega
field set), initialization must exit with a clear error message.

### 2.6 Requirement: Integration with global field registry

Dynamic fields are registered in `MetaData::AllFields` using the same `ArrayMetaData` and
`IOField` machinery as model state fields. Analysis operators that depend on dynamic fields
list them in `getInputFieldNames()` exactly as they would list any other model field. No
special handling is required in the Analysis orchestrator.

### 2.7 Requirement: Default to R8 storage

Integer fields (I4, I8) from the input file are promoted to R8 when stored in Omega. This
allows Analysis operators to treat weight fields uniformly with other floating-point fields.
R4 fields from the file are also stored as R8.

### 2.8 Requirement: Mesh dimension distributed, secondary dimension replicated

The mesh dimension of each dynamic field follows the standard Omega parallel decomposition:
each MPI task holds only the local portion of cells, edges, or vertices. The secondary
dimension (e.g. NMocBasins) is not distributed — every task holds all values across the
secondary dimension. SCORPIO decompositions for the mesh dimension are reused from the
existing IOEnv; 2D dynamic fields with a secondary dimension require dynamically created
decompositions (see Section 3).

### 2.9 Requirement: Re-read on every initialization

Dynamic streams are re-read on every model initialization, including restarts. Because these
fields are not written to the restart file, they must be reloaded from the original input
file each time the model starts.

### 2.10 Requirement: Initialization after mesh, before Analysis

Dynamic input streams must be read after Omega's mesh initialization (so that NCells, NEdges,
NVertices MetaDims and SCORPIO decompositions are available) and before the Analysis
orchestrator is constructed (so that all weight fields are in the registry when operators
resolve their dependencies).

### 2.12 Requirement: Automatic reading without hard-coded stream names

The model must provide a single initialization call that reads all streams marked
`DynamicFields: true` automatically, without the calling code knowing stream names.
This allows users to add new dynamic streams to `omega.yml` at runtime without
modifying any source code.

### 2.13 Desired: String name arrays (deferred)

Support for reading associated string name arrays (e.g., `regionNames(NMocBasins, StrLen)`,
`transectNames(NMocBasins, StrLen)`) is deferred to a future design iteration. When added,
name arrays would be registered as fields or as metadata attributes on the corresponding
numeric field.

## 3 Algorithmic Formulation

### 3.1 Dynamic Field Registration

The dynamic registration algorithm is invoked during `IOStream::Read` for each field in the
`contents` list when `dynamicFields` is true.

**Algorithm**: `IOStream::registerDynamicField`

**Input**: open file ID, field name string, IOEnv

**Output**: IOField registered in `MetaData::AllFields`; storage allocated; data read

1. **Query field info from file** using SCORPIO (`PIOc_inq_varid`, `PIOc_inq_varndims`,
   `PIOc_inq_vardimid`, `PIOc_inq_vartype`). Obtain dimension names and lengths, and the
   native type.

2. **Classify dimensions**. For each dimension of the field:
   - If the dimension name matches a known Omega mesh dimension (NCells, NEdges, NVertices):
     record it as the mesh dimension.
   - Otherwise: treat it as a secondary dimension.

3. **Validate structure**. Verify:
   - Exactly one mesh dimension is present.
   - At most one secondary dimension is present.
   - Exit with error if either condition is violated.

4. **Register secondary dimension** (if present):
   - Look up the dimension name in `MetaDim::AllDims`.
   - If absent: call `MetaDim::create(dimName, dimLength)`.
   - If present with matching length: reuse it.
   - If present with different length: exit with error.

5. **Create ArrayMetaData** for the field:
   - Name: the field name string from `contents`.
   - Description, units, standard name: read from file variable attributes if present;
     otherwise use empty strings.
   - Dimensions: mesh MetaDim (first) then secondary MetaDim (if present).
   - ValidMin, ValidMax, FillValue: read from file attributes if present; otherwise defaults.

6. **Check for name collision** in `MetaData::AllFields`. If the name already exists, exit
   with error.

7. **Allocate storage**. Allocate a Kokkos array of type R8 with shape
   `(nLocalMeshDim, nSecondaryDim)` (or `(nLocalMeshDim)` for 1D fields). The local mesh
   dimension extent comes from the existing Omega decomposition for that mesh location.

8. **Create IOField** combining the ArrayMetaData shared pointer and the data array pointer.
   Register in `MetaData::AllFields`.

9. **Build or reuse SCORPIO decomposition** for reading:
   - For a 1D mesh field: use the existing 1D decomposition in IOEnv for that mesh location
     and data type.
   - For a 2D field `(nMesh, nSecondary)`: check `IOEnv::dynamicDecomps` for a cached entry
     keyed on `(meshLocation, nSecondary, R8)`. If absent, create a new PIO decomposition.
     The global offsets for a task owning local cells
     $\{c_0, c_1, \ldots, c_{k-1}\}$ (0-based global indices) are:

$$
\text{offset}(c_j, r) = c_j \cdot N_{\text{secondary}} + r, \quad r \in [0, N_{\text{secondary}})
$$

     Cache the new decomposition for reuse by subsequent fields or restarts.

10. **Read data** from file using the SCORPIO decomposition. Promote from native file type to
    R8 as needed (integer or R4 → R8 conversion applied during or immediately after reading).

## 4 Design

### 4.1 Data types and parameters

#### 4.1.1 Parameters

No new global configuration parameters are introduced. The `DynamicFields` option is
per-stream and specified inside the `IOStreams` section of the Omega configuration file.

#### 4.1.2 Class and struct changes

The `IOEnv` class described in earlier drafts of this design was not created.
Instead, the dynamic decomposition cache is stored as a static member of `IOStream`,
which is simpler and sufficient for the current scope.

##### IOStream — new members

```c++
class IOStream {
   private:
      // ... existing members ...

      /// If true, fields in this stream are not required to be pre-registered;
      /// their metadata and dimensions are discovered from the input file at read time.
      bool DynamicFields;

      /// Cache for dynamically-created 2D SCORPIO decompositions keyed on
      /// (mesh dimension name, secondary dimension size). All dynamic decompositions
      /// use IOTypeR8. Freed in IOStream::finalize().
      static std::map<std::tuple<std::string, I4>, int> DynamicDecomps;

      // ... rest of existing class ...
};
```

### 4.2 Methods

#### 4.2.1 IOStream::create — parsing `DynamicFields`

`IOStream::create()` reads the optional `DynamicFields` boolean from the stream's
YAML config block. If absent, it defaults to `false`. No change to the constructor
signature was required.

#### 4.2.2 IOStream::validate — bypass for dynamic streams

```c++
bool IOStream::validate() {
   // Dynamic streams skip field existence checks: fields are registered at read time.
   if (DynamicFields) {
      Validated = true;
      return true;
   }
   // ... existing validation logic ...
}
```

#### 4.2.3 IOStream::readStream — dynamic branch

The field-processing loop in `IOStream::readStream()` branches on `DynamicFields`:

```c++
for (auto IFld = Contents.begin(); IFld != Contents.end(); ++IFld) {
   std::string FieldName = *IFld;
   if (DynamicFields) {
      Error DynErr = registerAndReadDynamicField(InFileID, FieldName);
      CHECK_ERROR_ABORT(DynErr,
          "IOStream::readStream: Failed to register/read dynamic field {} "
          "in stream {}", FieldName, Name);
   } else {
      std::shared_ptr<Field> ThisField = Field::get(FieldName);
      int FieldID;
      Err += readFieldData(ThisField, InFileID, AllDimIDs, FieldID);
   }
}
```

Any error returned by `registerAndReadDynamicField` causes an immediate abort via
`CHECK_ERROR_ABORT`.

#### 4.2.4 IOStream::registerAndReadDynamicField

New private method implementing the algorithm from Section 3.1:

```c++
/// Discovers field metadata from an open file, registers Dimension and Field,
/// allocates R8 storage, builds a SCORPIO decomposition, reads data with
/// type promotion to R8, and fills the field's Kokkos host array.
/// Returns a non-zero Error on any failure.
Error IOStream::registerAndReadDynamicField(
      int FileID,                  ///< open SCORPIO file ID
      const std::string &FieldName ///< variable name in file and Omega
);
```

#### 4.2.5 IOStream::getOrCreateDynamicDecomp

New private method managing the static `DynamicDecomps` cache:

```c++
/// Returns a cached R8 SCORPIO decomposition for a 2D dynamic field.
/// Creates and caches a new decomposition on first use.
int IOStream::getOrCreateDynamicDecomp(
      const std::string &MeshDimName, ///< name of the mesh dimension
      I4 NGlobalMesh,                 ///< global size of mesh dimension
      I4 NSecondary                   ///< size of secondary dimension
);
```

The global offset for local mesh index `j` (0-based global index `globalJ`) at
secondary index `r` is `globalJ * NSecondary + r`, matching the row-major layout
of `HostArray2DR8`.

#### 4.2.6 IOStream::readAllDynamic — automatic reading of all dynamic streams

New public static method satisfying requirement 2.12:

```c++
/// Reads every stream with DynamicFields=true. Call once after HorzMesh::init()
/// and before any code that depends on the resulting fields.
static Error IOStream::readAllDynamic(const Clock *ModelClock);
```

The method iterates `AllStreams` and calls `readStream()` for each dynamic stream.
It is called in `initOmegaModulesImpl()` between `HorzMesh::init()` and
`VertCoord::init()`.  Users can add any number of dynamic streams to `omega.yml`
without modifying source code.

#### 4.2.7 IO::getVarInfo — new SCORPIO inquiry wrapper

New function in `IO.h` / `IO.cpp`:

```c++
/// Queries a variable's dimension names, global lengths, and native data type.
Error IO::getVarInfo(
    int FileID,
    const std::string &VarName,
    int &NVarDims,
    std::vector<std::string> &DimNames,
    std::vector<I4> &DimLengths,
    IO::IODataType &NativeType
);
```

Uses `PIOc_inq_varndims`, `PIOc_inq_vartype`, `PIOc_inq_vardimid`,
`PIOc_inq_dimname`, and `PIOc_inq_dimlen`.

### 4.3 Configuration

Dynamic streams use `DynamicFields: true` (case-sensitive key, boolean value).
The `Precision` field is not meaningful for dynamic read streams and may be omitted;
a pre-existing bug in `IOStream::create()` that caused an abort when `Precision`
was absent was fixed as part of this implementation (see Section 4.6).

```yaml
Omega:
  IOStreams:
    MocMasksAndTransects:
      UsePointerFile: false
      Filename: '/path/to/oQU240_mocBasinsAndTransects.nc'
      Mode: read
      Freq: 1
      FreqUnits: OnStartup
      DynamicFields: true
      Contents:
        - MocCellMasks
        - MocEdgeSigns
```

The field names in `Contents` must exactly match the variable names in the netCDF
file. These names become the Omega-internal names used to retrieve the fields.

### 4.4 Initialization ordering

Dynamic input streams must be read after mesh initialization and before any
operators that depend on the dynamic fields.  This is accomplished automatically
in `initOmegaModulesImpl()` via `IOStream::readAllDynamic()`:

1. Machine environment and MPI setup
2. Configuration file parsing
3. Decomposition initialization (Decomp)
4. Mesh initialization (HorzMesh) — registers NCells, NEdges, NVertices dimensions
5. **`IOStream::readAllDynamic(ModelClock)`** — reads all dynamic streams, registers dynamic fields and secondary dimensions
6. Vertical coordinate initialization (VertCoord)
7. Analysis orchestrator construction — resolves field dependencies

Note: the `IOEnv` initialization step in earlier drafts of this design is not
present; no separate `IOEnv` class was created.

### 4.5 Analysis operator access

This section is unchanged from the original design.  Analysis operators access
dynamic fields by name exactly as they access model state fields.  The
AnalysisOrchestrator integration is deferred to the branch where the
AnalysisOrchestrator itself is being implemented.

### 4.6 Bug fix: `Precision` error accumulation in `IOStream::create()`

`IOStream::create()` uses a single accumulated `Err` variable for required config
fields.  The `Precision` option is optional and falls back to `"double"` when
absent, but the original code did not call `Err.reset()` after the fallback.
This caused the subsequent `CHECK_ERROR_ABORT` for `Freq` to abort even when
`Freq` was correctly present.  The fix adds `Err.reset()` to the Precision
fallback, matching the existing pattern for `UseStartEnd`.  Any stream that omits
`Precision` from its YAML config now works correctly.

## 5 Verification and Testing

Tests are implemented as CTest executables in
`components/omega/test/infra/`.  Each test writes its own synthetic netCDF
input files using SCORPIO write functions and does not require pre-committed
binary data files.

### 5.1 Test: Basic dynamic field read (DYNAMIC_IOSTREAM_TEST)

Implemented in `DynamicInputStreamTest.cpp`.  Writes `DynTest.nc` containing
`MocCellMasks(NCells, NRegions=3)` and `MocEdgeSigns(NEdges, NRegions=3)` as
`I4`, and `DynTest2.nc` containing `MocCellMasks2(NCells, NRegions=3)`.
A single call to `IOStream::readAllDynamic()` reads both streams and verifies:
- All three fields are present in `Field::AllFields`.
- `Dimension::AllDims` contains `NRegions` with global length 3.
- Data values match expected `I4→R8` promoted values for every owned cell/edge.
- The dimension count increases by exactly 1 (NRegions registered once, not twice).

Tests requirements 2.1, 2.2, 2.3, 2.7, 2.12.

### 5.2 Test: Field name collision (DYNAMIC_IOSTREAM_COLLISION_TEST — WILL_FAIL)

Implemented in `DynamicInputStreamCollisionTest.cpp`.  Two streams both list
`MocCellMasks`.  `IOStream::readAllDynamic()` reads them in alphabetical order:
CollisionStream1 succeeds, CollisionStream2 aborts via `CHECK_ERROR_ABORT`
because `registerAndReadDynamicField` detects `Field::exists("MocCellMasks") == true`.
CTest marks this test `WILL_FAIL` so the non-zero exit counts as a pass.

Tests requirement 2.5.

### 5.3 Test: Dimension name conflict — different sizes (DYNAMIC_IOSTREAM_DIM_CONFLICT_TEST — WILL_FAIL)

Implemented in `DynamicInputStreamDimConflictTest.cpp`.  `IOStream::readAllDynamic()`
reads DimStream1 (registers `NRegions=3`) then DimStream2 (file has `NRegions=5`).
The abort fires in `IOStream::readAllDims()` (which checks all registered dimensions
against the file before the dynamic field loop runs) rather than in
`registerAndReadDynamicField` step 4, but the overall behaviour is correct:
initialization aborts with a dimension-conflict error.

Tests requirement 2.4.

### 5.4 Test: Dimension deduplication — same size (DYNAMIC_IOSTREAM_TEST)

Covered in the same `DynamicInputStreamTest.cpp` test as §5.1.  The single
`IOStream::readAllDynamic()` call reads DynTest (registers `NRegions=3`) then
DynTest2 (same `NRegions=3`), and verifies the dimension count increases by
exactly 1.

Tests requirement 2.4.

### 5.5 Test: Analysis operator dependency resolution

Deferred.  The AnalysisOrchestrator is being developed on a separate branch.
When that work is merged, a test exercising end-to-end operator dependency
resolution through dynamic fields should be added.

Tests requirement 2.6.

### 5.6 Test: Restart re-read

Not implemented as a unit test.  Requirement 2.9 is satisfied by architecture:
`ocnInit()` always calls `IOStream::init()` before `IOStream::readAllDynamic()`,
and `IOStream::init()` recreates all stream objects with `OnStartup=true`.
Consequently, `readAllDynamic()` naturally re-reads all dynamic streams on every
model start, including restarts. Verification is available at the
driver/integration test level.

Tests requirement 2.9 (by architecture).

### 5.8 Test: Non-mesh field rejection (DYNAMIC_IOSTREAM_BAD_FIELD_TEST — WILL_FAIL)

Implemented in `DynamicInputStreamBadFieldTest.cpp`.  A file contains
`RegionAreas(NRegions=3)` with no mesh dimension.  The read aborts in
`registerAndReadDynamicField` when `NumMeshDims == 0` is detected.

Tests requirement 2.3.
