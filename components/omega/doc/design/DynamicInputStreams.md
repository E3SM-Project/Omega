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

### 2.11 Desired: String name arrays (deferred)

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

No new global configuration parameters are introduced. The `dynamicFields` option is
per-stream and specified inside the `IOStreams` section of the Omega configuration file.

#### 4.1.2 Class and struct changes

##### IOStream — new member

```c++
class IOStream {
   private:
      // ... existing members (name, filename, mode, precision, sAlarm, ...) ...

      /// If true, fields in this stream are not required to be pre-registered;
      /// their metadata and dimensions are discovered from the input file at read time.
      bool dynamicFields = false;

      // ... rest of existing class ...
};
```

##### IOEnv — dynamic decomposition cache

```c++
class IOEnv {
   private:
      // ... existing decompositions (decompCell1DR8, etc.) ...

      /// Cache for dynamically-created 2D decompositions keyed on
      /// (mesh location, secondary dimension size, data type).
      /// Created on demand during dynamic stream reads; reused for restarts.
      std::map<std::tuple<MeshLocation, I4, IODataType>, int*> dynamicDecomps;

      // ... existing friend declaration ...
   friend class IOStreams;
};
```

### 4.2 Methods

#### 4.2.1 IOStream constructor — `dynamicFields` parameter

The existing IOStream constructor is extended with a `dynamicFields` boolean argument
(default `false`). When parsing the `IOStreams:` configuration section, the presence of
`dynamicFields: true` in a stream's YAML block sets this flag.

```c++
IOStream(int&  streamID,
         const std::string name,
         const std::string filename,
         const IOmode      mode,
         const IOPrecision precision,
         const IOIfExists  ifExists,
         const std::string freqUnits,
         const int         freq,
         const std::string pointerFile,
         const std::string startDate,
         const std::string endDate,
         const bool        dynamicFields = false   // <-- new parameter
         );
```

#### 4.2.2 IOStream::Read — dynamic branch

The existing `IOStreams::Read(streamName)` method gains a branch at the start of field
processing:

```c++
int IOStreams::Read(const std::string streamName) {
   // ... locate stream, open file ... (existing logic)

   for (const auto& fieldName : stream->contents) {
      if (stream->dynamicFields) {
         // New path: discover and register field from file, then read data
         Err = stream->registerAndReadDynamicField(fileID, fieldName, *IOEnvPtr);
      } else {
         // Existing path: look up pre-registered IOField and read data
         Err = IORead(fileID, IOEnvPtr->getField(fieldName));
      }
      if (Err != 0) return Err;
   }

   // ... close file, reset alarm ... (existing logic)
}
```

#### 4.2.3 IOStream::registerAndReadDynamicField

New private method implementing the algorithm from Section 3.1:

```c++
/// Discovers field metadata from an open file, registers MetaDim and ArrayMetaData,
/// allocates R8 storage, registers in MetaData::AllFields, builds a SCORPIO
/// decomposition if needed, reads data, and promotes to R8.
/// Returns 0 on success, non-zero error code on failure.
int IOStream::registerAndReadDynamicField(
      const int         fileID,    ///< open file ID from SCORPIO
      const std::string fieldName, ///< variable name in file and Omega
      IOEnv&            ioEnv      ///< I/O environment (decompositions)
);
```

#### 4.2.4 IOEnv::getOrCreateDynamicDecomp

New method on IOEnv for obtaining a 2D decomposition for a dynamic field:

```c++
/// Returns a SCORPIO decomposition descriptor for a 2D field on the given mesh
/// location with the given secondary dimension size. Creates and caches a new
/// decomposition if one does not already exist for this (location, nSecondary) pair.
int* IOEnv::getOrCreateDynamicDecomp(
      const MeshLocation location,    ///< NCells, NEdges, or NVertices
      const I4           nSecondary,  ///< size of secondary (non-mesh) dimension
      const IODataType   dataType     ///< data type (typically R8)
);
```

### 4.3 Configuration

Dynamic streams are configured inside the existing `IOStreams:` section of the Omega YAML
input file. They are distinguished only by the `dynamicFields: true` flag. All other stream
options (filename, freqUnits, freq, etc.) follow the standard IOStream conventions.

```yaml
IOStreams:

   mocMasksAndTransects:
      mode: read
      filename: '/path/to/oQU240_mocBasinsAndTransects.nc'
      freqUnits: initial
      freq: 1
      dynamicFields: true
      contents:
      - MocCellMasks        # (NCells, NMocBasins) in file
      - MocEdgeSigns        # (NEdges, NMocBasins) in file — implicitly paired
                            # with regionCellMasks via shared dim NMocBasins
```

The field names in `contents` must exactly match the variable names in the netCDF file. These
names become the Omega-internal names used by Analysis operators.

### 4.4 Initialization ordering

Dynamic input streams must be read in a dedicated phase:

1. Machine environment and MPI setup
2. Configuration file parsing
3. Decomposition initialization (Decomp)
4. Mesh initialization (HorzMesh) — registers NCells, NEdges, NVertices MetaDims
5. I/O environment initialization (IOEnv) — registers SCORPIO decompositions
6. **Dynamic input stream reading** — registers dynamic fields and secondary MetaDims
7. Vertical coordinate initialization (VertCoord)
8. Analysis orchestrator construction (AnalysisOrchestrator) — resolves field dependencies

### 4.5 Analysis operator access

Analysis operators access dynamic fields the same way they access any other Omega field:
by declaring the field name in `getInputFieldNames()`. The AnalysisOrchestrator resolves
dependencies from `MetaData::AllFields` regardless of whether a field was pre-registered
at compile time or registered dynamically.

**Example operator declaring a dependency on dynamic fields:**

```c++
class MOCOperator : public AnalysisOperator {
 public:
   const std::vector<std::string> getInputFieldNames() override {
      // These fields will be resolved from MetaData::AllFields.
      // They may come from dynamic streams or from the model state.
      return {"NormalVelocity",
              "PseudoThickness",
              "MocCellMasks",   // dynamic field: (NCells, NMocBasins)
              "MocEdgeSigns"};  // dynamic field: (NEdges, NMocBasins)
   }
   // ...
};
```

No changes to the AnalysisOrchestrator are required; it already resolves all inputs through
the same registry lookup.

## 5 Verification and Testing

### 5.1 Test: Basic dynamic field read

Create a small synthetic netCDF file containing `MocCellMasks (NCells, NRegions=3)` as
integer data and `MocEdgeSigns (NEdges, NRegions=3)` as integer data. Configure a
dynamic IOStream pointing to this file. Call `IOStreams::Read`. Verify:
- Both fields are present in `MetaData::AllFields`.
- Each field's ArrayMetaData has the correct dimensions (NCells/NEdges + NRegions).
- `MetaDim::AllDims` contains `NRegions` with length 3.
- Field data values match the file contents after I4→R8 promotion.

Tests requirements 2.1, 2.2, 2.3, 2.7, 2.9.

### 5.2 Test: Field name collision

Configure two dynamic streams each listing a field named `MocCellMasks`. Verify that
initialization exits with a non-zero error code and a descriptive error message.

Tests requirement 2.5.

### 5.3 Test: Dimension name conflict — different sizes

Configure two dynamic streams whose files both define a dimension named `NRegions` but with
different sizes (3 vs 5). Verify that reading the second stream exits with an error.

Tests requirement 2.4.

### 5.4 Test: Dimension deduplication — same size

Configure two dynamic streams whose files both define a dimension named `NRegions` with the
same size (3). Verify that both streams read successfully and that `MetaDim::AllDims` contains
exactly one `NRegions` entry.

Tests requirement 2.4.

### 5.5 Test: Analysis operator dependency resolution

Construct a minimal AnalysisOrchestrator with a mock MOC-style operator that lists
`MocCellMasks` and `MocEdgeSigns` in `getInputFieldNames()`. Read dynamic streams
before constructing the orchestrator. Verify that:
- Dependency resolution succeeds (no "field not found" errors).
- The operator's `compute()` call receives data arrays with the correct values.

Tests requirement 2.6.

### 5.6 Test: Restart re-read

Initialize the model, read dynamic streams, then simulate a restart by re-running
initialization. Verify that dynamic fields are re-read and their values remain consistent
across both initializations.

Tests requirement 2.9.

### 5.8 Test: Non-mesh field rejection

Configure a dynamic stream with a field that has no mesh dimension (e.g., a 1D array of
length NRegions only). Verify that initialization exits with an appropriate error.

Tests requirement 2.3.
