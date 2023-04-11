<!--- OMEGA metadata requirements and design ---------------------------------->

# OMEGA Requirements and Design:
## *Metadata*

## 1 Overview

All data files from OMEGA and E3SM must be in a self-describing format
to enable easier analysis, archiving and provenance. We describe here
requirements and design for managing metadata within the OMEGA code
itself to enable appropriate metadata in all model output. The
capabilities described here must interact with other parts of the
I/O within OMEGA, particularly the lower-level I/O interfaces that
perform the actual reading/writing to disk and the management of
files or streams and their contents.

## 2 Requirements

### 2.1 Requirement: Data types

Metadata can be in all supported OMEGA data types (I4, I8, R4, R8,
Boolean/logical, strings). Initially only scalar metadata
is required, but extending to vectors/arrays may be desireable.

### 2.2 Requirement: Global and variable metadata

Metadata refers not only to field/variable metadata but
also metadata associated with the model, code base, simulation
or data file. This capability must be able to track metadata
for these global/general metadata in addition to metadata
attached to particular fields or variables.

### 2.3 Requirement: Metadata conventions

Where possible, all metadata must conform to the climate/forecast (CF)
metadata conventions at https://cfconventions.org. 

### 2.4 Requirement: Required metadata

To comply with conventions and E3SM practices, there will be a minimum
set of required metadata and interfaces must enforce this minimum set.
This minimum will be defined later, but for variables, this would
typically include a name (short), units, long name/description,
standard CF name (if exists), and dimensions.

### 2.5 Requirement: Dimensions

For array variables, a definition of each dimension used is required.
In an unstructure model like OMEGA, the dimension describes the
global extent of the index space for each dimension that is then paired with
mesh fields for the full description of OMEGA coordinate locations.
For time dimensions, we require support for both a fixed length time
dimension (eg to support multiple time levels in a restart) as well
as an unlimited dimension to support time series.

### 2.6 Requirement: Available fields

To support other I/O components, we require a means for maintaining
a list of fields available for I/O based on the requested model
configuration. This enables the I/O streams to define contents of
a given file as part of model configuration input by supplying a
list of variable names to be included.

### 2.7 Desired: Model config??

In current MPAS-Ocean files, the full set of namelist input is also
often included as file metadata. This was thought to be useful for
provenance information. Do we still want this or make use of E3SM
provenance tools instead?  If needed, a simple
interface for replicating the model Config in metadata may be
desireable. Other provenance information may also be useful
(eg code github hash/tag, simulation info) but would require
mechanisms for making that info available within the code or input
configuration.

## 3 Algorithmic Formulation

No algorithms are used.

## 4 Design

The general philosophy in this design is to have modules that
"own" the variable define the field metadata and add it to a
list of available fields. This internally generated list ensures
that the list of available fields remains consistent with the
code base and specific model configuration. Note that this
requires most initialization to occur before processing output
streams so that the requested contents can be checked against
the list of defined fields.

Internally, the metadata class will make use of the C++ 
`std::map` that provides a useful capability for name, value
pairs of each data type.

### 4.1 Data types and parameters

#### 4.1.1 Parameters 

To enable general metadata associated with the code,
simulation or file, we define standard names for these
to be used as the "field" name during metadata definition.
For code and simulation metadata, we will use:

    static const std::string codeMeta = "code";
    static const std::string simMeta  = "simulation";

For file metadata (eg. time stamps, convention, etc.),
we will use the stream name as the field name for metadata
definition.

All defined metadata (available fields) is kept public so
that the module owning the data can define the metadata and the
I/O routines can access all defined metadata.

    static int nDefMeta; // number of fields defined
    static std::vector<Metadata> definedMeta;

Similarly, there will be a list of defined dimensions

    static int nDefDims; // number of defined dimensions
    static std::vector<MetadataDim> definedDims;

#### 4.1.2 Class/structs/data types

There are two classes associated with metadata. The first is
the actual metadata class to hold the collection of metadata
associated with a variable or general category (code, sim, file).
This class is a collection of maps (name, value pairs) for 
metadata of each supported data type. Note that there is none
associated to the Real type since the type will need to
correspond to the proper netCDF type. The Real type will
be assigned the appropriate R4 or R8 by the compiler through
the correct choice of aliased interfaces in the methods below.

    class Metadata{

       private:
           std::string name; // field name to which metadata attached
                             // use reserved names "code", "simulation"
                             // or the I/O stream name for general
                             // metadata associated with those

           I4 numI4;   // number of I4 metadata entries
           I4 numI8;   // number of I8 metadata entries
           I4 numR4;   // number of R4 metadata entries
           I4 numR8;   // number of R8 metadata entries
           I4 numBool; // number of Boolean metadata entries
           I4 numStr;  // number of string metadata entries

           // Metadata name,value pairs for each type
           std::map<std::string, I4> metaMapI4;
           std::map<std::string, I8> metaMapI8;
           std::map<std::string, R4> metaMapR4;
           std::map<std::string, R8> metaMapR8;
           std::map<std::string, bool> metaMapBool;
           std::map<std::string, std::string> metaMapStr;

           // Dimensions - only for defining fields/variables
           I4 nDims;   // number of dimensions
           std::vector<MetadataDim*> dims; // pointer to the defined
                                           // dimension for each dim

       public:
          [methods described below]
    };

A second class is needed to describe dimensions for arrays
and vectors. It is basically only a name, length pair to define
the length of each dimension.

    class MetadataDim{

       private:
          std::string name; // name of dimension
          I4 length;        // length of dimension

       public:
          [methods described below]
    };

### 4.2 Methods

#### 4.2.1 Metadata Constructors and Destructors

There will be two constructors for the metadata class. In addition
to initializing the metadata, these constructors will also add
the metadata to the list of defined fields so will also check that
the metadata has not already be defined and fail if already exists.

    // construct empty instance by name, primarily for global
    // metadata not assigned to a variable/field
    Metadata(const std::string name); // name to assign

    // construct metadata for a variable/field using
    // required metadata. Not that if stdName doesn't exist
    // an empty string can be supplied
    Metadata(const std::string name,        // name of var
             const std::string description, // long name or description
             const std::string units,       // units
             const std::string stdName,     // CF standard name
             const int         nDims,       // number of dimensions
             const std::vector<MetadataDims*> // vector of dim pointers
             );

Do we want to require any other metadata (valid min/max, fill value)?

A destructor will be provided to free metadata and remove from the
defined list, though we anticipate most metadata will need to be
persistent through a given simulation.

#### 4.2.2 Add/remove metadata

An interface for adding new metadata to a previously defined Metadata
instance will be provided. This can be used to add global metadata to
the code, simulation, file metadata, but also for variables if additional
metadata is desired. This will be an aliased function based on
the data type of the metadata.

    ierr = MetadataAdd(varName, name, value);

where varName is the name of the already defined variable, and
the other arguments are the usual name/value pair for the metadata
to be added. This function will return 0 if successful and an error
if either the varName has not been defined or if a metadata pair with
that name already exists (or maybe we want to allow overwrite?).

For symmetry, we will supply a remove function, though the
use case is likely rare.

    ierr = MetadataRemove(varName, name); 

This doesn't need a type, though it might make it easier
to search.

#### 4.2.3 Retrieve metadata

The most common use case will be retrieving metadata. For a single
metadata entry, there will be an explicit get function:

    myValue = MetadataGet(varName, name);

where varName is the field defined (or generic name for global
metadata) and name is the name associated with the metadata
to be retrieved.

During I/O stages, we will need a capability for retrieving
all metadata for a defined field or global category. We will
supply a function to retrieve the map associated with a given
data type. The `std::map` iterators and functionality can then
be used to extract all of the metadata for a given type. For
example, to retrieve all R4 metadata from a defined field:

   std::map<std::string, R4> r4meta = MetadataGetAll(varName);
   std::string metaName;
   R4 metaVal

   // use std::map functionality to loop through all metadata
   for (auto it=r4meta.begin(); it != r4.meta.end(), ++it){
      metaName = it->first;  // retrieve name part of meta pair
      metaVal  = it->second; // retrieve value part of meta pair

      // do whatever needed with metadata
   }

#### 4.2.4 Create dimension

For the dimension class, there will be a constructor that will
both construct the instance and add it to the list of defined
dimensions. It will generate an error if the dimension is already
defined.

    MetadataDim(std::string name, // name of dimension
                I4 length         // length of dimension
                );

A destructor will be supplied as well.

#### 4.2.5 Retrieve dim length

The only other function for a dimension will be a retrieval for
the dimension length:

    I4 myLength = myDim("name");

#### 4.2.6 Existence inquiry

For both metadata and dimensions, a function will be provided
to inquire whether metadata or dimensions with a given name
have already been defined:

    bool isMetadataDefined(const std::string name);
    bool isDimDefined     (const std::string name); 

## 5 Verification and Testing

Requirement 2.3 will be enforced as needed when defining available
fields so will not be tested as part of unit testing.

### 5.1 Test creation/retrieval

Create dimensions for a sample multi-dimensional field and then
create field metadata with metadata entries for all supported types.
Test by retrieving all meta data and dimension information and
ensuring it is identical.
  - tests requirement 2.1, 2.4, 2.5

### 5.2 Test global metadata

Add metadata to supported global metadata fields "code" and
"simulation" and test by retrieving the same and comparing.
  - tests requirement 2.2

### 5.3 Test inquiry functions

Test available field and available dimensions by inquiring if
above metadata and dimensions have been defined. Also test names
that don't exist to test failure modes.
  - tests requirement 2.6


