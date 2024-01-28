(omega-user-metadata)=

# What is Metadata in Omega?
It is important that OMEGA and E3SM data files be in a self-describing format
to streamline analysis, archiving, and provenance. This involves detailing the
integration of metadata management within OMEGA's code and its interaction with
the system's I/O components for effective data handling.

Metadata classes in Omega enable the writing and reading of metadata,
encompassing not only model variables but also broader instances, such as the
code or the model itself. Please refer to the Metadata Design Document for
further details.

# How to Use MetaData in Omega Source Code

o use the Omega metadata feature, include the `MetaData.h` header file.

```
#include "MetaData.h"
```

# Create/Destroy a MetaData Instance

To add or retrieve metadata for a variable or other entity such as code or
a model, users have to create a MetaData instance for the variable or entity.

There are three ways to create a MetaData Instance. The simplest way is to use
the create method with a metadata name string:

```
auto Data1 = MetaData::create("MyVar1");

```

The create method is a static function that returns a C++ `shared_ptr` object
pointing to the MetaData instance named "MyVar1".

If the specified name already exists, a runtime error is thrown.

Once it is created successfully, users can add metadata later, as explained
in the section below.

Another way to create a MetaData Instance is to use the create function with
`std::make_pair`, consisting of name-value pairs:

```
auto Data2 = MetaData::create(
                "MyVar2",
                {
                    std::make_pair("Meta1", 1),
                    std::make_pair("Meta2", true),
                    std::make_pair("Meta3", "MyValue3"),
                }
            );
```
This allows multiple key-value pairs to be added at once during creation.

In the case of an array-type variable, it is mandatory to add a certain set
of key-value pairs. Therefore, users should use the following method to create
array-type variables:

```
    auto Data3 = ArrayMetaData::create(
                    ArrayName,
                    "Description", /// long Name or description
                    "Units",       /// units
                    "StdName",     /// CF standard Name
                    std::numeric_limits<int>::min(), /// min valid value
                    std::numeric_limits<int>::max(), /// max valid value
                    FILL_VALUE,    /// scalar used for undefined entries
                    1,             /// number of dimensions
                    Dimensions     /// dim pointers
            );
```

The `destroy` function, used with the metadata name, removes a MetaData
instance so that the instance can no longer be accessed.

```
int ret;
ret = MetaData::destroy(ArrayMeta);
```
The return value of the `destroy` function is zero on success, -1 when there
is no metadata by that name, and -2 when the destroy action fails for
other reasons.

# Add/Remove a Metadata to/from a MetaData Instance

To add metadata to a `MetaData` instance, the `addMetaData` function is used:

```
    int ret;
    const R8 AValue = 2.0;
    ret = Data1->addMetaData("NewMeta", AValue);
```

The following data types are allowed as metadata values:

- I4
- I8
- R4
- R8
- std::string
- bool

`addMetaData` returns an integer value of zero upon success and returns -1
if the metadata name already exists.

To remove metadata from a MetaData instance, the `removeMetaData` function
is used.

```
    int ret;
    ret = Data1->removeMetaData("NewMeta");
```

`removeMetaData` returns an integer value of zero upon success, -1 when there
is no metadata by that name, and -2 when the removal action fails for other
reasons.

# Retreive a Metadata from a MetaData Instance

Users can retrieve metadata using the `getMetaData` function.

```
    R8 R8Value;
    ret = Data1->getMetaData("NewMeta", R8Value);
```

The first argument of `getMetaData` is the name of the metadata. The second
argument is passed by reference, allowing `getMetaData` to place
the metadata's value into it.

Note that `getMetaData` is overloaded with several methods, each having
different data types for the second argument. It is the user's responsibility
to match the metadata name with the correct data type of the value. Otherwise,
the function may throw a type-casting error exception.

# Retreive a MetaData Instance

To retrieve an instance of MetaData for a variable or another entity, such
as code, the `get` static function can be used.

```
    auto Data4 = MetaData::get("MyVar2");
```

The `get` method is a static function that returns a C++ `shared_ptr` object
in the same way as the create function.

# Create/Destroy a MetaDim Instance

In the above example for array-type variables, users need to add metadata for
the dimension information of the variable. In Omega, the `MetaDim` class is
used for representing dimensions.

The `create` function returns a C++ `shared_ptr` object to the MetaDim
instance.

```
const I4 DimLength{3};
auto Dim1 = MetaDim::create("MyDim1", dimension);

```

In case the specified name already exists, a runtime error is thrown.

# Retreive Dimension from a MetaDim Instance

Users can retrieve dimension information using `getLength`.

```
    I4 DimLength;
    Dim1->getLength(DimLength);
```

The argument is passed by reference, allowing `getLength` to assign the length
value to it.

As of now, the return value is always zero.

# Create/Destroy a MetaGroup Instance

Because lists of contents in the model configuration file can become lengthy,
it is useful to define and maintain a metadata group. This group is essentially
a shorthand for a set of common fields (e.g., `meshFields` or `prognosticVars`).
Common groups of fields can be defined as metadata groups (MetaGroup). Please
see the metadata design document for more details.

Similarly to the MetaData and MetaDim classes explained above, the create
function returns a C++ `shared_ptr` object to the MetaGroup instance.

```
auto Group1 = MetaGroup::create("MyGroup1");

```

In case the specified name already exists, a runtime error is thrown.

# Add, Retrieve, or Remove a MetaData Instance from a MetaGroup Instance

```
int ret;
const std::string FieldName{"MyField"};

auto Data1  = MetaData::create(FieldName);

ret         = Group1->addField(FieldName);
auto Data2  = Group1->getField(FieldName):
ret         = Group1->removeField(FieldName):
```

To add, retrieve, or remove a field (a MetaData instance) from
a `MetaGroup shared_ptr`, the `addField`, `getField`, and `removeField`
methods can be used, respectively. `addField` and `removeField` return zero
upon success, -1 if the group already contains the field name, and -2 if
another error occurs. `getField` throws an error if no field name is
specified in the argument.

# Utilities

All the classes described here have the following static utility functions:

## `has` static method

This static method returns a boolean value indicating whether an instance
specified by the argument string exists.

```
bool ret1, ret2, ret3;
ret1    = MetaData::has("MyVar");
ret2    = MetaDim::has("MyDim");
ret3    = MetaGroup::has("MyGroup");
```

## `get` static method

This static method returns an instance of the class specified by
the argument string.

```
auto Data   = MetaData::get("MyVar");
auto Dim    = MetaDim::get("MyDim");
auto Group  = MetaGroup::get("MyGroup");
```
