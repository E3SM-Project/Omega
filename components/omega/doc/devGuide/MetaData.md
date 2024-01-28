(omega-dev-metadata)=

# MetaData Code Structure

The MetaData feature in Ometa is constructed using four distinct C++ classes.
These classes are designed to facilitate the organization and retrieval of
metadata throughout the Omega source code.

## Class Descriptions

### 1. MetaData Class
- **Purpose**: Serves as the central class, holding key-value pair metadata
               for variables and other entities such as Ometa code or
               the model itself.
- **Structure**: Contains a static member that maintains references to all
               instances of the class, allowing for universal access
               throughout the codebase.

### 2. MetaDim Class
- **Purpose**: Provides a structured representation of dimension information,
               primarily for array-type entities.
- **Structure**: Defines dimension-related attributes and functionalities to
               aid in handling array data.

### 3. MetaGroup Class
- **Purpose**: Offers a mechanism to aggregate certain MetaData instances for
               streamlined management and retrieval.
- **Structure**: Facilitates grouping and accessing grouped metadata,
               enhancing organization and accessibility.

### 4. ArrayMetaData Class
- **Purpose**: Specializes the MetaData class to handle array-type variables
               specifically, ensuring that required metadata is consistently
               provided.
- **Structure**: Inherits from MetaData and introduces additional constraints
               and functionalities pertinent to array-type data.

# Function Naming Convention
The naming of functions within these classes adheres to the LLVM coding style
and is categorized based on the nature of the operation and the type of member
variables they interact with. Function names begin with a lowercase letter.

### Specific Conventions
- **Static Member Access**: Functions interacting with static members are
               typically short verbs. For instance, `get("Name");` is utilized
               to instantiate a MetaData class.
- **Non-Static Member Access**: Functions dealing with non-static members often
               conclude with a specific descriptor indicating the targeted
               information, e.g., `getMetaData("Name", Var);` retrieves
               metadata from a MetaData instance.
- **Instance Management**: The terms `create/destroy` are designated for
               functions that instantiate or delete an instance, respectively.
- **Instance Retrieval**: The `get` function is employed to obtain an instance
               from the static member holding all instances.
- **Existence Check**: The `has("Name");` function returns a boolean value,
               indicating whether an instance named "Name" exists.

# Return Type and Return Code

The create and get static functions are designed to return a `shared_ptr` of
the class instance. This design choice simplifies the usage and enhances
code readability by allowing the use of the auto keyword, eliminating the
need for repeatedly specifying the `shared_ptr` type.

Other Functions typically returns an integer code indicating the operation
status representing a successful operation with no errors. Negative Values
correspond to specific error cases, each value indicating a distinct error
type.

# Inheriting MetaData class

The MetaData class is a versatile construct intended for attaching metadata
not only to model variables but also to other elements related to the model
and code. It serves as a foundational class for metadata management.

The ArrayMetaData class is a specialized version of the MetaData class,
focusing specifically on array-type model variables. It defines a singular
method, `create`, which is an overridden version of the one in the MetaData
class. The create function in ArrayMetaData is tailored to ensure that all
necessary metadata for array-type variables is provided upon creation.
This enforcement guarantees the integrity and completeness of metadata
for array-type model variables.
