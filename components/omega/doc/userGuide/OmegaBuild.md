(omega-user-build)=

Omega's build system is built upon the CMake build tool, which provides
a robust foundation for managing the build process.

## Standalone Build

To perform a standalone build, you need to execute the cmake command
with the necessary CMake and Omega parameters. For instance, you can
specify the CC parameter as the C++ compiler and enable the ctest
option to guide the Omega build process. Here is an example of the
cmake command:

```sh
cmake \
  -DOMEGA_BUILD_TEST=ON \
  -DCMAKE_CXX_COMPILER=CC \
  ${E3SM_HOME}/components/omega
```

## E3SM Component Build

In the E3SM component build, the compilation is triggered by
the E3SM build system.

### Manual Preparation for E3SM Component Build

Although all the build configurations in the E3SM component build
should ideally be managed through the E3SM build configuration,
the current Omega build is not yet integrated into the E3SM build.
Therefore, users need to make specific modifications to enable
the Omega build within the E3SM build process.

Step 1: Modify `${E3SM_ROOT}/components/CMakeLists.txt`

Add the lines indicated as "added" in the `CMakeLists.txt` file.

```cmake
# Include function definitions
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_util.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_mpas_model.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_omega.cmake) # <= added
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_eamxx.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_model.cmake)
...
set(BUILDCONF ${CASEROOT}/Buildconf)

build_mpas_models()
build_omega() # <= added
```

Step 2: Create `${E3SM_ROOT}/components/cmake/build_omega.cmake`

Create the `cmake` file and copy the following content into it.

```cmake
function(build_omega)

  # Set CIME source path relative to components
  set(CIMESRC_PATH "../cime/src")

  add_subdirectory("omega")

endfunction(build_omega)
```

Step 3: Create an E3SM Case

At this stage, you can create any E3SM case as usual without
any specific configuration for Omega.
