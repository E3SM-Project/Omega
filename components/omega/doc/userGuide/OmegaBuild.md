(omega-user-build)=

Omega's build system is based on the CMake build tool.

## Standalone Build

When performing a standalone build, you need to execute the `cmake`
command with various CMake and Omega parameters. For example, you
can specify the `CC` parameter as the C++ compiler and enable the
`ctest` option to guide the Omega build process. Here's an example
`cmake` command:

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

Eventually, all the build configurations in the E3SM component build
should be done through the E3SM build configuration. However, as the
current Omega build is not merged into the E3SM build, users must
make the following modifications.

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

At this step, you can create any E3SM case as usual without any 
configuration for Omega.

