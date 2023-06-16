(omega-dev-cmake-build)=

# Omega Build with CMake

The Omega build system is based on CMake, a popular build tool.

![CMake-based Omega build process](../_static/cmakebuild.png)

The overall build process is specified in the CMakeLists.txt file
located at the top-level directory of Omega. The build process
consists of four consecutive steps.

## Step 1: Setup

In this step, the build-controlling CMake variables are configured.
The Omega build system supports two modes of builds: standalone and
E3SM Omega-component. During this step, the build mode is detected
and any mode-specific differences are resolved, ensuring a consistent
build process.

## Step 2: Update

This step configures CMake and other external variables based on
the settings defined in the previous Setup step. The integrity of
the build setup is verified at the end of this step.

## Step 3: Build

During this step, the actual build process takes place. It includes
building external libraries such as YAKL and utilizing the Omega
source files for the build. Additionally, an optional test build
can be performed.

## Step 4: Output

This final step is optional and involves copying a subset of build
artifacts to designated locations or generating dynamic outputs
as required.

Note: Until the Omega build is merged into the E3SM build, users
need to follow specific modifications to trigger the Omega build
from within the E3SM build process.


