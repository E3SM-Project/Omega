(omega-user-build-standalone)=

# Standalone Build

Omega's build system uses CMake. This page describes how to build Omega
standalone, outside Polaris or within other workflows. To build it as part of
an E3SM case, see {ref}`omega-user-build-e3sm`; every CMake option is listed
in {ref}`omega-user-build-options`.

## If You Are Using Polaris

If you are already using Polaris to set up tasks or suites, you typically do
not need to build Omega manually. In that workflow, use Polaris auto-build
during `polaris setup` or `polaris suite`.

- Polaris quick start (developers):
  <https://docs.e3sm.org/polaris/main/developers_guide/quick_start.html>
- Building E3SM components:
  <https://docs.e3sm.org/polaris/main/developers_guide/quick_start.html#building-e3sm-components>
- `omega_pr` suite:
  <https://docs.e3sm.org/polaris/main/users_guide/ocean/suites.html#omega-pr-suite>

Standalone/manual and Polaris-based build workflows are both supported.

## Requirements

A standalone build creates a generic E3SM case and reads the build
configuration for the machine from it. It then configures CMake variables and
environment variables for the machine, along with any options given on the
CMake command line.

The build needs a Python interpreter and CMake 3.21 or newer.

## Building

CMake and Omega prefer an out-of-source build. This enables a user to build
and maintain multiple executables from the same source directory.
The standard practice is for a user to create a separate directory where
the build should take place and the commands below should be launched from
that directory.

Choose the machine and compiler by their names in the E3SM machine
configuration (`cime_config/machines/config_machines.xml`), and give the path
to ParMETIS, which must be built separately. For example, to build for GPUs on
Perlmutter:

```sh
>> cmake \
  -DOMEGA_CIME_COMPILER=nvidiagpu \
  -DOMEGA_CIME_MACHINE=pm-gpu \
  -DOMEGA_ARCH=CUDA \
  -DOMEGA_PARMETIS_ROOT=/path/to/parmetis \
  ${E3SM_HOME}/components/omega
```

{ref}`omega-user-build-options` describes these and the other options, such
as `OMEGA_BUILD_TEST` to build the unit tests and `OMEGA_BUILD_TYPE` for a
debug build.

Once the command completes successfully, several scripts will be created
in the build directory.

* omega\_env.sh   : load specific modules and set env. variables read from CIME
* omega\_build.sh : run `make` command after sourcing `omega_env.sh`
* omega\_run.sh   : run `./src/omega.exe` after sourcing `omega_env.sh`
* omega\_ctest.sh : run `ctest` after sourcing `omega_env.sh`

Run omega\_build.sh in the build directory to build Omega.

Once the `cmake` command succeeds, the directory where the command was
executed will contain the following files and directories.

```sh
>> ls -1 ${WORKDIR}
CMakeCache.txt
CMakeFiles
cmake_install.cmake
CTestTestfile.cmake
external
Makefile
omega_build.sh
omega_ctest.sh
omega_env.sh
omega_run.sh
src
test
```

To build the Omega, execute the `./omega_build.sh` command in the build directory.

If the build succeeds, the Omega library and executable are created in the
`src` sub-directory of the build directory.

The `./omega_ctest.sh` command runs Omega unit tests. To run the tests, MPI
parallel job launcher such as SLURM srun should be available. You may first
get allocation of an interactive computing node or use batch system.
In addition, you must either copy (or soft link) a mesh file into the same
directory as the unit test executables or specify in a configuration file
(not implemented yet) the path to that mesh input file. For detailed instructions
on how to obtain appropiate test meshes see {ref}`omega-dev-quick-start-getting-meshes`.

For example, if you are on an interactive computing node, you can run
Omega unit test by running `omega_ctest.sh` as shown below.

```sh
>> ./omega_ctest.sh
Test project <cmake working directory>
    Start 1: DATA_TYPES_TEST
1/1 Test #1: DATA_TYPES_TEST ..................   Passed    0.03 sec

100% tests passed, 0 tests failed out of 1

Total Test time (real) =   0.04 sec
```
