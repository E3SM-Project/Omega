(omega-user-build-options)=

# Build Options

These CMake variables control a standalone build. Pass each one to `cmake`
as `-D<NAME>=<value>`; see {ref}`omega-user-build-standalone` for the full
command. In an E3SM case, CIME sets them from the case instead
({ref}`omega-user-build-e3sm`). The {ref}`omega-dev-cmake-build` page lists
every variable the build system uses, including internal ones.

## Machine and compilers

| Variable | Default | Meaning |
| -------- | ------- | ------- |
| `OMEGA_CIME_MACHINE` | CIME's default for the host | E3SM machine name, as defined in `cime_config/machines/config_machines.xml` |
| `OMEGA_CIME_COMPILER` | CIME's default for the machine | E3SM compiler name for that machine, as defined in the same file |
| `OMEGA_CIME_PROJECT` | CIME's default for the machine | Account to give CIME during the build, on machines that do not define a default `PROJECT` |
| `OMEGA_CXX_COMPILER` | from `OMEGA_CIME_COMPILER` | Path to a C++ compiler; overrides `OMEGA_CIME_COMPILER` |
| `OMEGA_ARCH` | detected | Programming model: `CUDA`, `HIP`, `SYCL`, `OPENMP` or `SERIAL`. If not given, it is chosen from the machine and compiler; `OPENMP` when the machine configuration builds threaded. Set `CUDA` or `HIP` to build for GPUs. |

## Libraries

METIS, ParMETIS and, optionally, GKlib must be built separately, as static
libraries; the other external libraries are built from the E3SM source tree.

| Variable | Default | Meaning |
| -------- | ------- | ------- |
| `OMEGA_PARMETIS_ROOT` | none (required) | ParMETIS installation directory |
| `OMEGA_METIS_ROOT` | `OMEGA_PARMETIS_ROOT` | METIS installation directory |
| `OMEGA_GKLIB_ROOT` | `OMEGA_METIS_ROOT` | GKlib installation directory |

## Build type and installation

| Variable | Default | Meaning |
| -------- | ------- | ------- |
| `OMEGA_BUILD_TYPE` | `Release` | `Release` or `Debug`. A debug build gives a stack trace for critical errors ({ref}`omega-user-run-errors`). |
| `OMEGA_BUILD_TEST` | off | Build the CTest unit tests |
| `OMEGA_INSTALL_PREFIX` | none | Directory that `make install` copies the library and executable into, under `lib` and `bin` |

## Precision and performance

| Variable | Default | Meaning |
| -------- | ------- | ------- |
| `OMEGA_SINGLE_PRECISION` | off | Make the default `Real` type single precision (4 bytes) instead of double (8 bytes). Use the default unless you are studying the performance or accuracy of single precision. |
| `OMEGA_VECTOR_LENGTH` | 1 | Vector length for CPU loops |
| `OMEGA_MPI_ON_DEVICE` | on | Let MPI send and receive halo buffers directly from device memory |
