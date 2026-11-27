(omega-dev-code-organization)=

# Code Organization

Omega lives in `components/omega` of the E3SM repository. Its top level holds
the CMake build (`CMakeLists.txt`, `OmegaBuild.cmake`), the development
environment (`pixi.toml`, with `dev-conda.txt` as the conda fallback; see
{ref}`omega-dev-env`), the linter settings and these directories:

| Directory | Contents |
| --------- | -------- |
| `src/` | the model source code, by layer (below) |
| `test/` | the CTest unit tests, in the same subdirectories as `src/` |
| `configs/` | `Default.yml`, the default value of every configuration option |
| `cime_config/` | the E3SM component interface: compsets, `buildlib_cmake`, `buildnml` and its `omega_buildnml` package, and the CIME test definitions |
| `ci/` | the scripts that GitHub Actions runs to build and test Omega in the E3SM container images |
| `external/` | submodules Omega builds itself: GSW-C, cpptrace, spdlog and yaml-cpp |
| `doc/` | this documentation |

## Source layers

Code in `src/` is arranged from the lowest-level infrastructure up to the
drivers:

| Directory | Contents | Described in |
| --------- | -------- | ------------ |
| `src/base/` | the lowest level: data types, the machine environment, decomposition, halos, broadcasts, reductions, parallel I/O, fill values and tridiagonal solvers | {ref}`omega-dev-framework` |
| `src/infra/` | infrastructure shared by many routines: configuration, logging, errors, fields and dimensions, I/O streams, CF units, the time manager and Kokkos wrappers | {ref}`omega-dev-framework` |
| `src/ocn/` | the ocean model: state, mesh and operators, vertical coordinate, auxiliary variables, tendencies, the equation of state, forcing and the parameterizations | {ref}`omega-dev-model` |
| `src/timeStepping/` | the time steppers | {ref}`omega-dev-time-stepping` |
| `src/analysis/` | in-situ analysis operators and groups | {ref}`omega-dev-analysis` |
| `src/drivers/` | the top-level drivers, `standalone/` and `coupled/` | {ref}`omega-dev-driver` |
