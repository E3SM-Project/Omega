(omega-user-build-e3sm)=

# Building in E3SM

Omega can be built as the ocean component (`ocn`) of a coupled E3SM simulation,
using a compset where Omega is selected as the ocean model (see
`components/omega/cime_config/config_compsets.xml` for supported compsets, such
as `COMEGA`, `GOMEGA`, or `WCYCL...-OMEGA`).

When an E3SM case is built using CIME's `case.build`, CIME automatically configures
and builds Omega via `components/omega/cime_config/buildlib_cmake`. In this mode,
CIME configures the build with `OMEGA_BUILD_MODE=E3SM`, which builds `OmegaLib`
together with the coupled driver bridge in `src/drivers/coupled/` (the Fortran/C++
bridge and MCT cap), creating the `ocn` library that CIME links into the coupled
`e3sm.exe` executable.

During case setup and execution, CIME runs Omega's configuration generator
(`components/omega/cime_config/buildnml`) in-process within Python to generate
`omega.yml` from YAML configuration files. Because `buildnml` requires `pyyaml`,
the Python environment running CIME commands must have `pyyaml` available
(see {ref}`omega-dev-cime-env`).

For comprehensive instructions on setting up, building, and running coupled cases
and CIME system tests, see {ref}`omega-dev-cime-testing`.
