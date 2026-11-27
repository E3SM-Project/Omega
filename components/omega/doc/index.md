(omega-home)=
# Omega

The Ocean Model for E3SM Global Applications (Omega) is a global ocean model
under development by the [E3SM](https://e3sm.org/) ocean team. It is
initially eddy-resolving, with a non-eddying configuration to follow.

The model is written in C++ using the [Kokkos](https://github.com/kokkos)
framework for performance portability. Omega is based on the
[TRSK formulation](https://doi.org/10.1016/j.jcp.2009.08.006) for geophysical
models on unstructured meshes. The first version of Omega is primarily a
direct port of the
[MPAS-Ocean](https://e3sm.org/model/e3sm-model-description/v1-description/v1-ocean-sea-ice-land-ice/)
component of E3SM, for comparison purposes.

Development takes place at <https://github.com/E3SM-Project/Omega>.

## How these docs are organized

Each part answers one kind of question:

- The {ref}`omega-user-guide` answers *how do I set up and run Omega, and
  what can I control?* It covers building, running, configuration, input and
  output, and the model options.
- The {ref}`omega-dev-guide` answers *how do I work on the code?* It covers
  the development environment, contributing, the framework classes, how the
  ocean model is put together, and how-to recipes.
- The {ref}`omega-tech-guide` answers *what does the model compute, and why?*
  It holds the equations, their discretization and the parameterizations.
- The {ref}`omega-design` record *what we were thinking when we built it*.
  They are frozen when written; the guides describe the code as it is.

Running Omega standalone is done through
[Polaris](https://docs.e3sm.org/polaris/main/), whose documentation covers
setting up and running tasks and suites.

```{toctree}
:hidden:
:maxdepth: 2

userGuide/index
devGuide/index
techGuide/index
design/index
```
