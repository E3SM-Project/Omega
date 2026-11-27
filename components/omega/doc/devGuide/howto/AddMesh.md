(omega-dev-howto-add-mesh)=

# Add a Supported Mesh

A supported mesh is one that a coupled E3SM case can use as its ocean grid:
`buildnml` finds its input files and any overrides from the case's `OCN_GRID`.
How `buildnml` layers the configuration is described in
{ref}`omega-dev-buildnml`.

1. Confirm the mesh's grid alias is already defined for E3SM in
   `cime_config/config_grids.xml` at the repository root; `buildnml` looks
   up input files by the case's `OCN_GRID` value.
2. Add an entry for the mesh to `cime_config/omega_buildnml/data/input_files.yaml`, listing the input
   file(s) that provide its `HorzMeshIn`, `InitialVertCoord`, and
   `InitialState` streams (split across multiple `inputs` entries if the
   initial condition is a separate file from the mesh).
3. If the mesh needs overrides beyond Omega's defaults and the coupled
   overrides, add a `meshes.<mesh name>` entry to `cime_config/omega_buildnml/data/config_overrides.yaml`.
4. Confirm the new entries are complete and consistent with `Default.yml`
   (see {ref}`omega-dev-buildnml`):

   ```bash
   cd components/omega/cime_config
   ./validate_config.py
   ```
