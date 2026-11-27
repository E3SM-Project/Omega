(omega-user-io-fields)=

# Fields and Metadata

Within OMEGA, each module registers fields that it owns together with the
metadata associated with that field. Metadata must be compliant with the
Climate and Forecast [(CF) metadata conventions.](http://cfconventions.org/)
This information can then be used by IOStreams or other modules that may
require metadata information. Two special fields are defined internally as
CodeMeta and FieldMeta to hold any global metadata not associated with a
data field, but may be required for provenance. Required metadata are enforced
by the Field create method but additional metadata can be added after field
creation. The data array holding field values can be attached/retrieved and
updated whenever needed.

Groups of fields can also be defined to provide a shortcut to a list of fields
commonly used together, like the model state or tracer groups. This reduces the
size of the input YAML configuration. The user can select which fields or field
groups they wish to include via the streams section of the input YAML
configuration file. This is described further in {ref}`omega-user-io-streams`.
There are no user-configurable options for the Field class itself since all
fields are defined and updated internally, but further details on the methods
and implementation can be found in the {ref}`omega-dev-field` section of the
Developer's Guide.

(omega-user-io-fields-inactive-layers)=

## Inactive and boundary layers in output

Each output field declares a fill value (matching the NetCDF-C standard `_FillValue`) that marks entries outside the active domain.
In a cell column only the active layers (`[MinLayerCell, MaxLayerCell]`) carry valid data; deeper, inactive layers contain the fill value.
Edge and vertex fields have an intermediate boundary range where some, but not all, of the surrounding cells are active.
For flux-type edge fields (such as `NormalVelocity`), the flux through a face bordering land vanishes, so those boundary layers carry `0`; layers fully outside the valid range carry the fill value.
Vertex fields, by contrast, keep their computed (generally non-zero) values throughout the valid range, since a boundary vertex still has active surrounding cells contributing to quantities like vorticity; only layers outside the valid range carry the fill value.
After an initial-condition or restart read, `VertCoord` enforces this pattern on the state and tracer fields so that analysis tools (ncview, Xarray, NCO, etc.) recognize and mask the inactive entries automatically.

(omega-user-io-fields-auxiliary)=

## Auxiliary variables

Omega computes a number of auxiliary variables from the state at every time
step, and they can be written to output like any other field. The following
auxiliary variables are currently available:
| Name | Description |
| ------------------- | ------- |
| KineticEnergyCell | kinetic energy of horizontal velocity on cells
| VelocityDivCell | divergence of horizontal velocity
| FluxPseudoThickEdge | pseudo-thickness used for fluxes through edges. May be centered, upwinded, or a combination of the two
| MeanPseudoThickEdge | pseudo-thickness averaged from cell center to edges
| RelVortVertex | curl of horizontal velocity, defined at vertices
| NormRelVortVertex | curl of horizontal velocity divided by pseudo-thickness
| NormPlanetVortVertex | earth's rotational rate (Coriolis parameter, f) divided by pseudo-thickness
| NormRelVortEdge | curl of horizontal velocity divided by pseudo-thickness, averaged from vertices to edges
| NormPlanetVortEdge | earth's rotational rate (Coriolis parameter, f) divided by pseudo-thickness, averaged from vertices to edges
| VelDel2Edge | laplacian of horizontal velocity on edges
| VelDel2DivCell | divergence of laplacian of horizontal velocity on cells
| VelDel2RelVortVertex | curl of laplacian of horizontal velocity on cells
| HTracersEdge | thickness-weighted tracers used for fluxes through edges. May be centered, upwinded or a combination of the two
| Del2TracersCell | laplacian of tracers on cells
| SurfTracerRestoringDiffsCell | surface tracer restoring differences on cells
| TracersMonthlySurfClimoCell | monthly climatology values to restore to for surface tracer on cells
| DenMixLayerDepth | density-threshold mixed-layer depth on cells
