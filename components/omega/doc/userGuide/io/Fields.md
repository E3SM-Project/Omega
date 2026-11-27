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
