(omega-user-model-tendencies)=

# Tendency Terms

Tendencies for updating state variables are computed using functors, which are
class objects that can be called like functions. These functors are used within
Kokkos parallel loops to compute the full tendency arrays. The following
tendency terms are currently implemented:
| Name | Description |
| ---------------- | ---------------- |
| PseudoThicknessFluxDivOnCell | divergence of pseudo-pseudo-thickness flux, defined at cell center
| PotentialVortHAdvOnEdge | horizontal advection of potential vorticity, defined on edges
| KEGradOnEdge | gradient of kinetic energy, defined on edges
| SSHGradOnEdge | gradient of sea-surface height, multiplied by gravitational acceleration, defined on edges
| VelocityDiffusionOnEdge | Laplacian horizontal mixing, defined on edges
| VelocityHyperDiffOnEdge | biharmonic horizontal mixing, defined on edges
| TracerHorzAdvOnCell | horizontal advection of thickness-weighted tracers
| TracerDiffOnCell | horizontal diffusion of thickness-weighted tracers
| TracerHyperDiffOnCell | biharmonic horizontal mixing of thickness-weighted tracers
| SfcStressForcingOnEdge | forcing by surface stress (e.g. wind), defined on edges
| BottomDragOnEdge | bottom drag, defined on edges
| SfcThicknessForcingOnCell | surface pseudo-thickness forcing from coupled freshwater and salt fluxes, defined on cells
| SfcTracerForcingOnCell | surface tracer forcing from coupled heat and salt fluxes, with direct heat always and mass-flux enthalpy terms gated by thickness forcing, defined on cells
| PenetratingShortwaveOnCell | depth-dependent three-band shortwave heating, defined on cells
| SurfaceTracerRestoringOnCell | surface tracer restoring, defined on cells

Among the internal data stored by each functor is a `bool` which can enable or
disable the contribution of that particular term to the tendency. These flags
need to be provided in the configuration file:
```yaml
Tendencies:
   ThicknessFluxTendencyEnable: true
```

Some functors have member variables that can be set by the user by the
configuration file. The following are the user-configurable parameters for the
currently available tendency terms; the coefficients of the horizontal mixing
terms are described in {ref}`omega-user-model-horizontal-mixing`, and the
advection schemes in {ref}`omega-user-model-advection`:
| Term | Parameter | Description
| ------------ | ------------ | ------------ |
| PseudoThicknessFluxDivOnCell | ThicknessFluxTendencyEnable | enable/disable term
| PotentialVortHAdvOnEdge | PVTendencyEnable | enable/disable term
| KEGradOnEdge | KETendencyEnable | enable/disable term
| SSHGradONEdge | SSHTendencyEnable | enable/disable term
| VelocityDiffusionOnEdge | VelDiffTendencyEnable | enable/disable term
| VelocityHyperDiffOnEdge | VelHyperDiffTendencyEnable | enable/disable term
| TracerHorzAdvOnCell | TracerHorzAdvTendencyEnable | enable/disable term
| TracerDiffOnCell | TracerDiffTendencyEnable | enable/disable term
| TracerHyperDiffOnCell | TracerHyperDiffTendencyEnable | enable/disable term
| SfcStressForcingOnEdge | SfcStressForcingTendencyEnable | enable/disable term
| BottomDragOnEdge | BottomDragTendency | enable/disable term
| | BottomDragTendency:Mode | bottom drag mode; `Implicit` or `Explicit`
| | BottomDragTendency:Type | bottom drag type; `Constant`
| | BottomDragTendency:BottomDragCoeff | bottom drag coefficient
| SfcThicknessForcingOnCell | SfcThicknessForcingTendencyEnable | enable/disable term
| SfcTracerForcingOnCell | SfcTracerForcingTendencyEnable | enable/disable term
| PenetratingShortwaveOnCell | PenetratingShortwaveTendency:Enable | enable/disable depth-dependent shortwave heating; its parameters are described in {ref}`omega-user-model-forcing`
| SurfaceTracerRestoringOnCell | SurfaceTracerRestoringEnable | enable/disable term

The vertical advection, pressure gradient and vertical mixing tendencies are not
implemented as the functors above, but they are enabled from the same
`Tendencies` section of the configuration file:

| Parameter | Description
| ------------ | ------------ |
| ThicknessVertAdvTendencyEnable | enable/disable vertical advection of thickness
| VelocityVertAdvTendencyEnable | enable/disable vertical advection of velocity
| TracerVertAdvTendencyEnable | enable/disable vertical advection of tracers
| PressureGradTendencyEnable | enable/disable the pressure gradient tendency
| VelVertMixTendencyEnable | enable/disable vertical mixing of velocity; required when bottom drag uses `Implicit` mode
| TracerVertMixTendencyEnable | enable/disable vertical mixing of tracers
| KPPNonLocalTracerFluxTendencyEnable | enable/disable the KPP non-local tracer flux; see [KPP Boundary Layer Mixing](#omega-user-model-kpp)
| KPPNonLocalTracerDiagnosticsEnable | enable/disable diagnostic output of the non-local flux tendency

The non-local options are only read when velocity or tracer vertical mixing is
enabled. If omitted, `KPPNonLocalTracerFluxTendencyEnable` defaults to disabled
and `KPPNonLocalTracerDiagnosticsEnable` defaults to enabled.
The non-local tracer tendency is applied only when
`TracerVertMixTendencyEnable` is also enabled; velocity-only vertical mixing
does not apply a tracer mixing tendency.

Two further flags control tendency sources rather than individual terms:
`UseCustomTendency` enables user-supplied tendencies and
`ManufacturedSolutionTendency` adds the forcing used by the manufactured
solution test case.

## See Also

Additional information on forcing, including surface stress forcing,
surface thickness and tracer flux forcing, and surface tracer restoring, is detailed in
[](omega-user-model-forcing).
