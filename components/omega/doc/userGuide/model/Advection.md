(omega-user-model-advection)=

# Advection

Omega advects pseudo-thickness, normal velocity and tracers horizontally and
vertically. The schemes are chosen in the `Advection` section of the
configuration file, and each advection term is enabled in the `Tendencies`
section. The schemes themselves are described in
{ref}`omega-tech-tracer-advection`.

## Horizontal advection

Horizontal advection of tracers is enabled with
`Tendencies: TracerHorzAdvTendencyEnable`. Within the `Advection` group:

| Configuration name | Description | Options |
| ------------------ | ----------- | ------- |
| HorzTracerFluxOrder | order of the horizontal tracer flux | 2 for standard linear advection, 3 for the second order advection algorithm |
| HorzTracerFluxLimiterEnable | enable/disable monotonic flux-corrected transport (FCT); requires `HorzTracerFluxOrder` greater than 2 | true / false |
| HorzTracerFluxLimiterBudgetsEnable | enable/disable budget diagnostics when FCT is enabled | true / false |
| HorzTracerFluxLimiterMonotonicityCheckEnable | enable/disable a check for non-monotone values when FCT is enabled | true / false |
| FluxThicknessType | pseudo-thickness used in advective fluxes through edges | `Center` or `Upwind` |

When monotonic flux limiting is used then there are two diagnostic options that can be enabled,
HorzTracerFluxLimiterBudgetsEnable and HorzTracerFluxLimiterMonotonicityCheckEnable. The
HorzTracerFluxLimiterMonotonicityCheckEnable flag enables a post-transport check that the
resulting values are indeed monotone as required by the algorithm and any discrepancies are
printed along with information on where they occur.

The HorzTracerFluxLimiterBudgetsEnable option enables the output of two diagnostic fields to the
output mesh file.  One is the "FCTActiveTracerHorizontalAdvectionEdgeFlux" variable which is the
edge flux across every side of every element. The other is the "FCTActiveTracerHorizontalAdvectionTendency"
variable which is the cell-centered value of the advection tendency.

## Vertical advection

In Omega, vertical advection represents the transport of mass, momentum, and
tracer quantities along the vertical coordinate due to vertical motion within a
water column. In the context of Omega's
[governing equations](omega-design-governing-eqns-omega1), vertical advection
contributes to the tendencies of pseudo-thickness, horizontal velocity, and
tracer fields. Vertical motion is inferred from the divergence of horizontal
flow through the continuity equation and is necessary for accurately
representing vertical transport in three-dimensional ocean simulations.

The `VertAdv` class contains variables and functions for computing the vertical
velocity and the tendencies of `PseudoThickness`, `NormalVelocity`, and
`Tracers`. The algorithms used are determined by options specified in the
configuration file.

### Configuration options

Within the `Tendencies` group of the configuration file, flags are used to
enable or disable the contribution of vertical advection to the tendencies of
thickness, velocity, and tracers:
| Configuration name | Options |
| ------------------ | ------- |
| ThicknessVertAdvTendencyEnable | true / false |
| VelocityVertAdvTendencyEnable | true / false |
| TracerVertAdvTendencyEnable | true / false |

Within the `Advection` group, options are provided to configure the algorithms
used to compute the tracer tendencies:
| Configuration name | Description | Options |
| ------------------ | ----------- | ------- |
| VerticalTracerFluxOrder | order of accuracy used for the tracer flux calculation | 2, 3, or 4 |
| VerticalTracerFluxLimiterEnable | selects the tracer tendency algorithm | true / false |
| Coef3rdOrder | coefficient used for blending 3rd- and 4th-order fluxes when `VerticalTracerFluxOrder == 3` | 0.0 - 1.0 |

For `VerticalTracerFluxLimiterEnable`, `true` enables flux-corrected transport,
and `false` selects the standard algorithm. For `Coef3rdOrder`, `1` corresponds
to purely 3rd-order, `0` purely 4th-order.

The `VertAdv` class provides four fields that can be written to output by
adding them to the contents of an output stream in the configuration file:
`VerticalPseudoVelocity`, `TotalVerticalPseudoVelocity`, `VertFlux`, and
`LowOrderVertFlux`. These fields collectively make up the `VertAdv` output group.
