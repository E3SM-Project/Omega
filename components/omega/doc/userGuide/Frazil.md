(omega-user-frazil)=

# Frazil

This page describes user-facing configuration for frazil tendencies in Omega.
Frazil physics is used to represent the formation and melt of frazil ice within
the ocean water column. It impacts the local layer pseudo-thickness (i.e. mass),
temperature, and salinity tendencies. The vertical sum of the frazil energy,
mass of water and mass of salt are passed to the coupler (if coupled) or
discarded (in ocean standalone mode).

## Configuration overview

Frazil behavior is controlled by one enable switch in `Tendencies` and one
`Frazil` configuration block:

```yaml
Omega:
  Tendencies:
    FrazilTendencyEnable: true

  Frazil:
    FrazilType: FixedProperty
    MassLimit: 0.1
    Phi: 0.75
    DepthLimit: -1.0
    ConservationCheck: false
```

- `Tendencies.FrazilTendencyEnable`
  - enables/disables application of frazil tendency contributions
- `Frazil.FrazilType`
  - selects frazil option
  - supported options in current code: `FixedProperty` and `Teos10`
  - `basic` is a deprecated compatibility alias for `FixedProperty`
  - `simple` exists as a placeholder name but is not supported
- `Frazil.MassLimit`
  - limits per-layer frazil mass/thickness tendency magnitude (applied to formation and melt)
- `Frazil.Phi`
  - liquid-mass fraction of frazil (used by teos frazil formation, or by fixed-property frazil when using porosity rather than constant salinity).
- `Frazil.DepthLimit`
  - limits the depth range over which frazil is computed
  - negative values mean no depth limit
- `Frazil.ConservationCheck`
  - enables a column-level diagnostic conservation check with logging

## Available frazil options

Omega currently includes two active frazil pathways:

- `FixedProperty`
  - freezing is based on the formation of fresh solid ice, to which salt is added (similar the mpas-ocean implementation).
- `Teos10`
  - teos-10-based option using Gibbs SeaWater thermodynamic routines.

Both pathways contribute to:

- pseudo-thickness tendency
- temperature tracer tendency
- salinity tracer tendency

## Notes

- Frazil tendencies are applied through the `Tendencies` tracer-step workflow.
- The frazil tendency hook assumes tracer names include `Temperature` and
  `Salinity`.
- For implementation and algorithm details, see
  [Developer Frazil Guide](../devGuide/Frazil.md).
