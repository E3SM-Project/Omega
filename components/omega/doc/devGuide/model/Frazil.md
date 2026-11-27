(omega-dev-frazil)=

# Frazil

This page describes frazil design and implementation details in Omega,
including both the `FixedProperty` and `Teos10` pathways.

## Purpose and coupling points

Frazil computes frazil-related tendencies that modify:

- pseudo-thickness tendency
- temperature tracer tendency
- salinity tracer tendency

The tendency hook-up is implemented through `FrazilOnCell` in the tracer
tendency phase, where frazil contributions are added to the accumulated
`PseudoThicknessTend` and `TracerTend` arrays.

It also accumulates the column-integrated frazil terms over an ocean (outer) timestep to store the fluxes to be passed through the coupler.

## Data flow and call sequence

1. `Tendencies::computeTracerTendenciesOnly` checks
   `Tendencies.FrazilTendencyEnable`.
2. If enabled, `FrazilOnCell::operator()` retrieves the default `Frazil`
   object and zeros frazil tendency/accumulator arrays.
3. `FrazilOnCell` extracts `Temperature` and `Salinity` tracer subviews,
   then calls `Frazil::computeFrazil(CT, SA, PressureMid, PseudoThickness)`.
4. `Frazil::computeFrazil` dispatches to `computeFrazilFixedPropertyImpl` or
   `computeFrazilTeosImpl` based on `FrazilType`.
5. Returned frazil tendencies are added into `PseudoThicknessTend` and
   `TracerTend` for all active cell layers.

## Configuration coupling

Frazil behavior is configured with:

- `Omega.Tendencies.FrazilTendencyEnable`
  - global switch for applying frazil tendency terms
- `Omega.Frazil.FrazilType`
  - implementation choice (`FixedProperty` or `Teos10`)
- `Omega.Frazil.LayerMassFracMax`
  - per-layer mass/thickness limiter used by frazil formation/melt pathways
- `Omega.Frazil.Phi`
  - teos pathway liquid-fraction parameter for new frazil partitioning
- `Omega.Frazil.DepthLimit`
  - optional depth cutoff for frazil activity; negative means no cutoff
- `Omega.Frazil.ConservationCheck`
  - optional post-compute column conservation diagnostic/logging

## Physics/algorithm summary

The physics and the two pathways (`FixedProperty` and `Teos10`) are described
in {ref}`omega-tech-frazil`.

## Existing ctest coverage

The existing frazil test driver is in
`components/omega/test/ocn/FrazilTest.cpp` and covers:

- teos frazil formation in cold and warm single-layer states
- fixed-property frazil formation in cold and warm single-layer states
- mixed warm/cold column behavior with sign checks for branch switching
- depth-limit behavior ensuring excluded deep layers have zero frazil tendency


## Related pages

- User-facing options: [User Frazil Guide](#omega-user-model-frazil)
- Tendency container/hook-up: [Tendencies](#omega-dev-tendencies)
