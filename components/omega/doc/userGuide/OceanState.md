(omega-user-state)=

## Ocean State

The `OceanState` class provides a container for the non-tracer prognostic variables in Omega:
`NormalVelocity` and `PseudoThickness`.
Upon creation of an `OceanState` instance, these variables are allocated and registered with the
IO infrastructure as part of the `State` and `Restart` field groups, so they are read from the
initial-condition file and written to restart files.
The class contains a method to update the time levels for the state variables between timesteps.
This involves a halo update, time level index update, and updating the `IOFields` data references.

The relative pressure at the top of the ocean column, `SurfacePressure`, is owned by the
[vertical coordinate](omega-user-vert-coord) rather than the `OceanState`, since it is the top
boundary condition of the pressure field. It is still registered in the `State` and `Restart`
field groups, so it continues to be written to restart files and read from the initial-condition
and restart files when present; if it is absent, it defaults to zero.

(omega-user-state-validation)=

### State validation

After each time step, Omega checks `PseudoThickness`, `KineticEnergyCell`, `Temperature` and
`Salinity` for NaNs and for values outside a plausible range. For each problem it finds, the rank
that owns the cells writes a critical message giving the field, the number of entries and the
global cell ID, latitude, longitude and layer of the first one, with its value if it is out of
range. Ranks that do not log write these messages to standard error (see
[logging](omega-user-logging)). What happens next is set under `State: Validation`:

```yaml
Omega:
  State:
    Validation:
      AbortOnOOB: false
      AbortOnNan: true
      WriteStreamOnAbort: History
```

`AbortOnNan` and `AbortOnOOB` decide whether NaNs and out-of-range values abort the run; the
decision is based on the counts summed over all ranks. Before aborting, every rank writes the
stream named by `WriteStreamOnAbort`, so the bad state can be inspected. The record has the time
at the start of the failing step. Set it to `''` to write nothing; a stream that is not defined
is skipped. When the run continues instead, the master rank logs the total counts.
