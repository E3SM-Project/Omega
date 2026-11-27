(omega-user-model-tracers)=

# Tracers

Tracers refer to either heat (temperature) or material carried by a fluid
parcel (e.g., salt, chemical, or biological constituents).

The tracers that can be selected are defined in the code; to add one, see
{ref}`omega-dev-howto-add-tracer`.

## Selecting tracers using YAML configuration file

Note that not all tracers defined in `TracerDefs.inc` will be used during
a simulation. To select tracers and groups of tracers, users will configure
YAML files in OMEGA, such as `Default.yml` in the `omega/configs` directory,
as shown below:

```yaml
omega:
  Tracers:
    Base: [Temp, Salt]
    Debug: [Debug1, Debug2, Debug3]
    [other individual tracers or groups as needed]
```

In the above example, two tracer groups (Base and Debug) are selected. The
Base group includes the `Temp` and `Salt` tracers, while the Debug group
includes `Debug1`, `Debug2`, and `Debug3`.
