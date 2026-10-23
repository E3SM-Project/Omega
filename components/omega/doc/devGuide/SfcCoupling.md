(omega-dev-sfc-coupling)=

# Surface Coupling

The `SfcCoupling` class manages the variables exchanged to (`o2x`) and from
(`x2o`) the coupler. It handles import/export of the raw coupler data,
application of imported fields to the [`Forcing`](#omega-dev-forcing) object,
and accumulation of export fields over the coupling interval. It is possible
to have multiple `SfcCoupling` instances; every instance has a name and is
tracked in a static C++ map called `AllSfcCoupling`.

## Initialization

The static method:
```c++
OMEGA::SfcCoupling::init(CouplingInitParams);
```
initializes the default `SfcCoupling` given a `CouplingInitParams` struct
(number of import/export fields, name-to-index maps, coupling time step, and
coupling `Layout`, i.e. `MCT` or `MOAB`) supplied by the coupler. A pointer to
the default instance can be retrieved at any time using:
```c++
OMEGA::SfcCoupling* DefSfcCoupling = OMEGA::SfcCoupling::getDefault();
```

## Creation of non-default surface coupling objects

A non-default instance can be created with the static `create` method,
which returns a pointer to the newly created object:
```c++
OMEGA::SfcCoupling* NewSfcCoupling = OMEGA::SfcCoupling::create(
    Name, Mesh, NImportFields, NExportFields, ImportIdxMap, ExportIdxMap,
    Stepper, CouplingTimeStep, Layout);
```
Given its name, a pointer to a named instance can be obtained at any time:
```c++
OMEGA::SfcCoupling* NewSfcCoupling = OMEGA::SfcCoupling::get(Name);
```

## Data exchange

Before import/export, unmanaged views must be attached to the raw coupler
data pointers:
```c++
SfcCoupling.attachData(CplToOcnData, OcnToCplData);
```
To copy data from the coupler into the `CplToOcn` fields, and from the
`OcnToCpl` fields back out to the coupler:
```c++
SfcCoupling.importFromCoupler();
SfcCoupling.exportToCoupler();
```
Imported fields are applied to a `Forcing` object with:
```c++
SfcCoupling.applyImportFields(ForcingPtr);
```
Export fields are accumulated (running average) each ocean time step with:
```c++
SfcCoupling.updateExportFields(State, TracerArray);
```

## Units

`OcnToCplFields`'s averaged temperature is Conservative Temperature (deg C)
on device (`AvgSfcTemperature`), invariant at all times. `copyToHost()`
converts it to in-situ temperature in Kelvin (via TEOS-10, when in use)
into the host mirror `AvgSfcTemperatureH`, once per coupling interval.
All other averaged fields keep the same units and value on host and device.

To keep this invariant, `OcnToCplFields`'s averaged device arrays are
private: they can only be written via `updateAverages()` (called each
ocean time step from `updateExportFields()`) and read via `copyToHost()`.

## Frazil Sign Conventions

The frazil exports use the following signs during formation:

| Coupler field | Units | Positive meaning |
| --- | --- | --- |
| `Fioo_frazilm` | kg m^-2 s^-1 | Mass removed from the ocean into frazil (counted positive out of the ocean)|
| `Fioo_frazils` | kg m^-2 s^-1 | Salt removed from the ocean with frazil (counted positive out of the ocean)|
| `Fioo_frazilh` | W m^-2 | Energy gain by the ocean from frazil formation (counted positive into the ocean)|

The current sign of `Fioo_frazilm` and `Fioo_frazils` goes against the
common coupling convention. The convention for these fluxes will be flipped
once the sea ice and coupler conventions are changed.

Newly formed frazil has negative enthalpy, so removing it increases ocean
energy. `FrazilEnergyFlux` already represents this positive-into-ocean
contribution and is exported without an additional sign change. Sea ice
uses the positive `Fioo_frazilh` to reconstruct negative ice enthalpy.

`OcnToCplFields::updateFrazilFields` negates the ocean-facing
`FrazilMassFlux` and `FrazilSaltFlux`, but preserves `FrazilEnergyFlux`,
and accumulates coupling-interval means. Mass and salt exports therefore
use a positive-out-of-ocean convention to match the sea-ice expectation,
whereas energy retains the CPL positive-into-ocean convention.

`Fioo_q` is separate, nonpositive melt potential, not formation heat or an
already-applied heat flux. It is zero when frazil production is exported;
otherwise it represents the surface ocean's available energy for melting.

## Removal of surface coupling objects

To erase a specific named instance use `erase`:
```c++
OMEGA::SfcCoupling::erase(Name);
```
To clear all instances do:
```c++
OMEGA::SfcCoupling::clear();
```
