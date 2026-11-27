(omega-dev-howto-add-tracer)=

# Add a Tracer

Tracers are defined in `src/ocn/TracerDefs.inc`. The file contains
Tracer index variables and `defineAllTracers` C++ function defintion that
contains the calls to `define` function per each tracer as shown below:

```c++
inline static I4 IndxTemp             = Tracers::IndxInvalid;
inline static I4 IndxSalt             = Tracers::IndxInvalid;
inline static I4 IndxMyBGCTracer      = Tracers::IndxInvalid;

// Tracer definitions packaged in a defineAllTracers function
static void defineAllTracers() {

   define("Temp",                            ///< [in] Name of tracer
          "Potential Temperature",           ///< [in] Long name or description
          "degree_C",                        ///< [in] Units
          "sea_water_potential_temperature", ///< [in] CF standard Name
          -273.15,                           ///< [in] min valid field value
          100.0,                             ///< [in] max valid field value
          FillValueReal,                     ///< [in] value for undef entries
          IndxTemp);                         ///< [out] (optional) static index

   define("Salt", "Salinity", "psu", "sea_water_salinity", 0.0, 50.0,
          FillValueReal, IndxSalt);
   define("Debug1", "Debug Tracer 1", "", "", 0.0, 100.0, FillValueReal);
   define("Debug2", "Debug Tracer 2", "", "", 0.0, 100.0, FillValueReal);
   define("Debug3", "Debug Tracer 3", "", "", 0.0, 100.0, FillValueReal);
}
```

The `FillValueReal` constant is defined in `FillValues.h` and matches the
NetCDF-C standard fill value (`NC_FILL_DOUBLE` or `NC_FILL_FLOAT` depending on
the build precision). It is automatically in scope in `TracerDefs.inc` via
`Field.h`. New tracers should always use `FillValueReal` (or `FillValueI4` /
`FillValueI8` for integer fields) rather than hardcoded literals.

To add a new tracer, simply call the `define` function with the appropriate
arguments. Index argument is optional one that allows to access the tracer
data using the given tracer index variable.

Once defined, a tracer can be selected in the configuration file as described
in {ref}`omega-user-model-tracers`. How tracers are stored and accessed is
described in {ref}`omega-dev-tracers`.
