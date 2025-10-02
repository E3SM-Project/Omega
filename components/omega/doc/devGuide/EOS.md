(omega-dev-eos) =

# Equation of State (EOS)

Omega includes an `Eos` class that provides functions that compute `SpecVol` and `SpecVolDisplaced`.
Current EOS options are a linear EOS or an EOS computed using the TEOS-10 75 term expansion from
[Roquet et al. 2015](https://www.sciencedirect.com/science/article/pii/S1463500315000566). If
`SpecVolDisplaced` is calculated with the linear EOS option, it will be equal to `SpecVol` as there
is no pressure/depth dependence for the linear EOS. `SpecVolDisplaced` computes specific volume
adiabatically displaced to `K + KDisp`.

## Eos type

An enumeration listing all implemented schemes is provided. It needs to be extended every time an
EOS is added. It is used to identify which EOS method is to be used at run time.

```c++
enum class EosType { LinearEos, Teos10Eos };
```

## Initialization

An instance of the `Eos` class requires a [`HorzMesh`](#omega-dev-horz-mesh), so the mesh class
and all of its dependencies need to be initialized before the `Eos` class can be. The static method:

```c++
OMEGA::Eos::init();
```

initializes the default `Eos`. A pointer to it can be retrieved at any time using:

```c++
OMEGA::Eos* DefEos = OMEGA::Eos::getInstance();
```

## Computation of Eos

To compute `SpecVol` for a particular set of temperature, salinity, and pressure arrays, do

```c++
Eos.computeSpecVol(ConsrvTemp, AbsSalinity, Pressure);
```

`SpecVolDisplaced` is calculated using local temperature and salinity values, but a pressure
value at `K + KDisp`. To compute `SpecVolDisplaced` for a particular set of temperature, salinity,
and pressure arrays and displaced vertical index level, do

```c++
Eos.computeSpecVolDisp(ConsrvTemp, AbsSalinity, Pressure, KDisp);
```

where `KDisp` is the number of `k` levels you want to displace each specific volume level to.
For example, to displace each level to one below, set `KDisp = 1`.

## Bounds check (and truncation) for the state variables (under TEOS-10)

The implemented 75-term polynomial for the calculation of the specific volume under TEOS-10 is considered valid for ocean states contained in the ''oceanographic funnel'' defined in [McDougall et al., 2003](https://journals.ametsoc.org/view/journals/atot/20/5/1520-0426_2003_20_730_aaceaf_2_0_co_2.xml). When using TEOS-10, the Eos uses member methods `calcSLimits(P)` and `calcTLimits(Sa, P)` to calculate the valid ranges of Sa and T given the pressure. The conservative temperature lower bound is set by the freezing temperature, using the member method `calcCtFreezing(Sa, P, SaturationFract)`. This method implements the polynomial approximation of the conservative freezing temperature (called `gsw_ct_freezing_poly` in the GSW package), which is known to produce erros in the (-5e-4 K, 6e-4 K) range. Once we calculate the upper and lower bounds of validity, the state variables are clipped to the valid range (if outside the bounds) before we run the specific volume calculation. The state fields themselves are not changed.

 ## Removal of Eos

To clear the Eos instance do:

```c++
OMEGA::Eos::destroyInstance();
