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

To compute `SpecVol` for a particular set of temperature `Ct`, salinity `Sa`, and pressure `P` arrays, do

```c++
Eos.computeSpecVol(Ct, Sa, P);
```

`SpecVolDisplaced` is calculated using local temperature and salinity values, but a pressure
value at `K + KDisp`. To compute `SpecVolDisplaced` for a particular set of temperature, salinity,
and pressure arrays and displaced vertical index level, do

```c++
Eos.computeSpecVolDisp(Ct, Sa, P, KDisp);
```

where `KDisp` is the number of `k` levels you want to displace each specific volume level to.
For example, to displace each level to one below, set `KDisp = 1`.

## Bounds check (and truncation) for the state variables (under TEOS-10)

The implemented 75-term polynomial for the calculation of the specific volume under TEOS-10 has been evaluated for ocean states in the ''cube" (-2-40 C; 0-42 g/kg; 0-10,000 dbar) and the ''oceanographic funnel'' defined in [McDougall et al., 2003](https://journals.ametsoc.org/view/journals/atot/20/5/1520-0426_2003_20_730_aaceaf_2_0_co_2.xml). When using TEOS-10, the Eos uses member methods `calcSLimits(P)` and `calcTLimits(Sa, P)` to calculate the valid ranges of Sa and Ct. When using the `Funnel` opion of `EosLimits`, the salinity limits are calculated as a function of pressure and the temperature limits as a function of pressure and salinity. The conservative temperature lower bound is set by the freezing temperature, using the member method `calcCtFreezing(Sa, P, SaturationFract)`. This method implements the polynomial approximation of the conservative freezing temperature (called `gsw_ct_freezing_poly` in the GSW package), which is known to produce errors in the (-5e-4 K, 6e-4 K) range. When using the `Cube` option of `EosLimits`, the bounds are constant. Once we calculate the upper and lower bounds of validity, warnings are issued when the state variables are outside the validity bounds. If `ClampingEnable` is true, the state variables are clipped to the valid range (if outside the bounds) before we run the specific volume calculation. The state fields themselves are not changed. If `ClampingEnable` is false, the warnings are issued but the specific volume is calculated based on the unchanged state variables.

## Underlying helper functions when using TEOS-10

The computation of the TEOS-10 specific volume relies on 3 underlying member functions:

```c++
calcPCoeffs(SpecVolPCoeffs, K, Ct, Sa, P);
```
which calculates the coefficents that will be applied to the pressure. These coefficients are dependent on the temperature and salinity variables but not the pressure. The pressure argument present in the function call is only used to set the bounds of Ct,Sa validity. This function is where the bulk of the polynomial calculation takes place and thus it is advised to reuse the `SpecVolPCoeffs` which are a class data member if possible to reduce computational expense.

```c++
calcRefProfile(Pressure);
```
which calculates the ocean reference profile as a function of pressure in the layers. This is a simple 6th order polynomial in P with constant coefficients. It could be reused within a timestep provided that the layer pressures do not change but is cheap to recalculate.

```c++
calcDelta(SpecVolPCoeffs, K, Pressure);
```
which applies the `SpecVolPCoeffs` calculated above to the pressure state variable. The resulting "delta" in specific volume is added to the reference profile calculated above to form the specific volume. This step is a simple 5th order polynomial in pressure, which combines the effects of T,S (pre-calculated in `calcPCoeffs`) and the pressure.


 ## Removal of Eos

To clear the Eos instance do:

```c++
OMEGA::Eos::destroyInstance();
