(omega-user-model-eos)=

# Equation of State

The equation of state gives the specific volume of seawater from Conservative
Temperature, Absolute Salinity and pressure. Omega provides three: a linear
EOS, a constant EOS and TEOS-10. The linear and constant options are for
simple idealized test cases only and are not appropriate for realistic ocean
simulations; TEOS-10 is the preferred EOS for real ocean simulations. Each is
described in {ref}`omega-tech-eos`.

The user-configurable option `EosType` can be set to `linear`, `teos10`, or `constant`. Parameters in the `Linear` subsection are used only when `EosType` is set to `linear`.

```yaml
Eos:
   EosType: teos10
   Linear:
      DRhoDT: -0.2
      DRhoDS: 0.8
      RhoT0S0: 1000.0
```

where `DRhoDT` is the thermal expansion coefficient ($\textrm{kg}/(\textrm{m}^3 \cdot ^{\circ}\textrm{C})$), `DRhoDS` is the saline contraction coefficient ($\textrm{kg}/\textrm{m}^3$), and `RhoT0S0` is the reference density at (T,S)=(0,0) (in $\textrm{kg}/\textrm{m}^3$).

## Output fields

In addition to `SpecVol`, the EOS computes the displaced specific volume
`SpecVolDisplaced`, the squared Brunt-Vaisala frequency `BruntVaisalaFreqSq`,
and the first derivatives of specific volume `SpecVolDCt`, `SpecVolDSa` and
`SpecVolDP`. All are in the `Eos` field group and can be requested in a
stream's contents.

There is no option for the derivatives. They are computed only when a part of
the model that needs them, such as the finite-volume pressure gradient, is
enabled, in the same pass that computes the specific volume, so selecting
`teos10` does not make the model slower unless such a scheme is enabled.
