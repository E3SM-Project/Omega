(omega-user-model-vertical-mixing)=

# Vertical Mixing

The vertical mixing module in Omega handles the parameterization of unresolved vertical mixing
processes in the ocean. It calculates vertical diffusivity and viscosity coefficients that
determine how properties (like momentum, heat, salt, and biogeochemical tracers) mix vertically
in the ocean model. Omega offers three linearly additive mixing processes within the water column:
(1) a constant background mixing value, (2) a convective instability mixing value, and (3) a
Richardson number-dependent shear-instability-driven mixing value. They are described in
{ref}`omega-tech-vertical-mixing`. Boundary-layer mixing from the K Profile Parameterization
is described in {ref}`omega-user-model-kpp`.

The user-configurable options are the following parameters in the yaml configuration file:

```yaml
VertMix:
  Background:
    Viscosity: 1.0e-4    # Background vertical viscosity (m²/s)
    Diffusivity: 1.0e-5  # Background vertical diffusivity (m²/s)
  Convective:
    Enable: true         # Enables the convective-induced mixing option
    Diffusivity: 1.0     # Convective mixing coefficient (m²/s)
    TriggerBVF: 0.0      # Squared Brunt-Vaisala frequency threshold
  Shear:
    Enable: true         # Enables the shear-instability driven mixing option
    BaseShearValue: 0.005 # Base viscosity/diffusivity value
    RiCrit: 0.7          # Critical Richardson number
    Exponent: 3.0        # Richardson number exponent
    RiSmoothLoops: 2     # Number of Richardson number smoothing loops
```

The parameters under `Convective` and `Shear` are the constants
$\kappa_{conv}$, $N^2_{crit}$, $\nu_o$, $Ri_{crit}$ and $p$ of
{ref}`omega-tech-vertical-mixing`.
