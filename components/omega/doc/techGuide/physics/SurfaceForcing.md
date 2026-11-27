(omega-tech-surface-forcing)=

# Surface Forcing

How Omega applies the shortwave heating it receives at the surface, and the
energy correction it reports to the coupler. The forcing options are described
in {ref}`omega-user-model-forcing`, and the implementation in
{ref}`omega-dev-forcing`; both terms are in `src/ocn/TendencyTerms.cpp`.

(omega-tech-surface-forcing-shortwave)=

## Penetrating shortwave radiation

Penetrating shortwave radiation distributes incoming solar heat through the
active water column instead of applying all of it to the surface layer. The
incoming flux $I_0$ (`ShortWaveHeatFluxCell`) is split into near-infrared, red
and blue/green bands, each attenuated exponentially with depth. For a cell with surface depth coordinate $z=0$, the downward flux at depth $z$
is defined following {cite:t}`Manizza2005` as

$$
I(z) = I_0 \left[
  f_{nir} e^{-k_{nir}z} +
  f_{red} e^{-k_r z} +
  f_{blue} e^{-k_b z}
\right].
$$

Here $f_{nir}$, $k_{nir}$, $f_{red}$, and $f_{blue}$ are configured with
`Omega.Tendencies.PenetratingShortwaveTendency.NearIrFraction`,
`Omega.Tendencies.PenetratingShortwaveTendency.NearIrCoeff`,
`Omega.Tendencies.PenetratingShortwaveTendency.RedFraction`, and
`Omega.Tendencies.PenetratingShortwaveTendency.BlueFraction`,
respectively. Their default values are $0.58$, $2.86\ \mathrm{m}^{-1}$, $0.21$,
and $0.21$, and their sum $f_{nir} + f_{red} + f_{blue}$ must equal 1.

For each active layer, the deposited heat is the flux at its upper interface
minus the flux at its lower interface. For the bottom active layer, the lower
interface flux is set to zero. This residual-flux rule makes the vertically
integrated deposited heat equal to `ShortWaveHeatFluxCell` even when the
analytic attenuation profile has not reached zero at the model bottom.

(omega-tech-surface-forcing-extra-energy)=

## Coupled energy correction

The extra energy flux associated with freshwater mass exchange is included
in the ocean's surface heat budget. The coupled system accounts for the
corresponding energy as an atmospheric adjustment. This correction is not an
additional ocean temperature tendency or a total heat flux.

$$
Q_{\mathrm{extra}} =
(F_{\mathrm{rain}} + F_{\mathrm{river}}) h_{\mathrm{fw,in}}
+ F_{\mathrm{evap}} h_{\mathrm{fw,out}}
+ (F_{\mathrm{snow}} + F_{\mathrm{iceRunoff}})
  (h_{\mathrm{ice}} + L_{\mathrm{ice}}).
$$

Here the mass fluxes retain their forcing-array signs, and the specific
potential enthalpies are in J kg$^{-1}$:

- $h_{\mathrm{fw,in}} = c^0_{p,sw}\max(C_{T,\mathrm{lim}}, C_{T,\mathrm{top}})$,
  with `Ct0Fw` as the limit for TEOS-10 and zero otherwise.
- $h_{\mathrm{fw,out}} = c^0_{p,sw} C_{T,\mathrm{top}}$.
- $h_{\mathrm{ice}} = -L_{\mathrm{ice}}$ in the current approximation, so the
  snow and ice-runoff contribution to this correction is currently zero.


Direct heat fluxes and sea-ice meltwater enthalpy are explicitly accounted for
in the coupling fluxes and are not added to this correction.
