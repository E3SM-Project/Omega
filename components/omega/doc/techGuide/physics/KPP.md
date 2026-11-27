(omega-tech-kpp)=

# KPP Boundary Layer Mixing

The K-Profile Parameterization (KPP) of {cite:t}`Large1994` gives the vertical
viscosity and diffusivity in the ocean surface boundary layer, and a
non-local tracer flux. Omega's implementation follows that of MPAS-Ocean and
is in `src/ocn/KPPMix.cpp`; its use is described in {ref}`omega-dev-kppmix`
and its options in {ref}`omega-user-model-kpp`.

## What KPP computes

KPP computes:

- Ocean surface boundary layer depth (`OSBLDepth`)
- Vertical viscosity (`VertVisc`)
- Vertical diffusivity (`VertDiff`)
- Non-local tracer flux profile (`VertNonLocalFlux`)

It uses a bulk Richardson depth search followed by profile-based coefficient
construction. The non-local flux is a standard part of the KPP formulation and
is included by default.

## Notation

KPP splits the water column at the ocean surface boundary layer (OSBL) depth `h`. Inside the OSBL, diffusivity and viscosity are prescribed via a cubic shape function
scaled by `h` and a turbulent velocity scale; below it, only interior mixing
 applies (e.g., shear instability driven mixing).

| Symbol | Code name | Units | Meaning |
| --- | --- | --- | --- |
| `h` | `HOSBL`, `OSBLDepth` | m | OSBL depth below the free surface |
| `d` | `ZDepth`, `ZCenter` | m | depth below the free surface |
| `sigma` | `Sigma` | - | normalized depth in `[-1,0]`, Omega sign convention: 0 at the surface, -1 at the OBL base |
| `sigma_mu` | `SigmaMu` | - | `-Sigma`, in `[0,1]`; the CVMix/Large et al. convention |
| `u*` | `UStar` | m/s | surface friction velocity from wind stress |
| `B_0` | `BuoyFlux` | m^2/s^3 | surface buoyancy flux; **negative is destabilizing** (convection) |
| `L` | `LMoninObukhov` | m | Monin-Obukhov length, `u*^3 / (kappa B_0)` |
| `zeta` | `Zeta` | - | stability coordinate `d/L`; negative is unstable |
| `kappa` | `VonKar` | - | von Karman constant |
| `epsilon` | `SurfaceLayerExtent` | - | surface layer as a fraction of `h` |
| `w_m`, `w_s` | `WMTurb`, `WSTurb` | m/s | turbulent velocity scales for momentum and scalars |
| `Vt^2` | `Vt2`, `UnresolvedShear` | m^2/s^2 | unresolved turbulent shear, {cite:t}`Large1994` Eq. 23 |
| `Ri_b` | `RiBulk`, `BulkRichardsonNumber` | - | bulk Richardson number, Eq. 21 |
| `G(sigma)` | `kppShape*` | - | non-dimensional profile shapes |
| `gamma` | `VertNonLocalFlux` | - | non-local (counter-gradient) tracer flux coefficient |

The OSBL depth is the shallowest `d` at which `Ri_b(d)` reaches the critical
value; the search is done per cell in `computeOSBLDepth`, and the crossing depth
is refined by a quadratic fit through the three nearest cell-center `Ri_b`
values.

(omega-tech-kpp-shape)=

## Shape and stability functions

The non-dimensional functions live in `src/ocn/KPPConstants.h` in namespace
`OMEGA::KPP`. All shape functions take `Sigma` in `[-1,0]` and convert to
`SigmaMu` internally, so callers never flip signs:

- `kppShapeMomentum`, `kppShapeScalar` -- gradient shapes for viscosity and
  diffusivity (`SimpleShapes`)
- `kppShapeMatched` -- gradient shape that additionally reaches a prescribed
  value at the OBL base, used by `MatchBoth`
- the non-local flux has no shape function of its own; it reuses whichever
  scalar shape is in effect, scaled by `C_s` instead of `h*w_s`
- `kppPhiInvMomentum`, `kppPhiInvScalar` -- **inverse** Monin-Obukhov stability
  functions `phi^-1(zeta)`
