(omega-user-model-horizontal-mixing)=

# Horizontal Mixing

Omega mixes normal velocity and tracers horizontally with Laplacian and
biharmonic terms. Each term is enabled in the `Tendencies` section (see
{ref}`omega-user-model-tendencies`), which also holds its coefficient:

| Term | Parameter | Description |
| ---- | --------- | ----------- |
| VelocityDiffusionOnEdge | ViscDel2 | horizontal viscosity |
| VelocityHyperDiffOnEdge | ViscDel4 | horizontal biharmonic mixing coefficient for normal velocity |
| | DivFactor | scale factor for the divergence term |
| TracerDiffOnCell | EddyDiff2 | horizontal diffusion coefficient |
| TracerHyperDiffOnCell | EddyDiff4 | biharmonic horizontal mixing coefficient for tracers |

## Mesh scaling

Omega scales horizontal viscosity and diffusivity with the local mesh size.
This scaling is enabled by default and is controlled by the following
configuration option:

```yaml
  HorzMesh:
    MeshScalingEnable: true
```

The coefficients are multiplied by a factor that is one at the global median
cell width and grows with the local cell width; see
{ref}`omega-tech-horizontal-mixing` for the definition. The parameters
`ViscDel2`, `ViscDel4`, `EddyDiff2` and `EddyDiff4` therefore specify the
mixing coefficients at the median resolution. Their values may need
adjustment when changing to a mesh with a different median cell width.

On a mesh with equal cell areas, both scaling factors are one. Setting
`MeshScalingEnable` to false also sets both factors to one, giving spatially
uniform mixing coefficients.

The initialization log reports whether mesh scaling is enabled, the median
cell width and the range of both scaling factors over interior edges.

## Converting from MPAS-Ocean

MPAS-Ocean in E3SM references its mixing coefficients to the finest cells of
the mesh, where its scaling factors are one, rather than to the median cell.
To reproduce the MPAS-Ocean mixing on a given mesh, multiply the MPAS-Ocean
Laplacian coefficients by $W_{med}/W_{ref}$ and the biharmonic coefficients by
$(W_{med}/W_{ref})^3$. Here $W_{med}$ is the median cell width and $W_{ref}$
is the nominal finest resolution of the mesh (for example, 30 km for
EC30to60E2r2). The table below lists the converted viscosities for some E3SM
meshes.

| Mesh | $W_{med}$ (km) | $W_{ref}$ (km) | Parameter | Omega | MPAS-Ocean |
| ---- | -------------: | -------------: | --------- | ----: | ---------: |
| EC30to60E2r2 | 35.92 | 30.00 | `ViscDel2` (m$^2$/s) | 1197 | 1000 |
| | | | `ViscDel4` (m$^4$/s) | 2.06e11 | 1.20e11 |
| RRS18to6v3 | 8.77 | 6.00 | `ViscDel4` (m$^4$/s) | 9.99e9 | 3.20e9 |
| RRSwISC6to18E3r5 | 8.48 | 6.00 | `ViscDel4` (m$^4$/s) | 9.03e9 | 3.20e9 |
