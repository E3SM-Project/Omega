(omega-tech-frazil)=

# Frazil

Frazil ice forms in the water column wherever the water becomes colder than its local freezing temperature (aka super-cooled water). The implementation in Omega checks each column for the presence of such super-cooled water. Where it is found, the excess heat deficit is converted into frazil ice: the layer gains heat (and so warms back towards the freezing point), loses mass to the ice, and becomes saltier because frazil has a lower average salinity than seawater. Because frazil ice is buoyant, it is assumed to rise, so an ocean column is processed from the bottom upwards and ice formed at depth is added to a reservoir to interact with the layers above. If one of those layers is warmer than its freezing point, some or all of the rising ice melts there, cooling and freshening that layer. Only the net frazil mass, salt, and energy (i.e. the frazil that survives to the top) are passed to the coupler (or discarded in standalone mode). The optional `DepthLimit` parameter confines this process to the upper part of the column, to avoid instability or unrealistically deep frazil formation. Although the frazil code is designed to bring super-cooled water to its freezing temperature instantly , its implementation inside multi-phase time steppers means super-cooled water is eliminated gradually over several time steps. The details of the salt, mass and energy changes with the frazil formation (and its melt) is dependent on the option chosen (`FixedPropery` vs `Teos`).

Frazil is implemented in `src/ocn/Frazil.cpp`, described in
{ref}`omega-dev-frazil`; its options are in {ref}`omega-user-model-frazil`.

## Common behavior

- Freezing-point checks are based on conservative temperature and absolute
  salinity with pressure-dependent freezing temperature.
- Vertical accumulation order is bottom-to-top within each active column.
- Frazil tendencies are not time-step scaled inside `Frazil`; they are
  accumulated as tendency contributions.

## Fixed-property pathway (`FrazilType: FixedProperty`)

- Uses simplified energetics: the energy of the super-cooled water sets the
amount of solid ice formed (used constant latent heat of fusion of fresh ice).
Salt is added based on a constant bulk salinity `IceRefSal` (default) or a
manual toggle (for now) using the local salinity and the frazil porosity. The salt contribution is included in the salinity tendency term, but is not accounted for in the pseudo-thickness tendency.
The fraction of existing frazil to be melted is set by the amount of pure ice that can be melted by the warm layer, and used to determine the energy, mass and salt added to the tendencies. There is also fractional thickness limit to the possible melt.
- Computes local layer tendencies (`HTend`, `TTend`, `STend`) and updates
  accumulated frazil stores.
- Converts accumulators to coupler units at the end of the column loop.

## TEOS-10 pathway (`FrazilType: Teos10`)

- Uses teos-10 Gibbs SeaWater routines for frazil formation and melt state
  transitions.
- Formation uses `Phi` and `LayerMassFracMax` to partition and limit newly formed
  frazil contributions.
- Melt computes fraction melted subject to available thermodynamic energy and
  mass-limit constraints.
