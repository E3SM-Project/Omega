(omega-design-frazil)=
# Frazil

## 1 Overview

The Frazil module computes the frazil tendencies in the ocean column from the local temperature, salinity, pressure, and layer-thickness state. It is used as a tendency term in the OMEGA model and accumulates completed-step totals for export to the coupler. The module is organized as a registry with a default instance and optional non-default instances, and it exposes per-column per-timestep totals for accumulated freshwater mass, salt mass, and energy.

The implementation is intended to be:
- local to a single ocean column,
- consistent with the thermodynamic assumptions in the EOS,
- conservative with respect to water mass, salt mass and energy
- compatible with the existing tendency infrastructure and enable flag.


## 2 Requirements

### 2.1 Requirement: Compute frazil tendencies from formation and melt

The module must compute the total pseudo-thickness (aka mass), salt and conservative temperature tendency contribution from frazil formation and melt as a function of the water column state.
The calculation must allow for tendency contributions from
- the formation of frazil ice from supercooled seawater (leading to the removal of the supercooled water)
- melt of existing frazil if frazil is formed below encounters a layer warm enough to melt it.
- salt and energy contributions from the above processes.


### 2.2 Requirement: Preserve conservation and physical consistency

The implementation must maintain column-level conservation and it must reject or flag situations where the computed mass, salt or energy budgets cannot be reconciled.

Incorrect mass or energy partitioning can corrupt the tracer and thickness budgets and lead to a mass, salt, or energy leak in the coupled model.

### 2.3 Requirement: Support configuration-based enable/disable

The frazil capability must be enabled by configuration. When disabled, the default frazil object must not be created and the tendency contribution must be skipped.

This requirement is needed because the model should be able to run without a frazil parameterization without introducing unresolved dependencies or unnecessary computation.

### 2.4 Requirement: Provide per-cell per-timestep accumulated diagnostics

The module must track column-accumulated totals for:
- frazil freshwater mass,
- salt mass,
- energy associated with the mass fluxes

These totals are used to compute completed-step contributions and are needed for coupling to the sea ice model.

### 2.5 Requirement: calculate frazil contributions consistent with the TEOS-10 thermodynamics

Omega uses ocean thermodynamics consistent with TEOS-10. We target to form (and melt) frazil in a way consistent with the choice of the TEOS-10 Equation Of State (EOS).

### 2.6 Requirement: Support IO field registration and output

The accumulated frazil totals should be exposed as fields so they can be selected in output streams. The fields must be attached to the correct arrays and registered under a group.

This is needed because we desire to test the frazil contributions in standalone mode before coupling. This will also be useful for diagnostics of production simulations.

### 2.6 Requirement: Include a thickness limit
For stability reasons, the total mass tendency removed (or added) to a layer should be limited to a maximum fraction of the layer thickness. Ten percent is a threshold used in previous implementations (and the intended default value). The requirement is for there to be a parameter that caps the mass tendency.

### 2.7 Desired: Include a depth limit
In other implementations, deep formation of frazil (e.g. below 1000m) was an issue. Thus we prefer to have a parameter restricting the frazil formation (and melt) to depths shallower than this parameter.

### 2.7 Desired: have a simpler frazil scheme
To limit the risk of instability, we desire to implement a frazil scheme to provide a simpler alternative to the TEOS-10 option. As much as possible, this simpler frazil scheme should rely on assumptions that past frazil implementations have made to increase the chances of providing a stable frazil scheme in a coupled model.

### 2.8 Desired: limit computational expense
We desire to limit unnecessary calls and computations. Thus, we add a check on the water column state (is there super cooled water in this column?) before running through the frazil calculation. This is inspired by MPAS-O implementation.

## 3 Algorithmic Formulation

The frazil calculation is implemented as a local column process. The model state provides:
- conservative temperature $T$,
- absolute salinity $S$,
- pressure $P$,
- layer thickness $H$.

For a given active-layer cell, the module computes the formation and melt contributions to:
- ice mass,
- liquid mass,
- salt mass,
- total energy.

The calculation is organized around a set of vertically-accumulated quantities:
$$
M_{ice},\quad M_{liq},\quad M_{salt},\quad E_{ice},\quad E_{liq}
$$
and the tendency terms
$$
\frac{\partial H}{\partial t},\quad
\frac{\partial T}{\partial t},\quad
\frac{\partial S}{\partial t}.
$$

For each stage of a timestep, the calculation resets the (vertical) accumulator arrays, applies the local frazil algorithm, and stores the resulting stage-local contributions in the tendency arrays. At the end of the completed outer-step update, the accumulated totals are normalized by the full ocean timestep:
$$
\bar{F}_{mass} = \frac{\sum_k w_k M_{mass,k}}{\Delta t},
\quad
\bar{F}_{salt} = \frac{\sum_k w_k M_{salt,k}}{\Delta t},
\quad
\bar{F}_{energy} = \frac{\sum_k w_k E_k}{\Delta t},
$$
where $w_k$ are the timestepper-specific update weights. This preserves the correct outer-step mean flux while remaining compatible with Runge-Kutta stage logic.
We prefer that the public output arrays captures mean fluxes rather than stage or timestep totals to make the Polaris testing independent of the timestep used.

## 4 Design

### 4.1 Data types and parameters

#### 4.1.1 Parameters

The Frazil configuration is stored under the `Frazil` config group. Relevant entries include:
- `FrazilType`: fixed-property or TEOS-based formulation
- `MassLimit`: maximum total frazil mass changed allowed in an ocean layer
- `Phi`: target liquid mass fraction for the formed frazil
- `DepthLimit`: optional depth limit below which frazil is disabled
- `ConservationCheck`: whether conservation must be validated

These parameters are read during object creation in [add-basic-frazil/components/omega/src/ocn/Frazil.cpp](add-basic-frazil/components/omega/src/ocn/Frazil.cpp).

#### 4.1.2 Class/structs/data types

The relevant public structures in the Frazil implementation are:

```c++
enum class FrazilType {
   FixedPropertyFrazil,
   TeosFrazil
};

struct FrazilFormation { ... };
struct FrazilMelt { ... };
struct FixedPropertyFrazilFormation { ... };
struct FixedPropertyFrazilMelt { ... };

class Frazil {
 public:
   static void init();
   static Frazil *getDefault();
   static Frazil *get(const std::string &Name);
   static void clear();
   static void erase(std::string InName);

   Array2DReal FrazilTTend;
   Array2DReal FrazilSTend;
   Array2DReal FrazilHTend;

   Array1DReal AccMIce;
   Array1DReal AccEIce;
   Array1DReal AccMLiq;
   Array1DReal AccELiq;
   Array1DReal AccMSalt;

   Array1DReal OcnDtFrazilMass;
   Array1DReal OcnDtFrazilSalt;
   Array1DReal OcnDtFrazilEnergy;

   void computeFrazil(...);
   void resetOcnStepRates();
   void accumulateOcnStepRates(Real FinalUpdateWeight, R8 TimeStepSeconds);
   void registerFields();
   void unregisterFields();
};
```

The design keeps the arrays and accumulation logic local to the Frazil object, which makes the lifetime management explicit and avoids leaking the data model into time-stepping logic.

### 4.2 Methods

#### 4.2.1 Initialization

The static `Frazil::init()` method:
- ensures the mesh and vertical coordinate exist,
- initializes EOS dependencies,
- checks the `Tendencies` config for `FrazilTendencyEnable`,
- creates the default Frazil instance when enabled.

This ensures the object exists only when the feature is active.

#### 4.2.2 Creation and retrieval

The creation pattern mirrors the rest of the Omega module structure:
- a map of all Frazil objects,
- a pointer to the default instance,
- a name-based retrieval method,
- methods for erase and clear.

This design allows future non-default instances while keeping the default case fast and simple.

#### 4.2.3 Computation

The `computeFrazil()` method is the core thermodynamic routine. It:
- reads the local tracer state,
- evaluates the local freezing point,
- computes mass and energy changes associated with frazil formation or melt,
- writes the results into the layer-by-layer tendency arrays.

The method is later invoked from the tendency operator, which is the integration point between the frazil module and the model tendency machinery.

#### 4.2.4 Completed-step accumulation

The `accumulateOcnStepRates()` method updates the final per-cell ocean-step sums. It:
- reads the stage-local accumulation arrays,
- applies the update weight,
- divides by the full ocean timestep seconds,
- writes the mean flux into the public output arrays.

This is a key design choice because it keeps the output/coupling values independent of the internal stage ordering of the timestepper.

#### 4.2.5 Field registration and output

The `registerFields()` method creates the three default output fields:
- `OcnDtFrazilMass`
- `OcnDtFrazilSalt`
- `OcnDtFrazilEnergy`

Each field is:
- 1D with `NCells` dimension,
- attached to the corresponding array,
- inserted into the `Frazil` field group.

The `unregisterFields()` method is guarded by a registration flag so only the object that successfully registered the fields destroys them. This avoids collisions when non-default instances are later introduced.

#### 4.2.6 Tendency hook

The `FrazilOnCell` object is the bridge between the Frazil module and the tendency infrastructure. It:
- checks whether frazil is enabled,
- resets the stage-local accumulators,
- calls `computeFrazil()`,
- accumulates the completed-step totals,
- writes the resultant tendency contributions into the thickness and tracer tendency arrays.

This is the integration point used by [add-basic-frazil/components/omega/src/ocn/TendencyTerms.cpp](add-basic-frazil/components/omega/src/ocn/TendencyTerms.cpp).

## 5 Verification and Testing

### 5.1 Test cold formation

The cold-case frazil formation test verifies that:
- accumulated ice mass is positive,
- accumulated salt mass is positive,
- accumulated energy is negative,
- temperature tendency is positive,
- salinity and thickness tendencies are negative.

This verifies the sign and directionality of the formation branch.

### 5.2 Test warm formation

The warm-case formation test verifies that the formation branch is zero when water is not supercooled. This checks that the algorithm respects the physical cutoff condition.

### 5.3 Test fixed-property formulation

The fixed-property frazil test exercises the simplified formulation and validates the same sign and physical expectations without requiring the full TEOS-based path.

### 5.4 Test column conservation

The column test verifies that the computation remains consistent across multiple active layers and that the integrated frazil contribution satisfies the expected conservation behavior with the enabled check.

### 5.5 Test depth limit

The depth-limit test ensures that any layers below the configured depth threshold do not generate frazil tendencies, which guards the parameterized cutoff behavior.

### 5.6 Test mass limit

A future test should check that the mass limit is implemented correctly.

### 5.6 Test per-tilmestep conservation

A future test should check that the ocean tendencies and the terms ready for coupling export are conserved. The above test 5.4 ensures conservation within the frazil call but does not check the full timestep totals, which will differ across time stepper choice.
