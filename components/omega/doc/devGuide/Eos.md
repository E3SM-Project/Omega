(omega - dev -
 eos) =

    ##Equation of State

    The Omega `Eos` class A state object is created by the `init` method,
 which assumes that `Decomp` and `HorzMesh` have already been initialized.
```c++ Err = OMEGA::Decomp::init();
Err         = OMEGA::HorzMesh::init();
Err         = OMEGA::Eos::init();
``` The create method :
```c++ Eos::create(const std::string &Name, ///< [in] Name for Eos
                    HorzMesh *Mesh,          ///< [in] Horizontal mesh
                    const int NVertLevels,   ///< [in] Number of vertical levels
);
```
allocates the `SpecVol` and `SpecVolDisplaced` arrays for the mesh and vertical level dimensions.

After initialization, the default state object can be retrieved via:
```
OMEGA::Eos *Eos = OMEGA::Eos::getDefault();
```

    The `Eos` is meant to be a container that calls either the linear or
    Teos - 10 functor dpending on the `EosType` chosen by the user.
