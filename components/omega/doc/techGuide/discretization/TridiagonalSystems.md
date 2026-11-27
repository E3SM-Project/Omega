(omega-tech-tridiagonal-systems)=

# Tridiagonal Systems

Implicit vertical mixing gives one tridiagonal system of equations per ocean
column. Omega solves them in batches with the solvers in
`src/base/TriDiagSolvers.h`, whose use is described in
{ref}`omega-dev-tridiagonal`. Systems of the diffusion-type form have their
own solver, with better stability properties than the general one.

## General tridiagonal system
A general tridiagonal system is represented by three vectors
`DL`, `D`, `DU`, representing the lower, main, and upper diagonal, respectively, and a right-hand-side vector `X`
$$
(DL)_i Y_{i - 1} + D_i Y_i + (DU)_i Y_{i + 1} = X_i.
$$

(omega-tech-tridiagonal-systems-diffusion)=

## Diffusion-type tridiagonal system

A diffusion-type system is represented by two vectors `G` and `H`, and a right-hand-side vector `X`. It has the form
$$
-G_{i - 1} Y_{i - 1} + (G_{i - 1} + G_i + H_i) Y_i - G_{i} Y_{i + 1} = X_i.
$$
