# src

Source code for `DeltaRCWA`

## Contents
- `DeltaRCWA.jl`: main module
- `DeltaRCWAProblem`: struct definitions for `DeltaRCWAProblem`s,
`DeltaRCWASolution`s, and the `solve` method for both
- `modes.jl`: definition of `PlanewaveModes`, which computes the discretization
of real space and momentum space
- `plotting.jl`: plotting recipes for `DeltaRCWASolution` types
- `sheets.jl`: methods to obtain scattering matrices from metasurface impedance
parameters in 1D unit cells and 2D unit cells (in progress)
- `slabs.jl`: methods to obtain scattering matrices from slabs of uniform materials
- `star_products.jl`: methods to compute star products