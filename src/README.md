# src

Source code for `DeltaRCWA`

## Contents
- `DeltaRCWA.jl`: main module, abstract type definitions, and fallback methods
- `DeltaRCWAProblem`: struct definitions for `DeltaRCWAProblem`s,
`DeltaRCWASolution`s, and the `solve` method for both
- `interfaces`: an unfinished attempt to implement Fresnel scattering at an interface
- `Luke_functions`: functions written by Luke for 1D unit cell scattering matrices,
which is now superseded by the contents of `sheets.jl`
- `modes.jl`: definition of `PlanewaveModes`, which computes the discretization
of real space and momentum space
- `plotting.jl`: plotting recipes for `DeltaRCWASolution` types
- `sheets.jl`: methods to obtain scattering matrices from metasurface impedance
parameters in 1D unit cells and 2D unit cells (in progress)
- `slabs.jl`: methods to obtain scattering matrices from slabs of uniform materials
- `star_products.jl`: methods to compute star products