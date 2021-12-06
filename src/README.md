# src

Source code for `DeltaRCWA`

## Contents
- `DeltaRCWA.jl`: main module
- `DeltaRCWAProblem.jl`: struct definitions for `DeltaRCWAProblem`s,
`DeltaRCWASolution`s, and the `solve` method for both
- `modes.jl`: definition of `PlaneWaves`, which computes the discretization
of real space and momentum space
- `media.jl`: definition of `UniformMedium` and helper functions for dispersion relation
- `plotting.jl`: plotting recipes for `DeltaRCWASolution` types
- `sheets.jl`: methods to obtain scattering matrices from metasurface impedance
parameters in 1D unit cells and 2D unit cells (in progress)
- `slabs.jl`: methods to obtain scattering matrices from slabs of uniform
  materials
- `gstcs.jl`: methods to help construct scattering matrices for sheets
- `stacks.jl`: Definition of `SheetStack`, the representation of the physical system
- `star_products.jl`: methods to compute star products