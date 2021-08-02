module DeltaRCWA

using Reexport: include
using Reexport

@reexport using LinearAlgebra

@reexport using BlockArrays
@reexport using FFTW

# Define the abstract types used in this module

export RCWAScatterer

"""
    RCWAScatterer{N, M}

Supertype for structures with `N` periodic dimensions off the propagation axis
that can be parametrized with `M` coordinate axes.
`N=0` is for a 1D photonic crystal, `N=1` is for 2D and so on.
`M=2` is for a sheet/surface/interface, `M=3` is for a slab/volume/layer.
Note that `M<=3, N<=2`.
"""
abstract type RCWAScatterer{N, M} end
abstract type RCWASheet{N} <: RCWAScatterer{N, 2} end
abstract type RCWASlab{N} <: RCWAScatterer{N, 3} end

const maxNdim = 2 # there are at most 2 periodic tangential dimensions
const minNdim = 0 # all tangential dimensions can be invariant

const maxMdim = 3 # Volumes are the largest dimensional manifolds in real space
const minMdim = 2 # curves and points do not have the symmetries for RCWA

# include code

include("UniformMedium.jl")
include("polarizations.jl")
include("PlanewaveModes.jl")
include("smatrix.jl")
include("RCWAProblem.jl")
include("star_products.jl")
include("concrete_scatterers/_index.jl")

include("Luke_functions.jl")

end
