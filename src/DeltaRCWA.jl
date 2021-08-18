"""
    DeltaRCWA

An approximate RCWA-based solver for the scattering of electromagnetic fields in
layered, periodic metasurfaces via a simplified model with infinitesimal sheets.

Developed by Luke Qi, Lorenzo Van Muñoz at Steven G. Johnson's group in 2021
"""
module DeltaRCWA

using LinearAlgebra: I, inv, Diagonal

using BlockArrays: mortar, Block, BlockVector, BlockMatrix, AbstractBlockMatrix
using FFTW: fftfreq, fft, ifft, bfft, Frequencies
using RecipesBase
using KrylovKit: linsolve

export RCWAScatterer, RCWASheet, RCWASlab

"""
    RCWAScatterer{T <: Number, N}

Supertype for structures with a unit cell having `N` periodic dimensions
(determines the dimensionality of the problem) whose geometry is parametrized
by parameters of type `T`.
`N=0` is for a 1D photonic crystal, `N=1` is for 2D and so on. Note that `N<=2`.
"""
abstract type RCWAScatterer{T <: Number, N} end

"""
    RCWASheet{T, N} <: RCWAScatterer{T, N}

An abstract type to dispatch methods for 2D structures whose scattering
can be modelled with Generalized Sheet Transition Conditions (GSTCs).
"""
abstract type RCWASheet{T, N} <: RCWAScatterer{T, N} end

"""
    RCWASlab{T, N} <: RCWAScatterer{T, N}

An abstract type to dispatch methods for 3D structures, such as uniform media.
General methods to calculate the eigenmodes of the wave equation in variable
impedance slabs could be implemented to make this software a full RCWA solver.
"""
abstract type RCWASlab{T, N} <: RCWAScatterer{T, N} end

include("UniformMedium.jl")
include("polarizations.jl")
include("modes.jl")
include("smatrix.jl")
include("smatrixfree.jl")
include("sheets.jl")
include("concrete_scatterers/_index.jl")
include("DeltaRCWAProblem.jl")
include("star_products.jl")
include("plotting.jl")

include("Luke_functions.jl")

end
