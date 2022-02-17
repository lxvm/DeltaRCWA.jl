"""
    DeltaRCWA

An approximate RCWA-based solver for the scattering of electromagnetic fields in
layered, periodic metasurfaces via a simplified model with infinitesimal sheets.

Developed by Luke Qi, Lorenzo Van Muñoz at Steven G. Johnson's group in 2021
"""
module DeltaRCWA

using LinearAlgebra: I, inv, Diagonal

using BlockArrays
using FFTW: fftfreq, fft, ifft, bfft, Frequencies
using RecipesBase
using KrylovKit: linsolve
using IterativeSolvers: gmres
using LinearMaps: LinearMap, BlockMap
using ChainRulesCore

export smatrix
"""
    smatrix([T=Matrix], scatterer, modes, pol)

Form the scattering matrix of a given return type `T` for the specified
scatterers, incident modes, and polarizations.
"""
smatrix(args...) = smatrix(Matrix, args...)
smatrix(::Type{Matrix}, args...) = _sMatrix(args...)
smatrix(::Type{BlockMatrix}, args...) = _sBlockMatrix(args...)
smatrix(::Type{LinearMap}, args...) = _sLinearMap(args...)
smatrix(::Type{BlockMap}, args...) = _sBlockMap(args...)

include("modes.jl")
include("media.jl")
include("slabs.jl")
include("gstcs.jl")
include("sheets.jl")
include("stacks.jl")
include("DeltaRCWAProblem.jl")
include("star_products.jl")
include("plotting.jl")
include("adjoints.jl")

end