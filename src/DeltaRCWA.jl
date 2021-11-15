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

export smatrix
smatrix(::Type{Matrix}, args...) = _sMatrix(args...)
smatrix(::Type{BlockMatrix}, args...) = _sBlockMatrix(args...)
smatrix(::Type{LinearMap}, args...) = _sLinearMap(args...)
smatrix(::Type{BlockMap}, args...) = _sBlockMap(args...)

export UniformMedium, Vacuum

"""
    UniformMedium{T <: Number}(ϵ::T, μ::T)

Represents a medium with a constant permitivitty and permeability
Units are such that ϵ, μ are relative to the vacuum value(ϵ=μ=1)
"""
struct UniformMedium{T <: Number}
    ϵ::T
    μ::T
end

Vacuum(x=Bool) = UniformMedium(one(x), one(x))

export TE, TM, Coupled

abstract type AbstractPolarization end
abstract type UncoupledPolarization <: AbstractPolarization end
abstract type CoupledPolarization <: AbstractPolarization end

struct TE <: UncoupledPolarization end
struct TM <: UncoupledPolarization end
struct Coupled <: CoupledPolarization end

"""
    enforce_N_pol(N::Integer, ::AbstractPolarization)

Check that a polarization is specified when there are symmetries which decouple
the TE and TM polarizations
"""
function enforce_N_pol(N::Integer, pol::AbstractPolarization)
    if 0 <= N <= 1
        @assert pol isa UncoupledPolarization "1D, 2D photonic crystals uncouple TE, TM"
    elseif N == 2
        @assert pol isa CoupledPolarization "3D photonic crystals couple TE, TM"
    else
        error("only 1D, 2D, 3D photonic crystals are possible")
    end
end

include("modes.jl")
include("slabs.jl")
include("sheets.jl")
include("stacks.jl")
include("DeltaRCWAProblem.jl")
include("star_products.jl")
include("plotting.jl")

end