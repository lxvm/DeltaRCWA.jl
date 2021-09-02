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
using IterativeSolvers: gmres
using LinearMaps: LinearMap, BlockMap

export RCWAScatterer, RCWASheet, RCWASlab, RCWAInterface, RCWAStack

"""
    RCWAScatterer{T <: Number, N}

Supertype for structures with a unit cell having `N` periodic dimensions
(determines the dimensionality of the problem) whose geometry is parametrized
by parameters of type `T`.
`N=0` is for a 1D photonic crystal, `N=1` is for 2D and so on. Note that `N<=2`.
"""
abstract type RCWAScatterer{N} end

# """
#     RCWAInterface{T, N} <: RCWAScatterer{T, N}

# An abstract type to dispatch methods for 2D structures at the interface of two
# homogenous media modelled with Generalized Sheet Transition Conditions (GSTCs).
# See:
# Kuester, Edward et al.
# "A generalized sheet transition condition model for a metafilm partially embedded in an interface"
# https://doi.org/10.1109/APS.2010.5562250
# """
# abstract type RCWAInterface{T, N} <: RCWAScatterer{T, N} end

"""
    RCWASheet{N} <: RCWAScatterer{N}

An abstract type to dispatch methods for 2D structures embedded in a homogenous
medium modelled with Generalized Sheet Transition Conditions (GSTCs).
See:
Kuester, Edward et al.
"Averaged transition conditions for electromagnetic fields at a metafilm"
https://doi.org/10.1109/TAP.2003.817560
"""
abstract type RCWASheet{N} <: RCWAScatterer{N} end

"""
    RCWASlab{N} <: RCWAScatterer{N}

An abstract type to dispatch methods for 3D structures, such as uniform media.
General methods to calculate the eigenmodes of the wave equation in variable
impedance slabs could be implemented to make this software a full RCWA solver.
See:
Liu, Victor et al.
"S4 : A free electromagnetic solver for layered periodic structures"
https://doi.org/10.1016/j.cpc.2012.04.026
"""
abstract type RCWASlab{N} <: RCWAScatterer{N} end

"""
    RCWAStack{N} <: RCWAScatterer{N}

An abstract type to dispatch methods for stacks of individual scatterers that
can be concatenated with a Redheffer star product.
See:
Rumpf, Raymond
"Improved Formulation of Scattering Matrices for Semi-Analytical Methods That Is Consistent with Convention"
https://doi.org/10.2528/PIERB11083107
"""
abstract type RCWAStack{N} <: RCWAScatterer{N} end

export smatrix, smatrixBlockMap, smatrixLinearMap

# Not to be confused with SMatrix from StaticArrays.jl!
# This is a function to get a scattering matrix from a scattering layers

"""
    smatrix(::RCWAScatterer, modes, pol)

Returns the scattering matrix of a sheet acting on the propagating plane wave modes.
The default method returns the identity scattering matrix.
"""
function smatrix(::RCWAScatterer, modes, pol)::BlockMatrix
    N = length(modes.kz)
    id = I(N)
    zo = (Bool(0) * I)(N)
    mortar(
        (zo, id),
        (id, zo),
    )
end

function smatrixBlockMap(::RCWAScatterer, modes, pol)::BlockMap
    N = length(modes.kz)
    id = LinearMap(I(N))
    zo = LinearMap((Bool(0) * I)(N))
    [
        zo  id;
        id  zo
    ]
end

function smatrixLinearMap(::RCWAScatterer, modes, pol)::LinearMap
    N = length(modes.kz)
    LinearMap(x::AbstractVector -> vcat(x[(N+1):2N], x[1:N]), 2N)
end

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
include("interfaces.jl")
include("stacks.jl")
include("DeltaRCWAProblem.jl")
include("star_products.jl")
include("plotting.jl")

include("Luke_functions.jl")

end