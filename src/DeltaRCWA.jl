module DeltaRCWA

using Base: modules_warned_for, valid_import_path
using Reexport

@reexport using LinearAlgebra

@reexport using BlockArrays
@reexport using FFTW

# Define the abstract types used in this module

export AbstractScatterer, DeltaScatterer

abstract type AbstractScatterer end
abstract type DeltaScatterer <: AbstractScatterer end

export AbstractScatteringProblem, PeriodicScatteringProblem

abstract type AbstractScatteringProblem end
abstract type PeriodicScatteringProblem <: AbstractScatteringProblem end

# include code

include("scatterers/scatterers.jl")

export TE, TM, Coupled
include("polarizations.jl")

export Modes, get_kz
include("modes.jl")

export RCWAProblem, solve
include("RCWAProblem.jl")

export smatrix
include("smatrix.jl")

export smat_star, rred_star
include("star_products.jl")

export get_all_modes, get_prop_modes, get_transmissivity_normalincident
include("Luke_functions.jl")

end
