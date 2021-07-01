module DeltaRCWA

using Reexport

@reexport using LinearAlgebra

@reexport using BlockArrays

export smat_star, rred_star, get_all_modes, get_prop_modes, get_transmissivity_normalincident

include("smat_star.jl")
include("rred_star.jl")
include("get_all_modes.jl")
include("get_prop_modes.jl")
include("get_transmissivity_normalincident.jl")

end
