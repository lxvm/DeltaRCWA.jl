export smatrix

# Not to be confused with SMatrix from StaticArrays.jl!
# This is a function to get a scattering matrix from a scattering layers

"""
    smatrix(::T where T <: RCWAScatterer, modes::PlanewaveModes)

Returns the scattering matrix of a scatterer acting on the propagating plane wave modes.
The default method returns the identity scattering matrix.
"""
function smatrix(::T where T <: RCWAScatterer, modes::PlanewaveModes)
    # total number of propagating modes
    N = sum(modes.is_propagating)
    mortar(reshape([(false*I)(N), I(N), I(N), (false*I)(N)], 2, 2))
end

"""
    smatrix(scatterer::RCWAScatterer, modes::PlanewaveModes)

Returns a 2x2 BlockMatrix for the scattering of modes

scatterer::AbstractScatterer
k::NTuple{N, Frequencies}, only implemented for N in [0, 1, 2]
ω::Float64
"""
function smatrix(scatterer::T where T <: RCWASheet{1}, modes::PlanewaveModes)

end