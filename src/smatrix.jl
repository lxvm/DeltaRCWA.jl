# Not to be confused with SMatrix from StaticArrays.jl!
# This is a function to get a scattering matrix from a scattering layers
"""
Returns a 2x2 BlockMatrix for the scattering of modes

scatterer::AbstractScatterer
k::NTuple{N, Frequencies}, only implemented for N in [0, 1, 2]
ω::Real
"""
function smatrix(scatterer::AbstractScatterer, k::Tuple{}, ω::Real)
    kz = get_kz(Air(0), k, ω)
    σᵉ = get_σₑ(scatterer)
    σᵐ = get_σₘ(scatterer)
    A = mortar(reshape([
        -1 + σₘ
    ], 2, 2))
    B = mortar(reshape([
        -1 + σₘ
    ], 2, 2))
    return A\B
end

function smatrix(scatterer::AbstractScatterer, k::Tuple{Frequencies}, ω::Real)
    kz = get_kz(Air(0), k, ω)

end


function smatrix(scatterer::AbstractScatterer, k::NTuple{2, Frequencies}, ω::Real)
    error("NotImplemented")
end