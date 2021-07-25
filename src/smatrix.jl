# Not to be confused with SMatrix from StaticArrays.jl!
# This is a function to get a scattering matrix from a scattering layers
"""
    DeltaSMatrix()

Internally build the dense scattering matrix 
"""
function DeltaSMatrix(vz)
    A = mortar(reshape([
        -I + σₘ * vz,
        -vz + σₑ,
        I - σₘ * vz,
        -vz + σₑ,
    ], 2, 2))
    B = mortar(reshape([
        -I + σₘ * vz,
        -vz + σₑ,
        I - σₘ * vz,
        -vz + σₑ,
    ], 2, 2))
    A\B
end

"""
    smatrix(scatterer::AbstractScatterer, k::NTuple{N, Frequencies}, ω::Float64)

Returns a 2x2 BlockMatrix for the scattering of modes

scatterer::AbstractScatterer
k::NTuple{N, Frequencies}, only implemented for N in [0, 1, 2]
ω::Float64
"""
function smatrix(scatterer::DeltaScatterer, modes::Modes)
    # 
    vᶻ = get_kz(Vacuum(), k, ω) ./ ω
    σᵉ = get_σₑ(scatterer)/2
    σᵐ = get_σₘ(scatterer)/2
    Deltasmatrix()
end

# function smatrix(scatterer::DeltaScatterer, k::Tuple{Frequencies}, ω::Float64)
#     kz = get_kz(Air(0), k, ω)

# end


# function smatrix(scatterer::DeltaScatterer, k::NTuple{2, Frequencies}, ω::Float64)
#     error("NotImplemented")
# end