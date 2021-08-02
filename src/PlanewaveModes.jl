### TODO: Design problem
# how do you specify the incident modes effectively?
# how do you unambiguously iterate over the product space of modes?
# just do it consistently with the users order, or use named tuples

export get_kz, PlanewaveModes

"""
    get_kz(k⃗::NTuple{N, Frequencies{T}}, ω::T, [medium::UniformMedium])::Array{Complex{T}, N} where {N, T}

Returns an array of kz from the dispersion relation k⋅k = ω²ϵμ where without
a medium, the vacuum values of ϵ = μ = 1 are used.
"""
get_kz(k⃗, ω) = _get_kz(k⃗, ω^2)
get_kz(k⃗, ω, M::UniformMedium{T} where T) = _get_kz(k⃗, (ω^2) * (M.ϵ * M.μ))

function _get_kz(k⃗::NTuple{N, Frequencies{T}}, k⃗²::T)::Array{Complex{T}, N} where {N, T}
    map(e⃗ -> sqrt(Complex(k⃗² - mapreduce(abs2, +, e⃗, init=zero(T)))), Iterators.product(k⃗...))
    # [sqrt(Complex(k⃗² - mapreduce(abs2, +, e⃗, init=zero(T)))) for e⃗ in Iterators.product(k⃗...)]
end

struct PlanewaveModes{N}
    dims::NTuple{N, Tuple{Int64, Float64}}
    ω::Float64
    k⃗::NTuple{N, Frequencies}
    M::UniformMedium
    kz::Array{ComplexF64, N}
    is_propagating::BitArray{N}
end

"""
    PlanewaveModes(
        ω::Float64,
        dims::NTuple{N, Tuple{Int64, Float64}},
        medium::UniformMedium
    ) where N

This constructor takes a frequency, and an NTuple of Tuple{Int64, Float64}.
Each tuple represents a periodic dimension, and the Int64 represents the number
of grid points per unit cell in that periodic dimension, and the Float64 gives
the lattice constant in that periodic dimension (length of unit cell).
"""
function PlanewaveModes(ω::Float64, dims::NTuple{N, Tuple{Int64, Float64}}, M::UniformMedium) where N
    k⃗ = Tuple(fftfreq(e[1], e[1] / e[2]) for e in dims)
    kz = get_kz(k⃗, ω, M)
    is_propagating = @. isreal(kz) & !(iszero(kz))
    PlanewaveModes{N}(dims, ω, k⃗, M, kz, is_propagating)
end