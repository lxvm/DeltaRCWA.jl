### TODO: Design problem
# how do you specify the incident modes effectively?
# how do you unambiguously iterate over the product space of modes?
# just do it consistently with the users order, or use named tuples

export get_kz, PlanewaveModes, WeightedModes

"""
    get_kz(k⃗, ω, [M::UniformMedium])

Returns an array of kz from the dispersion relation k⋅k = ω²ϵμ where without
a medium, the vacuum values of ϵ = μ = 1 are used.

Arguments:
k⃗ :: collection of collections (e.g. NTuple{N, Frequencies})
ω :: Number
"""
get_kz(k⃗, ω) = _get_kz(k⃗, ω^2)
get_kz(k⃗, ω, M::UniformMedium) = _get_kz(k⃗, (ω^2) * (M.ϵ * M.μ))

function _get_kz(k⃗, k⃗²)
    map(e⃗ -> sqrt(Complex(k⃗² - mapreduce(abs2, +, e⃗, init=zero(k⃗²)))), Iterators.product(k⃗...))
end

struct PlanewaveModes{T, N}
    ω::Float64
    M::UniformMedium{T}
    dims::NTuple{N, Tuple{Int64, Float64}}
    x⃗::NTuple{N, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}
    k⃗::NTuple{N, Frequencies{Float64}}
    kz::Array{ComplexF64, N}
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
    x⃗ = Tuple(range(0, step=e[2] / e[1], length=e[1]) for e in dims)
    k⃗ = Tuple(2π * fftfreq(e[1], e[1] / e[2]) for e in dims)
    kz = get_kz(k⃗, ω, M)
    # is_propagating = BitArray(@. isreal(kz) & !(iszero(kz)))
    PlanewaveModes(ω, M, dims, x⃗, k⃗, kz)
end


"""
Use this struct to define weights over the distribution of modes.
Make sure that the fastest-changing index of the weight array
(Julia is column-major) corresponds to the first set of Frequencies in modes.k⃗
"""
struct WeightedModes{T, N}
    weights::Array{ComplexF64, N}
    modes::PlanewaveModes{T, N}
    function WeightedModes(weights::Array{ComplexF64, N}, modes::PlanewaveModes{T, N}) where {T, N}
        @assert size(weights) == Tuple(e[1] for e in modes.dims)
        new{T, N}(weights, modes)
    end
end

function WeightedModes(weights::AbstractArray{<: Number, N}, modes::PlanewaveModes{T, N}) where {T, N}
    WeightedModes(convert(Array{ComplexF64}, weights), modes)
end