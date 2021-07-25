abstract type AbstractModes end

struct PlanewaveModes{N, P <: AbstractPolarization, T <: AbstractUniformMedium} <: AbstractModes
    ω::Float64
    k⃗::NTuple{N, Frequencies}
    medium::T
    kz::Array{ComplexF64, N}
    is_propagating::BitArray{N}
end

"""
    RCWAModes(
        ω::Float64,
        dims::NTuple{N, Tuple{Int64, Float64}},
        P::AbstractPolarization,
        [medium::AbstractUniformMedium]
    ) where N

This constructor takes a frequency, and an NTuple of Tuple{Int64, Float64}.
Each tuple represents a periodic dimension, and the Int64 represents the number
of grid points per unit cell in that periodic dimension, and the Float64 gives
the lattice constant in that periodic dimension (length of unit cell).

By default, the medium is Vacuum

The polarization should be <: UncoupledPolarization if 0 <= N <= 1
or <: CoupledPolarization if N == 2
"""
function RCWAModes(
    ω::Float64,
    dims::NTuple{N, Tuple{Int64, Float64}},
    ::P;
    medium::T=Vacuum(),
) where {N, P <: AbstractPolarization, T <: AbstractUniformMedium}
    if 0 <= N <= 1
        @assert P <: UncoupledPolarization "1D, 2D photonic crystals uncouple TE, TM"
    elseif N == 2
        @assert P <: CoupledPolarization "3D photonic crystals couple TE, TM"
    else
        error("only 1D, 2D, 3D photonic crystal are possible")
    end
    k⃗ = Tuple(fftfreq(e[1], e[1] / e[2]) for e in dims)
    kz = get_kz(medium, k⃗, ω)
    is_propagating = isreal.(kz)
    PlanewaveModes{N, P, T}(ω, k⃗, medium, kz, is_propagating)
end
### TODO: Design problem
# how do you specify the incident modes effectively?
# how do you unambiguously iterate over the product space of modes?
# just do it consistently with the users order, or use named tuples

"""
    get_kz(medium::AbstractUniformMedium, k⃗::NTuple{N, Frequencies{T}}, ω::T)::Array{Complex{T}, N} where {N, T}

Returns an array of kz from the dispersion relation k⋅k = ω²ϵμ
"""
function get_kz(medium::AbstractUniformMedium, k⃗::NTuple{N, Frequencies{T}}, ω::T)::Array{Complex{T}, N} where {N, T}
    k⃗² = ω²ϵμ(medium, ω)
	[sqrt(Complex(k⃗² - mapreduce(abs2, +, e⃗, init=0.0))) for e⃗ in Iterators.product(k⃗...)]
end