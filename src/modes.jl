export PlaneWaves, PolarizationStyle, TE, TM, Coupled

struct PlaneWaves{N, T}
    ω::Float64
    dims::NTuple{N, Tuple{Int64, Float64}}
    x⃗::NTuple{N, T}
    k⃗::NTuple{N, Frequencies{Float64}}
end

"""
    PlaneWaves(ω::Float64, dims::Tuple{Vararg{Tuple{Int64, Float64}}})

This constructor takes a frequency, and an NTuple of Tuple{Int64, Float64}.
Each tuple represents a periodic dimension, and the Int64 represents the number
of grid points per unit cell in that periodic dimension, and the Float64 gives
the lattice constant in that periodic dimension (length of unit cell).
"""
function PlaneWaves(ω::Float64, dims::Tuple{Vararg{Tuple{Int64, Float64}}})
    # v1, needs adjoint constructors for StepRangeLen and Base.TwicePrecision
    # x⃗ = Tuple(range(0, step=e[2] / e[1], length=e[1]) for e in dims)
    # v2, works with Zygote
    x⃗ = Tuple([i * e[2] / e[1] for i in 1:e[1]] for e in dims)
    k⃗ = Tuple(2π * fftfreq(e[1], e[1] / e[2]) for e in dims)
    PlaneWaves(ω, dims, x⃗, k⃗)
end

Base.size(pw::PlaneWaves) = Tuple(e[1] for e in pw.dims)
Base.length(pw::PlaneWaves) = prod(size(pw))

abstract type PolarizationStyle end
struct TE <: PolarizationStyle end
struct TM <: PolarizationStyle end
struct Coupled <: PolarizationStyle end

PolarizationStyle(pw::PlaneWaves) = PolarizationStyle(typeof(pw))
PolarizationStyle(::Type{<:PlaneWaves{1}}) = TM()
PolarizationStyle(::Type{<:PlaneWaves{2}}) = Coupled()