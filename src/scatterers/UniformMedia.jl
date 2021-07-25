"""
    AbstractUniformMedium

Types T <: AbstractUniformMedium should implement

    ω²ϵμ(::T, ω::Real)
"""
abstract type AbstractUniformMedium end

"""
    Vacuum()

A type to represent a vacuum, with units such that ϵ=μ=c=1
"""
struct Vacuum <: AbstractUniformMedium end

"""
    UniformMedium{T <: Number}(ϵ::T, μ::T)

Represents a medium with a constant permitivitty and permeability
Units are such that ϵ, μ are actually relative to vacuum (ϵ=μ=1)
"""
struct UniformMedium{T <: Number} <: AbstractUniformMedium
    ϵ::T
    μ::T
end

"""
    ω²ϵμ(medium::AbstractUniformMedium, ω::Float64)

Calculates the RHS of the dispersion relation k⃗² = k⃗⋅k⃗ = ω²ϵμ
"""
ω²ϵμ(::Vacuum, ω::Real) = ω^2
ω²ϵμ(medium::UniformMedium, ω::Real) = (ω^2) * medium.ϵ * medium.μ