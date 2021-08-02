export UniformMedium, Vacuum

"""
    UniformMedium{T <: Number}(ϵ::T, μ::T)

Represents a medium with a constant permitivitty and permeability
Units are such that ϵ, μ are actually relative to vacuum (ϵ=μ=1)
"""
struct UniformMedium{T <: Number}
    ϵ::T
    μ::T
end

Vacuum() = UniformMedium(true, true)