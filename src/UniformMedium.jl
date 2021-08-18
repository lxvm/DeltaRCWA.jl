export UniformMedium, Vacuum

"""
    UniformMedium{T <: Number}(ϵ::T, μ::T)

Represents a medium with a constant permitivitty and permeability
Units are such that ϵ, μ are relative to the vacuum value(ϵ=μ=1)
"""
struct UniformMedium{T <: Number}
    ϵ::T
    μ::T
end

Vacuum(x=Bool) = UniformMedium(one(x), one(x))