abstract type AbstractScatteringLayer end
abstract type PeriodicScatteringLayer <: AbstractScatteringLayer end

struct UniformLayer <: PeriodicScatteringLayer
    lz::Real
    ϵ::Number
    μ::Number
end

struct XYPeriodicDeltaLayer <: PeriodicScatteringLayer
    lx::Real
    ly::Real
    ϵ::Function
    μ::Function
end
