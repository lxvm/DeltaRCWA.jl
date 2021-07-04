abstract type AbstractScatteringStructure{T <: AbstractScatteringLayer} <: AbstractVector{T} end
abstract type PeriodicScatteringStructure{T <: PeriodicScatteringLayer} <: AbstractScatteringStructure{T} end

mutable struct XYPeriodicScatteringStructure{T} <: PeriodicScatteringStructure{T}
    lx::Real
    ly::Real

end