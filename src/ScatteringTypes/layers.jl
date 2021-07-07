abstract type AbstractScatteringLayer end
abstract type ScatteringLayer <: AbstractScatteringLayer end
abstract type DeltaScatteringLayer <: AbstractScatteringLayer end
abstract type PeriodicScatteringLayer <: AbstractScatteringLayer end

struct UniformLayer <: AbstractScatteringLayer
    depth::Real
    ϵ::Number
    μ::Number
end

Air(depth) = UniformLayer(depth, 1, 1)

struct GaussianDeltaLayer <: DeltaScatteringLayer
    μx::Real
    σx²::Real
    μy::Real
    σy²::Real
end

struct HCellDeltaLayer <: PeriodicScatteringLayer, DeltaScatteringLayer
    height::Real
    center_gap::Real
    center_bridge::Real
    leg_width::Real
end