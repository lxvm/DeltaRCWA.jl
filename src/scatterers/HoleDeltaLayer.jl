struct HoleDeltaLayer <: DeltaScatterer
    radius::Real
    ox::Real
    oy::Real
    ϵhole::Number
    ϵwall::Number
    μhole::Number
    μwall::Number
end