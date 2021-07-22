struct HCellDeltaLayer <: DeltaScatterer
    height::Real
    center_gap::Real
    center_bridge::Real
    leg_width::Real
    ϵhole::Number
    ϵwall::Number
    μhole::Number
    μwall::Number
end