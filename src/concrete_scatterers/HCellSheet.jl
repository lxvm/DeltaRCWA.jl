export HCellSheet

struct HCellSheet <: RCWASheet{2}
    height::Real
    center_gap::Real
    center_bridge::Real
    leg_width::Real
    ϵ::Number
    μ::Number
end