export HCellSheet

struct HCellSheet{T₁, T₂} <: RCWASheet{T₁, 2}
    height::T₁
    center_gap::T₁
    center_bridge::T₁
    leg_width::T₁
    ϵ::T₂
    μ::T₂
end