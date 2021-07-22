struct GaussianDeltaLayer <: DeltaScatterer
    ϵx::Real
    σϵx²::Real
    ϵy::Real
    σϵy²::Real
    μx::Real
    σμx²::Real
    μy::Real
    σμy²::Real
end