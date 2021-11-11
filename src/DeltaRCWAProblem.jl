export DeltaRCWAProblem, solve

"""
    DeltaRCWAProblem(
        structure::RCWAScatterer{N},
        modes::PlanewaveModes{T, N},
        pol::P where P <: AbstractPolarization,
        I₁::AbstractArray{<:Number, N},
        I₂::AbstractArray{<:Number, N},
    ) where {T, N}

This type of problem assumes that the only type of scatterer are sheets and that
the solution can be expressed entirely by one set of modes
"""
struct DeltaRCWAProblem{N, T₁<:RCWAScatterer{N}, T₂}
    structure::T₁
    modes::PlanewaveModes{T₂, N}
    pol::AbstractPolarization
    I₁::Array{ComplexF64, N}
    I₂::Array{ComplexF64, N}
    function DeltaRCWAProblem(
        structure::T₁,
        modes::PlanewaveModes{T₂, N},
        pol::AbstractPolarization,
        I₁::AbstractArray{<:Number, N},
        I₂::AbstractArray{<:Number, N},
    ) where {N, T₁, T₂}
        check_mode_size(N, modes.dims, I₁, I₂)
        enforce_N_pol(N, pol)
        new{N, T₁, T₂}(structure, modes, pol, I₁, I₂)
    end
end

function check_mode_size(N, dims, weights...)
    for w in weights
        @assert size(w) == Tuple(i == length(dims) ? N * dims[i][1] : dims[i][1] for i in eachindex(dims)) ""
    end
end

"""
    DeltaRCWAProblem(structure, dims, ω, pol, I₁, I₂, medium=Vacuum())

Constructor
"""
function DeltaRCWAProblem(structure, dims, ω, pol, I₁, I₂, medium=Vacuum())
    DeltaRCWAProblem(structure, PlanewaveModes(ω, dims, medium), pol, I₁, I₂)
end

struct DeltaRCWASolution{N, T}
    modes::PlanewaveModes{T, N}
    pol::AbstractPolarization
    I₁::Array{ComplexF64, N}
    I₂::Array{ComplexF64, N}
    O₁::Array{ComplexF64, N}
    O₂::Array{ComplexF64, N}
    function DeltaRCWASolution(
        modes::PlanewaveModes{T, N},
        pol::AbstractPolarization,
        I₁::AbstractArray{<:Number, N},
        I₂::AbstractArray{<:Number, N},
        O₁::AbstractArray{<:Number, N},
        O₂::AbstractArray{<:Number, N},
    ) where {T, N}
        check_mode_size(N, modes.dims, I₁, I₂, O₁, O₂)
        enforce_N_pol(N, pol)
        new{N, T}(modes, pol, I₁, I₂, O₁, O₂)
    end
end

"""
    solve(::DeltaRCWAProblem; T=Matrix)::DeltaRCWASolution

Evaluates the solution to the problem, with the solver method to be selected by
the user. `T` can be `Matrix`, `BlockMatrix`, `LinearMap`, or `BlockMap` and
selects the underlying representation of the scattering matrices.
"""
function solve(prob::DeltaRCWAProblem{N, T₁, T₂}; T=Matrix)::DeltaRCWASolution{N, T₂} where {N, T₁, T₂}
    structure, modes, pol, I₁, I₂ = prob.structure, prob.modes, prob.pol, prob.I₁, prob.I₂
    sol = smatrix(T, structure, modes, pol) * vcat(reshape.([I₁, I₂], :)...)
    O₁ = reshape(sol[1:length(I₁)], size(I₁))
    O₂ = reshape(sol[(length(I₁)+1):sum(length.([I₁, I₂]))], size(I₂))
    DeltaRCWASolution(modes, pol, I₁, I₂, O₁, O₂)
end