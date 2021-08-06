export DeltaRCWAProblem, DeltaRCWASolution, solve

"""
    DeltaRCWAProblem(
        structure::SheetStack{N, L, D} where {L, D},
        modes::PlanewaveModes{N},
        pol::P where P <: AbstractPolarization,
        I₁::AbstractArray{<: Number, N},
        I₂::AbstractArray{<: Number, N},
    ) where {T <: Number, N}

This type of problem assumes that the only type of scatterer are sheets and that
the solution can be expressed entirely by one set of modes
"""
struct DeltaRCWAProblem{N}
    structure::RCWAScatterer{N, M} where M
    modes::PlanewaveModes{N}
    pol::AbstractPolarization
    I₁::Array{ComplexF64, N}
    I₂::Array{ComplexF64, N}
    function DeltaRCWAProblem(
        structure::RCWAScatterer{N, M} where M,
        modes::PlanewaveModes{N},
        pol::P where P <: AbstractPolarization,
        I₁::T where T <: AbstractArray{<: Number, N},
        I₂::S where S <: AbstractArray{<: Number, N},
    ) where N
        @assert all(size(weights) == Tuple(e[1] for e in modes.dims) for weights in (I₁, I₂))
        enforce_N_pol(N, pol)
        new{N}(structure, modes, pol,
            convert(Array{ComplexF64}, I₁),
            convert(Array{ComplexF64}, I₂),
        )
    end
end

struct DeltaRCWASolution{N}
    modes::PlanewaveModes{N}
    pol::AbstractPolarization
    I₁::Array{ComplexF64, N}
    I₂::Array{ComplexF64, N}
    O₁::Array{ComplexF64, N}
    O₂::Array{ComplexF64, N}
    function DeltaRCWASolution(
        modes::PlanewaveModes{N},
        pol::AbstractPolarization,
        I₁::AbstractArray{<: Number, N},
        I₂::AbstractArray{<: Number, N},
        O₁::AbstractArray{<: Number, N},
        O₂::AbstractArray{<: Number, N},
    ) where N
        @assert all(size(weights) == Tuple(e[1] for e in modes.dims) for weights in (I₁, I₂, O₁, O₂))
        enforce_N_pol(N, pol)
        new{N}(modes, pol, 
            convert(Array{ComplexF64}, I₁),
            convert(Array{ComplexF64}, I₂),
            convert(Array{ComplexF64}, O₁),
            convert(Array{ComplexF64}, O₂),
        )
    end
end

function solve(prob::DeltaRCWAProblem{N})::DeltaRCWASolution{N} where N
    sol = smatrix(prob.structure, prob.modes, prob.pol) * BlockVector(
        vcat(reshape(prob.I₁, :), reshape(prob.I₂, :)),
        [length(prob.I₁), length(prob.I₂)]
    )
    DeltaRCWASolution(
        prob.modes, prob.pol, prob.I₁, prob.I₂,
        reshape(sol[Block(1)], size(prob.I₁)),
        reshape(sol[Block(2)], size(prob.I₂))
    )
end