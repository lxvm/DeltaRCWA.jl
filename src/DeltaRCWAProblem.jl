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
struct DeltaRCWAProblem{T₁, T₂, N, L}
    structure::SheetStack{T₁, N, L}
    modes::PlanewaveModes{T₂, N}
    pol::AbstractPolarization
    I₁::Array{ComplexF64, N}
    I₂::Array{ComplexF64, N}
    function DeltaRCWAProblem(
        structure::SheetStack{T₁, N, L},
        modes::PlanewaveModes{T₂, N},
        pol::AbstractPolarization,
        I₁::AbstractArray{<: Number, N},
        I₂::AbstractArray{<: Number, N},
    ) where {T₁, T₂, N, L}
        check_mode_size(N, modes.dims, I₁, I₂)
        # @assert all(size(weights) == Tuple(i == length(modes.dims) ? N * modes.dims[i][1] : modes.dims[i][1] for i in eachindex(modes.dims)) for weights in (I₁, I₂))
        enforce_N_pol(N, pol)
        new{T₁, T₂, N, L}(structure, modes, pol, convert.(Array{ComplexF64}, (I₁, I₂))...)
    end
end

function check_mode_size(N, dims, weights...)
    for w in weights
        @assert size(w) == Tuple(i == length(dims) ? N * dims[i][1] : dims[i][1] for i in eachindex(dims)) ""
    end
end

function DeltaRCWAProblem(sheet::RCWASheet, modes, pol, I₁, I₂)
    DeltaRCWAProblem(SheetStack((sheet, ), ()), modes, pol, I₁, I₂)
end

"""
    DeltaRCWAProblem(structure, dims, ω, pol, I₁, I₂, medium=Vacuum())

Constructor 
"""
function DeltaRCWAProblem(structure, dims, ω, pol, I₁, I₂, medium=Vacuum())
    DeltaRCWAProblem(structure, PlanewaveModes(ω, dims, medium), pol, I₁, I₂)
end

struct DeltaRCWASolution{T, N}
    modes::PlanewaveModes{T, N}
    pol::AbstractPolarization
    I₁::Array{ComplexF64, N}
    I₂::Array{ComplexF64, N}
    O₁::Array{ComplexF64, N}
    O₂::Array{ComplexF64, N}
    function DeltaRCWASolution(
        modes::PlanewaveModes{T, N},
        pol::AbstractPolarization,
        I₁::AbstractArray{<: Number, N},
        I₂::AbstractArray{<: Number, N},
        O₁::AbstractArray{<: Number, N},
        O₂::AbstractArray{<: Number, N},
    ) where {T, N}
        check_mode_size(N, modes.dims, I₁, I₂, O₁, O₂)
        # @assert all(size(weights) == Tuple(e[1] for e in modes.dims) for weights in (I₁, I₂, O₁, O₂))
        enforce_N_pol(N, pol)
        new{T, N}(modes, pol, convert.(Array{ComplexF64}, (I₁, I₂, O₁, O₂))...)
    end
end

"""
    solve(::DeltaRCWAProblem; method=smatrix)::DeltaRCWASolution

Evaluates the solution to the problem, with the solver method to be selected by
the user. `method`s currently implemented are `smatrix` for a dense solver and
`smatrixBlockMap` and `smatrixLinearMap` for matrix free methods using
`BlockMap`s and `LinearMap`s respectively.
"""
function solve(prob::DeltaRCWAProblem{T₁, T₂, N, L}; method=smatrix)::DeltaRCWASolution{T₂, N} where {T₁, T₂, N, L}
    structure, modes, pol, I₁, I₂ = prob.structure, prob.modes, prob.pol, prob.I₁, prob.I₂
    sol = method(structure, modes, pol) * vcat(reshape.([I₁, I₂], :)...)
    O₁ = reshape(sol[1:length(I₁)], size(I₁))
    O₂ = reshape(sol[(length(I₁)+1):sum(length.([I₁, I₂]))], size(I₂))
    DeltaRCWASolution(modes, pol, I₁, I₂, O₁, O₂)
end