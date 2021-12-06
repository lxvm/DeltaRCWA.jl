export DeltaRCWAProblem, solve

"""
    DeltaRCWAProblem(
        stack::T,
        pw::PlaneWaves{N},
        I₁::AbstractArray{<:Number, N},
        I₂::AbstractArray{<:Number, N},
    ) where {N, T<:SheetStack}

This type of problem assumes that the only type of scatterer are sheets and that
the solution can be expressed entirely by one set of pw
"""
struct DeltaRCWAProblem{N, T<:SheetStack}
    stack::T
    pw::PlaneWaves{N}
    I₁::Array{ComplexF64, N}
    I₂::Array{ComplexF64, N}
    function DeltaRCWAProblem(
        stack::T,
        pw::PlaneWaves{N},
        I₁::AbstractArray{<:Number, N},
        I₂::AbstractArray{<:Number, N},
    ) where {N, T<:SheetStack}
        check_mode_size(N, pw.dims, I₁, I₂)
        new{N, T}(stack, pw, I₁, I₂)
    end
end

function check_mode_size(N, dims, weights...)
    for w in weights
        @assert size(w) == Tuple(i == length(dims) ? N * dims[i][1] : dims[i][1] for i in eachindex(dims)) ""
    end
end

"""
    DeltaRCWAProblem(structure, dims, ω, I₁, I₂)

Outer constructor for convenience
"""
function DeltaRCWAProblem(structure, dims, ω, I₁, I₂)
    DeltaRCWAProblem(convert(SheetStack, structure), PlaneWaves(ω, dims), I₁, I₂)
end

FieldStyle(prob::DeltaRCWAProblem) = FieldStyle(typeof(prob))
FieldStyle(::Type{<:DeltaRCWAProblem}) = HField()

struct DeltaRCWASolution{N, T<:SheetStack}
    stack::T
    pw::PlaneWaves{N}
    I₁::Array{ComplexF64, N}
    I₂::Array{ComplexF64, N}
    O₁::Array{ComplexF64, N}
    O₂::Array{ComplexF64, N}
    function DeltaRCWASolution(
        stack::T,
        pw::PlaneWaves{N},
        I₁::AbstractArray{<:Number, N},
        I₂::AbstractArray{<:Number, N},
        O₁::AbstractArray{<:Number, N},
        O₂::AbstractArray{<:Number, N},
    ) where {N, T<:SheetStack}
        check_mode_size(N, pw.dims, I₁, I₂, O₁, O₂)
        new{N, T}(stack, pw, I₁, I₂, O₁, O₂)
    end
end

"""
    solve(::DeltaRCWAProblem; T=Matrix)::DeltaRCWASolution

Evaluates the solution to the problem, with the solver method to be selected by
the user. `T` can be `Matrix`, `BlockMatrix`, `LinearMap`, or `BlockMap` and
selects the underlying representation of the scattering matrices.
"""
function solve(prob::DeltaRCWAProblem{N}; T=Matrix)::DeltaRCWASolution{N} where N
    stack, pw, I₁, I₂ = prob.stack, prob.pw, prob.I₁, prob.I₂
    sol = smatrix(T, FieldStyle(prob), pw, stack) * vcat(reshape.([I₁, I₂], :)...)
    O₁ = reshape(sol[1:length(I₁)], size(I₁))
    O₂ = reshape(sol[(length(I₁)+1):sum(length.([I₁, I₂]))], size(I₂))
    DeltaRCWASolution(stack, pw, I₁, I₂, O₁, O₂)
end