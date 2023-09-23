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
    DeltaRCWAProblem(SheetStack(structure), PlaneWaves(ω, dims), I₁, I₂)
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
    solve(::DeltaRCWAProblem, T=Matrix)::DeltaRCWASolution

Evaluates the solution to the problem, with the solver method to be selected by
the user. `T` can be `Matrix`, `BlockMatrix`, `LinearMap`, or `BlockMap` and
selects the underlying representation of the scattering matrices.
"""
function solve(prob::DeltaRCWAProblem{N}, ::Type{T}=Matrix)::DeltaRCWASolution{N} where {N,T}
    stack, pw, I₁, I₂ = prob.stack, prob.pw, prob.I₁, prob.I₂
    sol = smatrix(T, FieldStyle(prob), pw, stack) * vcat(reshape.([I₁, I₂], :)...)
    O₁ = reshape(sol[1:length(I₁)], size(I₁))
    O₂ = reshape(sol[(length(I₁)+1):sum(length.([I₁, I₂]))], size(I₂))
    DeltaRCWASolution(stack, pw, I₁, I₂, O₁, O₂)
end

struct ScatteringProblem{N, T<:SheetStack}
    stack::T
    pw::PlaneWaves{N}
end

ScatteringProblem(structure, pw) = ScatteringProblem(SheetStack(structure), pw)

struct ScatteringSolution{N, T<:SheetStack, S, F}
    stack::T
    pw::PlaneWaves{N}
    smatrix::S
    field::F
end

FieldStyle(prob::ScatteringProblem) = FieldStyle(typeof(prob))
FieldStyle(::Type{<:ScatteringProblem}) = HField()

function solve(prob::ScatteringProblem{N}, ::Type{T}=Matrix)::ScatteringSolution{N} where {N,T}
    stack, pw = prob.stack, prob.pw
    fs = FieldStyle(prob)
    ScatteringSolution(stack, pw, smatrix(T, fs, pw, stack), fs)
end

function normalization(::EField, pw::PlaneWaves{1}, m, mask)
    # TE mode # mask needed to avoid complex k
    eps, mu = _get_ϵμ(m)
    kz = _get_kz(pw, m)
    inv.(sqrt.(kz[mask] ./ (2*mu*pw.ω)))
end
function normalization(::HField, pw::PlaneWaves{1}, m, mask)
    # TM mode # mask needed to avoid complex k
    eps, mu = _get_ϵμ(m)
    kz = _get_kz(pw, m)
    inv.(sqrt.(kz[mask] ./ (2*eps*pw.ω)))
end
function normalization(::HField, pw::PlaneWaves{2}, m)
    twiddle = vcat(
        inv.(sqrt.(kz[mask] .* (1 .+ abs2.(getindex.(collect(Iterators.product(pw.k⃗...))[mask], 1)) ./ abs2.(kz[mask])))),
        inv.(sqrt.(kz[mask] .* (1 .+ abs2.(getindex.(collect(Iterators.product(pw.k⃗...))[mask], 2)) ./ abs2.(kz[mask])))),
    )
end
function normalization(::EField, pw::PlaneWaves{2}, m)
    twiddle = vcat(
        inv.(sqrt.(kz[mask] .* (1 .+ abs2.(getindex.(collect(Iterators.product(pw.k⃗...))[mask], 1)) ./ abs2.(kz[mask])))),
        inv.(sqrt.(kz[mask] .* (1 .+ abs2.(getindex.(collect(Iterators.product(pw.k⃗...))[mask], 2)) ./ abs2.(kz[mask])))),
    )
end

function propagating_unitary_smatrix(sol::ScatteringSolution)
    S, pw, stack, fs = sol.smatrix, sol.pw, sol.stack, sol.field
    m1 = stack.media[begin]
    m2 = stack.media[end]
    mask1 = .!iszero.(isreal.(_get_kz(pw, m1)))
    mask2 = .!iszero.(isreal.(_get_kz(pw, m2)))
    mask = vcat(mask1, mask2)
    norm1 = normalization(fs, pw, m1, mask1)
    norm2 = normalization(fs, pw, m2, mask2)
    D = Diagonal(vcat(norm1, norm2))
    U = D \ S[mask,mask] * D    # should be unitary!
    return U, mask
end

# returns the normalized transmission coefficients between 
function transmission_coefficient(sol::ScatteringSolution{1})
    S, mask = propagating_unitary_smatrix(sol)
    n = length(Iterators.product(sol.pw.k⃗...))
    m = count(mask[1:n+1])
    S[m,1], S[1,m]
end