# Defines the trivial fall-back methods that should be implemented by `RCWASheet`s
export
    RCWASheet,
    ImpedanceStyle,
    Impedanceρₑσₘ,
    Impedanceσₑσₘ,
    σₑˣˣ, σₑˣʸ, σₑʸˣ, σₑʸʸ, σₘˣˣ, σₘˣʸ, σₘʸˣ, σₘʸʸ, ρₑˣˣ, ρₑˣʸ, ρₑʸˣ, ρₑʸʸ

"""
    RCWASheet

An abstract type to dispatch methods for 2D structures embedded in a homogenous
medium modelled with Generalized Sheet Transition Conditions (GSTCs).
See:
Kuester, Edward et al.
"Averaged transition conditions for electromagnetic fields at a metafilm"
https://doi.org/10.1109/TAP.2003.817560
"""
abstract type RCWASheet end

abstract type ImpedanceStyle end

"""
    Impedanceρₑσₘ()

Subtype of `ImpedanceStyle` to describe sheets that implement methods returning
the components of ρₑ and of σₘ, which performs better when there is a perfect
electric conductor and no electric insulator.
"""
struct Impedanceρₑσₘ <: ImpedanceStyle end

"""
    Impedanceσₑσₘ()

Subtype of `ImpedanceStyle` to describe sheets that implement methods returning
the components of σₑ and of σₘ, which performs better when there is a perfect
electric insulator and no electric conductor.
"""
struct Impedanceσₑσₘ <: ImpedanceStyle end

"""
    ImpedanceStyle(::RCWASheet)

Based on IndexStyle in Base, provides an interface for user-defined sheets to
choose between a formulation of the same problem in terms of their preferred
conductivity parameters so long as they extend this method to choose a concrete
impedance style for their type.
"""
ImpedanceStyle(::RCWASheet) = Impedanceρₑσₘ()

# Set default admittance/impedance parameters
# insulating
for σ in (:σₑˣˣ, :σₑˣʸ, :σₑʸˣ, :σₑʸʸ, :σₘˣˣ, :σₘˣʸ, :σₘʸˣ, :σₘʸʸ)
    @eval $(σ)(::RCWASheet, x⃗) = zero(eltype(x⃗))
end
# conducting
for ρ in (:ρₑˣˣ, :ρₑˣʸ, :ρₑʸˣ, :ρₑʸʸ)
    @eval $(ρ)(::RCWASheet, x⃗) = zero(eltype(x⃗))
end

"""
    smatrix(sheet::RCWASheet{1}, modes, ::UncoupledPolarization)

Returns a 2x2 BlockMatrix for the scattering of modes specific to the TE or TM 
polarization
"""
function _sMatrix(sheet::RCWASheet, modes::PlanewaveModes{T, 1} where T, pol::UncoupledPolarization)
    n = length(modes.kz)
    Z = ImpedanceStyle(sheet)
    Zˣˣ, Zʸʸ, ωη = get_1D_uncoupled_GSTC_params(Z, sheet, modes, pol)
    Z̃ˣˣ, Z̃ʸʸ = diagFT.((Zˣˣ, Zʸʸ))
    A, B = get_1D_uncoupled_GSTC_matrices(Z, Z̃ˣˣ, Z̃ʸʸ, Diagonal(modes.kz), ωη)
    A\B
end

function _sBlockMatrix(sheet::RCWASheet, modes::PlanewaveModes{T, N}, pol) where {T, N}
    n = length(modes.kz)
    BlockMatrix(_sMatrix(sheet, modes, pol), [N*n, N*n], [N*n, N*n])
end

@generated function get_1D_uncoupled_GSTC_params(Z::ImpedanceStyle, sheet, modes, pol)
    if Z === Impedanceρₑσₘ
        Zˣˣ, Zʸʸ = :ρ, :σ
    else # Z === Impedanceσₑσₘ
        Zˣˣ, Zʸʸ = :σ, :σ
    end
    if pol === TE
        η, xx, yy = :μ, :ₘ, :ₑ
    else # pol === TM
        η, xx, yy = :ϵ, :ₑ, :ₘ
    end
    quote
        Zx = $(Symbol(Zˣˣ, xx, :ˣˣ)).(Ref(sheet), Iterators.product(modes.x⃗...))
        Zy = $(Symbol(Zʸʸ, yy, :ʸʸ)).(Ref(sheet), Iterators.product(modes.x⃗...))
        ωη = modes.ω * modes.M.$(η)
        Zx, Zy, ωη
    end
end

"""
    diagFT(A::AbstractArray{<:Number, N}) where N

Accepts a diagonal endomorphism in N-D position space and returns the Fourier
transformed  rank-2N tensor reshaped to a matrix with flattened indices
"""
function diagFT(A::AbstractArray{<:Number, ndim}) where ndim
    n⃗ = size(A)
    n = prod(n⃗)
    Ã = Array(reshape(Diagonal(reshape(A, n)), (n⃗..., n⃗...)))
    Ã = fft(ifft(Ã, 1:ndim), (ndim+1):2ndim)
    reshape(Ã, (n, n))
end

function get_1D_uncoupled_GSTC_matrices(::Impedanceρₑσₘ, Z̃ˣˣ, Z̃ʸʸ, kz, ωη)
    _get_1D_uncoupled_GSTC_matrices(I, Z̃ˣˣ, I, Z̃ʸʸ, kz, ωη)
end
function get_1D_uncoupled_GSTC_matrices(::Impedanceσₑσₘ, Z̃ˣˣ, Z̃ʸʸ, kz, ωη)
    _get_1D_uncoupled_GSTC_matrices(Z̃ˣˣ, I, I, Z̃ʸʸ, kz, ωη)
end
function _get_1D_uncoupled_GSTC_matrices(xxwithkz, xxother, yywithkz, yyother, kz, ωη)
    A = [
        (xxwithkz * kz + 2ωη * xxother)   (-xxwithkz * kz - 2ωη * xxother);
        (yywithkz * kz + 0.5ωη * yyother) (yywithkz * kz + 0.5ωη * yyother)
    ]
    B = [
        (xxwithkz * kz - 2ωη * xxother)   (-xxwithkz * kz + 2ωη * xxother);
        (yywithkz * kz - 0.5ωη * yyother) (yywithkz * kz - 0.5ωη * yyother)
    ]
    A, B
end

"""
    smatrixLinearMap(sheet, modes, pol)

Return a function that will compute the matrix-vector product with the
scattering matrix of a RCWASheet in a matrix-free fashion.
"""
function _sLinearMap(sheet::RCWASheet, modes::PlanewaveModes{T, 1} where T, pol::UncoupledPolarization)
    Z = ImpedanceStyle(sheet)
    Zˣˣ, Zʸʸ, ωη = get_1D_uncoupled_GSTC_params(Z, sheet, modes, pol)
    Z̃ˣˣ = ifft ∘ (x -> Diagonal(Zˣˣ) * x) ∘ fft
    Z̃ʸʸ = ifft ∘ (x -> Diagonal(Zʸʸ) * x) ∘ fft
    A, B = get_1D_uncoupled_GSTC_LinearMaps(Z, Z̃ˣˣ, Z̃ʸʸ, Diagonal(modes.kz), ωη)
    LinearMap(x -> gmres(A, B(x)), 2length(modes.kz))
end

function get_1D_uncoupled_GSTC_LinearMaps(::Impedanceρₑσₘ, Z̃ˣˣ, Z̃ʸʸ, kz, ωη)
    _get_1D_uncoupled_GSTC_LinearMaps(x -> x, Z̃ˣˣ, x -> x, Z̃ʸʸ, kz, ωη)
end
function get_1D_uncoupled_GSTC_LinearMaps(::Impedanceσₑσₘ, Z̃ˣˣ, Z̃ʸʸ, kz, ωη)
    _get_1D_uncoupled_GSTC_LinearMaps(Z̃ˣˣ, x -> x, x -> x, Z̃ʸʸ, kz, ωη)
end     
function _get_1D_uncoupled_GSTC_LinearMaps(xxwithkz, xxother, yywithkz, yyother, kz, ωη)
    N = size(kz, 1)
    A = LinearMap(2N) do x
        I₁ = x[1:N]
        I₂ = x[(N+1):2N]
        sumI₁I₂ = I₁ + I₂
        diffI₁I₂ = I₁ - I₂
        vcat(
            xxwithkz(kz * diffI₁I₂) + xxother(2ωη * diffI₁I₂),
            yywithkz(kz * sumI₁I₂) + yyother(0.5ωη * sumI₁I₂)
        )
    end
    B = LinearMap(2N) do x
        I₁ = x[1:N]
        I₂ = x[(N+1):2N]
        sumI₁I₂ = I₁ + I₂
        diffI₁I₂ = I₁ - I₂
        vcat(
            xxwithkz(kz * diffI₁I₂) - xxother(2ωη * diffI₁I₂),
            yywithkz(kz * sumI₁I₂) - yyother(0.5ωη * sumI₁I₂)
        )
    end
    A, B
end

"""
    smatrixBlockMap(sheet, modes, pol)

Return a function that will compute the matrix-vector product with the
scattering matrix of a RCWASheet in a matrix-free fashion.
"""
function _sBlockMap(sheet::RCWASheet, modes::PlanewaveModes{T, 1} where T, pol::UncoupledPolarization)
    N = length(modes.kz)
    Z = ImpedanceStyle(sheet)
    Zˣˣ, Zʸʸ, ωη = get_1D_uncoupled_GSTC_params(Z, sheet, modes, pol)
    _fft = LinearMap{ComplexF64}(fft, ifft, N)
    Z̃ˣˣ = _fft' * LinearMap(Diagonal(Zˣˣ)) * _fft
    Z̃ʸʸ = _fft' * LinearMap(Diagonal(Zʸʸ)) * _fft
    kz = LinearMap(Diagonal(modes.kz))
    A, B = get_1D_uncoupled_GSTC_matrices(Z, Z̃ˣˣ, Z̃ʸʸ, kz, ωη)
    # Build the matrix inversion block-by-block
    B₁₁, B₁₂, B₂₁, B₂₂ = B.maps
    invA_B₁₁ = LinearMap(x -> gmres(A, vcat(B₁₁ * x, B₂₁ * x))[1:N], N)
    invA_B₁₂ = LinearMap(x -> gmres(A, vcat(B₁₂ * x, B₂₂ * x))[1:N], N)
    invA_B₂₁ = LinearMap(x -> gmres(A, vcat(B₁₁ * x, B₂₁ * x))[(N+1):2N], N)
    invA_B₂₂ = LinearMap(x -> gmres(A, vcat(B₁₂ * x, B₂₂ * x))[(N+1):2N], N)
    [
        invA_B₁₁    invA_B₁₂;
        invA_B₂₁    invA_B₂₂;
    ]
end

"""
    smatrix(sheet, modes::PlanewaveModes{T, 2} where T, ::CoupledPolarization)

Returns a 2x2 BlockMatrix for the scattering of modes specific to the TE or TM 
polarization
"""
function _sMatrix(sheet::RCWASheet, modes::PlanewaveModes{T, 2} where T, ::CoupledPolarization)
    n = length(modes.kz)
    A, B = _2Dsheetsmatrix(sheet, modes)
    A\B
end

"Return the variables needed to construct the dense scattering matrix"
function _get_params_2Dsheetsmatrix(sheet, modes)
    n = length(modes.kz)
    R = (
        xx = (0I)(n),
        xy = (-I)(n),
        yx = (1I)(n),
        yy = (0I)(n),
    )
    ωμ = modes.ω * modes.M.μ
    k⃗² = modes.ω * modes.M.ϵ * ωμ
    ρₑ = (
        xx = ρₑˣˣ.(Ref(sheet), Iterators.product(modes.x⃗...)),
        xy = ρₑˣʸ.(Ref(sheet), Iterators.product(modes.x⃗...)),
        yx = ρₑʸˣ.(Ref(sheet), Iterators.product(modes.x⃗...)),
        yy = ρₑʸʸ.(Ref(sheet), Iterators.product(modes.x⃗...)),
    )
    σₘ = (
        xx = σₘˣˣ.(Ref(sheet), Iterators.product(modes.x⃗...)),
        xy = σₘˣʸ.(Ref(sheet), Iterators.product(modes.x⃗...)),
        yx = σₘʸˣ.(Ref(sheet), Iterators.product(modes.x⃗...)),
        yy = σₘʸʸ.(Ref(sheet), Iterators.product(modes.x⃗...)),
    )
    K = (
        xx = [k⃗[1] .* k⃗[1] for k⃗ in Iterators.product(modes.k⃗...)],
        xy = [k⃗[1] .* k⃗[2] for k⃗ in Iterators.product(modes.k⃗...)],
        # yx = [k⃗[2] * k⃗[1] for k⃗ in Iterators.product(modes.k⃗...)],
        yy = [k⃗[2] .* k⃗[2] for k⃗ in Iterators.product(modes.k⃗...)],
    )
    ρₑ, σₘ, R, K, modes.kz, k⃗², ωμ
end

function _2Dsheetsmatrix(sheet, modes)
    ρₑ, σₘ, R, K, kz, k⃗², ωμ = _get_params_2Dsheetsmatrix(sheet, modes)
    K_term = mortar(
        Diagonal.(reshape.((K.xx .- k⃗², K.xy), :)),
        Diagonal.(reshape.((K.xy, K.yy .- k⃗²), :)),
    )
    R_term = mortar(
        (R.xx, R.xy),
        (R.yx, R.yy),
    )
    # need to do Fourier transform of each block
    ρ̃_term = mortar(
        diagFT.((ρₑ.xx, ρₑ.xy)),
        diagFT.((ρₑ.yx, ρₑ.yy)),
    )
    ρ̃_term = Matrix(K_term * transpose(R_term) * ρ̃_term * R_term)
    σ̃_term = mortar(
        diagFT.((σₘ.xx, σₘ.xy)),
        diagFT.((σₘ.yx, σₘ.yy)),
    )
    σ̃_term = Matrix(K_term * σ̃_term)
    # return ρ̃_term, ρₑ
    kz_term = ωμ * Diagonal(vcat(reshape.((kz, kz), :)...))
    _build_2D_GSTC_smatrix(ρ̃_term, σ̃_term, kz_term)
end

function _build_2D_GSTC_smatrix(ρ̃_term, σ̃_term, kz_term)
    A = [
        -ρ̃_term + 0.5*kz_term   ρ̃_term - 0.5*kz_term;
        -σ̃_term + 2*kz_term     -σ̃_term + 2*kz_term
    ]
    B = [
        ρ̃_term + 0.5*kz_term   -ρ̃_term - 0.5*kz_term;
        σ̃_term + 2*kz_term     σ̃_term - 2*kz_term
    ]
    A, B
end