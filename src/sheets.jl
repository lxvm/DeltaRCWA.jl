# Define default behavior of impedance/admittance parameters
for R in (:Z, :Y), s in (:ₑ, :ₘ), i in (:ˣ, :ʸ), j in (:ˣ, :ʸ)
    @eval begin
        export $(Symbol(R, s, i, j))
        $(Symbol(R, s, i, j))(::Sheet, x⃗) = zero(eltype(x⃗))
    end
end

"""
    _sMatrix(sheet::Sheet, modes::PlanewaveModes{T, 1} where T, ::UncoupledPolarization)

Returns a 2x2 Matrix for the scattering of modes specific to the TE or TM 
polarization
"""
function _sMatrix(sheet::Sheet, modes::PlanewaveModes{T, 1} where T, pol::UncoupledPolarization)
    n = length(modes.kz)
    Rₑ = ElectricResponseStyle(sheet)
    Rₘ = MagneticResponseStyle(sheet)
    Rˣˣ, Rʸʸ, ωη = _get_1D_uncoupled_GSTC_params(Rₑ, Rₘ, sheet, modes, pol)
    R̃ˣˣ, R̃ʸʸ = diagFT.((Rˣˣ, Rʸʸ))
    A, B = _get_1D_uncoupled_GSTC_matrices(Rₑ, Rₘ, R̃ˣˣ, R̃ʸʸ, Diagonal(modes.kz), ωη)
    A\B
end

function _sBlockMatrix(sheet::Sheet, modes::PlanewaveModes{T, N}, pol) where {T, N}
    n = length(modes.kz)
    BlockMatrix(_sMatrix(sheet, modes, pol), [N*n, N*n], [N*n, N*n])
end

@generated function _get_1D_uncoupled_GSTC_params(Rₑ::ResponseStyle, Rₘ::ResponseStyle, sheet, modes, pol)
    Rˣ = Rₑ === Impedance ? :Z : :Y
    Rʸ = Rₘ === Impedance ? :Z : :Y
    if pol === TE
        η, sˣ, sʸ = :μ, :ₘ, :ₑ
    else # pol === TM
        η, sˣ, sʸ = :ϵ, :ₑ, :ₘ
    end
    quote
        Rˣˣ = $(Symbol(Rˣ, sˣ, :ˣˣ)).(Ref(sheet), Iterators.product(modes.x⃗...))
        Rʸʸ = $(Symbol(Rʸ, sʸ, :ʸʸ)).(Ref(sheet), Iterators.product(modes.x⃗...))
        ωη = modes.ω * modes.M.$(η)
        Rˣˣ, Rʸʸ, ωη
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

@generated function _get_1D_uncoupled_GSTC_matrices(Rₑ::ResponseStyle, Rₘ::ResponseStyle, R̃ˣˣ, R̃ʸʸ, kz, ωη)
    (A, B) = Rₑ === Impedance ? (:I, :R̃ˣˣ) : (:R̃ˣˣ, :I)
    (C, D) = Rₘ === Impedance ? (:R̃ʸʸ, :I) : (:I, :R̃ʸʸ)
    :(__get_1D_uncoupled_GSTC_matrices($A, $B, $C, $D, kz, ωη))
end
function __get_1D_uncoupled_GSTC_matrices(xxwithkz, xxother, yywithkz, yyother, kz, ωη)
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
    _sLinearMap(sheet, modes, pol)

Return a function that will compute the matrix-vector product with the
scattering matrix of a Sheet in a matrix-free fashion.
"""
function _sLinearMap(sheet::Sheet, modes::PlanewaveModes{T, 1} where T, pol::UncoupledPolarization)
    Rₑ = ElectricResponseStyle(sheet)
    Rₘ = MagneticResponseStyle(sheet)
    Rˣˣ, Rʸʸ, ωη = _get_1D_uncoupled_GSTC_params(Rₑ, Rₘ, sheet, modes, pol)
    R̃ˣˣ, R̃ʸʸ = diagFT.((Rˣˣ, Rʸʸ))
    R̃ˣˣ = ifft ∘ (x -> Diagonal(Rˣˣ) * x) ∘ fft
    R̃ʸʸ = ifft ∘ (x -> Diagonal(Rʸʸ) * x) ∘ fft
    A, B = _get_1D_uncoupled_GSTC_LinearMaps(Rₑ, Rₘ, R̃ˣˣ, R̃ʸʸ, Diagonal(modes.kz), ωη)
    LinearMap(x -> gmres(A, B(x)), 2length(modes.kz))
end

@generated function _get_1D_uncoupled_GSTC_LinearMaps(Rₑ::ResponseStyle, Rₘ::ResponseStyle, R̃ˣˣ, R̃ʸʸ, kz, ωη)
    (A, B) = Rₑ === Impedance ? (:identity, :R̃ˣˣ) : (:R̃ˣˣ, :identity)
    (C, D) = Rₘ === Impedance ? (:R̃ʸʸ, :identity) : (:identity, :R̃ʸʸ)
    :(__get_1D_uncoupled_GSTC_LinearMaps($A, $B, $C, $D, kz, ωη))
end
function __get_1D_uncoupled_GSTC_LinearMaps(xxwithkz, xxother, yywithkz, yyother, kz, ωη)
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
scattering matrix of a Sheet in a matrix-free fashion.
"""
function _sBlockMap(sheet::Sheet, modes::PlanewaveModes{T, 1} where T, pol::UncoupledPolarization)
    N = length(modes.kz)
    Rₑ = ElectricResponseStyle(sheet)
    Rₘ = MagneticResponseStyle(sheet)
    Rˣˣ, Rʸʸ, ωη = _get_1D_uncoupled_GSTC_params(Rₑ, Rₘ, sheet, modes, pol)
    _fft = LinearMap{ComplexF64}(fft, ifft, N)
    R̃ˣˣ = _fft' * LinearMap(Diagonal(Rˣˣ)) * _fft
    R̃ʸʸ = _fft' * LinearMap(Diagonal(Rʸʸ)) * _fft
    kz = LinearMap(Diagonal(modes.kz))
    A, B = _get_1D_uncoupled_GSTC_matrices(Rₑ, Rₘ, R̃ˣˣ, R̃ʸʸ, kz, ωη)
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
function _sMatrix(sheet::Sheet, modes::PlanewaveModes{T, 2} where T, ::CoupledPolarization)
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
    Zₑ = (
        xx = Zₑˣˣ.(Ref(sheet), Iterators.product(modes.x⃗...)),
        xy = Zₑˣʸ.(Ref(sheet), Iterators.product(modes.x⃗...)),
        yx = Zₑʸˣ.(Ref(sheet), Iterators.product(modes.x⃗...)),
        yy = Zₑʸʸ.(Ref(sheet), Iterators.product(modes.x⃗...)),
    )
    Yₘ = (
        xx = Yₘˣˣ.(Ref(sheet), Iterators.product(modes.x⃗...)),
        xy = Yₘˣʸ.(Ref(sheet), Iterators.product(modes.x⃗...)),
        yx = Yₘʸˣ.(Ref(sheet), Iterators.product(modes.x⃗...)),
        yy = Yₘʸʸ.(Ref(sheet), Iterators.product(modes.x⃗...)),
    )
    K = (
        xx = [k⃗[1] .* k⃗[1] for k⃗ in Iterators.product(modes.k⃗...)],
        xy = [k⃗[1] .* k⃗[2] for k⃗ in Iterators.product(modes.k⃗...)],
        # yx = [k⃗[2] * k⃗[1] for k⃗ in Iterators.product(modes.k⃗...)],
        yy = [k⃗[2] .* k⃗[2] for k⃗ in Iterators.product(modes.k⃗...)],
    )
    Zₑ, Yₘ, R, K, modes.kz, k⃗², ωμ
end

function _2Dsheetsmatrix(sheet, modes)
    Zₑ, Yₘ, R, K, kz, k⃗², ωμ = _get_params_2Dsheetsmatrix(sheet, modes)
    K_term = mortar(
        Diagonal.(reshape.((K.xx .- k⃗², K.xy), :)),
        Diagonal.(reshape.((K.xy, K.yy .- k⃗²), :)),
    )
    R_term = mortar(
        (R.xx, R.xy),
        (R.yx, R.yy),
    )
    # need to do Fourier transform of each block
    Z̃ₑ_term = mortar(
        diagFT.((Zₑ.xx, Zₑ.xy)),
        diagFT.((Zₑ.yx, Zₑ.yy)),
    )
    Z̃ₑ_term = Matrix(K_term * transpose(R_term) * Z̃ₑ_term * R_term)
    Ỹₘ_term = mortar(
        diagFT.((Yₘ.xx, Yₘ.xy)),
        diagFT.((Yₘ.yx, Yₘ.yy)),
    )
    Ỹₘ_term = Matrix(K_term * Ỹₘ_term)
    # return Z̃ₑ_term, Zₑ
    kz_term = ωμ * Diagonal(vcat(reshape.((kz, kz), :)...))
    _build_2D_GSTC_smatrix(Z̃ₑ_term, Ỹₘ_term, kz_term)
end

function _build_2D_GSTC_smatrix(Z̃ₑ_term, Ỹₘ_term, kz_term)
    A = [
        -Z̃ₑ_term + 0.5*kz_term   Z̃ₑ_term - 0.5*kz_term;
        -Ỹₘ_term + 2*kz_term     -Ỹₘ_term + 2*kz_term
    ]
    B = [
        Z̃ₑ_term + 0.5*kz_term   -Z̃ₑ_term - 0.5*kz_term;
        Ỹₘ_term + 2*kz_term     Ỹₘ_term - 2*kz_term
    ]
    A, B
end