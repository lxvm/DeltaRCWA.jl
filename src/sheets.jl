export Sheet, ElectricResponseStyle, MagneticResponseStyle, Impedance, Admittance

"""
    Sheet

An abstract type to dispatch methods for 2D structures embedded in a homogenous
medium modelled with Generalized Sheet Transition Conditions (GSTCs).
See:
Kuester, Edward et al.
"Averaged transition conditions for electromagnetic fields at a metafilm"
https://doi.org/10.1109/TAP.2003.817560
"""
abstract type Sheet end

"""
    ResponseStyle

Based on `IndexStyle` in Base, provides an interface for user-defined types to
choose between a formulation of the same problem in terms of their preferred
material response parameters.
"""
abstract type ResponseStyle end

"""
    Impedance()

`ResponseStyle` to describe materials whose response is given by an impedance.
"""
struct Impedance <: ResponseStyle end

"""
    Admittance()

`ResponseStyle` to describe materials whose response is given by an admittance.
"""
struct Admittance <: ResponseStyle end

"""
    ElectricResponseStyle(sheet)
    ElectricResponseStyle(typeof(sheet))
"""
ElectricResponseStyle(sheet::Sheet) = ElectricResponseStyle(typeof(sheet))
ElectricResponseStyle(::Type{<:Sheet}) = Impedance()

"""
    MagneticResponseStyle(sheet)
    MagneticResponseStyle(typeof(sheet))
"""
MagneticResponseStyle(sheet::Sheet) = MagneticResponseStyle(typeof(sheet))
MagneticResponseStyle(::Type{<:Sheet}) = Admittance()

# Define default behavior of impedance/admittance parameters
for R in (:Z, :Y), s in (:ₑ, :ₘ), i in (:ˣ, :ʸ), j in (:ˣ, :ʸ)
    @eval begin
        export $(Symbol(R, s, i, j))
        $(Symbol(R, s, i, j))(::Sheet, x⃗) = zero(eltype(x⃗))
    end
end

"""
    _sMatrix(::FieldStyle, ::PlaneWaves, ::Sheet, ::UniformMedium, ::UniformMedium)
"""
_sMatrix(fs, pw::PlaneWaves{1}, sheet::Sheet, um₁, um₂) = _sMatrix(fs, PolarizationStyle(pw), pw, sheet, um₁, um₂)
_sMatrix(fs::EField, pol::TE, pw, sheet, um₁, um₂) = __sMatrix(fs, pol, pw, ElectricResponseStyle(sheet), MagneticResponseStyle(sheet), sheet, um₁, um₂)
_sMatrix(fs::HField, pol::TM, pw, sheet, um₁, um₂) = __sMatrix(fs, pol, pw, ElectricResponseStyle(sheet), MagneticResponseStyle(sheet), sheet, um₁, um₂)
function __sMatrix(fs, pol, pw, rsₑ, rsₘ, sheet, um₁, um₂)
    R̃ₑ = diagFT(_get_electric_response(pol, rsₑ, pw, sheet))
    R̃ₘ = diagFT(_get_magnetic_response(pol, rsₘ, pw, sheet))
    kz₁ = Diagonal(_get_kz(pw, um₁))
    kz₂ = Diagonal(_get_kz(pw, um₂))
    ωη₁ = pw.ω * _get_η(pol, um₁)
    ωη₂ = pw.ω * _get_η(pol, um₂)
    out, inc = _gstcMatrix(fs,
        _set_electric_terms(pol, rsₑ, R̃ₑ, kz₁, kz₂, ωη₁, ωη₂)...,
        _set_magnetic_terms(pol, rsₘ, R̃ₘ, kz₁, kz₂, ωη₁, ωη₂)...,
    )
    out\inc
end
function __sMatrix(fs, pol, pw, rsₑ, rsₘ, sheet, um::T, ::T) where T<:UniformMedium
    R̃ₑ = diagFT(_get_electric_response(pol, rsₑ, pw, sheet))
    R̃ₘ = diagFT(_get_magnetic_response(pol, rsₘ, pw, sheet))
    kz = Diagonal(_get_kz(pw, um))
    ωη = pw.ω * _get_η(pol, um)
    out, inc = _gstcMatrix(fs,
        _set_electric_terms(pol, rsₑ, R̃ₑ, kz, ωη)...,
        _set_magnetic_terms(pol, rsₘ, R̃ₘ, kz, ωη)...,
    )
    out\inc
end

# v1, mutates arrays
# _get_electric_response(::TM, ::Admittance, pw, sheet) = Yₑˣˣ.(Ref(sheet), Iterators.product(pw.x⃗...))
# _get_electric_response(::TM, ::Impedance,  pw, sheet) = Zₑˣˣ.(Ref(sheet), Iterators.product(pw.x⃗...))
# _get_magnetic_response(::TM, ::Admittance, pw, sheet) = Yₘʸʸ.(Ref(sheet), Iterators.product(pw.x⃗...))
# _get_magnetic_response(::TM, ::Impedance,  pw, sheet) = Zₘʸʸ.(Ref(sheet), Iterators.product(pw.x⃗...))

# _get_electric_response(::TE, ::Admittance, pw, sheet) = Yₑʸʸ.(Ref(sheet), Iterators.product(pw.x⃗...))
# _get_electric_response(::TE, ::Impedance,  pw, sheet) = Zₑʸʸ.(Ref(sheet), Iterators.product(pw.x⃗...))
# _get_magnetic_response(::TE, ::Admittance, pw, sheet) = Yₘˣˣ.(Ref(sheet), Iterators.product(pw.x⃗...))
# _get_magnetic_response(::TE, ::Impedance,  pw, sheet) = Zₘˣˣ.(Ref(sheet), Iterators.product(pw.x⃗...))
# v2, works with zygote
_get_electric_response(::TM, ::Admittance, pw, sheet) = map(Yₑˣˣ, fill(sheet, length.(pw.x⃗)), Iterators.product(pw.x⃗...))
_get_electric_response(::TM, ::Impedance,  pw, sheet) = map(Zₑˣˣ, fill(sheet, length.(pw.x⃗)), Iterators.product(pw.x⃗...))
_get_magnetic_response(::TM, ::Admittance, pw, sheet) = map(Yₘʸʸ, fill(sheet, length.(pw.x⃗)), Iterators.product(pw.x⃗...))
_get_magnetic_response(::TM, ::Impedance,  pw, sheet) = map(Zₘʸʸ, fill(sheet, length.(pw.x⃗)), Iterators.product(pw.x⃗...))

_get_electric_response(::TE, ::Admittance, pw, sheet) = map(Yₑʸʸ, fill(sheet, length.(pw.x⃗)), Iterators.product(pw.x⃗...))
_get_electric_response(::TE, ::Impedance,  pw, sheet) = map(Zₑʸʸ, fill(sheet, length.(pw.x⃗)), Iterators.product(pw.x⃗...))
_get_magnetic_response(::TE, ::Admittance, pw, sheet) = map(Yₘˣˣ, fill(sheet, length.(pw.x⃗)), Iterators.product(pw.x⃗...))
_get_magnetic_response(::TE, ::Impedance,  pw, sheet) = map(Zₘˣˣ, fill(sheet, length.(pw.x⃗)), Iterators.product(pw.x⃗...))

"""
    diagFT(A::AbstractArray{<:Number, ndim}) where ndim

Performs a Fourier transform of a diagonal operator/scalar field over an Nd grid
"""
function diagFT(A::AbstractArray{<:Number, ndim}) where ndim
    n⃗ = size(A)
    n = prod(n⃗)
    Ã = Array(reshape(Diagonal(reshape(A, n)), (n⃗..., n⃗...)))
    Ã = fft(ifft(Ã, 1:ndim), (ndim+1):2ndim)
    reshape(Ã, (n, n))
end

_get_η(::TE, ::UniformMedium{ϵ, μ} where ϵ) where μ = μ
_get_η(::TM, ::UniformMedium{ϵ, μ} where μ) where ϵ = ϵ

_set_electric_terms(::TM, ::Impedance,  R̃, kz₁, kz₂, ωη₁, ωη₂) = (-R̃, kz₁ / 2ωη₁, kz₂ / 2ωη₂)
_set_electric_terms(::TM, ::Admittance, R̃, kz₁, kz₂, ωη₁, ωη₂) = (-I, R̃*(kz₁ / 2ωη₁), R̃*(kz₂ / 2ωη₂))
_set_magnetic_terms(::TM, ::Impedance,  R̃, kz₁, kz₂, ωη₁, ωη₂) = (R̃*(kz₁ / 0.5ωη₁), R̃*(kz₂ / 0.5ωη₂), -I)
_set_magnetic_terms(::TM, ::Admittance, R̃, kz₁, kz₂, ωη₁, ωη₂) = (kz₁ / 0.5ωη₁, kz₂ / 0.5ωη₂, -R̃)

_set_electric_terms(::TM, ::Impedance,  R̃, kz, ωη) = (-R̃, kz / 2ωη)
_set_electric_terms(::TM, ::Admittance, R̃, kz, ωη) = (-I, R̃*(kz / 2ωη))
_set_magnetic_terms(::TM, ::Impedance,  R̃, kz, ωη) = (R̃*(kz / 0.5ωη), -I)
_set_magnetic_terms(::TM, ::Admittance, R̃, kz, ωη) = (kz / 0.5ωη, -R̃)

_set_electric_terms(::TE, args...) = _set_magnetic_terms(TM(), args...)
_set_magnetic_terms(::TE, args...) = _set_electric_terms(TM(), args...)

# compute E field for TM mode by using TM relations to convert H to E
function _sMatrix(::EField, pol::TM, pw, sheet, um₁, um₂)
    S = _sMatrix(HField(), pol, pw, sheet, um₁, um₂)
    η₁ = _get_η(pol, um₁)
    η₂ = _get_η(pol, um₂)
    kz₁ = _get_kz(pw, um₁)
    kz₂ = _get_kz(pw, um₂)
    HtoE = Diagonal(vcat(-kz₁ / pw.ω*η₁, kz₂ / pw.ω*η₂))
    HtoE * S * inv(-HtoE)
end
# compute H field for TE mode by using TE relations to convert E to H
function _sMatrix(::HField, pol::TE, pw, sheet, um₁, um₂)
    S = _sMatrix(EField(), pol, pw, sheet, um₁, um₂)
    η₁ = _get_η(pol, um₁)
    η₂ = _get_η(pol, um₂)
    kz₁ = _get_kz(pw, um₁)
    kz₂ = _get_kz(pw, um₂)
    EtoH = Diagonal(vcat(kz₁ / pw.ω*η₁, -kz₂ / pw.ω*η₂))
    EtoH * S * inv(-EtoH)
end

function _sBlockMatrix(fs, pw::PlaneWaves{N}, sheet::Sheet, um₁, um₂) where N
    n = length(pw)
    BlockMatrix(_sMatrix(fs, pw, sheet, um₁, um₂), [N*n, N*n], [N*n, N*n])
end

"""
    _sLinearMap(pw, sheet, pol)

Return a function that will compute the matrix-vector product with the
scattering matrix of a Sheet in a matrix-free fashion.
"""
function _sLinearMap(fs, pw::PlaneWaves{1}, sheet::Sheet, um₁, um₂)
    rsₑ = ElectricResponseStyle(sheet)
    rsₘ = MagneticResponseStyle(sheet)
    pol = PolarizationStyle(pw)
    n = length(pw)
    R̃ₘ = LinearMap(ifft ∘ (x -> Diagonal(_get_magnetic_response(pol, rsₘ, pw, sheet)) * x) ∘ fft, n)
    R̃ₑ = LinearMap(ifft ∘ (x -> Diagonal(_get_electric_response(pol, rsₑ, pw, sheet)) * x) ∘ fft, n)
    kz₁ = LinearMap(Diagonal(_get_kz(pw, um₁)))
    kz₂ = LinearMap(Diagonal(_get_kz(pw, um₂)))
    ωη₁ = pw.ω * _get_η(pol, um₁)
    ωη₂ = pw.ω * _get_η(pol, um₂)
    out, inc = _gstcLinearMap(fs,
        _set_electric_terms(pol, rsₑ, R̃ₑ, kz₁, kz₂, ωη₁, ωη₂)...,
        _set_magnetic_terms(pol, rsₘ, R̃ₘ, kz₁, kz₂, ωη₁, ωη₂)...,
    )
    LinearMap(x -> gmres(out, inc(x)), 2n)
end

"""
    smatrixBlockMap(pw, sheet, pol)

Return a function that will compute the matrix-vector product with the
scattering matrix of a Sheet in a matrix-free fashion.
"""
function _sBlockMap(fs, pw::PlaneWaves{1}, sheet::Sheet, um₁, um₂)
    N = length(pw)
    rsₑ = ElectricResponseStyle(sheet)
    rsₘ = MagneticResponseStyle(sheet)
    pol = PolarizationStyle(pw)
    _fft = LinearMap{ComplexF64}(fft, ifft, N)
    R̃ₘ = _fft' * LinearMap(Diagonal(_get_magnetic_response(pol, rsₘ, pw, sheet))) * _fft
    R̃ₑ = _fft' * LinearMap(Diagonal(_get_electric_response(pol, rsₑ, pw, sheet))) * _fft
    kz₁ = LinearMap(Diagonal(_get_kz(pw, um₁)))
    kz₂ = LinearMap(Diagonal(_get_kz(pw, um₂)))
    ωη₁ = pw.ω * _get_η(pol, um₁)
    ωη₂ = pw.ω * _get_η(pol, um₂)
    out, inc = _gstcMatrix(fs,
        _set_electric_terms(pol, rsₑ, R̃ₑ, kz₁, kz₂, ωη₁, ωη₂)...,
        _set_magnetic_terms(pol, rsₘ, R̃ₘ, kz₁, kz₂, ωη₁, ωη₂)...,
    )
    out\inc
end

function Base.:\(A, B::BlockMap)
    B₁₁, B₁₂, B₂₁, B₂₂ = B.maps
    N = size(B₁₁, 1)
    # Build the matrix inversion block-by-block
    invA_B₁₁ = LinearMap(x -> gmres(A, vcat(B₁₁ * x, B₂₁ * x))[1:N], N)
    invA_B₁₂ = LinearMap(x -> gmres(A, vcat(B₁₂ * x, B₂₂ * x))[1:N], N)
    invA_B₂₁ = LinearMap(x -> gmres(A, vcat(B₁₁ * x, B₂₁ * x))[(N+1):2N], N)
    invA_B₂₂ = LinearMap(x -> gmres(A, vcat(B₁₂ * x, B₂₂ * x))[(N+1):2N], N)
    [
        invA_B₁₁    invA_B₁₂;
        invA_B₂₁    invA_B₂₂;
    ]
end

_sMatrix(fs, pw::PlaneWaves{2}, sheet::Sheet, um₁, um₂) = __sMatrix(fs, pw, ElectricResponseStyle(sheet), MagneticResponseStyle(sheet), sheet, um₁, um₂)
function __sMatrix(fs, pw, rsₑ, rsₘ, sheet, um₁, um₂)
    Rₑ = _get_electric_response_components(rsₑ, pw, sheet)
    Rₘ = _get_magnetic_response_components(rsₘ, pw, sheet)
    R̃ₑ = [
        diagFT(vec(Rₑ.xx))   diagFT(vec(Rₑ.xy))
        diagFT(vec(Rₑ.yx))   diagFT(vec(Rₑ.yy))
    ]
    R̃ₘ = [
        diagFT(vec(Rₘ.xx))   diagFT(vec(Rₘ.xy))
        diagFT(vec(Rₘ.yx))   diagFT(vec(Rₘ.yy))
    ]
    ωϵ₁, ωμ₁ = pw.ω .* _get_ϵμ(um₁)
    ωϵ₂, ωμ₂ = pw.ω .* _get_ϵμ(um₂)
    kz₁ = Diagonal(repeat(vec(_get_kz(pw, um₁)), 2))
    kz₂ = Diagonal(repeat(vec(_get_kz(pw, um₂)), 2))
    _K = _get_K_components(pw)
    K = [
        Diagonal(vec(_K.xx)) Diagonal(vec(_K.xy))
        Diagonal(vec(_K.yx)) Diagonal(vec(_K.yy))
    ]
    Ω₁ = K - (ωϵ₁ * ωμ₁)*I
    Ω₂ = K - (ωϵ₂ * ωμ₂)*I
    n = length(pw)
    C = [
        (0I)(n) (-I)(n)
        (1I)(n) (0I)(n)
    ]
    out, inc = _gstcMatrix(fs,
        _set_electric_terms(fs, rsₑ, R̃ₑ, kz₁, kz₂, ωϵ₁, ωμ₁, ωϵ₂, ωμ₂, Ω₁, Ω₂, C)...,
        _set_magnetic_terms(fs, rsₘ, R̃ₘ, kz₁, kz₂, ωϵ₁, ωμ₁, ωϵ₂, ωμ₂, Ω₁, Ω₂, C)...,
    )
    out\inc
end
function __sMatrix(fs, pw, rsₑ, rsₘ, sheet, um::T, ::T) where T<:UniformMedium
    Rₑ = _get_electric_response_components(rsₑ, pw, sheet)
    Rₘ = _get_magnetic_response_components(rsₘ, pw, sheet)
    R̃ₑ = [
        diagFT(vec(Rₑ.xx))   diagFT(vec(Rₑ.xy))
        diagFT(vec(Rₑ.yx))   diagFT(vec(Rₑ.yy))
    ]
    R̃ₘ = [
        diagFT(vec(Rₘ.xx))   diagFT(vec(Rₘ.xy))
        diagFT(vec(Rₘ.yx))   diagFT(vec(Rₘ.yy))
    ]
    ωϵ, ωμ = pw.ω .* _get_ϵμ(um)
    kz = Diagonal(repeat(vec(_get_kz(pw, um)), 2))
    _K = _get_K_components(pw)
    K = [
        Diagonal(vec(_K.xx)) Diagonal(vec(_K.xy))
        Diagonal(vec(_K.yx)) Diagonal(vec(_K.yy))
    ]
    Ω = K - (ωϵ * ωμ)*I
    n = length(pw)
    C = [
        (0I)(n) (-I)(n)
        (1I)(n) (0I)(n)
    ]
    out, inc = _gstcMatrix(fs,
        _set_electric_terms(fs, rsₑ, R̃ₑ, kz, ωϵ, ωμ, Ω, C)...,
        _set_magnetic_terms(fs, rsₘ, R̃ₘ, kz, ωϵ, ωμ, Ω, C)...,
    )
    out\inc
end

for (kind, ksym) in ((:electric, :ₑ), (:magnetic, :ₘ)), (rs, rsym) in ((:Impedance, :Z), (:Admittance, :Y))
    @eval function $(Symbol(:_get_, kind, :_response_components))(::$rs, pw, sheet)
        (
            xx = $(Symbol(rsym, ksym, :ˣˣ)).(Ref(sheet), Iterators.product(pw.x⃗...)),
            xy = $(Symbol(rsym, ksym, :ˣʸ)).(Ref(sheet), Iterators.product(pw.x⃗...)),
            yx = $(Symbol(rsym, ksym, :ʸˣ)).(Ref(sheet), Iterators.product(pw.x⃗...)),
            yy = $(Symbol(rsym, ksym, :ʸʸ)).(Ref(sheet), Iterators.product(pw.x⃗...)),
        )
    end
end

function _get_K_components(pw::PlaneWaves{2})
    (
        xx = [k⃗[1] .* k⃗[1] for k⃗ in Iterators.product(pw.k⃗...)],
        xy = [k⃗[1] .* k⃗[2] for k⃗ in Iterators.product(pw.k⃗...)],
        yx = [k⃗[2] .* k⃗[1] for k⃗ in Iterators.product(pw.k⃗...)],
        yy = [k⃗[2] .* k⃗[2] for k⃗ in Iterators.product(pw.k⃗...)],
    )
end

_set_electric_terms(::HField, ::Impedance,  R̃, kz₁, kz₂, ωϵ₁, ωμ₁, ωϵ₂, ωμ₂, Ω₁, Ω₂, C) = (Ω₁ * Ω₂ * -C * R̃ * C, Ω₂ * (0.5ωμ₁ * kz₁), Ω₁ * (0.5ωμ₂ * kz₂))
_set_electric_terms(::HField, ::Admittance, R̃, kz₁, kz₂, ωϵ₁, ωμ₁, ωϵ₂, ωμ₂, Ω₁, Ω₂, C) = (I, -C * R̃ * C * Ω₁ * inv(2ωϵ₁ * kz₁), -C * R̃ * C * Ω₁ * inv(2ωϵ₂ * kz₂))
_set_magnetic_terms(::HField, ::Impedance,  R̃, kz₁, kz₂, ωϵ₁, ωμ₁, ωϵ₂, ωμ₂, Ω₁, Ω₂, C) = (R̃ * -C * Ω₁ * C * inv(2ωϵ₁ * kz₁), R̃ * -C * Ω₂ * C * inv(0.5ωϵ₂ * kz₂), I)
_set_magnetic_terms(::HField, ::Admittance, R̃, kz₁, kz₂, ωϵ₁, ωμ₁, ωϵ₂, ωμ₂, Ω₁, Ω₂, C) = (Ω₂ * (2ωμ₁ * kz₁), Ω₁ * (2ωμ₂ * kz₂), Ω₁ * Ω₂ * R̃)

_set_electric_terms(::HField, ::Impedance,  R̃, kz, ωϵ, ωμ, Ω, C) = (Ω * -C * R̃ * C, 0.5ωμ * kz)
_set_electric_terms(::HField, ::Admittance, R̃, kz, ωϵ, ωμ, Ω, C) = (I, -C * R̃ * C * Ω * inv(2ωϵ * kz))
_set_magnetic_terms(::HField, ::Impedance,  R̃, kz, ωϵ, ωμ, Ω, C) = (R̃ * -C * Ω * C * inv(0.5ωϵ * kz), I)
_set_magnetic_terms(::HField, ::Admittance, R̃, kz, ωϵ, ωμ, Ω, C) = (2ωμ * kz, Ω * R̃)

_set_electric_terms(::EField, rs, R̃, kz, ωϵ₁, ωμ₁, ωϵ₂, ωμ₂, Ω₁, Ω₂, C) = _set_magnetic_terms(HField(), rs, R̃, kz, ωμ₁, ωϵ₁, ωμ₂, ωϵ₂, Ω₁, Ω₂, C)
_set_magnetic_terms(::EField, rs, R̃, kz, ωϵ₁, ωμ₁, ωϵ₂, ωμ₂, Ω₁, Ω₂, C) = _set_electric_terms(HField(), rs, R̃, kz, ωμ₁, ωϵ₁, ωμ₂, ωϵ₂, Ω₁, Ω₂, C)

_set_electric_terms(::EField, rs, R̃, kz, ωϵ, ωμ, Ω, C) = _set_magnetic_terms(HField(), rs, R̃, kz, ωμ, ωϵ, Ω, C)
_set_magnetic_terms(::EField, rs, R̃, kz, ωϵ, ωμ, Ω, C) = _set_electric_terms(HField(), rs, R̃, kz, ωμ, ωϵ, Ω, C)