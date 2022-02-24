export FieldStyle, PField, HField, EField

"""
    FieldStyle

Supertype for traits that describe kinds of electromagnetic fields. `EField` and
other types are subtypes of this.
"""
abstract type FieldStyle end

"""
    PField

Subtype of `FieldStyle` that describes fields normalized to unit power.
"""
struct PField <: FieldStyle end

"""
    HField

Subtype of `FieldStyle` that describes the auxiliary magnetic field.
"""
struct HField <: FieldStyle end

"""
    EField

Subtype of `FieldStyle` that describes the electric field.
"""
struct EField <: FieldStyle end

"""
    _gstcMatrix(::Union{HField, EField}, A, B, C, D, )

Converts gstcs `A(H² - H¹) = B(E² + E¹)` and `C(E² - E¹) = D(H² + H¹)` into a
corresponding scattering matrix for either the E or H field
"""
function _gstcMatrix(::HField, A, B₁, B₂, C₁, C₂, D)
    out = [
        -A + B₁    A - B₂
        C₁ - D     C₂ - D
    ]
    inc = [
        A + B₁     -A - B₂
        C₁ + D     C₂ + D
    ]
    out, inc
end
_gstcMatrix(fs::HField, A, B, C, D) = _gstcMatrix(fs, A, B, B, C, C, D)

function _gstcMatrix(::EField, A₁, A₂, B, C, D₁, D₂)
    out, inc = _gstcMatrix(HField(), C, D₁, D₂, A₁, A₂, B)
    n = Int(size(out, 1)//2)
    prop_sign = Diagonal(vcat(ones(n), -ones(n)))
    out * -prop_sign, inc * prop_sign
end
_gstcMatrix(fs::EField, A, B, C, D) = _gstcMatrix(fs, A, A, B, C, D, D)

function _get_halves(x::AbstractVector)
    n = Int(length(x)//2)
    I₁ = x[1:n]
    I₂ = x[(n+1):2n]
    I₁, I₂
end

function _gstcLinearMap(::HField, A, B₁, B₂, C₁, C₂, D)
    n = 2size(A, 1)
    out = LinearMap(n) do x
        I₁, I₂ = _get_halves(x)
        vcat(
            A*(-I₁ + I₂) + B₁*I₁ - B₂*I₂,
            C₁*I₁ + C₂*I₂ - D*(I₁ + I₂),
        )
    end
    inc = LinearMap(n) do x
        I₁, I₂ = _get_halves(x)
        vcat(
            A*(I₁ - I₂) + B₁*I₁ - B₂*I₂,
            C₁*I₁ + C₂*I₂ + D*(I₁ + I₂),
        )
    end
    out, inc
end

function _gstcLinearMap(::HField, A, B, C, D)
    n = 2size(A, 1)
    out = LinearMap(n) do x
        I₁, I₂ = _get_halves(x)
        sumI₁I₂ = I₁ + I₂
        diffI₁I₂ = I₁ - I₂
        vcat(
            -A*diffI₁I₂ + B*diffI₁I₂,
            C*sumI₁I₂ - D*sumI₁I₂,
        )
    end
    inc = LinearMap(n) do x
        I₁, I₂ = _get_halves(x)
        sumI₁I₂ = I₁ + I₂
        diffI₁I₂ = I₁ - I₂
        vcat(
            A*diffI₁I₂ + B*diffI₁I₂,
            C*sumI₁I₂ + D*sumI₁I₂,
        )
    end
    out, inc
end

function _gstcLinearMap(::EField, A₁, A₂, B, C, D₁, D₂)
    out, inc = _gstcLinearMap(HField(), C, D₁, D₂, A₁, A₂, B)
    n = Int(size(out, 1)//2)
    prop_sign = Diagonal(vcat(ones(n), -ones(n)))
    out * -prop_sign, inc * prop_sign
end

function _gstcLinearMap(::EField, A, B, C, D)
    out, inc = _gstcLinearMap(HField(), C, D, A, B)
    n = Int(size(out, 1)//2)
    prop_sign = Diagonal(vcat(ones(n), -ones(n)))
    out * -prop_sign, inc * prop_sign
end