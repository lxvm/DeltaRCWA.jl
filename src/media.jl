export UniformMedium, Vacuum

"""
    UniformMedium{ϵ, μ}()

Represents a medium with a constant permitivitty and permeability
Units are such that ϵ, μ are relative to the vacuum value(ϵ=μ=1)
"""
struct UniformMedium{ϵ, μ}
    function UniformMedium{ϵ, μ}() where {ϵ, μ}
        @assert ϵ isa Number
        @assert μ isa Number
        new{ϵ, μ}()
    end
end

const Vacuum = UniformMedium{1, 1}()

_get_ϵμ(::UniformMedium{ϵ, μ}) where {ϵ, μ} = (ϵ, μ)

"""
    _get_kz(k⃗::Tuple{Vararg{Number}}, k⃗²)
    _get_kz(k⃗::Tuple{Vararg{AbstractVector}}, k⃗²)
    _get_kz(k⃗::Tuple{Vararg{AbstractVector}}, ω, ϵ, μ)
    _get_kz(::PlaneWaves, ::UniformMedium)

Returns an array of kz from the dispersion relation k⃗² = k⃗⋅k⃗ = ω²ϵμ.
"""
_get_kz(pw::PlaneWaves, um::UniformMedium) = _get_kz(pw.k⃗, pw.ω, _get_ϵμ(um)...)
_get_kz(k⃗::Tuple{Vararg{AbstractVector}}, ω, ϵ, μ) = _get_kz(k⃗, ω^2 * ϵ * μ)
_get_kz(k⃗::Tuple{Vararg{AbstractVector}}, k⃗²) = _get_kz.(Iterators.product(k⃗...), k⃗²)
_get_kz(k⃗::Tuple{Vararg{Number}}, k⃗²) = sqrt(Complex(k⃗² - mapreduce(abs2, +, k⃗, init=zero(k⃗²))))