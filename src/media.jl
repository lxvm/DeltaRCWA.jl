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
# v1: ERROR: LoadError: Mutating arrays is not supported -- called copyto!(::Vector{Tuple{Float64}}, _...)
# _get_kz(k⃗::Tuple{Vararg{AbstractVector}}, k⃗²) = _get_kz.(Iterators.product(k⃗...), k⃗²)
# _get_kz(k⃗::Tuple{Vararg{Number}}, k⃗²) = sqrt(Complex(k⃗² - mapreduce(abs2, +, k⃗, init=zero(k⃗²))))
# v2, works with Zygote!
# _get_kz(k⃗::Tuple{Vararg{AbstractVector}}, k⃗²) = reshape(_get_kz.(k⃗², map(x -> repeat(x, Int(prod(length.(k⃗))//length(x))), k⃗)...), length.(k⃗))
# _get_kz(k⃗², k⃗...) = sqrt(Complex(k⃗² - mapreduce(abs2, +, k⃗, init=zero(k⃗²))))
# v3, works with Zygote!
_get_kz(k⃗::Tuple{Vararg{AbstractVector}}, k⃗²) = map(_get_kz, Iterators.product(k⃗...), fill(k⃗², length.(k⃗)))
_get_kz(k⃗::Tuple{Vararg{Number}}, k⃗²) = sqrt(Complex(k⃗² - mapreduce(abs2, +, k⃗, init=zero(k⃗²))))
# v4: ERROR: LoadError: MethodError: no method matching (::ChainRulesCore.ProjectTo{Float64, NamedTuple{(), Tuple{}}})(::Array{Float64, 0})
# _get_kz(k⃗::Tuple{Vararg{AbstractVector}}, k⃗²) = map(_get_kz, Iterators.product(k⃗², k⃗...))
# _get_kz(k⃗) = sqrt(Complex(k⃗[1] - mapreduce(abs2, +, k⃗[2:end], init=zero(k⃗[1]))))
# v5: ERROR: LoadError: MethodError: no method matching (::ChainRulesCore.ProjectTo{Float64, NamedTuple{(), Tuple{}}})(::Array{Float64, 0})
# _get_kz(k⃗::Tuple{Vararg{AbstractVector}}, k⃗²) = map(Base.splat(_get_kz), Iterators.product(k⃗², k⃗...))
# _get_kz(k⃗², k⃗...) = sqrt(Complex(k⃗² - mapreduce(abs2, +, k⃗, init=zero(k⃗²))))
# v6: ERROR: LoadError: MethodError: no method matching (::ChainRulesCore.ProjectTo{Float64, NamedTuple{(), Tuple{}}})(::Array{Float64, 0})
# _get_kz(k⃗::Tuple{Vararg{AbstractVector}}, k⃗²) = map(Base.splat(_get_kz), Iterators.product(Iterators.product(k⃗...), k⃗²))
# _get_kz(k⃗::Tuple{Vararg{Number}}, k⃗²) = sqrt(Complex(k⃗² - mapreduce(abs2, +, k⃗, init=zero(k⃗²))))
