function ChainRulesCore.rrule(::Type{<:Frequencies}, n_nonnegative, n, multiplier)
  y = Frequencies(n_nonnegative, n, multiplier)
  function Frequencies_pullback(ȳ)
    x̄ = transpose(Frequencies(n_nonnegative, n, 1)) * ȳ
    NoTangent(), NoTangent(), NoTangent(), x̄
  end
  return y, Frequencies_pullback
end