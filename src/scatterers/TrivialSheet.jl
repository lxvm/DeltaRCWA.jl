struct TrivialSheet <: DeltaScatterer end

get_σᵉ(::TrivialSheet, kwargs...) = zeros(2, 2)
get_σᵐ(::TrivialSheet, kwargs...) = zeros(2, 2)