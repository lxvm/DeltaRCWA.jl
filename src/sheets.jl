# Defines the trivial fall-back methods that should be implemented by `RCWASheet`s
export σₑˣˣ, σₑˣʸ, σₑʸˣ, σₑʸʸ, σₘˣˣ, σₘˣʸ, σₘʸˣ, σₘʸʸ

nonconducting(x⃗) = zeros(Bool, length.(x⃗))

σₑˣˣ(::RCWASheet, x⃗) = nonconducting(x⃗)
σₑˣʸ(::RCWASheet, x⃗) = nonconducting(x⃗)
σₑʸˣ(::RCWASheet, x⃗) = nonconducting(x⃗)
σₑʸʸ(::RCWASheet, x⃗) = nonconducting(x⃗)
σₘˣˣ(::RCWASheet, x⃗) = nonconducting(x⃗)
σₘˣʸ(::RCWASheet, x⃗) = nonconducting(x⃗)
σₘʸˣ(::RCWASheet, x⃗) = nonconducting(x⃗)
σₘʸʸ(::RCWASheet, x⃗) = nonconducting(x⃗)