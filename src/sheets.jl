# Defines the trivial fall-back methods that should be implemented by `RCWASheet`s
export σₑˣˣ, σₑˣʸ, σₑʸˣ, σₑʸʸ, σₘˣˣ, σₘˣʸ, σₘʸˣ, σₘʸʸ

function nonconducting(::RCWASheet{N}, x⃗::NTuple{N, StepRangeLen})::BitArray{N} where N
    BitArray(false for e⃗ in Iterators.product(x⃗...))
end

σₑˣˣ(sheet, x⃗) = nonconducting(sheet, x⃗)
σₑˣʸ(sheet, x⃗) = nonconducting(sheet, x⃗)
σₑʸˣ(sheet, x⃗) = nonconducting(sheet, x⃗)
σₑʸʸ(sheet, x⃗) = nonconducting(sheet, x⃗)
σₘˣˣ(sheet, x⃗) = nonconducting(sheet, x⃗)
σₘˣʸ(sheet, x⃗) = nonconducting(sheet, x⃗)
σₘʸˣ(sheet, x⃗) = nonconducting(sheet, x⃗)
σₘʸʸ(sheet, x⃗) = nonconducting(sheet, x⃗)