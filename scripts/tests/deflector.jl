# run this script in the "scripts" environment
using JLD2

include("../Carlos_method.jl")
include("../Luke_method.jl")
include("../Lorenzo_method.jl")

include("../deflector_def.jl")

jldopen("scripts/deflector.jld2", "a+") do file
    for e in exponents
        N = 2^e
        prefix = "deflector/k$(k)/n$(N)/i$(i)/"
        file[prefix * "BIE"] = compute_BIE_method(η, μ, L, N, k, i)
        file[prefix * "RCWA"] = compute_RCWA_method(η, μ, L, N, k, i)
        file[prefix * "DeltaRCWA"] = compute_DeltaRCWA_method(η, μ, L, N, k, i)
    end
end