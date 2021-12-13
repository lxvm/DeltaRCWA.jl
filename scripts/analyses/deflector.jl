using JLD2
using Plots

include("show.jl")
include("errors.jl")
include("interpolate.jl")

include("../deflector_def.jl")

paramplots = Vector{Plots.Plot}(undef, 0)
fieldBIEplots = Vector{Plots.Plot}(undef, 0)
fieldRCWAplots = Vector{Plots.Plot}(undef, 0)
fieldDeltaRCWAplots = Vector{Plots.Plot}(undef, 0)
fielderrorplots = Vector{Plots.Plot}(undef, 0)
errorbydistanceplots = Vector{Plots.Plot}(undef, 0)
L2errors = Vector{Any}(undef, 0)
RCWAfielderrorplots = Vector{Plots.Plot}(undef, 0)
RCWAerrorbydistanceplots = Vector{Plots.Plot}(undef, 0)
RCWAL2errors = Vector{Any}(undef, 0)

jldopen("scripts/deflector.jld2", "r") do file
    for e in exponents
        N = 2^e
        prefix = "deflector/k$(k)/n$(N)/i$(i)/"
        BIEsol = file[prefix * "BIE"]
        _RCWAsol = file[prefix * "RCWA"]
        _DeltaRCWAsol = file[prefix * "DeltaRCWA"]
        x, y₁, y₂ = BIEsol.x, BIEsol.y₁, BIEsol.y₂
        RCWAsol = interpolate(x, y₁, y₂, _RCWAsol)
        DeltaRCWAsol = interpolate(x, y₁, y₂, _DeltaRCWAsol)
        push!(paramplots, show_params(η, μ, RCWAsol.x))
        push!(fieldBIEplots, show_fields(BIEsol))
        push!(fieldRCWAplots, show_fields(RCWAsol))
        push!(fieldDeltaRCWAplots, show_fields(DeltaRCWAsol))
        BIEerrors = get_errors(BIEsol, DeltaRCWAsol)
        RCWAerrors = get_errors(RCWAsol, DeltaRCWAsol)
        push!(fielderrorplots, show_fields(BIEerrors))
        push!(errorbydistanceplots, show_errorbydistance(BIEerrors))
        push!(L2errors, errorL2(BIEerrors))
        push!(RCWAfielderrorplots, show_fields(RCWAerrors))
        push!(RCWAerrorbydistanceplots, show_errorbydistance(RCWAerrors))
        push!(RCWAL2errors, errorL2(RCWAerrors))
    end
end

for (j, p) in enumerate(paramplots)
    savefig(p, "scripts/analyses/paramplot$(j).png")
end
for (j, p) in enumerate(fieldBIEplots)
    savefig(p, "scripts/analyses/fieldBIEplot$(j).png")
end
for (j, p) in enumerate(fieldRCWAplots)
    savefig(p, "scripts/analyses/fieldRCWAplot$(j).png")
end
for (j, p) in enumerate(fieldDeltaRCWAplots)
    savefig(p, "scripts/analyses/fieldDeltaRCWAplot$(j).png")
end
for (j, p) in enumerate(fielderrorplots)
    savefig(p, "scripts/analyses/fielderrorplot$(j).png")
end
for (j, p) in enumerate(errorbydistanceplots)
    savefig(p, "scripts/analyses/errorbydistanceplot$(j).png")
end
for (j, p) in enumerate(RCWAfielderrorplots)
    savefig(p, "scripts/analyses/RCWAfielderrorplot$(j).png")
end
for (j, p) in enumerate(RCWAerrorbydistanceplots)
    savefig(p, "scripts/analyses/RCWAerrorbydistanceplot$(j).png")
end
savefig(show_errorL2(exponents, L2errors), "scripts/analyses/L2error.png")
savefig(show_errorL2(exponents, RCWAL2errors), "scripts/analyses/RCWAL2error.png")