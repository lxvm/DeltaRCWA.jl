using Plots

function show_params(sim)
    (; sheet, param) = sim
    x = range(0.0, sheet.L, length=(param.n+1))
    electric_data = sheet.μ.(x)
    magnetic_data = sheet.η.(x)
    out = plot(title="Sheet parameters: " * sheet.name, xguide="Unit cell position", yguide="Sheet response")
    plot!(x, real.(electric_data), label="Re(Zₑˣˣ)")
    plot!(x, imag.(electric_data), label="Im(Zₑˣˣ)", ls=:dash)
    plot!(x, real.(magnetic_data), label="Re(Yₘʸʸ)")
    plot!(x, imag.(magnetic_data), label="Im(Yₘʸʸ)", ls=:dash)
    out
end

function show_fields(sim)
    (; sol) = sim
    heatmap(sol.x, vcat(sol.y₁, sol.y₂),
        real.(cat(sol.I₁ .+ sol.O₁, sol.I₂ .+ sol.O₂; dims=1)),
        title=string("Total Hₓ field ", name_sim(sim)), xguide="x", yguide="y",
        aspect_ratio=1, color=:RdBu,
    )
end

function show_errorbydistance(solerrors)
    err = errorbydistance(solerrors)
    out = plot(title="L2 error of fields by distance from sheet", xguide="y", yguide="L2 norm")
    plot!(err.y₁, vec(err.I₁), label="I₁")
    plot!(err.y₁, vec(err.O₁), label="O₁")
    plot!(err.y₂, vec(err.I₂), label="I₂")
    plot!(err.y₂, vec(err.O₂), label="O₂")
    out
end

function show_errorL2(exponents, L2errors)
    out = plot(title="convergence of L2 errors of total fields", xguide="log2(grid points)", yguide="log2(L2 norm)")
    plot!(exponents, log2.(getproperty.(L2errors, :I₁)), label="I₁")
    plot!(exponents, log2.(getproperty.(L2errors, :O₁)), label="O₁")
    plot!(exponents, log2.(getproperty.(L2errors, :I₂)), label="I₂")
    plot!(exponents, log2.(getproperty.(L2errors, :O₂)), label="O₂")
    out
end