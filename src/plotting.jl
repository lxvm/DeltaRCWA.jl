# examples: https://github.com/JuliaPlots/ExamplePlots.jl/blob/master/notebooks/usertype_recipes.ipynb

@recipe function f(sol::DeltaRCWASolution{1}; part=:real)
    # Nz = sol.modes.dims[1][1]
    # dz = sol.modes.dims[1][2] / Nz
    # z⃗ = range(0, step=dz, length=Nz)
    x⃗ = sol.modes.x⃗[1]
    z⃗₂ = x⃗
    z⃗₁ = -reverse(z⃗₂)
    AO₂ = zeros(length.((x⃗, z⃗₂)))
    AI₁ = zeros(length.((x⃗, z⃗₁)))
    AO₁ = zeros(length.((x⃗, z⃗₁)))
    AI₂ = zeros(length.((x⃗, z⃗₂)))
    for i in eachindex(z⃗₂)
        # AO₂ = eval(part).(ifft(exp.(kron(im * z⃗₂, sol.modes.kz)), dims=1))
        AO₂[:, i] = eval(part).(ifft(exp.((z⃗₂[i] * im) .* sol.modes.kz) .* sol.O₂))
        AI₂[:, i] = eval(part).(ifft(exp.((z⃗₂[i] * im) .* sol.modes.kz) .* sol.I₂))
        AI₁[:, i] = eval(part).(ifft(exp.((z⃗₁[i] * im) .* sol.modes.kz) .* sol.I₁))
        AO₁[:, i] = eval(part).(ifft(exp.((z⃗₁[i] * im) .* sol.modes.kz) .* sol.O₁))
    end
    layout := 4
    # plot_title := String(part) * " part of $(typeof(sol.pol)) field"
    # subplots are in row-major order
    @series begin
        seriestype := :heatmap
        subplot := 1
        title := "Scattered, port 2"
        x⃗, z⃗₂, transpose(AO₂)
    end
    @series begin
        seriestype := :heatmap
        subplot := 2
        title := "Incident, port 2"
        x⃗, z⃗₂, transpose(AI₂)
    end
    @series begin
        seriestype := :heatmap
        subplot := 3
        title := "Incident, port 1"
        x⃗, z⃗₁, transpose(AI₁)
    end
    @series begin
        seriestype := :heatmap
        subplot := 4
        title := "Scattered, port 1"
        x⃗, z⃗₁, transpose(AO₁)
    end
end