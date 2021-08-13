export RCWAplot

"""
    plot(::DeltaRCWASolution{T₁, T₂, 1})

Visualize the RCWA solution in terms of the incident and scattered fields

Keyword arguments:
- combine{Bool} :: whether or not to plot the total fields or their components
- method{Symbol} :: choose between :series and :fft (default, fallback)
- part{Function} :: A function that gives a real from a complex (default: real)
"""
RCWAplot

@recipe function f(sol::DeltaRCWASolution{T, 1} where T; combine=false, part=real, method=:fft)
    # Nz = sol.modes.dims[1][1]
    # dz = sol.modes.dims[1][2] / Nz
    # z⃗ = range(0, step=dz, length=Nz)
    x⃗ = sol.modes.x⃗[1]
    z⃗₂ = x⃗
    z⃗₁ = -reverse(z⃗₂)
    kz = sol.modes.kz
    if method == :series
        AO₂ = zeros(length.((x⃗, z⃗₂)))
        AI₁ = zeros(length.((x⃗, z⃗₁)))
        AO₁ = zeros(length.((x⃗, z⃗₁)))
        AI₂ = zeros(length.((x⃗, z⃗₂)))
        k⃗ = sol.modes.k⃗[1]
        for j in eachindex(x⃗)
            for i in eachindex(z⃗₂)
                AO₂[j, i] = part.(sum(exp.((z⃗₂[i] * im) .*  kz .+ (x⃗[j] * im) .* k⃗) .* sol.O₂))
                AI₂[j, i] = part.(sum(exp.((z⃗₂[i] * im) .* -kz .+ (x⃗[j] * im) .* k⃗) .* sol.I₂))
                AI₁[j, i] = part.(sum(exp.((z⃗₁[i] * im) .*  kz .+ (x⃗[j] * im) .* k⃗) .* sol.I₁))
                AO₁[j, i] = part.(sum(exp.((z⃗₁[i] * im) .* -kz .+ (x⃗[j] * im) .* k⃗) .* sol.O₁))
            end
        end
    else
        # using bfft here because it matches the summation
        # but in reality ifft (resp. bfft) are off by a factor of 
        # sqrt(length(x⃗)) (resp. 1/sqrt(length(x⃗))) when it comes to preserving
        # the norm/energy of the frequency domain representation.
        # Note that fft doesn't preserve the norm (increases by sqrt(length(x⃗)))
        # Ultimately this is just a choice of convention that can be modified.
        AO₂ = part.(bfft(exp.( kz * transpose(im * z⃗₂)) .* sol.O₂, 1))
        AI₂ = part.(bfft(exp.(-kz * transpose(im * z⃗₂)) .* sol.I₂, 1))
        AI₁ = part.(bfft(exp.( kz * transpose(im * z⃗₁)) .* sol.I₁, 1))
        AO₁ = part.(bfft(exp.(-kz * transpose(im * z⃗₁)) .* sol.O₁, 1))
    end
    if combine
       @series begin
            seriestype := :heatmap
            title := "port 1 | port 2"
            xguide := "z"
            yguide := "x"
            vcat(z⃗₁, z⃗₂), x⃗, cat(AO₁ .+ AI₁, AO₂ .+ AI₂; dims=2)
        end
    else
        layout := 4
        # plot_title := String(part) * " part of $(typeof(sol.pol)) field"
        # subplots are in row-major order
        @series begin
            seriestype := :heatmap
            subplot := 1
            title := "Scattered, port 1"
            xaxis := nothing
            yguide := "x"
            z⃗₁, x⃗, AO₁
        end
        @series begin
            seriestype := :heatmap
            subplot := 2
            title := "Incident, port 2"
            xaxis := nothing
            yaxis := nothing
            z⃗₂, x⃗, AI₂
        end
        @series begin
            seriestype := :heatmap
            subplot := 3
            title := "Incident, port 1"
            xguide := "z"
            yguide := "x"
            z⃗₁, x⃗, AI₁
        end
        @series begin
            seriestype := :heatmap
            subplot := 4
            title := "Scattered, port 2"
            yaxis := nothing
            xguide := "z"
            z⃗₂, x⃗, AO₂
        end
    end
end