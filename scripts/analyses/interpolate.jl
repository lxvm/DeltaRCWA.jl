import FFTW: bfft, fftfreq

to_position_basis(sim) = to_position_basis(sim.method, sim)
to_position_basis(::typeof(compute_BIE_method), sim) = sim
function to_position_basis(::Union{typeof(compute_DeltaRCWA_method), typeof(compute_RCWA_method)}, sim)
    (;method, param, sheet, sol) = sim
    x = range(0.0, step=sheet.L/param.n, length=param.n)
    y = range(-sheet.L, sheet.L, length=2param.n)
    y₁ = y[y .< 0] # y points in domain below sheet
    y₂ = y[y .> 0] # y points in domain above sheet
    (method=method, param=param, sheet=sheet, sol=interpolate(x, y₁, y₂, sol))
end

function interpolate(x, y₁, y₂, RCWAsol)
    sol = RCWAsol
    @assert x == sol.x "grid doesn't match, so likely incorrect"
    # O₂ = zeros(ComplexF64, length.((y₂, x)))
    # I₁ = zeros(ComplexF64, length.((y₁, x)))
    # O₁ = zeros(ComplexF64, length.((y₁, x)))
    # I₂ = zeros(ComplexF64, length.((y₂, x)))
    # kx = 2π * fftfreq(length(x), 1/(x[2] - x[1]))
    # for j in eachindex(x)
    #     for i in eachindex(y₂)
    #         O₁[i, j] = sum(exp.((y₁[i] * im) .* -sol.kz₁ .+ (x[j] * im) .* kx) .* sol.Õ₁)
    #         I₂[i, j] = sum(exp.((y₂[i] * im) .* -sol.kz₂ .+ (x[j] * im) .* kx) .* sol.Ĩ₂)
    #         I₁[i, j] = sum(exp.((y₁[i] * im) .*  sol.kz₁ .+ (x[j] * im) .* kx) .* sol.Ĩ₁)
    #         O₂[i, j] = sum(exp.((y₂[i] * im) .*  sol.kz₂ .+ (x[j] * im) .* kx) .* sol.Õ₂)
    #     end
    # end
    O₁ = rotr90(bfft(exp.(-sol.kz₁ * transpose(im * y₁)) .* sol.Õ₁, 1))
    I₂ = rotr90(bfft(exp.(-sol.kz₂ * transpose(im * y₂)) .* sol.Ĩ₂, 1))
    I₁ = rotr90(bfft(exp.( sol.kz₁ * transpose(im * y₁)) .* sol.Ĩ₁, 1))
    O₂ = rotr90(bfft(exp.( sol.kz₂ * transpose(im * y₂)) .* sol.Õ₂, 1))
    (x=x, y₁=y₁, y₂=y₂, I₁=I₁, O₁=O₁, I₂=I₂, O₂=O₂)
end