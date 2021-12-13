import FFTW: bfft

function interpolate(x, y₁, y₂, RCWAsol)
    sol = RCWAsol
    @assert x == sol.x "grid doesn't match, so likely incorrect"
    O₁ = rotr90(bfft(exp.(-sol.kz₁ * transpose(im * y₁)) .* sol.Õ₁, 1))
    I₂ = rotr90(bfft(exp.(-sol.kz₂ * transpose(im * y₂)) .* sol.Ĩ₂, 1))
    I₁ = rotr90(bfft(exp.( sol.kz₁ * transpose(im * y₁)) .* sol.Ĩ₁, 1))
    O₂ = rotr90(bfft(exp.( sol.kz₂ * transpose(im * y₂)) .* sol.Õ₂, 1))
    (x=x, y₁=y₁, y₂=y₂, I₁=I₁, O₁=O₁, I₂=I₂, O₂=O₂)
end