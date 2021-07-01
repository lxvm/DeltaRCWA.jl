"""
pass in array of M and N, L, lambda -> outputs S matrix and array of betas

assumes epsilon and mu are 1 
L is the width of the cell
λ is in units of L
M, N are surface permitvities, expressed as a diagonal matrix
"""
function get_all_modes(M, N, L, λ)
    ω₀ = 2*pi/λ
    k = ω₀
    nvec = 0:length(M)^0.5-1 # 100 modes
    kₓ = 2*pi*fftfreq(length(nvec), length(nvec)/L) 
    β = @. sqrt(Complex(ω₀^2 - kₓ^2))
    
    A = [Diagonal(β)+k*ifft(fft(M, 2), 1)    Diagonal(β)+k*ifft(fft(M, 2), 1);
        Diagonal(β)+k*ifft(fft(N, 2), 1)     -Diagonal(β)-k*ifft(fft(N, 2), 1)]
    B = [Diagonal(β)-k*ifft(fft(M, 2), 1)    Diagonal(β)-k*ifft(fft(M, 2), 1);
        Diagonal(β)-k*ifft(fft(N, 2), 1)      -Diagonal(β)+k*ifft(fft(N, 2), 1)]

    # TODO: replace inversion with iterative solver
    S = A\B
    return S, β
end
