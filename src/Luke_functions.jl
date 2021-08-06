export get_all_modes, get_prop_modes, get_transmissivity_normalincident

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

"pass in S and beta, -> outputs propagating modes only"
function get_prop_modes(S, β)
    i_real = findall(isreal, β)
    ip = [i_real; i_real.+length(β)]
    Sp = S[ip,ip]
    βp = β[i_real]
#    print(norm(Sp'Sp-I)/norm(Sp'Sp))
    
    return Sp, βp
end

"""
compute transmission for kx = 0 (normal incident transitivity)

assumes epsilon and mu are 1 
L is the width of the cell
λ is in units of L
M, N are surface permitvities, expressed as a diagonal matrix
"""
function get_transmissivity_normalincident(M, N, L, λ, u_in₁, u_in₂)
    
    S, β = get_all_modes(M, N, L, λ)
    u_in = [u_in₁; u_in₂]
    u_out = S*u_in
    u_out₁ = u_out[1:length(u_in₁)]
    u_out₂ = u_out[length(u_in₁)+1:length(u_out)]
    
    T = (abs(u_out₂[1])/abs(u_in₁[1]))^2
    phase = angle(u_out₂[1]) - angle(u_in₁[1])
    
    return T, phase
end