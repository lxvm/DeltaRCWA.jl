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
