"pass in S and beta, -> outputs propagating modes only"
function get_prop_modes(S, β)
    i_real = findall(isreal, β)
    ip = [i_real; i_real.+length(β)]
    Sp = S[ip,ip]
    βp = β[i_real]
#    print(norm(Sp'Sp-I)/norm(Sp'Sp))
    
    return Sp, βp
end

