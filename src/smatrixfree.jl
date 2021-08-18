export smatrixfree

"""
    smatrixfree(sheet, modes, pol)

Return a function that will compute the matrix-vector product with the
scattering matrix of a RCWASheet in a matrix-free fashion.
"""
function smatrixfree(sheet::RCWASheet{T, 1} where T, modes, pol::UncoupledPolarization)
    _1Dsheetsmatrixfree(_get_params_smatrixfree(sheet, modes, pol)...)
end

function _get_params_smatrixfree(sheet, modes, pol)
    pol_ω, σˣˣ, σʸʸ = _getσωbypol(sheet, modes, pol)
    halfσˣˣ = 0.5Diagonal(reshape(σˣˣ, :))
    halfσʸʸ = 0.5Diagonal(reshape(σʸʸ, :))
    pol_kz = Diagonal(reshape(modes.kz, :)) / pol_ω
    halfσˣˣ, halfσʸʸ, pol_kz
end

function _1Dsheetsmatrixfree(halfσˣˣ, halfσʸʸ, pol_kz)
    # closure: should return a function that will solve for the output
    scattered_side = _get_matrixfree_scattered_side(halfσˣˣ, halfσʸʸ, pol_kz)
    incident_side = _get_matrixfree_incident_side(halfσˣˣ, halfσʸʸ, pol_kz)
    return function smatrix_gstc(input::Vector{ComplexF64})
        x, info = linsolve(scattered_side, incident_side(input))
        return x
        if info.converged == 1
            throw(info)
        else
            return x
        end
    end
end

function _get_matrixfree_scattered_side(halfσˣˣ, halfσʸʸ, pol_kz)
    N = size(pol_kz, 1)
    return function (input::Vector{ComplexF64})
        I₁ = input[1:N]
        I₂ = input[(N+1):2N]
        sumI₁I₂ = I₁ + I₂
        diffI₁I₂ = I₁ - I₂
        vcat(
            diffI₁I₂ + ifft(halfσˣˣ * fft(pol_kz * diffI₁I₂)),
            (pol_kz * sumI₁I₂) + ifft(halfσʸʸ * fft(sumI₁I₂))
        )
    end
end

function _get_matrixfree_incident_side(halfσˣˣ, halfσʸʸ, pol_kz)
    N = size(pol_kz, 1)
    return function (input::Vector{ComplexF64})
        I₁ = input[1:N]
        I₂ = input[(N+1):2N]
        sumI₁I₂ = I₁ + I₂
        diffI₁I₂ = I₁ - I₂
        vcat(
            -diffI₁I₂ + ifft(halfσˣˣ * fft(pol_kz * diffI₁I₂)),
            (pol_kz * sumI₁I₂) - ifft(halfσʸʸ * fft(sumI₁I₂))
        )
    end
end