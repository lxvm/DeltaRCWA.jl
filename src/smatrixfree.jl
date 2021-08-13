export smatrixfree

"""
    smatrixfree(sheet, modes, pol)

Return a function that will compute the matrix-vector product with the
scattering matrix of a RCWASheet in a matrix-free fashion.
"""
function smatrixfree(sheet::RCWASheet{T, 1} where T, modes, pol::UncoupledPolarization)
    ω, σˣˣ, σʸʸ = _getσωbypol(sheet, modes, pol)
    halfσˣˣ = 0.5Diagonal(reshape(σˣˣ, :))
    halfσʸʸ = 0.5Diagonal(reshape(σʸʸ, :))
    pol_kz = Diagonal(reshape(modes.kz, :)) / ω
    _1Dsheetsmatrixfree(halfσˣˣ, halfσʸʸ, pol_kz)
end

function _1Dsheetsmatrixfree(halfσˣˣ, halfσʸʸ, pol_kz)
    # should return a function that will solve for the output
    N = size(pol_kz, 1)
    function Lgstc(input::Vector{ComplexF64})
        I₁ = input[1:N]
        I₂ = input[(N+1):2N]
        sumI₁I₂ = I₁ + I₂
        diffI₁I₂ = I₁ - I₂
        vcat(
            diffI₁I₂ + fft(halfσˣˣ * ifft(pol_kz * diffI₁I₂)),
            pol_kz * (sumI₁I₂) - fft(halfσʸʸ * ifft(sumI₁I₂))
        )
    end
    # return Lgstc
    function Rgstc(input::Vector{ComplexF64})
        I₁ = input[1:N]
        I₂ = input[(N+1):2N]
        sumI₁I₂ = I₁ + I₂
        diffI₂I₁ = I₂ - I₁
        vcat(
            diffI₂I₁ + fft(halfσˣˣ * ifft(pol_kz * diffI₂I₁)),
            pol_kz * (sumI₁I₂) + fft(halfσʸʸ * ifft(sumI₁I₂))
        )
    end
    # return Rgstc
    return function Sgstc(input::Vector{ComplexF64})
        x, info = linsolve(Lgstc, Rgstc(input), rand(ComplexF64, 2N))
        if info.converged == 1
            throw(info)
        else
            return x
        end
    end
end