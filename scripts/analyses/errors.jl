import LinearAlgebra: norm

function get_errors(BIEsol, RCWAsol)
    x, y₁, y₂ = BIEsol.x, BIEsol.y₁, BIEsol.y₂
    err = [
        getproperty(BIEsol, sym) - getproperty(RCWAsol, sym)
        for sym in (:I₁, :O₁, :I₂, :O₂)
    ]
    (x=x, y₁=y₁, y₂=y₂, I₁=err[1], O₁=err[2], I₂=err[3], O₂=err[4])
end

function L2cubature(sol; dims)
    margL2 = [
        mapslices(norm, e; dims=dims) ./ sqrt(prod(size(e)[dims]))
        for e in (sol.I₁, sol.O₁, sol.I₂, sol.O₂)
    ]
    (I₁=margL2[1], O₁=margL2[2], I₂=margL2[3], O₂=margL2[4])
end

function errorL2(solerrors)
    err = L2cubature(solerrors; dims=1:2)
    (I₁=err.I₁[1], O₁=err.O₁[1], I₂=err.I₂[1], O₂=err.O₂[1])
end
function errorbydistance(solerrors)
    err = L2cubature(solerrors; dims=1)
    (y₁=solerrors.y₁, y₂=solerrors.y₂, I₁=vec(err.I₁), O₁=vec(err.O₁), I₂=vec(err.I₂), O₂=vec(err.O₂))
end