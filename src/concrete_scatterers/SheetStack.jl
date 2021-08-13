export SheetStack

"""
A simple tuple of sheets with a complementary tuple of the spacings between
sheets in a UniformMedium
"""
struct SheetStack{T, N, L} <: RCWAScatterer{T, N}
    sheets::Tuple{RCWASheet{T, N}, Vararg{RCWASheet{T, N}, L}}
    depths::Tuple{Vararg{Float64, L}}
end

function smatrix(stack::SheetStack, modes, pol)
    S = smatrix(stack.sheets[1], modes, pol)
    for i in eachindex(stack.depths)
        S = smat_star(S, smat_star(
            smatrix(UniformSlab(stack.depths[i]), modes, pol),
            smatrix(stack.sheets[i+1], modes, pol)
        ))
    end
    S
end