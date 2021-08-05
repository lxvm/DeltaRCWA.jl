export SheetStack

"""
A simple tuple of sheets with a complementary tuple of the spacings between
sheets in a UniformMedium
"""
struct SheetStack{N, L, D} <: RCWAScatterer{N, 3}
    sheets::NTuple{L, RCWASheet{N}} # this is a tuple of pointers
    depths::NTuple{D, <: Real} # this is a tuple of a concrete subtype of real
    function SheetStack(sheets::NTuple{L, RCWASheet{N}}, depths::NTuple{D, <: Real}) where {N, L, D}
        L-1==D || throw(ArgumentError("there must be one gap less than the number of sheets"))
        new{N, L, D}(sheets, depths)
    end
end

function smatrix(stack::SheetStack, modes::PlanewaveModes, pol::AbstractPolarization)::BlockMatrix
    S = smatrix(stack.sheets[1], modes, pol)
    for i in eachindex(stack.depths)
        S = smat_star(S, smat_star(
            smatrix(UniformSlab(stack.depths[i]), modes, pol),
            smatrix(stack.sheets[i+1], modes, pol)
        ))
    end
    S
end