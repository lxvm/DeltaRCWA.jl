export SheetStack

"""
A simple tuple of sheets with a complementary tuple of the spacings between
sheets in a UniformMedium
"""
struct SheetStack{N, L} <: RCWAStack{N}
    sheets::Tuple{RCWASheet{N}, Vararg{RCWASheet{N}, L}}
    depths::Tuple{Vararg{Float64, L}}
end

for method in (:smatrix, :smatrixBlockMap, :smatrixLinearMap)
    @eval function $(method)(stack::SheetStack, modes, pol)
        S = $(method)(stack.sheets[1], modes, pol)
        for i in eachindex(stack.depths)
            propagator = $(method)(UniformSlab(stack.depths[i]), modes, pol)
            S = S ⋆ (propagator ⋆ $(method)(stack.sheets[i+1], modes, pol))
        end
        S
    end
end