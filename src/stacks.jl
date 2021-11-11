export SheetStack

"""
A simple tuple of sheets with a complementary tuple of the spacings between
sheets in a UniformMedium
"""
struct SheetStack{N, L, T<:Tuple{RCWASheet{N}, Vararg{RCWASheet{N}, L}}} <: RCWAStack{N}
    sheets::T
    depths::Tuple{Vararg{Float64, L}}
end

for method in (:_sMatrix, :_sBlockMatrix, :_sLinearMap, :_sBlockMap)
    @eval function $method(stack::SheetStack, modes, pol)
        S = $method(stack.sheets[1], modes, pol)
        for i in eachindex(stack.depths)
            propagator = $method(UniformSlab(stack.depths[i]), modes, pol)
            S = S ⋆ₛ (propagator ⋆ₛ $method(stack.sheets[i+1], modes, pol))
        end
        S
    end
end