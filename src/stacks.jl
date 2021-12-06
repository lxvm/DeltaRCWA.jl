export SheetStack

"""
    SheetStack(::Sheet)
    SheetStack(::Sheet, ::NTuple{2, UniformMedium})
    SheetStack(::Tuple{Vararg{Float64, L}}, ::Tuple{Sheet, Vararg{Sheet, L}}) where L
    SheetStack(::Tuple{Vararg{Float64, L}}, ::Tuple{Sheet, Vararg{Sheet, L}}, ::Tuple{UniformMedium, UniformMedium, Vararg{UniformMedium, L}}) where L

A stack of `Sheet`s with distances between sheets and a `UniformMedium` for each
side of each `Sheet`
"""
struct SheetStack{L, S<:Tuple{Sheet, Vararg{Sheet, L}}, M<:Tuple{UniformMedium, UniformMedium, Vararg{UniformMedium, L}}}
    depths::Tuple{Vararg{Float64, L}}
    sheets::S
    media::M
end

SheetStack(sheet::Sheet) = SheetStack((), (sheet,))
SheetStack(sheet::Sheet, um::NTuple{2, UniformMedium}) = SheetStack((), (sheet,), um)
function SheetStack(depths::Tuple{Vararg{Float64, L}}, sheets::Tuple{Sheet, Vararg{Sheet, L}}) where L
    SheetStack(depths, sheets, Tuple(Vacuum for _ in 1:(length(depths)+2)))
end

Base.convert(::Type{T}, S::SheetStack) where T<:SheetStack = S
Base.convert(::Type{T}, S::Sheet) where T<:SheetStack = SheetStack(S)

for method in (:_sMatrix, :_sBlockMatrix, :_sLinearMap, :_sBlockMap)
    @eval function $method(fs, pw, stack::SheetStack)
        depths, sheets, media = stack.depths, stack.sheets, stack.media
        S = $method(fs, pw, sheets[1], media[1:2]...)
        for i in eachindex(depths)
            propagator = $method(pw, UniformSlab(depths[i]), media[i+1])
            S = S ⋆ₛ (propagator ⋆ₛ $method(fs, pw, sheets[i+1], media[(i+1):(i+2)]...))
        end
        S
    end
end