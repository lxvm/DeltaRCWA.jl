module deltaRCWA

    using Reexport

    @reexport using LinearAlgebra

    @reexport using BlockArrays

    export smat_star, rred_star

    include("smat_star.jl")
    include("rred_star.jl")

end
