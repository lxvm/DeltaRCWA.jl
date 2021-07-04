module ScatteringTypes

export AbstractScatteringLayer, PeriodicScatteringLayer, UniformLayer,
    XYPeriodicDeltaLayer

include("layers.jl")

export AbstractScatteringStructure, PeriodicScatteringStructure,
    XYPeriodicScatteringStructure

include("structures.jl")

export AbstractScatteringProblem, LinearScatteringProblem, RCWAProblem

include("problems.jl")

end