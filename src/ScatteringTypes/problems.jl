abstract type AbstractScatteringProblem end
abstract type LinearScatteringProblem <: AbstractScatteringProblem end

struct RCWAProblem <: LinearScatteringProblem
    structure::AbstractScatteringStructure
    modes::AbstractVector
end