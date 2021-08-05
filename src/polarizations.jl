export TE, TM, Coupled

abstract type AbstractPolarization end
abstract type UncoupledPolarization <: AbstractPolarization end
abstract type CoupledPolarization <: AbstractPolarization end

struct TE <: UncoupledPolarization end
struct TM <: UncoupledPolarization end
struct Coupled <: CoupledPolarization end

"""
    enforce_N_pol(N::Int64, ::P where P <: AbstractPolarization)

Check that a polarization is specified when there are symmetries which decouple
the TE and TM polarizations
"""
function enforce_N_pol(N::Int64, ::P) where P <: AbstractPolarization
    if 0 <= N <= 1
        @assert P <: UncoupledPolarization "1D, 2D photonic crystals uncouple TE, TM"
    elseif N == 2
        @assert P <: CoupledPolarization "3D photonic crystals couple TE, TM"
    else
        error("only 1D, 2D, 3D photonic crystals are possible")
    end
end