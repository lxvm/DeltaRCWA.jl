export TE, TM, Coupled

abstract type AbstractPolarization end
abstract type UncoupledPolarization <: AbstractPolarization end
abstract type CoupledPolarization <: AbstractPolarization end

struct TE <: UncoupledPolarization end
struct TM <: UncoupledPolarization end
struct Coupled <: CoupledPolarization end