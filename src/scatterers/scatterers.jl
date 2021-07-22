# Use this file to include (and export) all of your scatterers
# Each scatterer defines its own parametrization and functions to obtain
# scattering matrices, or generically conductivity matrices

export ScatteringStack
include("ScatteringStack.jl")

export UniformLayer
include("UniformLayer.jl")

export UniformInterface
include("UniformInterface.jl")

export GaussianDeltaLayer
include("GaussianDeltaLayer.jl")

export HCellDeltaLayer
include("HCellDeltaLayer.jl")

export HoleDeltaLayer
include("HoleDeltaLayer.jl")