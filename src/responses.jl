export ElectricResponseStyle, MagneticResponseStyle, Impedance, Admittance

"""
Based on `IndexStyle` in Base, provides an interface for user-defined types to
choose between a formulation of the same problem in terms of their preferred
material response parameters.
"""
abstract type ResponseStyle end

"""
    Impedance()

`ResponseStyle` to describe materials whose response is given by an impedance.
"""
struct Impedance <: ResponseStyle end

"""
    Admittance()

`ResponseStyle` to describe materials whose response is given by an admittance.
"""
struct Admittance <: ResponseStyle end

"""
    ElectricResponseStyle(sheet)
    ElectricResponseStyle(typeof(sheet))
"""
ElectricResponseStyle(sheet::RCWASheet) = ElectricResponseStyle(typeof(sheet))
ElectricResponseStyle(::Type{<:RCWASheet}) = Impedance()

"""
    MagneticResponseStyle(sheet)
    MagneticResponseStyle(typeof(sheet))
"""
MagneticResponseStyle(sheet::RCWASheet) = MagneticResponseStyle(typeof(sheet))
MagneticResponseStyle(::Type{<:RCWASheet}) = Admittance()