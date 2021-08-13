export SlitSheet

"""
A geometry for a diffraction grating with slits
"""
struct SlitSheet{T} <: RCWASheet{T, 1}
    gap_center::T
    gap_width::T
end