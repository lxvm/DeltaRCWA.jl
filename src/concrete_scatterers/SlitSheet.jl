export SlitSheet

"""
A geometry for a diffraction grating with slits
"""
struct SlitSheet{T <: Real} <: RCWASheet{1}
    gap_center::T
    gap_width::T
end