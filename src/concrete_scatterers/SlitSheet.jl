export SlitSheet

"""
A geometry for a diffraction grating with slits
"""
struct SlitSheet <: RCWASheet{1}
    gap_center::Float64
    gap_width::Float64
end