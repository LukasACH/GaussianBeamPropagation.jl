"""
    ThickLens(front_surface::Curvature, thickness::Real, back_surface::Curvature, n::Real)
    ThickLens(front_surface::Curvature, thickness::Real, back_surface::Curvature, n::Real, incidence::IncidenceAngle)

# Matrix representing a thick lens

This is a matrix builder that constructs the matrix for a thick lens.

The arguments are the following:

  - `front_surface::Curvature`: curvature of the front lens surface given by [`Curvature`](@ref), where a positive curvature denotes that the center of curvature lies after the lens, i.e., the front surface is convex.
  - `thickness::Real`: thickness of the lens where the ray or beam passes through.
  - `back_surface::Curvature`: curvature of the back given by [`Curvature`](@ref), where a positive curvature denotes that the center of curvature lies after the lens, i.e., the back surface is concave.
  - `n::Real`: ratio of the refractive index of the lens material to the refractive index of the surrounding material.
  - `incidence::IncidenceAngle` (optional): if given, defines the direction of the incident ray or beam using [`IncidenceAngle`](@ref). If absent, it is assumed that the ray is parallel to the surface normal in the center (disregarding any surface curvature).

!!! note "Coordinate system"

    The coordinate system used by the matrix is right handed, where the ray or beam travels in the +z direction.
    Both the orientation angle of the surface and the orientation of the incidence plane are rotations in the x-y-plane relative to the x̂ direction.
"""
function ThickLens(
    front_surface::Curvature,
    thickness::Real,
    back_surface::Curvature,
    n::Real,
    incidence::IncidenceAngle=IncidenceAngle(),
)
    incidence_2 = IncidenceAngle(; ι=asin(sin(incidence.ι) / n), θ=incidence.θ)

    front_interface = OpticalInterface(front_surface, n, incidence)

    space = FreeSpace(thickness)

    back_interface = OpticalInterface(back_surface, 1/n, incidence_2)

    return back_interface * space * front_interface
end
