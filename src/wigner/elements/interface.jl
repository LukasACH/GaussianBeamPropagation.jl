"""
    OpticalInterface(surface::Curvature, n::Real)
    OpticalInterface(surface::Curvature, n::Real, incidence::IncidenceAngle)

# Matrix representing an optical interface

This is a matrix builder that constructs the matrix for a refractive interface between two materials.

The arguments are the following:

  - `surface::Curvature`: curvature of the interface given by [`Curvature`](@ref), where a positive curvature denotes that the center of curvature lies after the interface.
  - `n::Real`: ratio of the refractive index of the material after the interface to the refractive index of the material before the interface.
  - `incidence::IncidenceAngle` (optional): if given, defines the direction of the incident ray or beam using [`IncidenceAngle`](@ref). If absent, it is assumed that the ray is parallel to the surface normal in the center (disregarding any surface curvature).

!!! note "Coordinate system"

    The coordinate system used by the matrix is right handed, where the ray or beam travels in the +z direction.
    Both the orientation angle of the surface and the orientation of the incidence plane are rotations in the x-y-plane relative to the x̂ direction.
"""
struct OpticalInterface{T} <: AbstractElement
    n::T
    surface::Curvature{T}
    incidence::IncidenceAngle{T}

    rtm::SMatrix{4,4,T,16}

    function OpticalInterface(surface::Curvature, n::Real, incidence=IncidenceAngle())
        κxx, κxy, κyy = rotate3(surface.xx, surface.xy, surface.yy, -incidence.θ)

        V = sqrt(n^2 - sin(incidence.ι)^2)

        Axx, Axy, Ayy = rotate2(V / cos(incidence.ι) / n, 1.0f0, incidence.θ)
        Dxx, Dxy, Dyy = rotate2(cos(incidence.ι) / V, 1 / n, incidence.θ)
        Mxx, Mxy, Myx, Myy = rotate4(
            κxx * (V - cos(incidence.ι)) / V / cos(incidence.ι),
            κxy * (V - cos(incidence.ι)) / V,
            κxy * (V - cos(incidence.ι)) / cos(incidence.ι) / n,
            κyy * (V - cos(incidence.ι)) / n,
            incidence.θ,
        )

        Axx, Ayy, Axy, Dxx, Dyy, Dxy, Mxx, Mxy, Myx, Myy =
            promote(Axx, Ayy, Axy, Dxx, Dyy, Dxy, Mxx, Mxy, Myx, Myy)

        R = typeof(Axx)

        rtm = SA{R}[
            Axx Axy 0.0 0.0
            Axy Ayy 0.0 0.0
            Mxx Mxy Dxx Dxy
            Myx Myy Dxy Dyy
        ]

        return new{R}(n, surface, incidence, rtm)
    end
end
