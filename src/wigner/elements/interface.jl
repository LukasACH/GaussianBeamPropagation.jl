"""
    OpticalInterface(surface::Curvature, n::Real)
    OpticalInterface(surface::Curvature, n::Real, incidence::IncidenceAngle)

# Matrix representing an optical interface

This is a matrix builder that constructs the matrix for a refractive interface between two materials.

The arguments are the following:

  - `surface::Curvature`: curvature of the interface given by [`Curvature`](@ref), where a positive curvature denotes that the center of curvature lies after the interface, i.e., the surface is convex towards the incoming beam.
  - `n::Real`: ratio of the refractive index of the material after the interface to the refractive index of the material before the interface.
  - `incidence::IncidenceAngle` (optional): if given, defines the direction of the incident ray or beam using [`IncidenceAngle`](@ref). If absent, it is assumed that the ray is parallel to the surface normal in the center (disregarding any surface curvature).

!!! note "Coordinate system"

    The coordinate system used by the matrix is right handed, where the ray or beam travels in the +z direction.
    Both the orientation angle of the surface and the orientation of the incidence plane are rotations in the x-y-plane relative to the x̂ direction.
"""
struct OpticalInterface{T<:AbstractFloat} <: AbstractElement
    n::T
    surface::Curvature{T}
    incidence::IncidenceAngle{T}

    function OpticalInterface(
        surface::Curvature{A},
        n::B,
        incidence::IncidenceAngle{C}=IncidenceAngle(),
    ) where {A<:Real,B<:Real,C<:Real}
        return new{float(promote_type(A, B, C))}(n, surface, incidence)
    end
end

function get_transfer_matrix(e::OpticalInterface{T})::SMatrix{4,4,T,16} where {T}
    κxx, κxy, κyy = rotate_symmetric(e.surface.xx, e.surface.xy, e.surface.yy, -e.incidence.θ)

    V = sqrt(e.n^2 - sin(e.incidence.ι)^2)

    Axx, Axy, Ayy = tensor_from_principal(V / cos(e.incidence.ι) / e.n, 1, e.incidence.θ)
    Dxx, Dxy, Dyy = tensor_from_principal(cos(e.incidence.ι) / V, 1 / e.n, e.incidence.θ)
    # A positive curvature (convex towards the incoming beam) focuses for n > 1
    Mxx, Mxy, Myx, Myy = rotate_general(
        κxx * (cos(e.incidence.ι) - V) / V / cos(e.incidence.ι),
        κxy * (cos(e.incidence.ι) - V) / V,
        κxy * (cos(e.incidence.ι) - V) / cos(e.incidence.ι) / e.n,
        κyy * (cos(e.incidence.ι) - V) / e.n,
        e.incidence.θ,
    )

    Axx, Ayy, Axy, Dxx, Dyy, Dxy, Mxx, Mxy, Myx, Myy =
        promote(Axx, Ayy, Axy, Dxx, Dyy, Dxy, Mxx, Mxy, Myx, Myy)

    return SA{T}[
        Axx Axy 0.0 0.0
        Axy Ayy 0.0 0.0
        Mxx Mxy Dxx Dxy
        Myx Myy Dxy Dyy
    ]
end
