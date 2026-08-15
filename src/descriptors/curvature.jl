"""
    Curvature(κ::Real)
    Curvature(κx::Real, κy::Real; φ::Real)
    Curvature(κ_xx::Real, κ_xy::Real, κ_yy::Real)

# Spherical curvature

Defines the curvature of a spherical surface from a single curvature _κ_, which is defined as the reciprocal radius of curvature _R_, such that _κ_ = _R_⁻¹.
It has a single argument:

  - `κ::Real`: defines the spherical curvature of the surface.

# Toroidal curvature

To define a toroidal curvature (including a cylindrical surface, two radii of curvature in two principal axes are required:

  - `κx::Real`: defines the surface curvature in the x' direction.
  - `κy::Real`: defines the surface curvature in the y' direction.
  - `φ::Real` (keyword, optional): defines the rotation of the x' and y' directions relative to the global x and y directions.

# Projected curvature

Instead of providing the curvature in the principal directions with a rotation angle, you can provide the components of the vectors already rotated:

  - `κ_xx::Real`: projected curvature in the x direction.
  - `κ_xy::Real`: (κx - κy) * cos(φ) * sin(φ)
  - `κ_yy::Real`: projected curvature in the y direction

!!! note "Coordinate system"

    The coordinate system used by the matrix is right handed, where the ray or beam travels in the +z direction.
    Both the orientation angle of the surface and the orientation of the incidence plane are rotations in the x-y-plane relative to the x̂ direction.

The curvature _κ_ of a surface is the inverse of its radius of curvature _R_, such that _R_ ≡ _κ_⁻¹. It is given in units of inverse metre.
"""
struct Curvature{T<:Real}
    xx::T
    xy::T
    yy::T
end

Base.convert(::Type{Curvature{T}}, κ::Curvature) where {T<:Real} =
    Curvature{T}(κ)

Curvature{T}(κ::Curvature) where {T<:Real} =
    Curvature{T}(T(κ.xx), T(κ.xy), T(κ.yy))

Curvature(κ_xx::Real, κ_xy::Real, κ_yy::Real) =
    Curvature(promote(κ_xx, κ_xy, κ_yy)...)

function Curvature(κx::Real, κy::Real; φ::Real=0.0f0)
    return Curvature(rotate2(κx, κy, convert(promote_type(typeof(κx), typeof(κy)), φ))...)
end

Curvature(κ::Real) = Curvature(κ, κ)
Curvature() = Curvature(0.0f0)

Base.:-(κ::Curvature) = Curvature(-κ.xx, -κ.xy, -κ.yy)

function Base.:-(lhs::Curvature, rhs::Curvature)
    return Curvature(lhs.xx - rhs.xx, lhs.xy - rhs.xy, lhs.yy - rhs.yy)
end

function Base.:+(lhs::Curvature, rhs::Curvature)
    return Curvature(lhs.xx + rhs.xx, lhs.xy + rhs.xy, lhs.yy + rhs.yy)
end

Base.:*(a::Real, κ::Curvature) =
    Curvature(a * κ.xx, a * κ.xy, a * κ.yy)

Base.:*(κ::Curvature, a::Real) =
    Curvature(a * κ.xx, a * κ.xy, a * κ.yy)

Base.:/(κ::Curvature, a::Real) =
    Curvature(κ.xx / a, κ.xy / a, κ.yy/a)

function Base.show(io::IO, ::MIME"text/plain", κ::Curvature{T}) where {T}
    x, y, φ = get_x_y_φ(κ.xx, κ.xy, κ.yy)

    return print(
        io,
        Base.repr(κ),
        " with κx = ",
        x,
        " (Rx = ",
        1/x,
        "), κy = ",
        y,
        " (Ry = ",
        1/y,
        "), φ = ",
        rad2deg(φ),
        "°",
    )
end
