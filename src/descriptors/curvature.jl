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
struct Curvature{T<:Real} <: SymmetricDescriptor{T}
    xx::T
    xy::T
    yy::T
end

@symmetric_descriptor Curvature "κ" "R"
