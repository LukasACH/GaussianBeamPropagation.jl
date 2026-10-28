"""
    OpticalPower(P)
    OpticalPower(Px, Py; φ)
    OpticalPower(P_xx, P_xy, P_yy)

The optical power _P_ of a lens or mirror describes the change in convergence or divergence of light, equal to the reciprocal focal length _f_ such that _P_ ≡ _f_⁻¹. It is given in units of inverse metre.
"""
struct OpticalPower{T<:Real} <: SymmetricDescriptor{T}
    xx::T
    xy::T
    yy::T
end

@symmetric_descriptor OpticalPower "P" "f"

"""
    lensmaker(front::Curvature, back::Curvature, n::Real)

Optical power of a thin lens with the surface curvatures `front` and `back` and the relative
refractive index `n` (lensmaker's equation).
"""
function lensmaker(front::Curvature, back::Curvature, n::Real)
    P = (n - 1) * (front - back)
    return OpticalPower(P.xx, P.xy, P.yy)
end
