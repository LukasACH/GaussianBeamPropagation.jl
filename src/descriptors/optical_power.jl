"""
    OpticalPower(P)
    OpticalPower(Px, Py; φ)
    OpticalPower(P_xx, P_xy, P_yy)

The optical power _P_ of a lens or mirror describes the change in convergence or divergence of light, equal to the reciprocal focal length _f_ such that _P_ ≡ _f_⁻¹. It is given in units of inverse metre.
"""
struct OpticalPower{T<:Real}
    xx::T
    xy::T
    yy::T
end

Base.convert(::Type{OpticalPower{T}}, P::OpticalPower) where {T<:Real} = OpticalPower{T}(P)

function OpticalPower{T}(P::OpticalPower) where {T<:Real}
    return OpticalPower{T}(T(P.xx), T(P.xy), T(P.yy))
end
function OpticalPower(P_xx::Real, P_xy::Real, P_yy::Real)
    return OpticalPower(promote(P_xx, P_xy, P_yy)...)
end

function OpticalPower(Px::Real, Py::Real; φ::Real=0.0f0)
    return OpticalPower(rotate2(
        Px,
        Py,
        convert(promote_type(typeof(Px), typeof(Py)), φ),
    )...)
end

OpticalPower(P::Real) = OpticalPower(P, P)
OpticalPower() = OpticalPower(0.0f0)

Base.:-(P::OpticalPower) = OpticalPower(-P.xx, -P.xy, -P.yy)

function Base.:-(lhs::OpticalPower, rhs::OpticalPower)
    return OpticalPower(lhs.xx - rhs.xx, lhs.xy - rhs.xy, lhs.yy - rhs.yy)
end

function Base.:+(lhs::OpticalPower, rhs::OpticalPower)
    return OpticalPower(lhs.xx + rhs.xx, lhs.xy + rhs.xy, lhs.yy + rhs.yy)
end

Base.:*(a::Real, P::OpticalPower) =
    OpticalPower(a * P.xx, a * P.xy, a * P.yy)

Base.:*(P::OpticalPower, a::Real) =
    OpticalPower(a * P.xx, a * P.xy, a * P.yy)

Base.:/(P::OpticalPower, a::Real) =
    OpticalPower(P.xx / a, P.xy / a, P.yy/a)

function Base.show(io::IO, ::MIME"text/plain", P::OpticalPower{T}) where {T}
    x, y, φ = get_x_y_φ(P.xx, P.xy, P.yy)

    return print(
        io,
        Base.repr(P),
        " with Px = ",
        x,
        " (fx = ",
        1/x,
        "), Py = ",
        y,
        " (fy = ",
        1/y,
        "), φ = ",
        rad2deg(φ),
        "°",
    )
end
