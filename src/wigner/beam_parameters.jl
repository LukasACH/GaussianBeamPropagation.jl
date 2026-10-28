"""
    radius(::Beam [, φ::Real])
    radius(::SecondOrderMoments [, φ::Real])
    radius(xx, xy, yy [, φ::Real])

Calculates the 1/e² radius from a `Beam`, `SecondOrderMoments` or measured second order moments (`xx`, `xy`, and `yy`).
The return value is a `SVector` with three elements, the radius in both _x_ and _y_, and the angle _φ_.
The values _x_ and _y_ are the given in metres, and the angle _φ_ is given in radians.
It follows the definition of ISO 11146 for `x`, `y`, and the angle `φ`: _x_ is the principal axis closest to the
_x_ direction, and _φ_ ∈ [−π/4, π/4] is its angle to the _x_ direction.

If the optional argument `φ` is provided, the radius is calculated along the direction of the angle `φ` (in radians).
In that case, the return value is a `SVector` with only two elements, the radius in the _x_ and _y_ directions, both given in metres.

It uses a right-handed coordinate system.
"""
radius(beam::Beam) = radius(beam.inner)
radius(beam::Beam, φ::Real) = radius(beam.inner, φ)

"""
    divergence(::Beam)
    divergence(::SecondOrderMoments)

Calculates the divergence angle from a `Beam` or `SecondOrderMoments`.
The return value is a `SVector` with three elements, the divergence in both _x_ and _y_ (in radians), and the angle _φ_ (in radians).
The axes and the angle follow the same convention as [`radius`](@ref): _φ_ ∈ [−π/4, π/4].

It uses a right-handed coordinate system.
"""
divergence(beam::Beam) = divergence(beam.inner)

"""
    phase_curvature(::Beam)
    phase_curvature(::SecondOrderMoments)

Calculates the phase-front curvature from a `Beam` or `SecondOrderMoments`.
The return value is a `SVector` with three elements, the phase curvature in both _x_ and _y_ (in metres/radians), and the angle _φ_ (in radians).
The axes and the angle follow the same convention as [`radius`](@ref): _φ_ ∈ [−π/4, π/4].

It uses a right-handed coordinate system.
"""
phase_curvature(beam::Beam) = phase_curvature(beam.inner)

"""
    twist(::Beam)
    twist(::SecondOrderMoments)

Calculates the twist from a `Beam` or `SecondOrderMoments`.
The return value is a single `Float64` in units of metres²/radians².
"""
twist(beam::Beam) = twist(beam.inner)

function radius(som::SecondOrderMoments)
    return radius(som.rxrx, som.rxry, som.ryry)
end

function radius(som::SecondOrderMoments, φ::Real)
    return radius(som.rxrx, som.rxry, som.ryry, φ)
end

function radius(xx::Real, xy::Real, yy::Real, φ0::Real)
    x, y, φ = radius(xx, xy, yy)

    s, c = sincos(φ0 - φ)

    xn = sqrt((x * c)^2 + (y * s)^2)
    yn = sqrt((x * s)^2 + (y * c)^2)

    return SA[xn, yn]
end

function radius(xx::Real, xy::Real, yy::Real)
    x, y, φ = principal_axes(xx, xy, yy)
    return SA[2sqrt(x), 2sqrt(y), φ]
end

function divergence(som::SecondOrderMoments)
    x, y, φ = principal_axes(som.θxθx, som.θxθy, som.θyθy)
    return SA[2sqrt(x), 2sqrt(y), φ]
end

function phase_curvature(som::SecondOrderMoments)
    W = SA[som.rxrx som.rxry; som.rxry som.ryry]
    M = SA[som.rxθx som.rxθy; som.ryθx som.ryθy]

    A = SA[0 1; -1 0] / W
    B = M / W
    C = (A - A') \ (B - B') * A - B

    x, y, φ = principal_axes(C[1, 1], (C[1, 2] + C[2, 1]) / 2, C[2, 2])
    return SA[-1 / x, -1 / y, φ]
end

twist(som::SecondOrderMoments) = som.rxθy - som.ryθx
