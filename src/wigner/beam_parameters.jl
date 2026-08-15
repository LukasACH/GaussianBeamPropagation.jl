"""
    radius(::Beam)
    radius(::SecondOrderMoments)
    radius(xx, xy, yy)

Calculates the 1/e² radius from a `Beam`, `SecondOrderMoments` or measured second order moments (`xx`, `xy`, and `yy`).
The return value is a `NamedTuple` with entries `(:x, :y, :φ)`, where both _x_ and _y_ are given in metres, and the angle _φ_ is given in radians.
It follows the definition of ISO 11146 for `x`, `y`, and the angle `φ`.

It uses a right-handed coordinate system.
"""
radius(beam::Beam{F,SecondOrderMoments{F}}) where {F} = radius(beam.inner)
radius(beam::Beam{F}) where {F} = radius(Beam{F,SecondOrderMoments{F}}(beam).inner)

"""
    divergence(::Beam)
    divergence(::SecondOrderMoments)

Calculates the divergence angle from a `Beam` or `SecondOrderMoments`.
The return value is a `NamedTuple` with entries `(:x, :y, :φ)`, where both _x_ and _y_ are given in radians, and the angle _φ_ is given in radians as well.

It uses a right-handed coordinate system.
"""
divergence(beam::Beam{F,SecondOrderMoments{F}}) where {F} = divergence(beam.inner)
divergence(beam::Beam{F}) where {F} = divergence(Beam{F,SecondOrderMoments{F}}(beam).inner)

"""
    phase_curvature(::Beam)
    phase_curvature(::SecondOrderMoments)

Calculates the phase-front curvature from a `Beam` or `SecondOrderMoments`.
The return value is a `NamedTuple` with entries `(:x, :y, :φ)`, where both _x_ and _y_ are given in metres/radians, and the angle _φ_ is given in radians.

It uses a right-handed coordinate system.
"""
phase_curvature(beam::Beam{F,SecondOrderMoments{F}}) where {F} = phase_curvature(beam.inner)
function phase_curvature(beam::Beam{F}) where {F}
    return phase_curvature(Beam{F,SecondOrderMoments{F}}(beam).inner)
end

"""
    twist(::Beam)
    twist(::SecondOrderMoments)

Calculates the twist from a `Beam` or `SecondOrderMoments`.
The return value is a single `Float64` in units of metres²/radians².
"""
twist(beam::Beam{F,SecondOrderMoments{F}}) where {F} = twist(beam.inner)
twist(beam::Beam{F}) where {F} = twist(Beam{F,SecondOrderMoments{F}}(beam).inner)

function radius(som::SecondOrderMoments)::NamedTuple{(:x, :y, :φ),NTuple{3,Float64}}
    return radius(som.rxrx, som.rxry, som.ryry)
end

function radius(xx, xy, yy)::NamedTuple{(:x, :y, :φ),NTuple{3,Float64}}
    diff_xx_yy = ifelse(xx ≈ yy, 0.0, xx-yy)

    a = xx + yy
    b = flipsign(sqrt((xx - yy)^2 + 4 * xy^2), diff_xx_yy)

    x = sqrt(2 * (a + b))
    y = sqrt(2 * (a - b))

    φ = atan(flipsign(2 * xy, diff_xx_yy), abs(xx - yy)) / 2

    return (; x, y, φ)
end

function divergence(som::SecondOrderMoments)
    a = (som.θxθx + som.θyθy)
    b = sqrt((som.θxθx - som.θyθy)^2 + 4 * som.θxθy^2)
    τ = 1 # sign(som.θxθx - som.θyθy)

    x = sqrt(2 * (a + τ * b))
    y = sqrt(2 * (a - τ * b))
    φ = atan(2 * som.θxθy, (som.θxθx - som.θyθy)) / 2

    φ = mod2pi(2φ) / 2
    return (; x, y, φ)
end

function phase_curvature(som::SecondOrderMoments)
    # rtm = generate_wigner_matrix(som)
    # W = SMatrix{2,2}(@view rtm[1:2, 1:2]) # Extract position-position submatrix
    # U = SMatrix{2,2}(@view rtm[3:4, 3:4]) # Extract angle-angle submatrix
    # M = SMatrix{2,2}(@view rtm[1:2, 3:4]) # Extract position-angle submatrix
    W = SA_F64[som.rxrx som.rxry; som.rxry som.ryry] # .* u"m^2"
    # U = SA_F64[som.θxθx som.θxθy; som.θxθy som.θyθy]
    M = SA_F64[som.rxθx som.rxθy; som.ryθx som.ryθy] # .* u"m*rad"

    A = SA[0 1; -1 0] / W
    B = M / W
    C = (A - A') \ (B - B') * A - B

    Cxx = C[1, 1]
    Cxy = (C[1, 2] + C[2, 1]) / 2
    Cyy = C[2, 2]

    a = (Cxx + Cyy)
    b = sqrt((Cxx - Cyy)^2 + 4 * Cxy^2)
    μ = 1 # sign(Cxx - Cyy)

    x = -2 / (a + μ * b)
    y = -2 / (a - μ * b)
    φ = atan(2 * Cxy, (Cxx - Cyy)) / 2

    φ = mod2pi(2φ) / 2
    return (; x, y, φ)
end

twist(som::SecondOrderMoments) = som.rxθy - som.ryθx

#=
function phase_curvature(som::SecondOrderMoments, λ::Real)
    @info "Calculating phase from Wigner Second Order Moments"
    a = (som.rxθx + som.ryθy)
    b = sqrt((som.rxθx - som.ryθy)^2 + 4 * som.rxθy * som.ryθx)
    μ = 1 # sign(som.rxθx - som.ryθy)

    x = sqrt(2 * (a + υ * b))
    y = sqrt(2 * (a - υ * b))
    φ = atan(2 * som.rxθy, (som.rxθx - som.ryθy)) / 2

    return (; x, y, φ)
end
=#