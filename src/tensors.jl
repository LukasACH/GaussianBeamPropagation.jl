"""
    principal_axes(xx, xy, yy)

Principal values and orientation `(x, y, φ)` of the symmetric tensor `[xx xy; xy yy]`, such that
`tensor_from_principal(x, y, φ) == (xx, xy, yy)`.

`x` belongs to the principal axis closest to the _x_ direction, and `φ ∈ [−π/4, π/4]` is the
angle between the two. All quantities given along _x_ and _y_ use this convention.
"""
principal_axes(params::NTuple{3,<:Real}) = principal_axes(params...)

function principal_axes(xx::Real, xy::Real, yy::Real)
    diff_xx_yy = ifelse(xx ≈ yy, zero(xx - yy), xx - yy)

    a = xx + yy
    b = flipsign(sqrt((xx - yy)^2 + 4 * xy^2), diff_xx_yy)

    x = (a + b) / 2
    y = (a - b) / 2

    φ = atan(flipsign(2 * xy, diff_xx_yy), abs(xx - yy)) / 2

    return x, y, φ
end

"""
    tensor_from_principal(x, y, α)

Components `(xx, xy, yy)` of the symmetric tensor with the principal values `x` and `y`, whose
principal axes are rotated by `α`. Inverse of [`principal_axes`](@ref).
"""
function tensor_from_principal(x::Real, y::Real, α::Real)
    s, c = sincos(α)
    s2 = s^2
    c2 = c^2
    cs = c * s

    xx_ = x * c2 + y * s2
    xy_ = (x - y) * cs
    yy_ = y * c2 + x * s2

    return xx_, xy_, yy_
end

"""
    rotate_symmetric(xx, xy, yy, α)

Rotate the symmetric tensor `[xx xy; xy yy]` by `α`.
"""
function rotate_symmetric(xx::Real, xy::Real, yy::Real, α::Real)
    s, c = sincos(α)
    s2 = s^2
    c2 = c^2
    cs = c * s

    xx_ = xx * c2 - 2xy * cs + yy * s2
    xy_ = xy * c2 + (xx - yy) * cs - xy * s2
    yy_ = yy * c2 + 2xy * cs + xx * s2

    return xx_, xy_, yy_
end

"""
    rotate_general(xx, xy, yx, yy, α)

Rotate the (not necessarily symmetric) tensor `[xx xy; yx yy]` by `α`.
"""
function rotate_general(xx::Real, xy::Real, yx::Real, yy::Real, α::Real)
    s, c = sincos(α)
    s2 = s^2
    c2 = c^2
    cs = c * s

    xx_ = xx * c2 - (xy + yx) * cs + yy * s2
    xy_ = xy * c2 + (xx - yy) * cs - yx * s2
    yx_ = yx * c2 + (xx - yy) * cs - xy * s2
    yy_ = yy * c2 + (xy + yx) * cs + xx * s2

    return xx_, xy_, yx_, yy_
end

@deprecate rotate(x, y, α) tensor_from_principal(x, y, α) false
@deprecate rotate_around(A, B, C, α) rotate_symmetric(A, C, B, α) false
@deprecate rotate_around(A, B, α) tensor_from_principal(A, B, α) false
@deprecate rotate_around(XX, XY, YX, YY, α) rotate_general(XX, XY, YX, YY, α) false
@deprecate rotate_around(α; xx, xy, yx, yy) rotate_general(xx, xy, yx, yy, α) false
