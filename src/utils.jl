mod_pi(x) = mod2pi(2x) / 2
mod_pi2(x) = mod_pi(2x) / 2
mod_mpi2_pi2(x) = mod_pi(x + π / 2) - π / 2
mod_mpi4_pi4(x) = mod_mpi2_pi2(2x) / 2

get_x_y_φ(params::NTuple{3,<:Real}) = get_x_y_φ(params...)

function get_x_y_φ(xx::Real, xy::Real, yy::Real)
    diff_xx_yy = ifelse(xx ≈ yy, 0.0, xx-yy)

    a = xx + yy
    b = flipsign(sqrt((xx - yy)^2 + 4 * xy^2), diff_xx_yy)

    x = (a + b) / 2
    y = (a - b) / 2

    φ = atan(flipsign(2 * xy, diff_xx_yy), abs(xx - yy)) / 2

    return x, y, φ
end

function rotate2(x::Real, y::Real, α::Real)
    s, c = sincos(α)
    s2 = s^2
    c2 = c^2
    cs = c * s

    xx_ = x * c2 + y * s2
    xy_ = (x - y) * cs
    yy_ = y * c2 + x * s2

    return xx_, xy_, yy_
end

function rotate3(xx::Real, xy::Real, yy::Real, α::Real)
    s, c = sincos(α)
    s2 = s^2
    c2 = c^2
    cs = c * s

    xx_ = xx * c2 - 2xy * cs + yy * s2
    xy_ = xy * c2 + (xx - yy) * cs - xy * s2
    yy_ = yy * c2 + 2xy * cs + xx * s2

    return xx_, xy_, yy_
end

function rotate4(xx::Real, xy::Real, yx::Real, yy::Real, α::Real)
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

@deprecate rotate(x, y, α) rotate2(x, y, α)
@deprecate rotate_around(A, B, C, α) rotate3(A, C, B, α)
@deprecate rotate_around(A, B, α) rotate2(A, B, α)
@deprecate rotate_around(XX, XY, YX, YY, α) rotate4(XX, XY, YX, YY, α)
@deprecate rotate_around(α; xx, xy, yx, yy) rotate4(xx, xy, yx, yy, α)
