"""
    SecondOrderMoments(; 
        rxrx, rxry, ryry,
        θxθx, θxθy, θyθy,
        rxθx, rxθy, ryθx, ryθy,
    )
"""
struct SecondOrderMoments{T}
    rxrx::T
    rxry::T
    ryry::T
    θxθx::T
    θxθy::T
    θyθy::T
    rxθx::T
    rxθy::T
    ryθx::T
    ryθy::T
end

function SecondOrderMoments(; rx, ry, φ_r, θx, θy, φ_θ, rθx, rθy, φ_rθ, twist)
    rxrx, rxry, ryry = rotate2(rx, ry, φ_r)
    θxθx, θxθy, θyθy = rotate2(θx, θy, φ_θ)
    rxθx, rθavg, ryθy = rotate2(rθx, rθy, φ_rθ)
    rxθy = rθavg + twist / 2
    ryθx = rθavg - twist / 2
    return SecondOrderMoments(rxrx, rxry, ryry, θxθx, θxθy, θyθy, rxθx, rxθy, ryθx, ryθy)
end

function SecondOrderMoments(rtm::SMatrix{4,4,F,16}) where {F}
    rxrx = rtm[1, 1]
    rxry = (rtm[1, 2] + rtm[2, 1]) / 2
    ryry = rtm[2, 2]

    θxθx = rtm[3, 3]
    θxθy = (rtm[3, 4] + rtm[4, 3]) / 2
    θyθy = rtm[4, 4]

    rxθx = (rtm[1, 3] + rtm[3, 1]) / 2
    rxθy = (rtm[1, 4] + rtm[4, 1]) / 2
    ryθx = (rtm[2, 3] + rtm[3, 2]) / 2
    ryθy = (rtm[2, 4] + rtm[4, 2]) / 2

    return SecondOrderMoments(rxrx, rxry, ryry, θxθx, θxθy, θyθy, rxθx, rxθy, ryθx, ryθy)
end

function WSOM_from_Q(λ::Real, q::Union{Complex,Missing})
    ismissing(q) && return missing

    ε = λ / 4π

    rr = ε * abs2(q) / imag(q)
    rθ = ε * real(q) / imag(q)
    θθ = ε / imag(q)

    return rr, rθ, θθ
end

function WSOM_from_z_zR(λ::Real, z::Union{Real,Missing}, zR::Union{Real,Missing})
    (ismissing(z) || ismissing(zR)) && return missing

    ε = λ / 4π

    rr = ε * (z^2 / zR + zR)
    rθ = ε * z / zR
    θθ = ε / zR

    return rr, rθ, θθ
end

function WSOM_from_w_r(λ::Real, w::Union{Real,Missing}, r::Union{Real,Missing})
    (ismissing(w) || ismissing(r)) && return missing

    ε = λ / 4π

    w2 = (w / 2)^2

    rr = w2
    rθ = w2 / r
    θθ = w2 / r^2 + ε^2 / w2

    return rr, rθ, θθ
end

function SecondOrderMoments(
    λ::Real;
    zx::Union{Real,Missing}=missing,
    zRx::Union{Real,Missing}=missing,
    zy::Union{Real,Missing}=missing,
    zRy::Union{Real,Missing}=missing,
    wx::Union{Real,Missing}=missing,
    wy::Union{Real,Missing}=missing,
    rx::Union{Real,Missing}=missing,
    ry::Union{Real,Missing}=missing,
    qx::Union{Complex,Missing}=missing,
    qy::Union{Complex,Missing}=missing,
    φ::Real=0.0,
)
    rxrx_, rxθx_, θxθx_ = ifelse(
        !ismissing(zx) && !ismissing(zRx),
        WSOM_from_z_zR(λ, zx, zRx),
        ifelse(
            !ismissing(wx) && !ismissing(rx),
            WSOM_from_w_r(λ, wx, rx),
            ifelse(!ismissing(qx), WSOM_from_Q(λ, qx), (missing, missing, missing)),
        ),
    )
    ryry_, ryθy_, θyθy_ = ifelse(
        !ismissing(zy) && !ismissing(zRy),
        WSOM_from_z_zR(λ, zy, zRy),
        ifelse(
            !ismissing(wy) && !ismissing(ry),
            WSOM_from_w_r(λ, wy, ry),
            ifelse(!ismissing(qy), WSOM_from_Q(λ, qy), (missing, missing, missing)),
        ),
    )

    if ismissing(rxrx_) || ismissing(ryry_) || ismissing(θxθx_) || ismissing(θyθy_)
        throw(ArgumentError("insufficient parameters: each axis needs (z, zR), (w, r) or q"))
    end

    rxrx, rxry, ryry = rotate2(rxrx_, ryry_, φ)
    θxθx, θxθy, θyθy = rotate2(θxθx_, θyθy_, φ)
    rxθx, rxθy, ryθy = rotate2(rxθx_, ryθy_, φ)
    ryθx = rxθy

    return SecondOrderMoments(rxrx, rxry, ryry, θxθx, θxθy, θyθy, rxθx, rxθy, ryθx, ryθy)
end
