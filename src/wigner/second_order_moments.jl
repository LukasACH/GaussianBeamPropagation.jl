struct SymmetricMoments{T}
    xx::T
    xy::T
    yy::T
end

function SymmetricMoments(xx::T, xy::T, yy::T) where {T<:AbstractFloat}
    return SymmetricMoments{T}(xx, xy, yy)
end

function SymmetricMoments(xx::T, yy::T) where {T<:AbstractFloat}
    return SymmetricMoments{T}(xx, zero(T), yy)
end

calculate_matrix(m::SymmetricMoments{T}) where {T} = SA{T}[
    m.xx m.xy
    m.xy m.yy
]

struct AsymmetricMoments{T}
    xx::T
    xy::T
    yx::T
    yy::T
end

function AsymmetricMoments(xx::T, xy::T, yx::T, yy::T) where {T<:AbstractFloat}
    return AsymmetricMoments{T}(xx, xy, yx, yy)
end

function AsymmetricMoments(xx::T, xy::T, yy::T) where {T<:AbstractFloat}
    return AsymmetricMoments{T}(xx, xy, yx, yy)
end

function AsymmetricMoments(xx::T, yy::T; sum=zero(T), dif=zero(T)) where {T<:AbstractFloat}
    return AsymmetricMoments{T}(xx, (sum+dif)/2, (sum-dif)/2, yy)
end

calculate_matrix(m::AsymmetricMoments{T}) where {T} = SA{T}[
    m.xx m.xy
    m.yx m.yy
]

struct SecondOrderBeamMoments{T}
    spatial::SymmetricMoments{T}
    angluar::SymmetricMoments{T}
    mixed::AsymmetricMoments{T}
end

function SecondOrderMoments(rtm::SMatrix{4,4,T,16}) where {T}
    spatial = SymmetricMoments(rtm[1, 1], (rtm[1, 2] + rtm[2, 1]) / 2, rtm[2, 2])

    angular = SymmetricMoments(rtm[3, 3], (rtm[3, 4] + rtm[4, 3]) / 2, rtm[4, 4])

    mixed = AsymmetricMoments(
        (rtm[1, 3] + rtm[3, 1]) / 2,
        (rtm[1, 4] + rtm[4, 1]) / 2,
        (rtm[2, 3] + rtm[3, 2]) / 2,
        (rtm[2, 4] + rtm[4, 2]) / 2,
    )

    return SecondOrderBeamMoments{T}(spatial, angular, mixed)
end

function calculate_matrix(som::SecondOrderBeamMoments{T}) where {T}
    return SA{T}[
        som.spatial.xx som.spatial.xy som.mixed.xx som.mixed.xy
        som.spatial.xy som.spatial.yy som.mixed.yx som.mixed.yy
        som.mixed.xx som.mixed.yx som.angular.xx som.angular.xy
        som.mixed.xy som.mixed.yy som.angular.xy som.angular.yy
    ]
end

# struct SpatialBeamMoments{T}
#     moments::SymmetricMoments{T}
# end
# Base.getproperty(m::SpatialBeamMoments, f::Symbol) = getfield(m.moments, f)

# function SpatialBeamMoments(xx::T, xy::T, yy::T) where {T<:AbstractFloat}
#     return SpatialBeamMoments(SymmetricMoments{T}(xx, xy, yy))
# end

# function SpatialBeamMoments(xx::T, yy::T) where {T<:AbstractFloat}
#     return SpatialBeamMoments(SymmetricMoments{T}(xx, zero(T), yy))
# end

"""
    SecondOrderMoments(; 
        rxrx, rxry, ryry,
        θxθx, θxθy, θyθy,
        rxθx, rxθy, ryθx, ryθy,
    )
"""
struct SecondOrderMoments{T<:AbstractFloat}
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

function SecondOrderMoments(rtm::SMatrix{4,4,Float64,16})
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
        @error "Insufficient parameters to determine Wigner Second Order Moments"
        return missing
    end

    rxrx, rxry, ryry = rotate2(rxrx_, ryry_, φ)
    θxθx, θxθy, θyθy = rotate2(θxθx_, θyθy_, φ)
    rxθx, rxθy, ryθy = rotate2(rxθx_, ryθy_, φ)
    ryθx = rxθy

    return SecondOrderMoments(rxrx, rxry, ryry, θxθx, θxθy, θyθy, rxθx, rxθy, ryθx, ryθy)
end
