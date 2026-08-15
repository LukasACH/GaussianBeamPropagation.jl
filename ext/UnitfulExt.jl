module UnitfulExt

using GaussianBeamPropagation
using Unitful

function GaussianBeamPropagation.rotate2(
    x::Quantity{<:Real,D}, y::Quantity{<:Real,D}, α::Real
) where {D}
    s, c = sincos(α)
    s2 = s^2
    c2 = c^2
    cs = c * s

    xx_ = x * c2 + y * s2
    xy_ = (x - y) * cs
    yy_ = y * c2 + x * s2

    return xx_, xy_, yy_
end

function GaussianBeamPropagation.SecondOrderMoments(
    rxrx::Quantity{Float64,Unitful.𝐋^2},
    rxry::Quantity{Float64,Unitful.𝐋^2},
    ryry::Quantity{Float64,Unitful.𝐋^2},
    θxθx::Union{Float64,Quantity{Float64,NoDims}},
    θxθy::Union{Float64,Quantity{Float64,NoDims}},
    θyθy::Union{Float64,Quantity{Float64,NoDims}},
    rxθx::Quantity{Float64,Unitful.𝐋},
    rxθy::Quantity{Float64,Unitful.𝐋},
    ryθx::Quantity{Float64,Unitful.𝐋},
    ryθy::Quantity{Float64,Unitful.𝐋},
)
    return SecondOrderMoments(
        ustrip(u"m^2", rxrx),
        ustrip(u"m^2", rxry),
        ustrip(u"m^2", ryry),
        ustrip(NoUnits, θxθx),
        ustrip(NoUnits, θxθy),
        ustrip(NoUnits, θyθy),
        ustrip(u"m", rxθx),
        ustrip(u"m", rxθy),
        ustrip(u"m", ryθx),
        ustrip(u"m", ryθy),
    )
end

function GaussianBeamPropagation.WSOM_from_Q(
    λ::Quantity{<:Real}, q::Union{Quantity{<:Complex},Missing}
)
    ismissing(q) && return missing

    ε = λ / 4π

    rr = ε * abs2(q) / imag(q)
    rθ = ε * real(q) / imag(q)
    θθ = ε / imag(q)

    return rr, rθ, θθ
end

function GaussianBeamPropagation.WSOM_from_z_zR(
    λ::Quantity{<:Real},
    z::Union{Quantity{<:Real},Missing},
    zR::Union{Quantity{<:Real},Missing},
)
    (ismissing(z) || ismissing(zR)) && return missing

    ε = λ / 4π

    rr = ε * (z^2 / zR + zR)
    rθ = ε * z / zR
    θθ = ε / zR

    return rr, rθ, θθ
end

function GaussianBeamPropagation.WSOM_from_w_r(
    λ::Quantity{<:Real},
    w::Union{Quantity{<:Real},Missing},
    r::Union{Quantity{<:Real},Missing},
)
    (ismissing(w) || ismissing(r)) && return missing

    ε = λ / 4π

    w2 = (w / 2)^2

    rr = w2
    rθ = w2 / r
    θθ = w2 / r^2 + ε^2 / w2

    return rr, rθ, θθ
end

function GaussianBeamPropagation.Beam(λ::Quantity{<:Real,Unitful.𝐋}; kwargs...)
    return Beam(ustrip(u"m", λ), SecondOrderMoments(λ; kwargs...))
end

function GaussianBeamPropagation.SecondOrderMoments(
    λ::Quantity;
    zx::Union{Quantity,Missing}=missing,
    zRx::Union{Quantity,Missing}=missing,
    zy::Union{Quantity,Missing}=missing,
    zRy::Union{Quantity,Missing}=missing,
    wx::Union{Quantity,Missing}=missing,
    wy::Union{Quantity,Missing}=missing,
    rx::Union{Quantity,Missing}=missing,
    ry::Union{Quantity,Missing}=missing,
    qx::Union{Quantity{<:Complex},Missing}=missing,
    qy::Union{Quantity{<:Complex},Missing}=missing,
    φ::Real=0.0,
)
    rxrx_, rxθx_, θxθx_ = ifelse(
        !ismissing(zx) && !ismissing(zRx),
        GaussianBeamPropagation.WSOM_from_z_zR(λ, zx, zRx),
        ifelse(
            !ismissing(wx) && !ismissing(rx),
            GaussianBeamPropagation.WSOM_from_w_r(λ, wx, rx),
            ifelse(
                !ismissing(qx),
                GaussianBeamPropagation.WSOM_from_Q(λ, qx),
                (missing, missing, missing),
            ),
        ),
    )
    ryry_, ryθy_, θyθy_ = ifelse(
        !ismissing(zy) && !ismissing(zRy),
        GaussianBeamPropagation.WSOM_from_z_zR(λ, zy, zRy),
        ifelse(
            !ismissing(wy) && !ismissing(ry),
            GaussianBeamPropagation.WSOM_from_w_r(λ, wy, ry),
            ifelse(
                !ismissing(qy),
                GaussianBeamPropagation.WSOM_from_Q(λ, qy),
                (missing, missing, missing),
            ),
        ),
    )

    if ismissing(rxrx_) || ismissing(ryry_) || ismissing(θxθx_) || ismissing(θyθy_)
        @error "Insufficient parameters to determine Wigner Second Order Moments"
        return missing
    end

    rxrx, rxry, ryry = GaussianBeamPropagation.rotate2(rxrx_, ryry_, φ)
    θxθx, θxθy, θyθy = GaussianBeamPropagation.rotate2(θxθx_, θyθy_, φ)
    rxθx, rxθy, ryθy = GaussianBeamPropagation.rotate2(rxθx_, ryθy_, φ)
    ryθx = rxθy

    @debug "" rxrx rxry ryry θxθx θxθy θyθy rxθx rxθy ryθy

    return SecondOrderMoments(rxrx, rxry, ryry, θxθx, θxθy, θyθy, rxθx, rxθy, ryθx, ryθy)
end

function GaussianBeamPropagation.ThinLens(
    Dx::Quantity{<:Real,Unitful.𝐋^-1}, Dy::Quantity{<:Real,Unitful.𝐋^-1}; kwargs...
)
    return ThinLens(ustrip(u"m^-1", Dx), ustrip(u"m^-1", Dy); kwargs...)
end

function GaussianBeamPropagation.ThinLens(
    fx::Quantity{<:Real,Unitful.𝐋}, fy::Quantity{<:Real,Unitful.𝐋}; kwargs...
)
    return ThinLens(1 / ustrip(u"m", fx), 1 / ustrip(u"m", fy); kwargs...)
end

GaussianBeamPropagation.FreeSpace(d::Quantity) = FreeSpace(ustrip(u"m", d))

function GaussianBeamPropagation.ThinLens(D::Quantity{<:Real}; ι::Real=0.0, θ::Real=0.0)
    return ThinLens(D, D; ι, θ)
end

end
