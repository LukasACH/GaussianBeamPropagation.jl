"""
    ThinLens(power::OpticalPower)
    ThinLens(power::OpticalPower, n::Real, incidence::IncidenceAngle)
    ThinLens(front::Curvature, back::Curvature, n::Real)
    ThinLens(front::Curvature, back::Curvature, n::Real, incidence::IncidenceAngle)

The zero-thickness approximation of the thick lens. The oblique astigmatism depends on the refractive indices of the
lens material and the surrounding medium, which are given as a single ratio of internal / external, and the incidence angles are given by `IncidenceAngle`. The focussing of the lens can either be given as `RadiusOfCurvature`,
`FocalLength`, or `FocalPower`.
"""
struct ThinLens{T<:AbstractFloat} <: AbstractElement
    n::T
    power::OpticalPower{T}
    incidence::IncidenceAngle{T}

    function ThinLens(
        power::OpticalPower{A},
        n::B,
        incidence::IncidenceAngle{C},
    ) where {A<:Real,B<:Real,C<:Real}
        return new{promote_type(A, B, C)}(n, power, incidence)
    end
end

ThinLens(power::OpticalPower{T}) where {T} = ThinLens(power, T(3//2), IncidenceAngle{T}())

function ThinLens(
    front::Curvature,
    back::Curvature,
    n::Real,
    incidence::IncidenceAngle=IncidenceAngle(),
)
    P = (n - 1) * (front - back)

    return ThinLens(OpticalPower(P.xx, P.xy, P.yy), n, incidence)
end

function get_transfer_matrix(e::ThinLens{T})::SMatrix{4,4,T,16} where {T}
    Pxx, Pxy, Pyy = rotate3(e.power.xx, e.power.xy, e.power.yy, -e.incidence.θ)

    V = sqrt(e.n^2 - sin(e.incidence.ι)^2)
    W = (V * sec(e.incidence.ι) - 1) / (e.n - 1)

    Mxx, Mxy, Myy = rotate3(
        -Pxx * W * sec(e.incidence.ι),
        -Pxy * W,
        -Pyy * W * cos(e.incidence.ι),
        e.incidence.θ,
    )

    return SA{T}[
        1.0 0.0 0.0 0.0;
        0.0 1.0 0.0 0.0;
        Mxx Mxy 1.0 0.0;
        Mxy Myy 0.0 1.0;
    ]
end

function Base.getproperty(e::ThinLens, s::Symbol, args...)
    if s==:rtm
        Base.depwarn(
            "Accessing field rtm of ThinLens directly is deprecated. Use get_transfer_matrix instead!",
            :getproperty,
        )
        return get_transfer_matrix(e)
    else
        return getfield(e, s, args...)
        # throw(FieldError(typeof(e), s))
    end
end
