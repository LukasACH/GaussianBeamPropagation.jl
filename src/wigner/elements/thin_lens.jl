"""
    ThinLens(focussing::Union{FocalLength,FocalPower} [, n::Real] [, ia::IncidenceAngle])
    ThinLens(front::RadiusOfCurvature, back::RadiusOfCurvature [, n::Real] [, ia::IncidenceAngle])

The zero-thickness approximation of the thick lens. The oblique astigmatism depends on the refractive indices of the
lens material and the surrounding medium, which are given as a single ratio of internal / external, and the incidence angles are given by `IncidenceAngle`. The focussing of the lens can either be given as `RadiusOfCurvature`,
`FocalLength`, or `FocalPower`.
"""
struct ThinLens{T} <: AbstractElement
    n::T
    power::OpticalPower{T}
    incidence::IncidenceAngle{T}

    rtm::SMatrix{4,4,T,16}

    function ThinLens(
        power::OpticalPower;
        n::Real=1.5f0,
        incidence::IncidenceAngle=IncidenceAngle(),
    )
        # R=promote_type(A, B, typeof(n_internal), typeof(n_external), D)

        Pxx, Pxy, Pyy = rotate3(power.xx, power.xy, power.yy, -incidence.θ)

        V = sqrt(n^2 - sin(incidence.ι)^2)
        W = (V * sec(incidence.ι) - 1) / (n - 1)

        Mxx, Mxy, Myy = rotate3(
            -Pxx * W * sec(incidence.ι),
            -Pxy * W,
            -Pyy * W * cos(incidence.ι),
            incidence.θ,
        )

        R = typeof(Mxx)

        rtm = SA{R}[
            1.0 0.0 0.0 0.0;
            0.0 1.0 0.0 0.0;
            Mxx Mxy 1.0 0.0;
            Mxy Myy 0.0 1.0;
        ]

        return new{R}(n, power, incidence, rtm)
    end
end

function ThinLens(front::Curvature, back::Curvature; n::Real=1.5f0, kwargs...)
    P = (n - 1) * (front - back)

    return ThinLens(OpticalPower(P.xx, P.xy, P.yy); n, kwargs...)
end
