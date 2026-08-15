"""
    Mirror(surface::Curvature; incidence::IncidenceAngle)
    Mirror(power::OpticalPower; incidence::IncidenceAngle)

Matrix representing a mirror with surface curvature `surface` or an optical power of `power`. The incident direction of
the beam is given by the optional parameter `incidence`.
"""
struct Mirror{T} <: AbstractElement
    surface::Curvature{T}
    incidence::IncidenceAngle{T}

    rtm::SMatrix{4,4,T,16}

    function Mirror(surface::Curvature; incidence::IncidenceAngle=IncidenceAngle())
        κxx, κxy, κyy = rotate3(surface.xx, surface.xy, surface.yy, -incidence.θ)

        Mxx, Mxy, Myy = rotate3(
            κxx / 2 * sec(incidence.ι), κxy / 2, κyy / 2 * cos(incidence.ι), incidence.θ,
        )

        R = typeof(Mxx)

        rtm = SA{R}[
            1.0 0.0 0.0 0.0;
            0.0 1.0 0.0 0.0;
            Mxx Mxy 1.0 0.0;
            Mxy Myy 0.0 1.0;
        ]

        return new{R}(surface, incidence, rtm)
    end
end

function Mirror(power::OpticalPower; kwargs...)
    return Mirror(Curvature(-2power.xx, -2power.xy, -2power.yy); kwargs...)
end

# """
#     Mirror(D::Real; ι::Real=0.0, θ::Real=0.0)

# Matrix representing a mirror with focal power `D` in both principal axes. The
# beam is incident on the surface by `ι`, and the incidence plane is rotated by
# `θ` w.r.t. the beam coordinate system.
# """
# function Mirror(D::Real; ι::Real=0.0, θ::Real=0.0)
#     return Mirror(D, D; ι, θ)
# end

# function Mirror(Dx::Real, Dy::Real; φ::Real=0.0, ι::Real=0.0, θ::Real=0.0)
#     Mirror(FocalPower(Dx, Dy; φ); ia=IncidenceAngle(ι, θ))
# end

# function Mirror(
#     surface::Union{FocalPower,FocalLength,RadiusOfCurvature}, ia::IncidenceAngle
# )
#     Mirror(surface; ia)
# end

# function Mirror(surface::FocalLength; ia::IncidenceAngle=IncidenceAngle())
#     Mirror(FocalPower(surface); ia)
# end

# function Mirror(
#     surface::RadiusOfCurvature{T}; ia::IncidenceAngle=IncidenceAngle()
# ) where {T}
#     Mirror(FocalPower{T}(2*surface.Ra_inv, 2*surface.Ra_inv, 2*surface.Ra_inv); ia)
# end
