"""
    Mirror(surface::Curvature)
    Mirror(surface::Curvature, incidence::IncidenceAngle)
    Mirror(power::OpticalPower)
    Mirror(power::OpticalPower, incidence::IncidenceAngle)

Matrix representing a mirror with surface curvature `surface` or an optical power of `power`. The incident direction of
the beam is given by the optional parameter `incidence`.

The curvature follows the sign convention of [`Curvature`](@ref): a concave mirror, whose center of curvature lies
before the mirror, has a negative curvature _κ_ and focuses with the optical power _P_ = −2_κ_.
"""
struct Mirror{T<:AbstractFloat} <: AbstractElement
    surface::Curvature{T}
    incidence::IncidenceAngle{T}

    function Mirror(
        surface::Curvature{A},
        incidence::IncidenceAngle{B}=IncidenceAngle(),
    ) where {A<:Real,B<:Real}
        return new{float(promote_type(A, B))}(surface, incidence)
    end
end

function Mirror(power::OpticalPower, args...)
    return Mirror(Curvature(-power.xx / 2, -power.xy / 2, -power.yy / 2), args...)
end

function get_transfer_matrix(e::Mirror{T})::SMatrix{4,4,T,16} where {T}
    κxx, κxy, κyy = rotate_symmetric(e.surface.xx, e.surface.xy, e.surface.yy, -e.incidence.θ)

    Mxx, Mxy, Myy = rotate_symmetric(
        2κxx * sec(e.incidence.ι),
        2κxy,
        2κyy * cos(e.incidence.ι),
        e.incidence.θ,
    )

    return SA{T}[
        1.0 0.0 0.0 0.0;
        0.0 1.0 0.0 0.0;
        Mxx Mxy 1.0 0.0;
        Mxy Myy 0.0 1.0;
    ]
end
