"""
    Mirror(surface::Curvature)
    Mirror(surface::Curvature, incidence::IncidenceAngle)
    Mirror(power::OpticalPower)
    Mirror(power::OpticalPower, incidence::IncidenceAngle)

Matrix representing a mirror with surface curvature `surface` or an optical power of `power`. The incident direction of
the beam is given by the optional parameter `incidence`.
"""
struct Mirror{T<:AbstractFloat} <: AbstractElement
    surface::Curvature{T}
    incidence::IncidenceAngle{T}

    function Mirror(
        surface::Curvature{A},
        incidence::IncidenceAngle{B}=IncidenceAngle(),
    ) where {A<:Real,B<:Real}
        return new{promote_type(A, B)}(surface, incidence)
    end
end

function Mirror(power::OpticalPower, args...)
    return Mirror(Curvature(-2power.xx, -2power.xy, -2power.yy), args...)
end

function get_transfer_matrix(e::Mirror{T})::SMatrix{4,4,T,16} where {T}
    κxx, κxy, κyy = rotate3(e.surface.xx, e.surface.xy, e.surface.yy, -e.incidence.θ)

    Mxx, Mxy, Myy = rotate3(
        κxx / 2 * sec(e.incidence.ι),
        κxy / 2,
        κyy / 2 * cos(e.incidence.ι),
        e.incidence.θ,
    )

    return SA{T}[
        1.0 0.0 0.0 0.0;
        0.0 1.0 0.0 0.0;
        Mxx Mxy 1.0 0.0;
        Mxy Myy 0.0 1.0;
    ]
end

function Base.getproperty(e::Mirror, s::Symbol, args...)
    if s==:rtm
        Base.depwarn(
            "Accessing field rtm of Mirror directly is deprecated. Use get_transfer_matrix instead!",
            :getproperty,
        )
        return get_transfer_matrix(e)
    else
        return getfield(e, s, args...)
        # throw(FieldError(typeof(e), s))
    end
end
