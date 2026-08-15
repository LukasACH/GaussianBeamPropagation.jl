"""
    IncidenceAngle(ι, θ)

Description of an incidence direction different from the surface normal of the mirror, lens or interface. It is given
by an incidence angle _ι_, the angle between the incident direction and the surface normal, and the rotation angle _θ_
giving the orientation of the incidence plane relative to the _x_-direction in a right-handed coordinate system.
"""
struct IncidenceAngle{T}
    ι::T
    θ::T
end

# function IncidenceAngle{T}(; ι::Real=0, θ::Real=0) where {T}
#     IncidenceAngle{T}(convert(T, ι), convert(T, θ))
# end

function IncidenceAngle(ι::A, θ::B) where {A<:Real,B<:Real}
    T = promote_type(A, B)
    return IncidenceAngle{T}(T(ι), T(θ))
end

IncidenceAngle(; ι::Real=0, θ::Real=0) = IncidenceAngle(promote(ι, θ)...)

function Base.convert(::Type{IncidenceAngle{T}}, x::IncidenceAngle) where {T<:Real}
    return IncidenceAngle{T}(x)
end

IncidenceAngle{T}(x::IncidenceAngle) where {T} =
    IncidenceAngle{T}(x.ι, x.θ)
