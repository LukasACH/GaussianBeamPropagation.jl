"""
    FreeSpace(d::Real)

# Matrix representing a propagation of the beam through free space.

The arguments are the following:

  - `d::Real`: propagation distance
"""
struct FreeSpace{T<:Real} <: AbstractElement
    d::T

    FreeSpace(d::T) where {T<:Real} = new{float(T)}(d)
end

function get_transfer_matrix(e::FreeSpace{T})::SMatrix{4,4,T,16} where {T}
    return SA{T}[
        1.0 0.0 e.d 0.0
        0.0 1.0 0.0 e.d
        0.0 0.0 1.0 0.0
        0.0 0.0 0.0 1.0
    ]
end

function Base.getproperty(e::FreeSpace, s::Symbol, args...)
    if s==:rtm
        Base.depwarn(
            "Accessing field rtm of FreeSpace directly is deprecated. Use get_transfer_matrix instead!",
            :getproperty,
        )
        return get_transfer_matrix(e)
    else
        return getfield(e, s, args...)
        # throw(FieldError(typeof(e), s))
    end
end
