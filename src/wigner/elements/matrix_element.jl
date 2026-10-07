"""
    MatrixElement(matrix::SMatrix{4,4})

# Element that is defined by an arbitrary transfer matrix.
"""
struct MatrixElement{T} <: AbstractElement
    matrix::SMatrix{4,4,T,16}
end

function get_transfer_matrix(e::MatrixElement{T})::SMatrix{4,4,T,16} where {T}
    return e.matrix
end

function Base.getproperty(e::MatrixElement, s::Symbol, args...)
    if s==:rtm
        Base.depwarn(
            "Accessing field rtm of MatrixElement directly is deprecated. Use get_transfer_matrix instead!",
            :getproperty,
        )
        return get_transfer_matrix(e)
    else
        return getfield(e, s, args...)
        # throw(FieldError(typeof(e), s))
    end
end
