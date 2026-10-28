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
