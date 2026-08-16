"""
    MatrixElement

TODO
"""
struct MatrixElement <: AbstractElement
    rtm::SMatrix{4,4,Float64,16}
end

# function get_transfer_matrix(e::MatrixElement)
#     return e.rtm
# end
