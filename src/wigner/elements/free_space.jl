"""
    FreeSpace

TODO
"""
struct FreeSpace <: AbstractElement
    rtm::SMatrix{4,4,Float64,16}

    function FreeSpace(d::Real)
        rtm = SA_F64[
            1.0 0.0 d 0.0;
            0.0 1.0 0.0 d;
            0.0 0.0 1.0 0.0;
            0.0 0.0 0.0 1.0;
        ]
        return new(rtm)
    end
end

# struct FreeSpace{T<:AbstractFloat} <: AbstractElement
#     d::T
# end

# function get_transfer_matrixs(e::FreeSpace{T}) where {T}
#     return SA{T}[
#         1 0 e.d 0
#         0 1 0 e.d
#         0 0 1 0
#         0 0 0 1
#     ]
# end
