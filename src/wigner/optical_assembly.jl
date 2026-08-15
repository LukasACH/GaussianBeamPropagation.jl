struct OpticalAssembly{T<:AbstractFloat}
    rtm::TMatrix{T}
    elements::Vector{AbstractElement}
    positions::Vector{T}
    length::T
end

function Base.:*(
    lhs::OpticalAssembly{L}, rhs::OpticalAssembly{R},
) where {L<::AbstractFloat,R<:AbstractFloat}
    T = promote_type(L, R)
    return OpticalAssembly{T}(
        lhs.rtm * rhs.rtm,
        vcat(rhs.elements, lhs.elements),
        vcat(rhs.positions, length .+ lhs.positions),
        lhs.length+rhs.length,
    )
end

OpticalAssembly(e::AbstractElement) = OpticalAssembly(get_transfer_matrix(e), [e], [0.0])

function Base.:*(lhs::AbstractElement, rhs::AbstractElement) end