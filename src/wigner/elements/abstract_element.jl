"""
    AbstractElement

Supertype of all elements that can be used in the Wigner formalism. This includes optical interfaces, free space propagation, and other elements that can be represented by a transfer matrix.
"""
abstract type AbstractElement end

function Base.getproperty(e::AbstractElement, s::Symbol)
    if s === :rtm
        Base.depwarn(
            "Accessing field rtm of $(nameof(typeof(e))) directly is deprecated. Use get_transfer_matrix instead!",
            :getproperty,
        )
        return get_transfer_matrix(e)
    end
    return getfield(e, s)
end
