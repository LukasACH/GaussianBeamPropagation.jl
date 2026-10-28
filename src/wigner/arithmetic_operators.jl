function Base.:*(M::AbstractElement, beam::Beam)
    T = get_transfer_matrix(M)
    return Beam(beam.λ, T * moment_matrix(beam) * T')
end

function Base.:*(Mlhs::AbstractElement, Mrhs::AbstractElement)
    return MatrixElement(get_transfer_matrix(Mlhs) * get_transfer_matrix(Mrhs))
end

"""
    element^p

Pass `p` times through `element`. Negative powers give the inverse element, and `p = 0` the
identity.
"""
Base.:^(e::AbstractElement, p::Integer) = MatrixElement(get_transfer_matrix(e)^p)

"""
    inv(element)

The element that undoes `element`, e.g. `FreeSpace(-d)` for `FreeSpace(d)`.
"""
Base.inv(e::AbstractElement) = MatrixElement(inv(get_transfer_matrix(e)))
