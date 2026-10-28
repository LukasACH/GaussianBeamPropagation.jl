function Base.:*(M::AbstractElement, beam::Beam)
    M_rtm=get_transfer_matrix(M)
    beam_matrix = generate_wigner_matrix(beam.inner)
    return Beam(beam.λ, M_rtm * beam_matrix * M_rtm')
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
