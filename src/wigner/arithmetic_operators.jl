function Base.:*(M::AbstractElement, beam::Beam{F,TMatrix{F}}) where {F}
    return Beam(beam.λ, M.rtm * beam.inner * M.rtm')
end

Base.:*(M::AbstractElement, beam::Beam{F}) where {F} = M * Beam{F,TMatrix{F}}(beam)


Base.:*(Mlhs::AbstractElement, Mrhs::AbstractElement) = MatrixElement(Mlhs.rtm * Mrhs.rtm)

Base.:^(e::AbstractElement, p::Integer) =
    if p < zero(p)
        error("Not Implemented")
    elseif iszero(p)
        return FreeSpace(0) * FreeSpace(0)
    elseif isone(p)
        return e
    elseif iseven(p) && p > 2
        return Base.:^(e^2, p >> 1)
    else
        return e * Base.:^(e, p - 1)
    end
