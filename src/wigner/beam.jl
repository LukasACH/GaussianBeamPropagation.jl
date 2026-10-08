const TMatrix{F} = SMatrix{4,4,F,16} where {F<:AbstractFloat}

struct Beam{F<:Real}
    λ::F
    inner::SecondOrderMoments{F}
    # rtm::SMatrix{4,4,Float64,16}
    # som::SecondOrderMoments

    # Beam(λ, rtm::SMatrix{4,4,Float64,16}) = new(λ, rtm, SecondOrderMoments(rtm))
    # Beam(λ, som::SecondOrderMoments) = new(λ, generate_wigner_matrix(som), som)
end

Beam(λ::Real; kwargs...) = Beam(λ, SecondOrderMoments(λ; kwargs...))

Beam(λ::Real, rtm::SMatrix{4,4,<:Real,16}) = Beam(λ, SecondOrderMoments(rtm))

function generate_wigner_matrix(beam::Beam{F}) where {F}
    return generate_wigner_matrix(beam.inner)
end

function generate_wigner_matrix(som::SecondOrderMoments{F}) where {F}
    return SA{F}[
        som.rxrx som.rxry som.rxθx som.rxθy;
        som.rxry som.ryry som.ryθx som.ryθy;
        som.rxθx som.ryθx som.θxθx som.θxθy;
        som.rxθy som.ryθy som.θxθy som.θyθy;
    ]
end

SecondOrderMoments(beam::Beam) = beam.inner