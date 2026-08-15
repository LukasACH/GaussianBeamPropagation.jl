const TMatrix{F} = SMatrix{4,4,F,16} where {F<:AbstractFloat}

struct Beam{F<:AbstractFloat,T<:Union{SecondOrderMoments{F},TMatrix{F}}}
    λ::F
    inner::T
    # rtm::SMatrix{4,4,Float64,16}
    # som::SecondOrderMoments

    # Beam(λ, rtm::SMatrix{4,4,Float64,16}) = new(λ, rtm, SecondOrderMoments(rtm))
    # Beam(λ, som::SecondOrderMoments) = new(λ, generate_wigner_matrix(som), som)
end

Beam(λ, matrix::TMatrix{F}) where {F} = Beam{TMatrix{F}}(λ, matrix)
Beam(λ, som::SecondOrderMoments{F}) where {F} = Beam{F,SecondOrderMoments{F}}(λ, som)

function Beam{F,SecondOrderMoments{F}}(beam::Beam{F,TMatrix{F}}) where {F}
    return Beam{F,SecondOrderMoments{F}}(beam.λ, SecondOrderMoments(beam.inner))
end

function Beam{F,TMatrix{F}}(beam::Beam{F,SecondOrderMoments{F}}) where {F}
    return Beam{F,TMatrix{F}}(beam.λ, generate_wigner_matrix(beam.inner))
end

Beam(λ::Real; kwargs...) = Beam(λ, SecondOrderMoments(λ; kwargs...))

function get_transfer_matrix(beam::Beam{F,SecondOrderMoments{F}}) where {F}
    som = beam
    return SA{F}[
        som.rxrx som.rxry som.rxθx som.rxθy;
        som.rxry som.ryry som.ryθx som.ryθy;
        som.rxθx som.ryθx som.θxθx som.θxθy;
        som.rxθy som.ryθy som.θxθy som.θyθy;
    ]
end

function generate_wigner_matrix(som::SecondOrderMoments{F}) where {F}
    return SA{F}[
        som.rxrx som.rxry som.rxθx som.rxθy;
        som.rxry som.ryry som.ryθx som.ryθy;
        som.rxθx som.ryθx som.θxθx som.θxθy;
        som.rxθy som.ryθy som.θxθy som.θyθy;
    ]
end

SecondOrderMoments(beam::Beam) = beam.som