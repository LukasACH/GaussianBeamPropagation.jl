"""
    PrincipalMoments(x, y=x, φ=0)

Second moments `x` and `y` along two principal axes, which are rotated by `φ` relative to the
_x_ and _y_ directions. Used to build [`SecondOrderMoments`](@ref).
"""
struct PrincipalMoments{T<:Real}
    x::T
    y::T
    φ::T
end

PrincipalMoments(x::Real, y::Real=x, φ::Real=0) = PrincipalMoments(promote(x, y, φ)...)

"""
    SecondOrderMoments(; spatial::PrincipalMoments, angular::PrincipalMoments,
                         mixed::PrincipalMoments=PrincipalMoments(0.0), twist=0.0)
    SecondOrderMoments(rxrx, rxry, ryry, θxθx, θxθy, θyθy, rxθx, rxθy, ryθx, ryθy)
    SecondOrderMoments(matrix::SMatrix{4,4})
    SecondOrderMoments(beam::Beam)

The second-order moments of a beam: the spatial moments ``⟨x²⟩, ⟨xy⟩, ⟨y²⟩`` (`rxrx`, `rxry`,
`ryry`), the angular moments ``⟨θx²⟩, ⟨θxθy⟩, ⟨θy²⟩`` (`θxθx`, `θxθy`, `θyθy`) and the mixed
moments ``⟨xθx⟩, ⟨xθy⟩, ⟨yθx⟩, ⟨yθy⟩`` (`rxθx`, `rxθy`, `ryθx`, `ryθy`).

They can be given directly, as the symmetric 4×4 matrix of moments in the order
`(x, y, θx, θy)`, or by the principal axes of the spatial, angular and mixed moments together
with the `twist`, the antisymmetric part of the mixed moments.

The moments do not depend on the wavelength. A beam with a wavelength is a [`Beam`](@ref).
"""
struct SecondOrderMoments{T}
    rxrx::T
    rxry::T
    ryry::T
    θxθx::T
    θxθy::T
    θyθy::T
    rxθx::T
    rxθy::T
    ryθx::T
    ryθy::T
end

function SecondOrderMoments(;
    spatial::Union{PrincipalMoments,Nothing}=nothing,
    angular::Union{PrincipalMoments,Nothing}=nothing,
    mixed::PrincipalMoments=PrincipalMoments(0.0),
    twist::Real=0.0,
    kwargs...,
)
    if !isempty(kwargs)
        Base.depwarn(
            "SecondOrderMoments(; rx, ry, φ_r, θx, θy, φ_θ, rθx, rθy, φ_rθ, twist) is deprecated, " *
            "use SecondOrderMoments(; spatial=PrincipalMoments(rx, ry, φ_r), " *
            "angular=PrincipalMoments(θx, θy, φ_θ), mixed=PrincipalMoments(rθx, rθy, φ_rθ), twist) instead.",
            :SecondOrderMoments,
        )
        (; rx, ry, φ_r, θx, θy, φ_θ, rθx, rθy, φ_rθ) = (; kwargs...)
        spatial = PrincipalMoments(rx, ry, φ_r)
        angular = PrincipalMoments(θx, θy, φ_θ)
        mixed = PrincipalMoments(rθx, rθy, φ_rθ)
    end
    (isnothing(spatial) || isnothing(angular)) &&
        throw(ArgumentError("the spatial and angular moments are required"))

    rxrx, rxry, ryry = tensor_from_principal(spatial.x, spatial.y, spatial.φ)
    θxθx, θxθy, θyθy = tensor_from_principal(angular.x, angular.y, angular.φ)
    rxθx, rθavg, ryθy = tensor_from_principal(mixed.x, mixed.y, mixed.φ)
    rxθy = rθavg + twist / 2
    ryθx = rθavg - twist / 2
    return SecondOrderMoments(promote(rxrx, rxry, ryry, θxθx, θxθy, θyθy, rxθx, rxθy, ryθx, ryθy)...)
end

function SecondOrderMoments(matrix::SMatrix{4,4,F,16}) where {F}
    rxrx = matrix[1, 1]
    rxry = (matrix[1, 2] + matrix[2, 1]) / 2
    ryry = matrix[2, 2]

    θxθx = matrix[3, 3]
    θxθy = (matrix[3, 4] + matrix[4, 3]) / 2
    θyθy = matrix[4, 4]

    rxθx = (matrix[1, 3] + matrix[3, 1]) / 2
    rxθy = (matrix[1, 4] + matrix[4, 1]) / 2
    ryθx = (matrix[2, 3] + matrix[3, 2]) / 2
    ryθy = (matrix[2, 4] + matrix[4, 2]) / 2

    return SecondOrderMoments(rxrx, rxry, ryry, θxθx, θxθy, θyθy, rxθx, rxθy, ryθx, ryθy)
end

"""
    moment_matrix(moments::SecondOrderMoments)
    moment_matrix(beam::Beam)

The symmetric 4×4 matrix of second-order moments in the order `(x, y, θx, θy)`.
"""
function moment_matrix(som::SecondOrderMoments{F}) where {F}
    return SA{F}[
        som.rxrx som.rxry som.rxθx som.rxθy;
        som.rxry som.ryry som.ryθx som.ryθy;
        som.rxθx som.ryθx som.θxθx som.θxθy;
        som.rxθy som.ryθy som.θxθy som.θyθy;
    ]
end

function Base.isapprox(a::SecondOrderMoments, b::SecondOrderMoments; kwargs...)
    return isapprox(moment_matrix(a), moment_matrix(b); kwargs...)
end

function SecondOrderMoments(λ::Real; kwargs...)
    Base.depwarn(
        "SecondOrderMoments(λ; …) is deprecated, use Beam(λ, x, y; φ).moments instead.",
        :SecondOrderMoments,
    )
    return Beam(λ; kwargs...).moments
end
