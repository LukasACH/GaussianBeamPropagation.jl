"""
    Beam(λ, x::AxisParameters, y::AxisParameters=x; φ=0)
    Beam(λ, moments::SecondOrderMoments)
    Beam(λ, matrix::SMatrix{4,4})

A beam with a wavelength `λ`, described by its second-order moments.

A simple astigmatic Gaussian beam is given by its parameters along its two principal axes `x`
and `y` (see [`AxisParameters`](@ref)), which are rotated by `φ` relative to the _x_ and _y_
directions. A general beam is given by its [`SecondOrderMoments`](@ref), or by the symmetric 4×4
matrix of second-order moments in the order `(x, y, θx, θy)`.

The fields are the wavelength `λ` and the second-order `moments`.

A `Beam` is propagated by multiplying optical elements onto it from the left, and
its properties are obtained with, e.g., [`radius`](@ref), [`divergence`](@ref),
[`phase_curvature`](@ref) and [`twist`](@ref).

# Examples

```julia
# Elliptical beam at its waist, with the principal axes rotated by 30°
beam = Beam(1.03e-6, Waist(2e-3), Waist(1e-3); φ=deg2rad(30))
radius(FreeSpace(0.5) * beam)
```
"""
struct Beam{F<:Real}
    λ::F
    moments::SecondOrderMoments{F}
end

Beam(λ::Real, moments::SecondOrderMoments{T}) where {T} = Beam(T(λ), moments)

Beam(λ::Real, matrix::SMatrix{4,4,<:Real,16}) = Beam(λ, SecondOrderMoments(matrix))

function Beam(λ::Real, x::AxisParameters, y::AxisParameters=x; φ::Real=0)
    rxrx, rxθx, θxθx = axis_moments(λ, x)
    ryry, ryθy, θyθy = axis_moments(λ, y)
    moments = SecondOrderMoments(;
        spatial=PrincipalMoments(rxrx, ryry, φ),
        angular=PrincipalMoments(θxθx, θyθy, φ),
        mixed=PrincipalMoments(rxθx, ryθy, φ),
    )
    return Beam(λ, moments)
end

SecondOrderMoments(beam::Beam) = beam.moments

moment_matrix(beam::Beam) = moment_matrix(beam.moments)

function Base.isapprox(a::Beam, b::Beam; kwargs...)
    return isapprox(a.λ, b.λ; kwargs...) && isapprox(a.moments, b.moments; kwargs...)
end

function Base.getproperty(beam::Beam, s::Symbol)
    if s === :inner
        Base.depwarn("The field inner of Beam is deprecated, use moments instead.", :getproperty)
        return getfield(beam, :moments)
    end
    return getfield(beam, s)
end

function Beam(λ::Real; φ::Real=0, kwargs...)
    Base.depwarn(
        "Beam(λ; wx, wy, rx, ry, …) is deprecated, use Beam(λ, x, y; φ) with BeamRadius, " *
        "Waist, RayleighRange or ComplexBeamParameter for the axes instead.",
        :Beam,
    )
    return Beam(λ, legacy_axis_parameters(kwargs, :x), legacy_axis_parameters(kwargs, :y); φ)
end

# The axis parameters from the keywords of the deprecated keyword constructor
function legacy_axis_parameters(kwargs, axis::Symbol)
    value(name) = get(kwargs, Symbol(name, axis), missing)
    z, zR, w, r, q = value(:z), value(:zR), value(:w), value(:r), value(:q)

    !ismissing(z) && !ismissing(zR) && return RayleighRange(zR; z)
    !ismissing(w) && !ismissing(r) && return BeamRadius(w, r)
    !ismissing(q) && return ComplexBeamParameter(q)
    throw(ArgumentError("insufficient parameters: each axis needs (z, zR), (w, r) or q"))
end
