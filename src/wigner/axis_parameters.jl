"""
    AxisParameters

Supertype of the descriptions of a Gaussian beam along one of its principal axes, used to
construct a [`Beam`](@ref):

  - [`BeamRadius`](@ref): beam radius and radius of curvature of the phase front,
  - [`Waist`](@ref): waist radius and distance from the waist,
  - [`RayleighRange`](@ref): Rayleigh range and distance from the waist,
  - [`ComplexBeamParameter`](@ref): complex beam parameter _q_.

New descriptions are added by defining a subtype together with a method of
`GaussianBeamPropagation.axis_moments(λ, parameters)`, which returns the second moments
`(⟨r²⟩, ⟨rθ⟩, ⟨θ²⟩)` along that axis for the wavelength `λ`.
"""
abstract type AxisParameters end

"""
    axis_moments(λ, parameters::AxisParameters)

Second moments `(⟨r²⟩, ⟨rθ⟩, ⟨θ²⟩)` of a Gaussian beam with the wavelength `λ` along one axis.
"""
function axis_moments end

"""
    BeamRadius(w, R=Inf)

A beam with the 1/e² radius `w` and the radius of curvature `R` of its phase front. A positive
`R` belongs to a diverging beam, and `R = Inf` to a beam at its waist.
"""
struct BeamRadius{T<:Real} <: AxisParameters
    w::T
    R::T
end

BeamRadius(w::Real, R::Real=Inf) = BeamRadius(promote(w, R)...)

function axis_moments(λ::Real, p::BeamRadius)
    ε = λ / 4π
    w2 = (p.w / 2)^2
    return w2, w2 / p.R, w2 / p.R^2 + ε^2 / w2
end

"""
    Waist(w0; z=0)

A beam with the waist radius `w0`, at the distance `z` behind its waist. A positive `z` belongs
to a diverging beam.
"""
struct Waist{T<:Real} <: AxisParameters
    w0::T
    z::T
end

Waist(w0::Real; z::Real=0) = Waist(promote(w0, z)...)

axis_moments(λ::Real, p::Waist) = axis_moments(λ, RayleighRange(π * p.w0^2 / λ; z=p.z))

"""
    RayleighRange(zR; z=0)

A beam with the Rayleigh range `zR`, at the distance `z` behind its waist. A positive `z`
belongs to a diverging beam.
"""
struct RayleighRange{T<:Real} <: AxisParameters
    zR::T
    z::T
end

RayleighRange(zR::Real; z::Real=0) = RayleighRange(promote(zR, z)...)

function axis_moments(λ::Real, p::RayleighRange)
    ε = λ / 4π
    return ε * (p.z^2 / p.zR + p.zR), ε * p.z / p.zR, ε / p.zR
end

"""
    ComplexBeamParameter(q)

A beam with the complex beam parameter `q = z + i zR`.
"""
struct ComplexBeamParameter{T<:Real} <: AxisParameters
    q::Complex{T}
end

function axis_moments(λ::Real, p::ComplexBeamParameter)
    ε = λ / 4π
    return ε * abs2(p.q) / imag(p.q), ε * real(p.q) / imag(p.q), ε / imag(p.q)
end
