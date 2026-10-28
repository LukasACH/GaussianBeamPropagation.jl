"""
    Beam(λ, som::SecondOrderMoments)
    Beam(λ, rtm::SMatrix{4,4})

A beam with a wavelength `λ`, described by its second-order moments.

The moments can be given as a [`SecondOrderMoments`](@ref), or as the symmetric 4×4 matrix of
second-order moments in the order `(x, y, θx, θy)`. The second-order moments can be obtained
from a `Beam` with the function `SecondOrderMoments(beam)`.

A `Beam` is propagated by multiplying optical elements onto it from the left, and
its properties are obtained with, e.g., [`radius`](@ref), [`divergence`](@ref),
[`phase_curvature`](@ref) and [`twist`](@ref).

# Examples

```julia
beam = Beam(1.03e-6, SecondOrderMoments(;
    rx=(1e-3)^2, ry=(1e-3)^2, φ_r=0.0,
    θx=1e-8, θy=1e-8, φ_θ=0.0,
    rθx=0.0, rθy=0.0, φ_rθ=0.0, twist=0.0,
))
radius(FreeSpace(0.5) * beam)
```
"""
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