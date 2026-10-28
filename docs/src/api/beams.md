# Beams

A [`Beam`](@ref) consists of a wavelength and the [`SecondOrderMoments`](@ref) of the beam.

## Simple astigmatic beams

Most beams are described by their parameters along two perpendicular principal axes.
Each axis is given by one of the subtypes of [`AxisParameters`](@ref), and both axes can be described differently:

```julia
λ = 1.03e-6
Beam(λ, Waist(1e-3))                                    # round beam at its waist
Beam(λ, Waist(2e-3), Waist(1e-3); φ=deg2rad(30))        # elliptical beam, rotated by 30°
Beam(λ, BeamRadius(1e-3, 0.5), RayleighRange(2.0; z=1.0))
Beam(λ, ComplexBeamParameter(0.1 + 2.0im))
```

```@docs; canonical=false
AxisParameters
BeamRadius
Waist
RayleighRange
ComplexBeamParameter
```

New descriptions can be added by defining a subtype of `AxisParameters` and a method of `GaussianBeamPropagation.axis_moments`.

## General beams

Beams with general astigmatism are described by their second-order moments, which can be given by the principal axes of the spatial, angular and mixed moments, and the twist:

```julia
moments = SecondOrderMoments(;
    spatial=PrincipalMoments(1e-6, 4e-6, 0.1),
    angular=PrincipalMoments(1e-8, 2e-8, 0.4),
    twist=1e-10,
)
Beam(λ, moments)
```

```@docs; canonical=false
Beam
SecondOrderMoments
PrincipalMoments
```
