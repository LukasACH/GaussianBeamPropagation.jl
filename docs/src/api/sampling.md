# Sampling beams

To plot a beam quantity along a scan, e.g. along the propagation distance or a varied parameter, it has to be evaluated for many beams.
The quantities along _x_ and _y_ jump whenever the principal axes of the beam rotate past ±45° (see [`radius`](@ref)).
`sample_beam_continuous` follows each principal axis instead, which gives continuous curves.

```@example sampling
using GaussianBeamPropagation

λ = 1.03e-6
ε2 = (λ / 4π)^2
# The principal axes rotate by 180° over the scan
beams = [
    Beam(λ, SecondOrderMoments(;
        rx=4e-6, ry=1e-6, φ_r=α, θx=ε2 / 4e-6, θy=ε2 / 1e-6, φ_θ=α,
        rθx=0.0, rθy=0.0, φ_rθ=0.0, twist=0.0,
    )) for α in range(0, π; length=7)
]

sample_beam(beams).x
```

```@example sampling
sample_beam_continuous(beams).x
```

```@docs; canonical=false
sample_beam
sample_beam_continuous
sample_beam_at_angle
```
