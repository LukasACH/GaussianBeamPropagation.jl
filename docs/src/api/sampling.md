# Sampling beams

To plot a beam quantity along a scan, e.g. along the propagation distance or a varied parameter, it has to be evaluated for many beams.
The quantities along _x_ and _y_ jump whenever the principal axes of the beam rotate past ±45° (see [`radius`](@ref)).
`sample_beam_continuous` follows each principal axis instead, which gives continuous curves.

```@example sampling
using GaussianBeamPropagation

λ = 1.03e-6
# The principal axes rotate by 180° over the scan
beams = [
    Beam(λ, Waist(4e-3), Waist(2e-3); φ=α) for α in range(0, π; length=7)
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
