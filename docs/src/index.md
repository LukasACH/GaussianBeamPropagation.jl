# Example.jl Documentation

```@repl
using GaussianBeamPropagation
using Unitful # optional, but enables the use of the @u_str macros
λ = 1030u"nm"
w0_in = 2.75u"mm"
beam_in = Beam(λ;
    wx=w0_in,
    wy=w0_in,
    rx=Inf * u"m",
    ry=Inf * u"m",
)
M = FreeSpace(0.5) *
        ThinLens(OpticalPower(1 / 0.5)) *
        FreeSpace(0.5)
beam_out = M * beam_in
radius(beam_out)
phase_curvature(beam_out)
```

```@docs
AbstractElement
```

```@docs
FreeSpace
```

```@docs
MatrixElement
```

```@docs
Mirror
```

```@docs
OpticalInterface
```

```@docs
ThinLens
```

```@docs
ThickLens
```

```@docs
Curvature
```

```@docs
IncidenceAngle
```
