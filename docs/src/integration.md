# Integration with other packages

## Unitful.jl

By default, all values are given in units of metres.
When the [Unitful.jl](https://juliaphysics.github.io/Unitful.jl/stable) package is loaded, the following constructors also accept quantities with units.
The units are converted to metres when the object is constructed, so all results are plain numbers in metres.

| Constructor                                         | Quantities                     |
|:--------------------------------------------------- |:------------------------------ |
| [`Beam`](@ref)                                      | wavelength                     |
| [`BeamRadius`](@ref), [`Waist`](@ref), [`RayleighRange`](@ref), [`ComplexBeamParameter`](@ref) | lengths |
| [`Curvature`](@ref), [`OpticalPower`](@ref)         | inverse lengths                |
| [`SecondOrderMoments`](@ref) (ten moments)          | lengths², lengths and angles²  |

```@example unitful
using GaussianBeamPropagation, Unitful

beam = Beam(1030u"nm", Waist(2.75u"mm"), BeamRadius(2u"mm", 1u"m"))
radius(ThinLens(OpticalPower(2u"m^-1")) * beam)
```

!!! note

    Elements such as [`FreeSpace`](@ref) and the results of [`radius`](@ref) and the other beam quantities do not support units yet.
