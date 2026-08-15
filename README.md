# GaussianBeamPropagation.jl

> [!WARNING]
> This package is very much a work in progress.
> The final API is not stable yet, and may change drastically over time until I am happy with at and stabilise it with a 1.0 release.

The goal of this package is to provide a framework to propagate Gaussian beams, e.g., Gaussian laser beams, through an optical system supporting general astigmatism.
It focusses on a fast propagation over perfect accuracy to enable numerical optimisation strategies to, e.g., minimise astigmatism in optical systems.

It implements extended ABCD matrices used with the [ray transfer matrix analysis](https://en.wikipedia.org/wiki/Ray_transfer_matrix_analysis) that supports both [non-rotationally symmetric elements and third-order/oblique astigmatism](<https://en.wikipedia.org/wiki/Astigmatism_(optical_systems)#Forms_of_astigmatism>).

## Installation

> [!NOTE]
> This package is not yet registered in the official repository.
> The installation instructions will be updated as soon as it is.

Enter the package manager REPL by pressing <kbd>]</kbd> from the Julia REPL.
As the package is not yet officially registered, you can install the package using:

```julia-repl
pkg> add https://github.com/LukasACH/GaussianBeamPropagation.jl
```

## Integration with other packages

### Unitful.jl

By default, all values are given in units of metres.
However, when loading the [Unitful.jl](https://juliaphysics.github.io/Unitful.jl/stable) package, all relevant functions get a unit-supporting version.

> [!IMPORTANT]
> Most of these are not yet done, and will be added over time.

## Example

See the [Pluto.jl](https://github.com/fonsp/Pluto.jl) examples [here](examples/).
A simple 4f re-imaging system with a magnification of 2 can be expressed as.

```julia-repl
julia> using GaussianBeamPropagation

julia> using Untiful # optional, but enables the use of the @u_str macros

julia> λ = 1030u"nm"
1030 nm

julia> w0_in = 2.75u"mm"
2.75 mm

julia> beam_in = GaussianBeam(λ;
            wx=w0_in,
            wy=w0_in,
            rx=Inf * u"m",
            ry=Inf * u"m",
        )
Beam{Float64, SecondOrderMoments{Float64}}(1.03e-6, SecondOrderMoments{Float64}(1.890625e-6, 0.0, 1.890625e-6, 3.553442767806815e-9, 0.0, 3.553442767806815e-9, 0.0, 0.0, 0.0, 0.0))

julia> M = FreeSpace(0.5) *
               ThinLens(OpticalPower(1 / 0.5)) *
               FreeSpace(0.5)
MatrixElement([0.0 0.0 0.5 0.0; 0.0 0.0 0.0 0.5; -2.0 0.0 0.0 0.0; 0.0 -2.0 0.0 0.0])

julia> beam_out = M * beam_in
Beam{Float64, StaticArraysCore.SMatrix{4, 4, Float64, 16}}(1.03e-6, [8.883606919517038e-10 0.0 0.0 0.0; 0.0 8.883606919517038e-10 0.0 0.0; 0.0 0.0 7.5625e-6 0.0; 0.0 0.0 0.0 7.5625e-6])

julia> radius(beam_out)
(x = 5.961076050350989e-5, y = 5.961076050350989e-5, φ = 0.0)

julia> phase_curvature(beam_out)
(x = -Inf, y = -Inf, φ = 0.0)
```

## To-Dos

- [ ] Write list of To-Dos

## Contributing

For any discussion and issues, please open an issue here. If your proposed changes are small, you can directly create a PR.
If the proposed changes are larger, it might save time on all sides if we discuss the issue beforehand!

## Credits

The matrices were largely copied from [Dupraz, K., Cassou, K., Martens, A., et al. 2019, Optics Communications, 443, 172](https://doi.org/10.1016/j.optcom.2019.03.041).
