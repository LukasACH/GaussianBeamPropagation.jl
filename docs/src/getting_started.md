# How to get started

## Installation

!!! note

    This package is not yet registered in the official repository.
    The installation instructions will be updated as soon as it is.

Enter the package manager REPL by typing `]` from the Julia REPL.
As the package is not yet officially registered, you can install the package using:

```julia-repl
pkg> add https://github.com/LukasACH/GaussianBeamPropagation.jl
```

## Usage

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

## Pluto.jl notebook

An example notebook can be found in the repository under [`examples/`](https://github.com/LukasACH/GaussianBeamPropagation.jl/tree/main/examples).
