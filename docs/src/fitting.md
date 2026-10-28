# Fitting beams to measurements

A beam is described by its [`SecondOrderMoments`](@ref), and these moments can be determined from measured beam sizes.
This page describes how to fit the input beam, and optionally the parameters of an optical system, to such measurements.

## What is measured, and what is fitted

Each measurement consists of the beam size and orientation in some plane: the radii `wx` and `wy` along the principal axes of the beam profile and the angle `φ` of these axes.
These are exactly the values returned by [`radius`](@ref), which means that data from a beam profiler and results from this package can be mixed freely.
[`moments_from_radius`](@ref) converts the measurements into the second-order moments ``\langle x^2\rangle``, ``\langle xy\rangle`` and ``\langle y^2\rangle`` that the fit compares against.

The fit has ten parameters for the input beam, in the order

| Parameters                | Meaning                                                |
|:------------------------- |:------------------------------------------------------ |
| `rx`, `ry`, `φ_r`         | spatial moments along the principal axes, and rotation |
| `θx`, `θy`, `φ_θ`         | angular moments along the principal axes, and rotation |
| `rθx`, `rθy`, `φ_rθ`      | mixed moments along the principal axes, and rotation   |
| `twist`                   | twist of the beam                                      |

(see [`SecondOrderMoments`](@ref) for details).
The sum of squared differences between the measured and the predicted moments is minimised with the Nelder-Mead algorithm from [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl), and the logarithm of this sum is used as objective.

!!! warning "Twist is not fitted"

    Beam sizes after free-space propagation do not depend on the twist, so it cannot be determined from them.
    `beam_fit_parameters` starts it at zero, but since the bounds are not enforced (see below), the optimiser may change it without any effect on the objective.
    Ignore the fitted twist, and only use the fit for beams without twist.

!!! warning "Bounds are not enforced"

    The lower and upper bounds are passed to Optim.jl, but Nelder-Mead does not enforce them.
    Fitted values can leave the bounds; the radii are made positive internally.
    Treat the bounds as a hint for the scale of the problem, not as a guarantee.

## Fitting a propagation measurement

The common case is to measure the beam size at several distances and to fit the beam that produced them.
[`fit_propagation`](@ref) takes the wavelength and a matrix with one row per measurement and the columns `wx`, `wy`, `φ` and the propagation distance.

To demonstrate this without a real measurement, we create the data from a known beam:

```@example fitting
using GaussianBeamPropagation
using Logging # hide
Logging.disable_logging(Logging.Info) # hide

λ = 1.03e-6
ε2 = (λ / 4π)^2 # smallest possible product of spatial and angular moments

# An astigmatic beam with a waist of 3 mm × 2 mm, rotated by 0.5 rad
beam = Beam(
    λ,
    SecondOrderMoments(;
        spatial=PrincipalMoments((3e-3 / 2)^2, (2e-3 / 2)^2, 0.5),
        angular=PrincipalMoments(4ε2 / (3e-3 / 2)^2, 4ε2 / (2e-3 / 2)^2, 0.5),
    ),
)

# "Measure" the beam size at eight positions
distances = range(0.1, 2.0; length=8)
data = reduce(vcat, [(w = radius(FreeSpace(z) * beam); [w[1] w[2] w[3] z]) for z in distances])
```

The beam is then recovered with

```@example fitting
result = fit_propagation(λ, data)
radius(result.beam)
```

```@example fitting
radius(beam)
```

The returned [`SystemFit`](@ref) holds the fitted `beam` as an ordinary [`Beam`](@ref), so it can be propagated like any other beam:

```@example fitting
radius(FreeSpace(3.0) * result.beam)
```

```@example fitting
radius(FreeSpace(3.0) * beam)
```

The final value of the objective is available as `result.minimum`, and the full result of the optimiser as `result.optim`.
A good fit to noiseless data reaches very negative values, here

```@example fitting
result.minimum
```

### Elements between the measurements

By default the beam propagates through free space.
If other elements are placed between the beam and the measurement plane, pass the `system` keyword.
It is called with one row of the distance columns and returns the element that transforms the beam:

```@example fitting
lens(z) = FreeSpace(z) * ThinLens(OpticalPower(1 / 0.5))

data_lens = reduce(vcat, [(w = radius(lens(z) * beam); [w[1] w[2] w[3] z]) for z in distances])

result_lens = fit_propagation(λ, data_lens; system=(row, _) -> lens(row[1]))
radius(result_lens.beam)
```

## The general case: fitting an optical system

[`fit_propagation`](@ref) is a thin layer on top of [`fit_system`](@ref), which can fit the input beam together with parameters of an arbitrary optical system.
The system is a function `system(row, params)` that returns the transfer element for one measurement.

- `row` is one row of the matrix `vars`, which holds the settings of each measurement, for example a distance or the number of passes through a system.
- `params` are the current values of the system parameters that are being fitted.

The measured beam sizes are passed as the moments returned by [`moments_from_radius`](@ref), and the bounds and starting values of the system parameters as [`FitParameters`](@ref).
As an example, we fit the optical power of a lens together with the beam.
The settings of each measurement are the distance `z` behind the lens and the position `d` of the lens, which are varied independently:

```@example fitting
system(row, params) = FreeSpace(row[1]) * ThinLens(OpticalPower(params[1])) * FreeSpace(row[2])

vars = reduce(vcat, [[z d] for z in (0.5, 1.0, 1.5) for d in (0.2, 0.5, 0.8)])

# Measurements for a lens with an optical power of 2 dioptres
sizes = [radius(system(row, [2.0]) * beam) for row in eachrow(vars)]
observed = moments_from_radius(first.(sizes), getindex.(sizes, 2), getindex.(sizes, 3))

result_system = fit_system(
    system, λ, vars, observed;
    system_params=FitParameters([0.5 1.5 4.0]), # one row per parameter: lower, initial, upper
)
result_system.system_params
```

The bounds and starting values for the input beam default to [`beam_fit_parameters`](@ref), which takes a range of the beam radius for the lower bound, the starting value and the upper bound.
Choosing it close to the expected size helps the fit to converge:

```julia
beam_params = beam_fit_parameters(λ; beam_radius=(1e-3, 3e-3, 10e-3))
result = fit_system(system, λ, vars, observed; beam_params)
```

!!! note "Choose measurements that separate beam and system"

    The input beam and the system parameters are fitted at the same time.
    If the measurements cannot tell their effects apart, the fit is not unique.
    For example, a free distance in front of the measurement is indistinguishable from a different position of the beam waist.
    Vary the settings of the system between the measurements (a column of `vars`) to break this degeneracy.

## Reference

```@docs; canonical=false
fit_propagation
fit_system
FitParameters
beam_fit_parameters
moments_from_radius
SystemFit
```
