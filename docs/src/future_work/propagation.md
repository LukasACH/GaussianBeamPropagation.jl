# Propagation along a system

!!! note "Not part of the package"

    This code is not loaded by the package. It is kept here as a starting point for future work.

## Idea

To plot the beam size along an optical system, the beam has to be known at every position `z`, not only at the output.
Two approaches were started:

- **Sampling between elements.** `generate_propagator` walks through the elements of a system, stores the beam after each element together with its position, and returns a callable `BeamPropagator`.
  Calling it with a position `z` takes the last stored beam before `z` and propagates it through the remaining free space.
- **Analytic free-space propagation.** In free space, the spatial moments are quadratic polynomials in `z`,
  ``\langle x^2\rangle(z) = \langle x^2\rangle + 2z\langle xθ_x\rangle + z^2\langle θ_x^2\rangle``,
  so the beam radius along `z` can be evaluated without building a `Beam` at each position.
  The same idea is used by the fitting code.

## State

The code was written for an older version of the package and does not run anymore:

- `generate_propagator` uses the removed 2D types (`GaussianBeam2D`, `FreeSpace2D`, `BeamPropagator2D`), a `components` field and a `distance(element)` function that do not exist.
- `propagation` uses `SecondOrderBeamMoments` from the [moment types](moment_types.md) draft.
- `Propagator` is unused.

The natural place for this functionality is the planned [`OpticalAssembly`](optical_assembly.md), which knows the positions of its elements.

## Last state of the code

```julia
struct BeamPropagator
    distances::Any
    beams::Any
end

function (prop::BeamPropagator)(x)
    i = findlast(prop.distances .<= x)
    last_beam = prop.beams[i]
    rem_dist = x - prop.distances[i]
    return FreeSpace(rem_dist) * last_beam
end

function generate_propagator(assembly::Vector{Union{FreeSpace,ThinLens}}, beam::Beam)
    distances = Vector{Float64}([0.0])
    beams = Vector{GaussianBeam2D}([beam])

    current_beam = beam
    current_distance = 0.0

    for elem in reverse(assembly.components)
        current_distance += distance(elem)
        current_beam = elem * current_beam
        push!(distances, current_distance)
        push!(beams, current_beam)
    end

    # @info "Propagator values" distances beams

    function propagator(x)
        i = findlast(distances .<= x)
        last_beam = beams[i]
        rem_dist = x - distances[i]
        return FreeSpace2D(rem_dist) * last_beam
    end

    return BeamPropagator2D(distances, beams)
end

struct Propagator
    xx::Function
    xy::Function
    yy::Function
end

function propagation(som::SecondOrderBeamMoments, z::Vector{T}) where {T<:Real}
    xx(z) = som.spatial.xx + 2z * som.mixed.xx + z^2 * som.angluar.xx
    xy(z) = som.spatial.xy + z * (som.mixed.xy + som.mixed.yx) + z^2 * som.angluar.xy
    yy(z) = som.spatial.yy + 2z * som.mixed.yy + z^2 * som.angluar.yy

    return radius.(xx.(z), xy.(z), yy.(z))
end
```
