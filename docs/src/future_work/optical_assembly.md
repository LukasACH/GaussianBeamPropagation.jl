# Optical assemblies

!!! note "Not part of the package"

    This code is not loaded by the package. It is kept here as a starting point for future work.

## Idea

At the moment, multiplying two elements collapses them into a single `MatrixElement`: the transfer matrix is all that remains, and the individual elements and their positions are lost.

An `OpticalAssembly` should instead keep the structure of the system:

- the combined transfer matrix (so propagating a beam through the whole assembly stays a single matrix product),
- the list of elements, in the order the beam passes them,
- the position of each element along the optical axis, and the total length.

With this, a beam could be [sampled along the system](propagation.md), elements could be inspected or replaced, and plots of the beam along a system become possible.

The groundwork for this has been done: the elements (`FreeSpace`, `ThinLens`, `Mirror`, …) only store their parameters, and their matrix is computed on demand with `get_transfer_matrix`.
This keeps the elements small enough to be stored inside an assembly.

## State of the draft

The draft below is incomplete and was never included in the package:

- `L<::AbstractFloat` is a syntax error (should be `L<:AbstractFloat`).
- The multiplication uses `length` as if it were the length of the left-hand assembly; it should be `rhs.length`.
- `OpticalAssembly(e::AbstractElement)` does not pass a length, and the element length (e.g. the distance of a `FreeSpace`) is not defined anywhere.
- The empty method `Base.:*(lhs::AbstractElement, rhs::AbstractElement) end` would replace the existing element multiplication, which returns a `MatrixElement`. It has to be designed together with the existing method.
- The `elements` vector is abstractly typed, so a type-stable alternative (e.g. a tuple) might be preferable for performance.

## Last state of the code

```julia
struct OpticalAssembly{T<:AbstractFloat}
    rtm::TMatrix{T}
    elements::Vector{AbstractElement}
    positions::Vector{T}
    length::T
end

function Base.:*(
    lhs::OpticalAssembly{L},
    rhs::OpticalAssembly{R},
) where {L<::AbstractFloat,R<:AbstractFloat}
    T = promote_type(L, R)
    return OpticalAssembly{T}(
        lhs.rtm * rhs.rtm,
        vcat(rhs.elements, lhs.elements),
        vcat(rhs.positions, length .+ lhs.positions),
        lhs.length+rhs.length,
    )
end

OpticalAssembly(e::AbstractElement) = OpticalAssembly(get_transfer_matrix(e), [e], [0.0])

function Base.:*(lhs::AbstractElement, rhs::AbstractElement) end
```
