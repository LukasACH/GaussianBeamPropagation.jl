# Moment types

!!! note "Not part of the package"

    This code is not loaded by the package. It is kept here as a starting point for future work.

## Idea

The second-order moments of a beam consist of three 2×2 blocks:
the spatial moments ``\langle x^2\rangle, \langle xy\rangle, \langle y^2\rangle``,
the angular moments ``\langle θ_x^2\rangle, \langle θ_xθ_y\rangle, \langle θ_y^2\rangle``,
and the mixed moments ``\langle xθ_x\rangle, \langle xθ_y\rangle, \langle yθ_x\rangle, \langle yθ_y\rangle``.
The spatial and angular blocks are symmetric, the mixed block is not.

Instead of storing ten loose numbers in `SecondOrderMoments`, these blocks could be represented by their own types:

- `SymmetricMoments` for the spatial and angular blocks (three values),
- `AsymmetricMoments` for the mixed block (four values),
- `SecondOrderBeamMoments` combining the three blocks.

This would make it possible to reuse code through the type system:

- The rotation helpers `rotate3` and `rotate4` act on exactly these blocks, and could become methods of a single `rotate(m::AbstractMoments, α)`.
- If the types behave like 2×2 matrices (e.g. by subtyping `StaticArrays.StaticMatrix(2, 2)` or by implementing the `AbstractMatrix` interface), the beam quantities could be computed directly from them: [`radius`](@ref) and [`divergence`](@ref) are the eigen-decompositions of the spatial and angular blocks, and the [`twist`](@ref) is the antisymmetric part of the mixed block.
- The same principal-axis decomposition is used by `Curvature` and `OpticalPower`, which could share it.

The open question is whether this makes the package simpler or just more complicated.
The implementation below was started, but is not finished: nothing in the package uses it, and the rotation helpers do not act on it yet.

## Last state of the code

```julia
abstract type AbstractMoments end

struct SymmetricMoments{T} <: AbstractMoments
    xx::T
    xy::T
    yy::T
end

function SymmetricMoments(xx::T, xy::T, yy::T) where {T<:AbstractFloat}
    return SymmetricMoments{T}(xx, xy, yy)
end

function SymmetricMoments(xx::T, yy::T) where {T<:AbstractFloat}
    return SymmetricMoments{T}(xx, zero(T), yy)
end

calculate_matrix(m::SymmetricMoments{T}) where {T} = SA{T}[
    m.xx m.xy
    m.xy m.yy
]

struct AsymmetricMoments{T} <: AbstractMoments
    xx::T
    xy::T
    yx::T
    yy::T
end

function AsymmetricMoments(xx::T, xy::T, yx::T, yy::T) where {T<:AbstractFloat}
    return AsymmetricMoments{T}(xx, xy, yx, yy)
end

function AsymmetricMoments(xx::T, yy::T; sum=zero(T), dif=zero(T)) where {T<:AbstractFloat}
    return AsymmetricMoments{T}(xx, (sum+dif)/2, (sum-dif)/2, yy)
end

calculate_matrix(m::AsymmetricMoments{T}) where {T} = SA{T}[
    m.xx m.xy
    m.yx m.yy
]

struct SecondOrderBeamMoments{T}
    spatial::SymmetricMoments{T}
    angular::SymmetricMoments{T}
    mixed::AsymmetricMoments{T}
end

function SecondOrderBeamMoments(rtm::SMatrix{4,4,T,16}) where {T}
    spatial = SymmetricMoments(rtm[1, 1], (rtm[1, 2] + rtm[2, 1]) / 2, rtm[2, 2])

    angular = SymmetricMoments(rtm[3, 3], (rtm[3, 4] + rtm[4, 3]) / 2, rtm[4, 4])

    mixed = AsymmetricMoments(
        (rtm[1, 3] + rtm[3, 1]) / 2,
        (rtm[1, 4] + rtm[4, 1]) / 2,
        (rtm[2, 3] + rtm[3, 2]) / 2,
        (rtm[2, 4] + rtm[4, 2]) / 2,
    )

    return SecondOrderBeamMoments{T}(spatial, angular, mixed)
end

function calculate_matrix(som::SecondOrderBeamMoments{T}) where {T}
    return SA{T}[
        som.spatial.xx som.spatial.xy som.mixed.xx som.mixed.xy
        som.spatial.xy som.spatial.yy som.mixed.yx som.mixed.yy
        som.mixed.xx som.mixed.yx som.angular.xx som.angular.xy
        som.mixed.xy som.mixed.yy som.angular.xy som.angular.yy
    ]
end

# struct SpatialBeamMoments{T}
#     moments::SymmetricMoments{T}
# end
# Base.getproperty(m::SpatialBeamMoments, f::Symbol) = getfield(m.moments, f)

# function SpatialBeamMoments(xx::T, xy::T, yy::T) where {T<:AbstractFloat}
#     return SpatialBeamMoments(SymmetricMoments{T}(xx, xy, yy))
# end

# function SpatialBeamMoments(xx::T, yy::T) where {T<:AbstractFloat}
#     return SpatialBeamMoments(SymmetricMoments{T}(xx, zero(T), yy))
# end
```
