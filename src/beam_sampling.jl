"""
    unwrap(v::AbstractVector; range=2π)

Remove jumps of (multiples of) `range` between consecutive elements of `v`, such that the
difference between neighbours is always within `±range/2`.
"""
function unwrap(v::AbstractVector; range=2π)
    out = float.(collect(v))
    for i in Iterators.drop(eachindex(out), 1)
        out[i] = out[i - 1] + rem(out[i] - out[i - 1], range, RoundNearest)
    end
    return out
end

"""
    sample_beam([f=radius,] beams)

Evaluate `f` (e.g. [`radius`](@ref), [`divergence`](@ref) or [`phase_curvature`](@ref)) for each
beam in `beams`.

Returns a named tuple `(; x, y, φx, φy)` of vectors: the values along both principal axes and
the orientation of these axes, with `φx` ∈ [0, π) and `φy = φx + π/2` (modulo π).
"""
sample_beam(beams::AbstractVector{<:Beam}) = sample_beam(radius, beams)
function sample_beam(f, beams::AbstractVector{<:Beam})
    r = map(f, beams)
    x = getindex.(r, 1)
    y = getindex.(r, 2)
    φ = getindex.(r, 3)

    return (; x, y, φx=mod_pi.(φ), φy=mod_pi.(φ .+ π / 2))
end

"""
    sample_beam_continuous([f=radius,] beams)

Like [`sample_beam`](@ref), but follows each principal axis continuously through `beams`, for
example along a propagation distance or a scanned parameter.

The quantities along _x_ and _y_ (see [`radius`](@ref)) jump whenever the principal axes rotate
past ±45°, because _x_ is always the axis closest to the _x_ direction. Here, the axes are
swapped back whenever that happens, and the angles are unwrapped, such that `x`, `y`, `φx` and
`φy` are continuous for a continuous sequence of beams.
"""
function sample_beam_continuous(beams::AbstractVector{<:Beam})
    return sample_beam_continuous(radius, beams)
end
function sample_beam_continuous(f, beams::AbstractVector{<:Beam})
    (; x, y, φx, φy) = sample_beam(f, beams)

    swapped = isodd.(count_quarter_turns(φx))

    x_out = @. ifelse(swapped, y, x)
    y_out = @. ifelse(swapped, x, y)

    φx_out = unwrap(@.(mod_pi(ifelse(swapped, φy, φx))); range=π)
    φy_out = unwrap(@.(mod_pi(ifelse(swapped, φx, φy))); range=π)

    return (; x=x_out, y=y_out, φx=φx_out, φy=φy_out)
end

"""
    sample_beam_at_angle([f=radius,] beams, angle)

Like [`sample_beam`](@ref), but projects the values onto two fixed perpendicular directions,
at `angle` and `angle + π/2`, instead of the principal axes of each beam.
"""
function sample_beam_at_angle(beams::AbstractVector{<:Beam}, angle::Real)
    return sample_beam_at_angle(radius, beams, angle)
end
function sample_beam_at_angle(f, beams::AbstractVector{<:Beam}, angle::Real)
    (; x, y, φx) = sample_beam(f, beams)

    Δφ = φx .- angle

    x_out = @. sqrt((x * cos(Δφ))^2 + (y * sin(Δφ))^2)
    y_out = @. sqrt((x * sin(Δφ))^2 + (y * cos(Δφ))^2)

    φx_out = fill(mod_mpi2_pi2(angle), size(x_out))
    φy_out = fill(mod_mpi2_pi2(angle + π / 2), size(y_out))

    return (; x=x_out, y=y_out, φx=φx_out, φy=φy_out)
end

# Number of quarter turns by which the axis orientations in `φ` have rotated relative to a
# continuous orientation, starting from the first element. Odd numbers mean that x and y are
# swapped.
function count_quarter_turns(φ::AbstractVector{<:Real})
    turns_clamped = @. mod(2φ / π + 0.5, 1) - 0.5
    turns = @. mod(2φ / π + 1, 2) - 1
    turns_unwrapped = unwrap(turns_clamped; range=1)
    return @. Int(round(turns - turns_unwrapped))
end
