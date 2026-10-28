# Angle utilities

!!! note "Not part of the package"

    This code is not loaded by the package. It is kept here as a starting point for future work.

The functions [`sample_beam`](@ref), [`sample_beam_continuous`](@ref) and [`sample_beam_at_angle`](@ref) are part of the package.
They follow the principal axes of a beam through a sequence of beams, for example along a scan of a parameter.
For this, they need to detect when the axes swap, i.e. when the orientation of the principal axes rotates past ±45°.

Several approaches to this detection were tried. The one that works is used in the package (`count_quarter_turns`).
The others are collected here.

## Earlier variants of the swap detection

The first version of `find_continous_angle` computed an unwrapped angle in steps of π/2 before switching to the current approach, and contained a second, unreachable approach after its `return`:

```julia
function find_continous_angle(ang::AbstractVector{T}) where {T<:AbstractFloat}
    angle_90 = @. mod2pi(4ang + π) / 4 - π / 4
    angle_90_unwrapped = unwrap(angle_90; range=π / 2)
    ints = @. Int(round((ang - angle_90_unwrapped) / (π / 2)))
    # return ints

    angle_clamped = @. mod(2ang / π + 0.5, 1) - 0.5
    angle_clamped_2 = @. mod(2ang / π + 1, 2) - 1
    angle_unwrapped = unwrap(angle_clamped; range=1)
    ints = @. Int(round(angle_clamped_2 - angle_unwrapped))
    return ints

    # The angle is between 0 and 2π. We only care about 0 and π. We want to detect if the angle changed by π/2.
    # When we take the difference between two consecutive angles, if they are offset by more than π/4, we consider this a change.
    # For that, we want to first calculate the differences, and then check if they are between π/4 and 3π/4 or 5π/4 and 7π/4.
    angle_differece = @views ang[begin:(end - 1)] - ang[(begin + 1):end]
    # Now we map everything onto 0 and π.
    angle_differece_mod = @. mod2pi(2angle_differece) / 2

    bang = @. (1π / 4 < angle_differece_mod < 3π / 4)

    flip = Vector{Int}(undef, length(ang))
    fill!(flip, 0)
    # flip[begin] = 0
    # flip_view = @view flip[begin+1:end]
    # cumsum!(flip_view, bang)
    return flip
end
```

## Outlier removal and the remaining variants

`EliminateOutlier` was an attempt to remove single outliers from a sequence of angles, by estimating the "jolt" (the fourth derivative) of the angle with a five-point stencil and skipping values where it is large.
An iterator version with a predictor based on the angle and its first three derivatives was also tried, but did not work well (see the comment in the code).
`eliminate_outliers` was started with the same idea, but is unfinished and does not run (it indexes before the start of the vector, and uses undefined variables).

`find_continous_angle1`, `find_continous_angle2` and `find_continous_angle3` detect jumps between consecutive angles by thresholding their differences, or by comparing the second derivatives of both possible orientations.

```julia
struct EliminateOutlier{T}
    func::Any
    angles::Any
    threshold::Any

    function EliminateOutlier(func, angles::AbstractArray{T}; threshold=0.1) where {T}
        return new{T}(func, angles, threshold)
    end
end
EliminateOutlier(angles; threshold=0.1) = EliminateOutlier(identity, angles; threshold)

import Base: length, eltype
Base.IteratorSize(::EliminateOutlier) = Base.HasLength()
Base.length(x::EliminateOutlier) = length(x.angles)
Base.eltype(::Type{EliminateOutlier{T}}) where {T} = Union{T,Missing}

function Base.isdone(x::EliminateOutlier, state=(; n=firstindex(x.angles)))
    return state.n > lastindex(x.angles)
end

function Base.iterate(x::EliminateOutlier, state=(; n=firstindex(x.angles)))
    if firstindex(x.angles) + 2 <= state.n <= lastindex(x.angles) - 2
        angles = x.angles[(state.n - 2):(state.n + 2)]
        factors = [1, -4, +6, -4, +1]

        values = @. exp(im * x.func(angles))

        jolt = sum(factors .* values)
        println("Jolt: $jolt ; Decision: ", ifelse(abs(jolt) > x.threshold, "skip", "keep"))
        return (ifelse(abs(jolt) > x.threshold, NaN, x.angles[state.n]), (; n=state.n + 1))
    elseif state.n > lastindex(x.angles)
        return nothing
    else
        return (x.angles[state.n], (; n=state.n + 1))
    end
end

#= This doesnt work very well
function Base.iterate(
    x::EliminateOutlier,
    state=(; n=firstindex(x.angles), θ=first(x.angles), ω=0.0, α=0.0, ι=0.0, dt=1),
)
    α_pred = muladd(state.ι, state.dt, state.α)
    ω_pred = muladd(α_pred, state.dt, state.ω)
    θ_pred = muladd(ω_pred, state.dt, state.θ)

    dt = state.dt

    if state.n > lastindex(x.angles)
        return nothing
    end

    change = abs(θ_pred - x.angles[state.n])
    if change < x.threshold
        θ = x.angles[state.n]
        ω = (θ - state.θ) / state.dt
        α = (ω - state.ω) / state.dt
        ι = (α - state.α) / state.dt

        dt = 1

        println("kept    $change")
        return (θ, (; n=state.n + 1, θ, ω, α, ι, dt))
    else
        println("skipped $change")
        return (
            missing,
            (; n=state.n + 1, θ=state.θ, ω=state.ω, α=state.α, ι=state.ι, dt=state.dt + 1),
        )
    end
end =#

function eliminate_outliers(θs)
    t = 1
    for n in eachindex(θs)
        θ = θs[n - 1]
        ω = θs[n - 1] - θs[n - 2]
        α = θs[n - 1] - 2θs[n - 2] + θs[n - 3]
        ι = θs[n - 1] - 3θs[n - 2] + 3θs[n - 3] - θs[n - 4]

        θ_pred = muladd(muladd(muladd(ι, t, α), t, ω), t, θ)

        θ_pred - θs[n]
    end

    ω = θ[(begin + 1):end] .- θ[begin:(end - 1)]
    α = ω[(begin + 1):end] .- ω[begin:(end - 1)]
    ι = α[(begin + 1):end] .- α[begin:(end - 1)]

    θ[1] - θ[0]

    return nothing
end

function find_continous_angle1(ang::AbstractVector{T}) where {T<:AbstractFloat}
    N = length(ang)

    fang = Vector{Float64}(undef, N)
    iang = BitArray(undef, N)
    bang = Vector{Bool}(undef, N)
    bang[begin] = false

    fang .= mod2pi.(ang)
    # @turbo for i in firstindex(fang):(lastindex(fang) - 1)
    for i in firstindex(fang):(lastindex(fang) - 1)
        delta = abs(fang[i] - fang[i + 1])
        a = (1π / 4) < delta
        b = delta < (3π / 4)
        c = (5π / 4) < delta
        d = delta < (7π / 4)
        ab = a & b
        cd = c & d
        abcd = ab | cd
        bang[i + 1] = abcd
    end
    accumulate!(xor, iang, bang)
    return iang
end

function find_continous_angle2(ang::AbstractVector{T}) where {T<:AbstractFloat}
    # The angle is between 0 and 2π. We only care about 0 and π. We want to detect if the angle changed by π/2.
    # When we take the difference between two consecutive angles, if they are offset by more than π/4, we consider this a change.
    # For that, we want to first calculate the differences, and then check if they are between π/4 and 3π/4 or 5π/4 and 7π/4.
    angle_differece = @views ang[begin:(end - 1)] - ang[(begin + 1):end]
    # Now we map everything onto 0 and π.
    angle_differece_mod = @. mod2pi(2angle_differece) / 2

    bang = @. (1π / 4 < angle_differece_mod < 3π / 4)

    flip = Vector{Int}(undef, length(ang))
    flip[begin] = 0
    flip_view = @view flip[(begin + 1):end]
    cumsum!(flip_view, bang)
    return flip
end

function find_continous_angle3(ang::AbstractVector{T}) where {T<:AbstractFloat}
    ang = @. mod2pi(2ang)
    angprim = @. mod2pi(2ang + π)
    toflip = Vector{Bool}(undef, length(ang))

    for i in (firstindex(toflip) + 2):(lastindex(toflip))
        toflip[i] =
            abs(angle_second_derivative(ang[i - 2], ang[i - 1], ang[i])) >
            abs(angle_second_derivative(ang[i - 2], ang[i - 1], angprim[i]))
        if toflip[i]
            @info "" angle_second_derivative(ang[i - 2], ang[i - 1], ang[i]) angle_second_derivative(
                ang[i - 2],
                ang[i - 1],
                angprim[i],
            )
        end
        ang[i], angprim[i] =
            (angprim[i], ang[i]) .* (toflip[i]) .+ (ang[i], angprim[i]) .* (!toflip[i])
    end

    #d1_ang = @views @. mod2pi(angle_first_derivative(ang[begin:end-1], ang[begin+1:end]))
    #d2_ang = @views @. mod2pi(angle_second_derivative(ang[begin:end-2], ang[begin+1:end-1], ang[begin+2:end]))
    #d2_angprim = @views @. mod2pi(angle_second_derivative(ang[begin:end-2], ang[begin+1:end-1], angprim[begin+2:end]))
    #@info d2_ang d2_angprim
    #toflip = @views @. (π / 2 < d1_ang < 3π / 2)
    flip = Vector{Int}(undef, length(ang))
    flip[begin] = 0
    flip_view = @view flip[begin:end]
    cumsum!(flip_view, toflip)
    #display(plot(x=range(1, length(d1_ang)),y=d1_ang, Geom.line()))
    return flip
end

angle_first_derivative(x0, x1; h=1) = (mod2pi(x1 - x0 + π) - π) / h

function angle_second_derivative(x0, x1, x2; h=1)
    return angle_first_derivative(
        angle_first_derivative(x2, x1; h),
        angle_first_derivative(x1, x0; h);
        h,
    )
end
```
