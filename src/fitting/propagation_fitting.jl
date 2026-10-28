"""
    fit_propagation(λ, data; system=free_space_system, kwargs...)

Fit the input beam to beam measurements taken after propagating through free space.
This is [`fit_system`](@ref) with a fixed system and no system parameters.

`data` has one row per measurement with the columns beam radius along x, beam radius
along y, orientation angle and propagation distance. Elements between the measurements can
be included by passing `system(row, params)` with `row = [L]`. The remaining `kwargs` are
passed on to `fit_system`.
"""
function fit_propagation(λ::Real, data::AbstractMatrix; system=free_space_system, kwargs...)
    observed = moments_from_radius(data[:, 1], data[:, 2], data[:, 3])

    return fit_system(system, λ, data[:, 4:4], observed; kwargs...)
end

free_space_system(row, _) = FreeSpace(row[1])
