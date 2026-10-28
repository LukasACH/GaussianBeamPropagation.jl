"""
    FitParameters(lower, initial, upper)

Box-bounded parameter vector for a fit: element-wise `lower ≤ initial ≤ upper`.
Vectors of parameters are concatenated with `vcat`.
"""
struct FitParameters
    lower::Vector{Float64}
    initial::Vector{Float64}
    upper::Vector{Float64}

    function FitParameters(lower, initial, upper)
        length(lower) == length(initial) == length(upper) ||
            throw(DimensionMismatch("lower, initial and upper must have the same length"))
        return new(lower, initial, upper)
    end
end

FitParameters() = FitParameters(Float64[], Float64[], Float64[])

"""
    FitParameters(table)

Build from a matrix with one row per parameter and the columns `lower`, `initial`, `upper`.
"""
FitParameters(table::AbstractMatrix) = FitParameters(table[:, 1], table[:, 2], table[:, 3])

Base.length(p::FitParameters) = length(p.initial)
Base.vcat(a::FitParameters, b::FitParameters) =
    FitParameters(vcat(a.lower, b.lower), vcat(a.initial, b.initial), vcat(a.upper, b.upper))

"""
    beam_fit_parameters(λ; radius=(1e-3, 3e-3, 10e-3))

[`FitParameters`](@ref) for the ten parameters of the input beam, in the order
`rx, ry, φ_r, θx, θy, φ_θ, rθx, rθy, φ_rθ, twist`. `radius` is the `(lower, initial, upper)`
beam diameter in meters; the divergence bounds follow from the minimum emittance at `λ`.
"""
function beam_fit_parameters(λ::Real; radius=(1.0e-3, 3.0e-3, 10.0e-3))
    ε2 = (λ / 4π)^2
    w2_lb, w2, w2_up = (radius ./ 2) .^ 2

    return FitParameters(
        [w2_lb, w2_lb, -2π, ε2 / w2_up, ε2 / w2_up, -2π, -Inf, -Inf, -2π, 0.0],
        [w2, w2, 0.0, ε2 / w2, ε2 / w2, 0.0, 0.0, 0.0, 0.0, 0.0],
        [w2_up, w2_up, 2π, ε2 / w2_lb, ε2 / w2_lb, 2π, Inf, Inf, 2π, 0.0],
    )
end

"""
    observation(wx, wy, φ)

Second-order moments `(rxrx, rxry, ryry)` of measured beams with radii `wx`, `wy` along the
principal axes, rotated by `φ`.
"""
function observation(wx, wy, φ)
    rx = @. (wx / 2)^2
    ry = @. (wy / 2)^2

    xy_all = tensor_from_principal.(rx, ry, φ)
    rxrx = map(x -> x[1], xy_all)
    rxry = map(x -> x[2], xy_all)
    ryry = map(x -> x[3], xy_all)

    return (rxrx, rxry, ryry)
end

beam_from_parameters(λ, p) = Beam(
    λ,
    SecondOrderMoments(;
        spatial=PrincipalMoments(abs(p[1]), abs(p[2]), p[3]),
        angular=PrincipalMoments(p[4], p[5], p[6]),
        mixed=PrincipalMoments(p[7], p[8], p[9]),
        twist=p[10],
    ),
)

"""
    predict_moments(system, λ, vars, beam_params, system_params)

Second-order moments `(rxrx, rxry, ryry)` at the output of `system` for every row of `vars`.
"""
function predict_moments(system, λ::Real, vars::AbstractMatrix, beam_params, system_params)
    beam = beam_from_parameters(λ, beam_params)
    som = [SecondOrderMoments(system(row, system_params) * beam) for row in eachrow(vars)]

    return (getproperty.(som, :rxrx), getproperty.(som, :rxry), getproperty.(som, :ryry))
end

"""
    SystemFit

Result of [`fit_system`](@ref): the fitted input `beam`, the fitted `system_params`, the
objective value `minimum` and the raw `optim` result.
"""
struct SystemFit{B,O}
    beam::B
    system_params::Vector{Float64}
    minimum::Float64
    optim::O
end

"""
    fit_system(system, λ, vars, observed; system_params=FitParameters(), beam_params=beam_fit_parameters(λ),
               optimizer=NelderMead(), options=Optim.Options(; g_tol=1e-10, iterations=10000))

Fit the input beam and the parameters of an optical system to measured beam moments.

- `system(row, params)` returns the transfer element (anything that can be multiplied onto a
  [`Beam`](@ref)) for one measurement: `row` is a row of `vars` and `params` the current
  system parameters.
- `vars` is a matrix with one row of measurement settings per observation (e.g. a distance).
- `observed` is the tuple `(rxrx, rxry, ryry)` of measured moments, see [`observation`](@ref).
- `system_params` holds the bounds and initial values of the system parameters, and
  `beam_params` those of the input beam (see [`beam_fit_parameters`](@ref)). To fit only
  the system, or only the beam, give a fixed `FitParameters` (lower == initial == upper)
  for the other.

Returns a [`SystemFit`](@ref).
"""
function fit_system(
    system,
    λ::Real,
    vars::AbstractMatrix,
    observed::NTuple{3,AbstractVector};
    system_params::FitParameters=FitParameters(),
    beam_params::FitParameters=beam_fit_parameters(λ),
    optimizer=NelderMead(),
    options=Optim.Options(; g_tol=1.0e-10, iterations=10000),
)
    all(length.(observed) .== size(vars, 1)) ||
        throw(DimensionMismatch("observed moments and vars must have one entry per row"))

    p = vcat(beam_params, system_params)
    nb = length(beam_params)

    function objective(params)
        predicted = predict_moments(system, λ, vars, params[1:nb], params[(nb+1):end])
        return log(sum(sum((o .- m) .^ 2) for (o, m) in zip(observed, predicted)))
    end

    # @info "Fitting initial parameters" string(p.upper) string(p.initial) string(p.lower)
    # @info "Fit start"
    optim = optimize(objective, p.lower, p.upper, p.initial, optimizer, options)
    # @info "Fit complete"

    best = Optim.minimizer(optim)
    beam = beam_from_parameters(λ, best[1:nb])
    # @info "Fitted parameters" best
    # @info "Fitted Beam SOM" radius(beam)

    return SystemFit(beam, best[(nb+1):end], Optim.minimum(optim), optim)
end
