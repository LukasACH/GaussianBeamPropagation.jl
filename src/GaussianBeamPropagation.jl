module GaussianBeamPropagation

using StaticArrays, Optim

include("tensors.jl")

include("wigner/second_order_moments.jl")
include("wigner/axis_parameters.jl")
include("wigner/beam.jl")

include("descriptors/symmetric_descriptor.jl")
include("descriptors/curvature.jl")
include("descriptors/optical_power.jl")
include("descriptors/incidence_angle.jl")
include("descriptors/deprecated.jl")

include("wigner/elements/abstract_element.jl")
include("wigner/elements/matrix_element.jl")
include("wigner/elements/free_space.jl")
include("wigner/elements/thin_lens.jl")
include("wigner/elements/mirror.jl")
include("wigner/elements/interface.jl")
include("wigner/elements/thick_lens.jl")

include("wigner/arithmetic_operators.jl")

include("wigner/beam_parameters.jl")
include("beam_sampling.jl")

include("fitting/system_fitting.jl")
include("fitting/propagation_fitting.jl")

export Beam, SecondOrderMoments, PrincipalMoments
export AxisParameters, BeamRadius, Waist, RayleighRange, ComplexBeamParameter
export ThinLens,
    FreeSpace, OpticalInterface, MatrixElement, AbstractElement, ThickLens, Mirror

export get_transfer_matrix

export Curvature, OpticalPower, IncidenceAngle

export radius, divergence, phase_curvature, twist
export sample_beam, sample_beam_continuous, sample_beam_at_angle

export fit_system, fit_propagation, FitParameters, beam_fit_parameters, moments_from_radius, SystemFit
# Deprecated:
@deprecate iso_radius radius
@deprecate propagation_fitting fit_propagation
@deprecate observation moments_from_radius

end # module
