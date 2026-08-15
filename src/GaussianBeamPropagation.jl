module GaussianBeamPropagation

import Base: convert, *, getproperty, angle

using StaticArrays


include("utils.jl")

include("wigner/second_order_moments.jl")
include("wigner/beam.jl")

include("descriptors/curvature.jl")
include("descriptors/optical_power.jl")
include("descriptors/incidence_angle.jl")
include("descriptors/deprecated.jl")

include("utils2.jl")

include("wigner/elements/abstract_element.jl")
include("wigner/elements/matrix_element.jl")
include("wigner/elements/free_space.jl")
include("wigner/elements/thin_lens.jl")
include("wigner/elements/mirror.jl")
include("wigner/elements/interface.jl")
include("wigner/elements/thick_lens.jl")

include("wigner/arithmetic_operators.jl")

include("wigner/wigner.jl")

include("wigner/beam_parameters.jl")

export Beam, SecondOrderMoments
export ThinLens,
    FreeSpace, OpticalInterface, MatrixElement, AbstractElement, ThickLens, Mirror

export Curvature, OpticalPower, IncidenceAngle

export radius, divergence, phase_curvature, twist
# Deprecated:
@deprecate iso_radius radius

end # module
