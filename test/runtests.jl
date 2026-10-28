using Test
using GaussianBeamPropagation
using Documenter
using StaticArrays
using Unitful

@testset verbose = true "GaussianBeamPropagation.jl" begin
    doctest(GaussianBeamPropagation)

    include("descriptors.jl")
    include("tensors.jl")
    include("beams.jl")
    include("regression.jl")
    include("integration.jl")
    include("sampling.jl")
    include("fitting.jl")
    include("unitful.jl")
    include("deprecation.jl")
end
