@testset "Beam sampling" begin
    using GaussianBeamPropagation: unwrap

    @testset "unwrap" begin
        @test unwrap([0.0, 3.0, -3.0, -0.5]) ≈ [0.0, 3.0, 2π - 3.0, 2π - 0.5]
        @test unwrap([0.4, 0.1, 0.4]; range=0.5) ≈ [0.4, 0.6, 0.4]
        @test unwrap(Float64[]) == Float64[]
    end

    λ = 632.8e-9
    ε2 = (λ / 4π)^2
    αs = range(0, π; length=37)
    beams = [
        Beam(λ, SecondOrderMoments(;
            spatial=PrincipalMoments(4e-6, 1e-6, α),
            angular=PrincipalMoments(ε2 / 4e-6, ε2 / 1e-6, α),
        )) for α in αs
    ]

    @testset "sample_beam" begin
        (; x, y, φx, φy) = sample_beam(beams)
        @test x ≈ [radius(b)[1] for b in beams]
        @test y ≈ [radius(b)[2] for b in beams]
        @test all(0 .<= φx .< π)
        @test mod.(φy .- φx, π) ≈ fill(π / 2, length(beams))
    end

    @testset "sample_beam_continuous" begin
        (; x, y, φx) = sample_beam_continuous(beams)
        @test x ≈ fill(4e-3, length(beams))
        @test y ≈ fill(2e-3, length(beams))
        @test φx ≈ αs
        @test sample_beam_continuous(divergence, beams).x ≈
              fill(2sqrt(ε2 / 4e-6), length(beams))
    end

    @testset "sample_beam_at_angle" begin
        (; x, y, φx, φy) = sample_beam_at_angle(beams, deg2rad(10))
        @test x ≈ [radius(b, deg2rad(10))[1] for b in beams]
        @test y ≈ [radius(b, deg2rad(10))[2] for b in beams]
        @test all(φx .≈ deg2rad(10))
        @test all(φy .≈ deg2rad(-80))
    end
end
