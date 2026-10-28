@testset "Deprecation" begin
    @testset "rtm property" begin
        fs = FreeSpace(1)
        @test (@test_deprecated fs.rtm) == get_transfer_matrix(fs)

        i = OpticalInterface(Curvature(), 1.5)
        @test (@test_deprecated i.rtm) == get_transfer_matrix(i)

        m = Mirror(OpticalPower(1.0, 1.0))
        @test (@test_deprecated m.rtm) == get_transfer_matrix(m)

        l = ThinLens(OpticalPower(1.0, 1.0))
        @test (@test_deprecated l.rtm) == get_transfer_matrix(l)

        l2 = ThickLens(Curvature(), 1.0, Curvature(), 1.5)
        @test (@test_deprecated l2.rtm) == get_transfer_matrix(l2)
    end

    @testset "Beam keyword constructor" begin
        λ = 632.8e-9
        @test (@test_deprecated Beam(λ; wx=1e-3, wy=2e-3, rx=Inf, ry=1.0, φ=0.3)) ==
              Beam(λ, BeamRadius(1e-3), BeamRadius(2e-3, 1.0); φ=0.3)
        @test (@test_deprecated Beam(λ; zx=0.5, zRx=1.0, qy=0.5 + 1.0im)) ==
              Beam(λ, RayleighRange(1.0; z=0.5), ComplexBeamParameter(0.5 + 1.0im))
        @test (@test_deprecated SecondOrderMoments(λ; wx=1e-3, wy=1e-3, rx=Inf, ry=Inf)) ==
              Beam(λ, BeamRadius(1e-3)).moments
        @test_throws ArgumentError Beam(λ; wx=1e-3, wy=1e-3)
    end

    @testset "SecondOrderMoments keywords" begin
        @test (@test_deprecated SecondOrderMoments(;
            rx=4.0, ry=1.0, φ_r=0.1, θx=2.0, θy=1.0, φ_θ=0.2, rθx=0.5, rθy=0.1, φ_rθ=0.3, twist=0.2,
        )) == SecondOrderMoments(;
            spatial=PrincipalMoments(4.0, 1.0, 0.1),
            angular=PrincipalMoments(2.0, 1.0, 0.2),
            mixed=PrincipalMoments(0.5, 0.1, 0.3),
            twist=0.2,
        )
    end

    @testset "Fitting names" begin
        @test (@test_deprecated observation([1e-3], [2e-3], [0.1])) ==
              moments_from_radius([1e-3], [2e-3], [0.1])
        @test (@test_deprecated beam_fit_parameters(1e-6; radius=(1e-3, 2e-3, 3e-3))).initial ==
              beam_fit_parameters(1e-6; beam_radius=(1e-3, 2e-3, 3e-3)).initial
    end

    @testset "Beam field inner" begin
        beam = Beam(632.8e-9, Waist(1e-3))
        @test (@test_deprecated beam.inner) === beam.moments
    end
end
