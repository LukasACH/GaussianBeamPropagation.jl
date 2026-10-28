@testset "Beams and elements" begin
    @testset "Axis parameters" begin
        using GaussianBeamPropagation: axis_moments

        λ = 632.8e-9
        moments = (1.0071324798855137e-7, 5.035662399427568e-8, 5.035662399427568e-8)

        zR = 1.0
        z = 1.0
        @test all(moments .≈ axis_moments(λ, RayleighRange(zR; z)))

        w0 = sqrt(zR * λ / π)
        @test all(moments .≈ axis_moments(λ, Waist(w0; z)))

        r = z + zR^2 / z
        w = sqrt(λ / π * (zR + z^2 / zR))
        @test all(moments .≈ axis_moments(λ, BeamRadius(w, r)))

        q = 1 / (1 / r - im * λ / π / w^2)
        @test all(moments .≈ axis_moments(λ, ComplexBeamParameter(q)))

        # At the waist
        @test radius(Beam(λ, Waist(1e-3), Waist(2e-3)))[1:2] ≈ [1e-3, 2e-3]
        @test radius(Beam(λ, BeamRadius(1e-3)))[1:2] ≈ [1e-3, 1e-3]
    end

    @testset "Angles consistent" begin
        λ = 632.8e-9

        @testset "Main" begin
            wx = 2e-3
            wy = 1e-3
            φ = deg2rad(35)

            beam = Beam(λ, BeamRadius(wx), BeamRadius(wy); φ)
            calc_radius = radius(beam)
            set_M, set_m = minmax(wx, wy)
            calc_M, calc_m = minmax(calc_radius[1], calc_radius[2])
            set_φ = mod2pi(2φ + π * (wx < wy)) / 2
            calc_φ = mod2pi(2calc_radius[3] + π * (calc_radius[1] < calc_radius[2])) / 2
            @test calc_M ≈ set_M
            @test calc_m ≈ set_m
            @test calc_φ ≈ set_φ
        end

        @testset "ISO radius" begin
            @testset "Normal 1" begin
                radius_ =
                    radius(Beam(λ, BeamRadius(2e-3), BeamRadius(4e-3); φ=deg2rad(30)),)
                @test radius_ ≈ [2e-3, 4e-3, deg2rad(30)]
            end

            @testset "Normal 2" begin
                radius_ =
                    radius(Beam(λ, BeamRadius(2e-3), BeamRadius(4e-3); φ=deg2rad(-15)),)
                @test radius_ ≈ [2e-3, 4e-3, deg2rad(-15)]
            end
            @testset "Normal 3" begin
                radius_ =
                    radius(Beam(λ, BeamRadius(4e-3), BeamRadius(2e-3); φ=deg2rad(15)),)
                @test radius_ ≈ [4e-3, 2e-3, deg2rad(15)]
            end

            @testset "Normal 4" begin
                radius_ =
                    radius(Beam(λ, BeamRadius(4e-3), BeamRadius(2e-3); φ=deg2rad(-30)),)
                @test radius_ ≈ [4e-3, 2e-3, deg2rad(-30)]
            end

            @testset "45°, swapping x and y" begin
                radius_ =
                    radius(Beam(λ, BeamRadius(2e-3), BeamRadius(4e-3); φ=deg2rad(45)),)
                @test radius_ ≈ [4e-3, 2e-3, deg2rad(-45)]
            end

            @testset "-45°, swapping x and y" begin
                radius_ =
                    radius(Beam(λ, BeamRadius(2e-3), BeamRadius(4e-3); φ=deg2rad(-45)),)
                @test radius_ ≈ [4e-3, 2e-3, deg2rad(45)]
            end

            @testset "45°, not swapping x and y" begin
                radius_ =
                    radius(Beam(λ, BeamRadius(4e-3), BeamRadius(2e-3); φ=deg2rad(45)),)
                @test radius_ ≈ [4e-3, 2e-3, deg2rad(45)]
            end

            @testset "-45°, not swapping x and y" begin
                radius_ =
                    radius(Beam(λ, BeamRadius(4e-3), BeamRadius(2e-3); φ=deg2rad(-45)),)
                @test radius_ ≈ [4e-3, 2e-3, deg2rad(-45)]
            end
        end
    end

    @testset "Angle convention" begin
        # Principal axes rotated by 60°: the axis closest to x is the second one, at -30°
        ε2 = (632.8e-9 / 4π)^2
        som = SecondOrderMoments(;
            spatial=PrincipalMoments(4e-6, 1e-6, deg2rad(60)),
            angular=PrincipalMoments(ε2 / 4e-6, ε2 / 1e-6, deg2rad(60)),
        )
        beam = FreeSpace(1.0) * Beam(632.8e-9, som)

        for quantity in (radius, divergence, phase_curvature)
            @test quantity(beam)[3] ≈ deg2rad(-30)
        end
        @test radius(Beam(632.8e-9, som))[1:2] ≈ [2e-3, 4e-3]
        @test divergence(beam)[1:2] ≈ 2 .* sqrt.([ε2 / 1e-6, ε2 / 4e-6])
    end

    @testset "Radius along a direction" begin
        λ = 632.8e-9
        beam = Beam(λ, BeamRadius(2e-3), BeamRadius(1e-3); φ=deg2rad(20))

        # Along the principal axes, the radii are the principal radii
        @test radius(beam, deg2rad(20)) ≈ [2e-3, 1e-3]
        @test radius(beam, deg2rad(110)) ≈ [1e-3, 2e-3]

        # In between, the second moment (not the radius) interpolates
        α = deg2rad(30)
        @test radius(beam, deg2rad(20) + α) ≈
              [sqrt((2e-3 * cos(α))^2 + (1e-3 * sin(α))^2),
               sqrt((2e-3 * sin(α))^2 + (1e-3 * cos(α))^2)]
    end

    @testset "Insufficient parameters" begin
        @test_throws ArgumentError SecondOrderMoments(; angular=PrincipalMoments(1.0))
    end

    @testset "Mirror" begin
        @testset "Simple mirror" begin
            λ = 632.8e-9
            power = OpticalPower(1.0, 1.0; φ=deg2rad(0))
            incidence = IncidenceAngle(ι=deg2rad(45), θ=deg2rad(0))
            beam_out =
                FreeSpace(1) *
                    Mirror(power, incidence) *
                    Beam(λ, BeamRadius(1e-3), BeamRadius(1e-3); φ=deg2rad(0))

            @test radius(beam_out)[3] ≈ deg2rad(0)
        end

        @testset "Incidence angle" begin
            λ = 632.8e-9
            power = OpticalPower(1.0, 1.0; φ=deg2rad(10))
            incidence = IncidenceAngle(; ι=deg2rad(45), θ=deg2rad(-20))
            beam_out =
                FreeSpace(1) *
                    Mirror(power, incidence) *
                    Beam(λ, BeamRadius(1e-3), BeamRadius(1e-3); φ=deg2rad(0))

            @test radius(beam_out)[3] ≈ deg2rad(-20)
        end
        @testset "Toric mirror" begin
            λ = 632.8e-9
            power = OpticalPower(2.0, 1.0; φ=deg2rad(10))
            incidence = IncidenceAngle(ι=0, θ=deg2rad(-20))

            beam_out =
                FreeSpace(1) *
                    Mirror(power, incidence) *
                    Beam(λ, BeamRadius(1e-3), BeamRadius(1e-3); φ=deg2rad(0))

            @test radius(beam_out)[3] ≈ deg2rad(10)
        end
    end

    @testset "Matrix elements" begin
        @testset "Mirror" begin
            Dx = 2.0
            Dy = 1.0
            power = OpticalPower(Dx, Dy)
            mirror = Mirror(power)
            @test get_transfer_matrix(mirror) ≈ SA[
                1.0 0.0 0.0 0.0
                0.0 1.0 0.0 0.0
                -Dx 0.0 1.0 0.0
                0.0 -Dy 0.0 1.0
            ]

            # A concave mirror with R = -1 m focuses at f = R/2 = 0.5 m
            @test get_transfer_matrix(Mirror(Curvature(-1.0))) ≈
                get_transfer_matrix(Mirror(OpticalPower(2.0)))
        end

        @testset "Powers" begin
            lens = ThinLens(OpticalPower(2.0, 1.0; φ=0.3))
            @test get_transfer_matrix(lens^0) ≈ one(SMatrix{4,4,Float64})
            @test get_transfer_matrix(lens^1) ≈ get_transfer_matrix(lens)
            @test get_transfer_matrix(lens^3) ≈ get_transfer_matrix(lens * lens * lens)
            @test get_transfer_matrix(FreeSpace(1.0)^-2) ≈ get_transfer_matrix(FreeSpace(-2.0))
        end

        @testset "Integer input" begin
            @test Mirror(Curvature(1, 0, 1)) isa Mirror{Float64}
            @test ThinLens(OpticalPower(1, 0, 1)) isa ThinLens{Float64}
            @test OpticalInterface(Curvature(1, 0, 1), 2) isa OpticalInterface{Float64}
            @test MatrixElement(SMatrix{4,4}(1:16)) isa MatrixElement{Float64}
            @test FreeSpace(1) isa FreeSpace{Float64}
        end

        @testset "Free Space" begin
            L = 1.5
            free_space = FreeSpace(L)
            @test get_transfer_matrix(free_space) == SA[
                1.0 0.0 L 0.0
                0.0 1.0 0.0 L
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0
            ]
        end
    end
end
