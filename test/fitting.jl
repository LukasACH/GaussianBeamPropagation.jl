@testset "Propagation fitting" begin
    # using Logging: with_logger, NullLogger
    λ = 1.03e-6
    ε2 = (λ / 4π)^2
    positions = collect(range(0.1, 2.0; length=8))

    # Beam sizes and orientation at the given positions, in the layout expected by `fit_propagation`
    function measure(beam, positions; system=z -> FreeSpace(z))
        return reduce(
            vcat, [(w=radius(system(z) * beam); [w[1] w[2] w[3] z]) for z in positions]
        )
    end

    # Spatial and angular moments for a waist of diameter `wx`×`wy`, `m2` times the diffraction limit
    function make_beam(; wx, wy, m2=3, φ_r=0.0, φ_θ=0.0, rθx=0.0, rθy=0.0, φ_rθ=0.0)
        return Beam(
            λ,
            SecondOrderMoments(;
                spatial=PrincipalMoments((wx / 2)^2, (wy / 2)^2, φ_r),
                angular=PrincipalMoments(m2 * ε2 / (wx / 2)^2, m2 * ε2 / (wy / 2)^2, φ_θ),
                mixed=PrincipalMoments(rθx, rθy, φ_rθ),
            ),
        )
    end

    fit(args...; kwargs...) =
    # with_logger(NullLogger()) do
        fit_propagation(args...; kwargs...)
    # end

    function test_fit(beam, positions)
        system(z) = FreeSpace(z)
        result = fit(λ, measure(beam, positions; system))

        # The beam sizes are reproduced at the measured positions...
        for z in positions
            expected = radius(system(z) * beam)
            fitted = radius(system(z) * result.beam)
            @test fitted[1:2] ≈ expected[1:2] rtol=1e-3
        end

        # ...and at positions that were not part of the fit
        for z in (-0.2, 0.0, 3.0)
            expected = radius(system(z) * beam)
            fitted = radius(system(z) * result.beam)
            @test fitted[1:2] ≈ expected[1:2] rtol=1e-3
        end

        return result
    end

    @testset "moments_from_radius" begin
        wx, wy, φ = [2.0e-3, 1.5e-3], [1.0e-3, 1.0e-3], [0.4, -0.2]
        xx, xy, yy = moments_from_radius(wx, wy, φ)
        for i in eachindex(wx)
            @test radius(xx[i], xy[i], yy[i]) ≈ [wx[i], wy[i], φ[i]]
        end
    end

    @testset "Round beam" begin
        result = test_fit(make_beam(wx=1.5e-3, wy=1.5e-3), positions)
        @test result.minimum < -30
    end

    @testset "Elliptical beam" begin
        test_fit(make_beam(wx=2e-3, wy=1e-3), positions)
    end

    @testset "Rotated elliptical beam" begin
        test_fit(make_beam(wx=2e-3, wy=1e-3, φ_r=0.4, φ_θ=0.4), positions)
    end

    @testset "Converging beam" begin
        test_fit(make_beam(wx=3e-3, wy=2e-3, rθx=-1e-9, rθy=-0.5e-9), positions)
    end

    @testset "General astigmatic beam" begin
        beam = make_beam(
            wx=3e-3, wy=2e-3, m2=4, φ_r=0.5, φ_θ=-0.3, rθx=-1e-9, rθy=-0.5e-9, φ_rθ=0.2
        )
        test_fit(beam, positions)
    end

    @testset "Element between beam and measurement" begin
        lens(z) = FreeSpace(z) * ThinLens(OpticalPower(1 / 0.5))
        beam = make_beam(wx=3e-3, wy=2e-3, φ_r=0.3, φ_θ=0.3)
        data = measure(beam, positions; system=lens)

        result = fit(λ, data; system=(row, _) -> lens(row[1]))

        for z in positions
            @test radius(lens(z) * result.beam)[1:2] ≈ radius(lens(z) * beam)[1:2] rtol=1e-3
        end
        @test radius(result.beam)[1:2] ≈ radius(beam)[1:2] rtol=1e-3
    end

    @testset "Result" begin
        beam = make_beam(wx=2e-3, wy=1e-3)
        result = fit(λ, measure(beam, positions))

        @test result isa SystemFit
        @test isempty(result.system_params)
        @test result.minimum < -30
    end
end
