@testset "Integration tests" begin
    # Fig. 2 of K. Dupraz et al., Opt. Commun. 443, 172 (2019). Both cylindrical lenses are tilted
    # by 10° in the same plane of incidence (θ = 0); the first lens is then rotated by 45° around
    # its surface normal, which rotates its curvature (φ) but not the plane of incidence.
    # The figure gives unlabelled radii read off the plot to about 0.01 mm, hence the atol.
    @testset "Dupraz et al. 2019" begin
        lens_thickness = 5.93e-3
        lens_curvature = 155.1e-3
        lens_refractive_index = 1.5108 # Schott N-BK7 at 800nm
        beam =
            FreeSpace(250e-3 - lens_thickness) *
                ThickLens(
                    Curvature(1/(-lens_curvature), 0; φ=deg2rad(0)),
                    lens_thickness,
                    Curvature(),
                    lens_refractive_index, # Schott N-BK7 at 800nm
                    IncidenceAngle(ι=deg2rad(10), θ=deg2rad(0)),
                ) *
                FreeSpace(100e-3 - lens_thickness) *
                ThickLens(
                    Curvature(1/(-lens_curvature), 0; φ=deg2rad(45)),
                    lens_thickness,
                    Curvature(),
                    lens_refractive_index, # Schott N-BK7 at 800nm
                    IncidenceAngle(ι=deg2rad(10), θ=deg2rad(0)),
                ) *
                FreeSpace(100e-3) *
                Beam(800e-9, BeamRadius(10e-3), BeamRadius(10e-3); φ=0.0)

        @test radius(beam)[2] ≈ 7.23e-3 rtol=0.005 atol=0.02e-3
        @test radius(beam)[1] ≈ 5.36e-3 rtol=0.005 atol=0.02e-3

        beam_p100 = FreeSpace(100e-3) * beam
        @test radius(beam_p100)[1] ≈ 10.97e-3 rtol=0.005 atol=0.02e-3
        @test radius(beam_p100)[2] ≈ 5.59e-3 rtol=0.005 atol=0.02e-3

        beam_p200 = FreeSpace(200e-3) * beam
        @test radius(beam_p200)[1] ≈ 16.04e-3 rtol=0.005 atol=0.02e-3
        @test radius(beam_p200)[2] ≈ 4.53e-3 rtol=0.005 atol=0.02e-3

        beam_m100 = FreeSpace(-100e-3) * beam
        @test radius(beam_m100)[2] ≈ 8.08e-3 rtol=0.005 atol=0.02e-3
        @test radius(beam_m100)[1] ≈ 0.55e-3 rtol=0.005 atol=0.02e-3

        beam_m200 = FreeSpace(-200e-3) * beam
        @test radius(beam_m200)[2] ≈ 9.29e-3 rtol=0.005 atol=0.02e-3
        @test radius(beam_m200)[1] ≈ 4.41e-3 rtol=0.005 atol=0.02e-3

        beam_waist = FreeSpace(-111.2e-3) * beam
        @test radius(beam_waist)[1] ≈ 0 atol=0.02e-3
    end

    # The same system with each lens replaced by a thin lens behind its reduced thickness
    @testset "Dupraz et al. 2019 simple" begin
        lens_thickness = 5.93e-3
        lens_curvature = 155.1e-3
        lens_refractive_index = 1.5108 # Schott N-BK7 at 800nm
        lens_power = (lens_refractive_index - 1) / lens_curvature
        beam =
            FreeSpace(250e-3 - lens_thickness) *
                FreeSpace(lens_thickness/lens_refractive_index) *
                ThinLens(
                    OpticalPower(lens_power, 0; φ=deg2rad(0)),
                    lens_refractive_index,
                    IncidenceAngle(; ι=deg2rad(10), θ=deg2rad(0)),
                ) *
                FreeSpace(100e-3 - lens_thickness) *
                FreeSpace(lens_thickness/lens_refractive_index) *
                ThinLens(
                    OpticalPower(lens_power, 0; φ=deg2rad(45)),
                    lens_refractive_index,
                    IncidenceAngle(; ι=deg2rad(10), θ=deg2rad(0)),
                ) *
                FreeSpace(100e-3) *
                Beam(800e-9, BeamRadius(10e-3), BeamRadius(10e-3); φ=0.0)

        @test radius(beam)[2] ≈ 7.23e-3 rtol=0.005 atol=0.02e-3
        @test radius(beam)[1] ≈ 5.36e-3 rtol=0.005 atol=0.02e-3

        beam_p100 = FreeSpace(100e-3) * beam
        @test radius(beam_p100)[1] ≈ 10.97e-3 rtol=0.005 atol=0.02e-3
        @test radius(beam_p100)[2] ≈ 5.59e-3 rtol=0.005 atol=0.02e-3

        beam_p200 = FreeSpace(200e-3) * beam
        @test radius(beam_p200)[1] ≈ 16.04e-3 rtol=0.005 atol=0.02e-3
        @test radius(beam_p200)[2] ≈ 4.53e-3 rtol=0.005 atol=0.02e-3

        beam_m100 = FreeSpace(-100e-3) * beam
        @test radius(beam_m100)[2] ≈ 8.08e-3 rtol=0.005 atol=0.02e-3
        @test radius(beam_m100)[1] ≈ 0.55e-3 rtol=0.005 atol=0.02e-3

        beam_m200 = FreeSpace(-200e-3) * beam
        @test radius(beam_m200)[2] ≈ 9.29e-3 rtol=0.005 atol=0.02e-3
        @test radius(beam_m200)[1] ≈ 4.41e-3 rtol=0.005 atol=0.02e-3

        beam_waist = FreeSpace(-111.2e-3) * beam
        @test radius(beam_waist)[1] ≈ 0 atol=0.02e-3
    end

    @testset "ThinLens" begin
        lens = ThinLens(Curvature(), Curvature(), 1.5)
        @test get_transfer_matrix(lens) ≈
            SA[1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0;]

        lens = ThinLens(Curvature(1.23), Curvature(1.23), 1.5)
        @test get_transfer_matrix(lens) ≈
            SA[1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0;]

        lens = ThinLens(Curvature(-1.23), Curvature(-1.23), 1.5)
        @test get_transfer_matrix(lens) ≈
            SA[1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0;]

        lens = ThinLens(OpticalPower())
        @test get_transfer_matrix(lens) ≈
            SA[1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0;]

        lens = ThinLens(OpticalPower(2.0, 3.0))
        @test get_transfer_matrix(lens) ≈
            SA[1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; -2.0 0.0 1.0 0.0; 0.0 -3.0 0.0 1.0;]

        roc = 1.234
        lens_a = ThinLens(OpticalPower((1.5 - 1) * (1/roc)))
        lens_b = ThinLens(Curvature(1/roc), Curvature(), 1.5)
        @test get_transfer_matrix(lens_a) ≈ get_transfer_matrix(lens_b)

        roc = 1.234
        lens_a = ThinLens(OpticalPower((2.0/1.5 - 1) * (1/roc)))
        lens_b = ThinLens(Curvature(1/roc), Curvature(), 2.0 / 1.5)
        @test get_transfer_matrix(lens_a) ≈ get_transfer_matrix(lens_b)

        lens_b = ThinLens(Curvature(), Curvature(-1/roc), 2.0 / 1.5)
        @test get_transfer_matrix(lens_a) ≈ get_transfer_matrix(lens_b)
    end

    @testset "ThickLens" begin
        unit_matrix = SA[
            1.0 0.0 0.0 0.0;
            0.0 1.0 0.0 0.0;
            0.0 0.0 1.0 0.0;
            0.0 0.0 0.0 1.0;
        ]
        n = 1.5

        lens = ThickLens(Curvature(), 0, Curvature(), n)
        @test get_transfer_matrix(lens) ≈ unit_matrix rtol = 1e-5

        lens = ThickLens(Curvature(1.23), 0, Curvature(1.23), n)
        @test get_transfer_matrix(lens) ≈ unit_matrix rtol = 1e-5

        lens = ThickLens(Curvature(-1.23), 0, Curvature(-1.23), n)
        @test get_transfer_matrix(lens) ≈ unit_matrix rtol = 1e-5

        lens = ThickLens(Curvature(), 1.23, Curvature(), n)
        dist = FreeSpace(1.23 / n)
        @test get_transfer_matrix(lens) ≈ get_transfer_matrix(dist) rtol = 1e-5
    end

    @testset "ThinLens vs Thicklens" begin
        lens_thickness = 5.93e-3
        lens_curvature = 155.1e-3
        lens_refractive_index = 1.5168
        lens_focal_length = 0.30012

        @testset "Different refractive index outside" begin
            medium_refractive_index = 1.5

            n = lens_refractive_index / medium_refractive_index

            thick_front = ThickLens(
                Curvature(-1/lens_curvature, 0; φ=0),
                lens_thickness,
                Curvature(),
                n,
            )
            thin_front =
                FreeSpace(lens_thickness/n) *
                    ThinLens(OpticalPower((n-1)/lens_curvature, 0))

            @test get_transfer_matrix(thin_front) ≈ get_transfer_matrix(thick_front) rtol =
                0.00001

            thick_back = ThickLens(
                Curvature(),
                lens_thickness,
                Curvature(1/lens_curvature, 0; φ=0),
                n,
            )
            thin_back =
                ThinLens(OpticalPower((n-1)/lens_curvature, 0.0)) *
                    FreeSpace(lens_thickness/n)

            @test get_transfer_matrix(thin_back) ≈ get_transfer_matrix(thick_back) rtol =
                0.00001
        end

        @testset "Zero indcidence angle" begin
            thick_front = ThickLens(
                Curvature(-1/lens_curvature, 0; φ=0),
                lens_thickness,
                Curvature(),
                lens_refractive_index, # Schott N-BK7 at 800nm
            )
            thin_front =
                FreeSpace(lens_thickness/lens_refractive_index) *
                    ThinLens(OpticalPower(1/lens_focal_length, 0))

            @test get_transfer_matrix(thin_front) ≈ get_transfer_matrix(thick_front) rtol =
                0.00002

            thick_back = ThickLens(
                Curvature(),
                lens_thickness,
                Curvature(1/lens_curvature, 0; φ=0),
                lens_refractive_index, # Schott N-BK7 at 800nm
            )
            thin_back =
                ThinLens(
                    OpticalPower(1/lens_focal_length, 0),
                    1.5,
                    IncidenceAngle(; ι=0, θ=0),
                ) * FreeSpace(lens_thickness/lens_refractive_index)

            @test get_transfer_matrix(thin_back) ≈ get_transfer_matrix(thick_back) rtol =
                0.00002
        end

        @testset "Non-zero indcidence angle" begin
            φ = deg2rad(0)

            ι = deg2rad(0)
            θ = deg2rad(0)

            thick_front = ThickLens(
                Curvature(-1/lens_curvature, 0; φ),
                lens_thickness,
                Curvature(),
                lens_refractive_index, # Schott N-BK7 at 800nm
                IncidenceAngle(; ι, θ),
            )
            thin_front =
                FreeSpace(lens_thickness/lens_refractive_index) * ThinLens(
                    OpticalPower(1/lens_focal_length, 0; φ),
                    lens_refractive_index, # Schott N-BK7 at 800nm
                    IncidenceAngle(; ι, θ),
                )

            @test get_transfer_matrix(thin_front) ≈ get_transfer_matrix(thick_front) rtol=0.001 atol=1e-10

            thick_back = ThickLens(
                Curvature(),
                lens_thickness,
                Curvature(1/lens_curvature, 0; φ=deg2rad(0)),
                lens_refractive_index, # Schott N-BK7 at 800nm
                IncidenceAngle(; ι, θ=deg2rad(45)),
            )
            thin_back =
                ThinLens(
                    OpticalPower(1/lens_focal_length, 0),
                    1.5,
                    IncidenceAngle(; ι, θ=deg2rad(45)),
                ) * FreeSpace(lens_thickness/lens_refractive_index)

            @test get_transfer_matrix(thin_back) ≈ get_transfer_matrix(thick_back) rtol=0.001 atol=1e-10
        end
    end
end
