using Test
using GaussianBeamPropagation
using StaticArrays

@testset verbose = true "GaussianBeamPropagation.jl" begin
    @testset "Utils" begin
        using GaussianBeamPropagation: rotate2, rotate3, rotate4, get_x_y_φ

        xx = 0.44793749401064553
        xy = 0.38922522287528594
        yx = 0.6485114409721432
        yy = 0.3645526281257251
        α = 1.6328107606377502

        @testset "rotate2" begin
            @test rotate2(1, 1, α)[1] ≈ 1.0
            @test rotate2(1, 1, α)[2] ≈ 0.0 atol=eps(Float64)
            @test rotate2(1, 1, α)[3] ≈ 1.0

            @test rotate2(xx, yy, 0)[1] ≈ xx
            @test rotate2(xx, yy, 0)[2] ≈ 0.0 atol=eps(Float64)
            @test rotate2(xx, yy, 0)[3] ≈ yy

            @test rotate2(xx, yy, π)[1] ≈ xx
            @test rotate2(xx, yy, π)[2] ≈ 0.0 atol=eps(Float64)
            @test rotate2(xx, yy, π)[3] ≈ yy

            @test rotate2(xx, yy, π/2)[1] ≈ yy
            @test rotate2(xx, yy, π/2)[2] ≈ 0.0 atol=eps(Float64)
            @test rotate2(xx, yy, π/2)[3] ≈ xx

            @test rotate2(xx, yy, π/4)[1] ≈ (xx+yy)/2
            @test rotate2(xx, yy, π/4)[2] ≈ (xx-yy)/2
            @test rotate2(xx, yy, π/4)[3] ≈ (xx+yy)/2
        end

        @testset "rotate3" begin
            @test rotate3(xx, 0.0, yy, α)[1] .≈ rotate2(xx, yy, α)[1]
            @test rotate3(xx, 0.0, yy, α)[2] .≈ rotate2(xx, yy, α)[2]
            @test rotate3(xx, 0.0, yy, α)[3] .≈ rotate2(xx, yy, α)[3]
        end

        @testset "rotate4" begin
            @test rotate4(xx, 0.0, 0.0, yy, α)[1] ≈ rotate2(xx, yy, α)[1]
            @test rotate4(xx, 0.0, 0.0, yy, α)[2] ≈ rotate2(xx, yy, α)[2]
            @test rotate4(xx, 0.0, 0.0, yy, α)[3] ≈ rotate2(xx, yy, α)[2]
            @test rotate4(xx, 0.0, 0.0, yy, α)[4] ≈ rotate2(xx, yy, α)[3]

            @test rotate4(xx, xy, xy, yy, α)[1] ≈ rotate3(xx, xy, yy, α)[1]
            @test rotate4(xx, xy, xy, yy, α)[2] ≈ rotate3(xx, xy, yy, α)[2]
            @test rotate4(xx, xy, xy, yy, α)[3] ≈ rotate3(xx, xy, yy, α)[2]
            @test rotate4(xx, xy, xy, yy, α)[4] ≈ rotate3(xx, xy, yy, α)[3]

            @test rotate4(xx, xy, yx, yy, α)[2] - rotate4(xx, xy, yx, yy, α)[3] ≈ xy - yx
        end

        @testset "get_x_y_φ" begin
            @test get_x_y_φ(rotate2(1, 2, 0))[1] ≈ 1.0
            @test get_x_y_φ(rotate2(1, 2, 0))[2] ≈ 2.0
            @test get_x_y_φ(rotate2(1, 2, 0))[3] ≈ 0.0 atol=eps(Float64)

            @test get_x_y_φ(rotate2(1, 2, π/2))[1] ≈ 2.0
            @test get_x_y_φ(rotate2(1, 2, π/2))[2] ≈ 1.0
            @test get_x_y_φ(rotate2(1, 2, π/2))[3] ≈ 0.0 atol=eps(Float64)

            @test get_x_y_φ(rotate2(1, 2, π/4))[1] ≈ 2.0
            @test get_x_y_φ(rotate2(1, 2, π/4))[2] ≈ 1.0
            @test get_x_y_φ(rotate2(1, 2, π/4))[3] ≈ -π/4

            @test get_x_y_φ(rotate2(1, 2, -π/4))[1] ≈ 2.0
            @test get_x_y_φ(rotate2(1, 2, -π/4))[2] ≈ 1.0
            @test get_x_y_φ(rotate2(1, 2, -π/4))[3] ≈ π/4

            @test get_x_y_φ(rotate2(2, 1, π/4))[1] ≈ 2.0
            @test get_x_y_φ(rotate2(2, 1, π/4))[2] ≈ 1.0
            @test get_x_y_φ(rotate2(2, 1, π/4))[3] ≈ π/4

            @test get_x_y_φ(rotate2(2, 1, -π/4))[1] ≈ 2.0
            @test get_x_y_φ(rotate2(2, 1, -π/4))[2] ≈ 1.0
            @test get_x_y_φ(rotate2(2, 1, -π/4))[3] ≈ -π/4
        end
    end

    @testset "Wigner Beam" begin
        @testset "SOM from variables" begin
            λ = 632.8e-9
            wsom = (1.0071324798855137e-7, 5.035662399427568e-8, 5.035662399427568e-8)

            zR = 1.0
            z = 1.0
            wsom_zzR = GaussianBeamPropagation.WSOM_from_z_zR(λ, z, zR)
            @test all(wsom .≈ wsom_zzR)

            r = z + zR^2 / z
            w = sqrt(λ / π * (zR + z^2 / zR))
            wsom_wr = GaussianBeamPropagation.WSOM_from_w_r(λ, w, r)
            @test all(wsom .≈ wsom_wr)

            q = 1 / (1 / r - im * λ / π / w^2)
            wsom_q = GaussianBeamPropagation.WSOM_from_Q(λ, q)
            @test all(wsom .≈ wsom_q)
        end

        @testset "Angles consistent" begin
            λ = 632.8e-9

            @testset "Main" begin
                wx = 2e-3
                wy = 1e-3
                φ = deg2rad(35)

                beam = Beam(λ; wx, wy, rx=Inf, ry=Inf, φ=φ)
                calc_radius = radius(beam)
                set_M, set_m = minmax(wx, wy)
                calc_M, calc_m = minmax(calc_radius.x, calc_radius.y)
                set_φ = mod2pi(2φ + π * (wx < wy)) / 2
                calc_φ = mod2pi(2calc_radius.φ + π * (calc_radius.x < calc_radius.y)) / 2
                @test calc_M ≈ set_M
                @test calc_m ≈ set_m
                @test calc_φ ≈ set_φ
            end

            @testset "ISO radius" begin
                @testset "Normal 1" begin
                    radius_ =
                        radius(Beam(λ; wx=2e-3, wy=4e-3, rx=Inf, ry=Inf, φ=deg2rad(30)),)
                    @test radius_.x ≈ 2e-3
                    @test radius_.y ≈ 4e-3
                    @test radius_.φ ≈ deg2rad(30)
                end

                @testset "Normal 2" begin
                    radius_ =
                        radius(Beam(λ; wx=2e-3, wy=4e-3, rx=Inf, ry=Inf, φ=deg2rad(-15)),)
                    @test radius_.x ≈ 2e-3
                    @test radius_.y ≈ 4e-3
                    @test radius_.φ ≈ deg2rad(-15)
                end
                @testset "Normal 3" begin
                    radius_ =
                        radius(Beam(λ; wx=4e-3, wy=2e-3, rx=Inf, ry=Inf, φ=deg2rad(15)),)
                    @test radius_.x ≈ 4e-3
                    @test radius_.y ≈ 2e-3
                    @test radius_.φ ≈ deg2rad(15)
                end

                @testset "Normal 4" begin
                    radius_ =
                        radius(Beam(λ; wx=4e-3, wy=2e-3, rx=Inf, ry=Inf, φ=deg2rad(-30)),)
                    @test radius_.x ≈ 4e-3
                    @test radius_.y ≈ 2e-3
                    @test radius_.φ ≈ deg2rad(-30)
                end

                @testset "45°, swapping x and y" begin
                    radius_ =
                        radius(Beam(λ; wx=2e-3, wy=4e-3, rx=Inf, ry=Inf, φ=deg2rad(45)),)
                    @test radius_.x ≈ 4e-3
                    @test radius_.y ≈ 2e-3
                    @test radius_.φ ≈ deg2rad(-45)
                end

                @testset "-45°, swapping x and y" begin
                    radius_ =
                        radius(Beam(λ; wx=2e-3, wy=4e-3, rx=Inf, ry=Inf, φ=deg2rad(-45)),)
                    @test radius_.x ≈ 4e-3
                    @test radius_.y ≈ 2e-3
                    @test radius_.φ ≈ deg2rad(45)
                end

                @testset "45°, not swapping x and y" begin
                    radius_ =
                        radius(Beam(λ; wx=4e-3, wy=2e-3, rx=Inf, ry=Inf, φ=deg2rad(45)),)
                    @test radius_.x ≈ 4e-3
                    @test radius_.y ≈ 2e-3
                    @test radius_.φ ≈ deg2rad(45)
                end

                @testset "-45°, not swapping x and y" begin
                    radius_ =
                        radius(Beam(λ; wx=4e-3, wy=2e-3, rx=Inf, ry=Inf, φ=deg2rad(-45)),)
                    @test radius_.x ≈ 4e-3
                    @test radius_.y ≈ 2e-3
                    @test radius_.φ ≈ deg2rad(-45)
                end
            end
        end

        @testset "Mirror" begin
            @testset "Simple mirror" begin
                λ = 632.8e-9
                power = OpticalPower(1.0, 1.0; φ=deg2rad(0))
                incidence = IncidenceAngle(ι=deg2rad(45), θ=deg2rad(0))
                beam_out =
                    FreeSpace(1) *
                    Mirror(power; incidence) *
                    Beam(λ; wx=1e-3, wy=1e-3, rx=Inf, ry=Inf, φ=deg2rad(0))

                @test radius(beam_out).φ ≈ deg2rad(0)
            end

            @testset "Incidence angle" begin
                λ = 632.8e-9
                power = OpticalPower(1.0, 1.0; φ=deg2rad(10))
                incidence = IncidenceAngle(; ι=deg2rad(45), θ=deg2rad(-20))
                beam_out =
                    FreeSpace(1) *
                    Mirror(power; incidence) *
                    Beam(λ; wx=1e-3, wy=1e-3, rx=Inf, ry=Inf, φ=deg2rad(0))

                @test radius(beam_out).φ ≈ deg2rad(-20)
            end
            @testset "Toric mirror" begin
                λ = 632.8e-9
                power = OpticalPower(2.0, 1.0; φ=deg2rad(10))
                incidence = IncidenceAngle(ι=0, θ=deg2rad(-20))

                beam_out =
                    FreeSpace(1) *
                    Mirror(power; incidence) *
                    Beam(λ; wx=1e-3, wy=1e-3, rx=Inf, ry=Inf, φ=deg2rad(0))

                @test radius(beam_out).φ ≈ deg2rad(10)
            end
        end

        @testset "Matrix elements" begin
            @testset "Mirror" begin
                Dx = 2.0
                Dy = 1.0
                power = OpticalPower(Dx, Dy)
                rtm = Mirror(power).rtm
                @test rtm ≈ SA[
                    1.0 0.0 0.0 0.0
                    0.0 1.0 0.0 0.0
                    -Dx 0.0 1.0 0.0
                    0.0 -Dy 0.0 1.0
                ]
            end

            @testset "Free Space" begin
                L = 1.5
                rtm = FreeSpace(L).rtm
                @test rtm == SA[
                    1.0 0.0 L 0.0
                    0.0 1.0 0.0 L
                    0.0 0.0 1.0 0.0
                    0.0 0.0 0.0 1.0
                ]
            end
        end
    end

    @testset "Matrix Regression tests" begin
        @testset "Interface" begin
            n = 1.5
            incidence = IncidenceAngle(deg2rad(10), π/3)

            @test OpticalInterface(Curvature(), n).rtm ≈ SA[
                1.0 0.0 0.0 0.0;
                0.0 1.0 0.0 0.0;
                0.0 0.0 2/3 0.0;
                0.0 0.0 0.0 2/3;
            ]

            @test OpticalInterface(Curvature(1/2.0, 1/3.0), n).rtm ≈ SA[
                1.0 0.0 0.0 0.0;
                0.0 1.0 0.0 0.0;
                1/6 0.0 2/3 0.0;
                0.0 1/9 0.0 2/3;
            ]

            @test OpticalInterface(Curvature(1/2.0, 1/3.0; φ=π/6), n).rtm ≈ SA[
                1.0 0.0 0.0 0.0;
                0.0 1.0 0.0 0.0;
                11/72 0.02405626121623441 2/3 0.0;
                0.02405626121623441 1/8 0.0 2/3;
            ]
            @test OpticalInterface(Curvature(), n, incidence).rtm ≈ SA[
                1.002149867536771 0.0037236798032304766 0.0 0.0;
                0.0037236798032305035 1.0064496026103134 0.0 0.0;
                0.0 0.0 0.6652456417038981 -0.002461287434338719;
                0.0 0.0 -0.002461287434338687 0.6624035917783612;
            ]

            @test OpticalInterface(Curvature(1/2.0, 1/3.0), n, incidence).rtm ≈ SA[
                1.002149867536771 0.0037236798032304766 0.0 0.0;
                0.0037236798032305035 1.0064496026103134 0.0 0.0;
                0.16930659506108206 0.0014593981271599752 0.6652456417038981 -0.002461287434338719;
                0.0012489993990113287 0.1141244360464875 -0.002461287434338687 0.6624035917783612;
            ]

            @test OpticalInterface(Curvature(1/2.0, 1/3.0; φ=π/6), n, incidence).rtm ≈ SA[
                1.002149867536771 0.0037236798032304766 0.0 0.0;
                0.0037236798032305035 1.0064496026103134 0.0 0.0;
                0.15543228455938418 0.026032342106687215 0.6652456417038981 -0.002461287434338719;
                0.025821943378538575 0.1286245171875129 -0.002461287434338687 0.6624035917783612;
            ]
        end
    end

    @testset "Integration tests" begin
        @testset "Dupraz et al. 2019" begin
            lens_thickness = 5.93e-3
            lens_curvature = 155.1e-3
            lens_refractive_index = 1.5108
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
                    IncidenceAngle(ι=deg2rad(10), θ=deg2rad(45)),
                ) *
                FreeSpace(100e-3) *
                Beam(800e-9; wx=10e-3, wy=10e-3, rx=Inf, ry=Inf, φ=0.0)
            @test radius(beam).y≈7.235348837209305e-3 rtol=0.02 skip=false
            @test radius(beam).x≈5.3600000000000065e-3 rtol=0.02 skip=false

            beam_p100 = FreeSpace(100e-3) * beam
            @test radius(beam_p100).x≈10.973023255813956e-3 rtol=0.02 skip=false
            @test radius(beam_p100).y≈5.5944186046511675e-3 rtol=0.02 skip=false

            beam_p200 = FreeSpace(200e-3) * beam
            @test radius(beam_p200).x≈16.039069767441863e-3 rtol=0.02 skip=false
            @test radius(beam_p200).y≈4.539534883720936e-3 rtol=0.02 skip=true

            beam_m100 = FreeSpace(-100e-3) * beam
            @test radius(beam_m100).y≈8.081860465116284e-3 rtol=0.02 skip=false
            @test radius(beam_m100).x≈0.5544186046511648e-3 rtol=0.02 skip=true

            beam_m200 = FreeSpace(-200e-3) * beam
            @test radius(beam_m200).y≈9.293023255813958e-3 rtol=0.02 skip=false
            @test radius(beam_m200).x≈4.4093023255814e-3 rtol=0.02 skip=false

            beam_waist = FreeSpace(-112e-3) * beam
            @test radius(beam_waist).x≈0.0 atol=0.000001 skip=true

            # @info "" rms=sqrt(
            #     mapreduce(
            #         x->x^2,
            #         +,
            #         [
            #             radius(beam).x - 7.235348837209305e-3,
            #             radius(beam).y - 5.3600000000000065e-3,
            #             radius(beam_p100).x - 10.973023255813956e-3,
            #             radius(beam_p100).y - 5.5944186046511675e-3,
            #             radius(beam_p200).x - 16.039069767441863e-3,
            #             radius(beam_p200).y - 4.539534883720936e-3,
            #             radius(beam_m100).x - 8.081860465116284e-3,
            #             radius(beam_m100).y - 0.5544186046511648e-3,
            #             radius(beam_m200).x - 9.293023255813958e-3,
            #             radius(beam_m200).y - 4.4093023255814e-3,
            #         ],
            #     ),
            # )
        end

        @testset "Dupraz et al. 2019 simple" begin
            lens_thickness = 5.93e-3
            lens_curvature = 155.1e-3
            lens_refractive_index = 1.5108
            beam =
                FreeSpace(250e-3 - lens_thickness) *
                FreeSpace(lens_thickness/1.5168) *
                ThinLens(
                    OpticalPower(1/0.30012, 0; φ=deg2rad(0));
                    n=lens_refractive_index,
                    incidence=IncidenceAngle(; ι=deg2rad(10), θ=deg2rad(0)),
                ) *
                FreeSpace(100e-3 - lens_thickness) *
                FreeSpace(lens_thickness/1.5168) *
                ThinLens(
                    OpticalPower(1/0.30012, 0; φ=deg2rad(45));
                    n=lens_refractive_index,
                    incidence=IncidenceAngle(; ι=deg2rad(10), θ=deg2rad(45)),
                ) *
                FreeSpace(100e-3) *
                Beam(800e-9; wx=10e-3, wy=10e-3, rx=Inf, ry=Inf, φ=0.0)
            @test radius(beam).x≈7.235348837209305e-3 rtol=0.01 skip=true
            @test radius(beam).y≈5.3600000000000065e-3 rtol=0.01 skip=true

            beam_p100 = FreeSpace(100e-3) * beam
            @test radius(beam_p100).x≈10.973023255813956e-3 rtol=0.01 skip=true
            @test radius(beam_p100).y≈5.5944186046511675e-3 rtol=0.01 skip=true

            beam_p200 = FreeSpace(200e-3) * beam
            @test radius(beam_p200).x≈16.039069767441863e-3 rtol=0.01 skip=true
            @test radius(beam_p200).y≈4.539534883720936e-3 rtol=0.01 skip=true

            beam_m100 = FreeSpace(-100e-3) * beam
            @test radius(beam_m100).x≈8.081860465116284e-3 rtol=0.01 skip=true
            @test radius(beam_m100).y≈0.5544186046511648e-3 rtol=0.01 skip=true

            beam_m200 = FreeSpace(-200e-3) * beam
            @test radius(beam_m200).x≈9.293023255813958e-3 rtol=0.01 skip=true
            @test radius(beam_m200).y≈4.4093023255814e-3 rtol=0.01 skip=true

            # beam_waist = FreeSpace(-109e-3) * beam
            # @test radius(beam_waist).y ≈ 0.0 atol=0.000001

            # @info "" rms=sqrt(
            #     mapreduce(
            #         x->x^2,
            #         +,
            #         [
            #             radius(beam).x - 7.235348837209305e-3,
            #             radius(beam).y - 5.3600000000000065e-3,
            #             radius(beam_p100).x - 10.973023255813956e-3,
            #             radius(beam_p100).y - 5.5944186046511675e-3,
            #             radius(beam_p200).x - 16.039069767441863e-3,
            #             radius(beam_p200).y - 4.539534883720936e-3,
            #             radius(beam_m100).x - 8.081860465116284e-3,
            #             radius(beam_m100).y - 0.5544186046511648e-3,
            #             radius(beam_m200).x - 9.293023255813958e-3,
            #             radius(beam_m200).y - 4.4093023255814e-3,
            #         ],
            #     ),
            # )
        end

        @testset "ThinLens" begin
            lens = ThinLens(Curvature(), Curvature()).rtm
            @test lens ≈
                  SA[1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0;]

            lens = ThinLens(Curvature(1.23), Curvature(1.23)).rtm
            @test lens ≈
                  SA[1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0;]

            lens = ThinLens(Curvature(-1.23), Curvature(-1.23)).rtm
            @test lens ≈
                  SA[1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0;]

            lens = ThinLens(OpticalPower()).rtm
            @test lens ≈
                  SA[1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0;]

            lens = ThinLens(OpticalPower(2.0, 3.0)).rtm
            @test lens ≈
                  SA[1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; -2.0 0.0 1.0 0.0; 0.0 -3.0 0.0 1.0;]

            roc = 1.234
            lens_a = ThinLens(OpticalPower((1.5 - 1) * (1/roc))).rtm
            lens_b = ThinLens(Curvature(1/roc), Curvature()).rtm
            @test lens_a ≈ lens_b

            roc = 1.234
            lens_a = ThinLens(OpticalPower((2.0/1.5 - 1) * (1/roc)); n=2.0 / 1.5).rtm
            lens_b = ThinLens(Curvature(1/roc), Curvature(); n=2.0 / 1.5).rtm
            @test lens_a ≈ lens_b

            lens_b = ThinLens(Curvature(), Curvature(-1/roc); n=2.0 / 1.5).rtm
            @test lens_a ≈ lens_b
        end

        @testset "ThickLens" begin
            unit_matrix = SA[
                1.0 0.0 0.0 0.0;
                0.0 1.0 0.0 0.0;
                0.0 0.0 1.0 0.0;
                0.0 0.0 0.0 1.0;
            ]
            n = 1.5

            lens = ThickLens(Curvature(), 0, Curvature(), n).rtm
            @test lens ≈ unit_matrix rtol = 1e-5

            lens = ThickLens(Curvature(1.23), 0, Curvature(1.23), n).rtm
            @test lens ≈ unit_matrix rtol = 1e-5

            lens = ThickLens(Curvature(-1.23), 0, Curvature(-1.23), n).rtm
            @test lens ≈ unit_matrix rtol = 1e-5

            lens = ThickLens(Curvature(), 1.23, Curvature(), n).rtm
            dist = FreeSpace(1.23 / n).rtm
            @test lens ≈ dist rtol = 1e-5
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
                    ThinLens(OpticalPower((n-1)/lens_curvature, 0); n)

                @test thin_front.rtm ≈ thick_front.rtm rtol = 0.00001

                thick_back = ThickLens(
                    Curvature(),
                    lens_thickness,
                    Curvature(1/lens_curvature, 0; φ=0),
                    n,
                )
                thin_back =
                    ThinLens(OpticalPower((n-1)/lens_curvature, 0.0); n) *
                    FreeSpace(lens_thickness/n)

                @test thin_back.rtm ≈ thick_back.rtm rtol = 0.00001
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

                @test thin_front.rtm ≈ thick_front.rtm rtol = 0.00002

                thick_back = ThickLens(
                    Curvature(),
                    lens_thickness,
                    Curvature(1/lens_curvature, 0; φ=0),
                    lens_refractive_index, # Schott N-BK7 at 800nm
                )
                thin_back =
                    ThinLens(
                        OpticalPower(1/lens_focal_length, 0);
                        incidence=IncidenceAngle(; ι=0, θ=0),
                    ) * FreeSpace(lens_thickness/lens_refractive_index)

                @test thin_back.rtm ≈ thick_back.rtm rtol = 0.00002
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
                        OpticalPower(1/lens_focal_length, 0; φ);
                        n=lens_refractive_index, # Schott N-BK7 at 800nm
                        incidence=IncidenceAngle(; ι, θ),
                    )

                @test thin_front.rtm≈thick_front.rtm rtol=0.001 atol=1e-10

                thick_back = ThickLens(
                    Curvature(),
                    lens_thickness,
                    Curvature(1/lens_curvature, 0; φ=deg2rad(0)),
                    lens_refractive_index, # Schott N-BK7 at 800nm
                    IncidenceAngle(; ι, θ=deg2rad(45)),
                )
                thin_back =
                    ThinLens(
                        OpticalPower(1/lens_focal_length, 0);
                        incidence=IncidenceAngle(; ι, θ=deg2rad(45)),
                    ) * FreeSpace(lens_thickness/lens_refractive_index)

                @test thin_back.rtm≈thick_back.rtm rtol=0.001 atol=1e-10
            end
        end
    end
end
#=
@testset "2D matrices" begin
    @test GaussianBeamPropagation.RTM(GaussianBeamPropagation.ThinLens2D(1)) ≈
        GaussianBeamPropagation.RTM(GaussianBeamPropagation.ThinLens2D(Inf, 1) * GaussianBeamPropagation.ThinLens2D(1, Inf))

    @test GaussianBeamPropagation.RTM(GaussianBeamPropagation.ThinLens2D(Inf)) ≈
        GaussianBeamPropagation.RTM(GaussianBeamPropagation.ThinLens2D(1) * GaussianBeamPropagation.ThinLens2D(-1))

    @test GaussianBeamPropagation.RTM(GaussianBeamPropagation.ThinLens2D(1, 2)) ≈ GaussianBeamPropagation.RTM(GaussianBeamPropagation.ThinLens2D(2, 1; β=π / 2))
    @test GaussianBeamPropagation.RTM(GaussianBeamPropagation.ThinLens2D(1, 2; α=π / 12)) ≈
        GaussianBeamPropagation.RTM(GaussianBeamPropagation.ThinLens2D(2, 1; α=π / 12, β=π / 2))

    @test GaussianBeamPropagation.RTM(GaussianBeamPropagation.ThinLens2D(1, 2)) ≈ GaussianBeamPropagation.RTM(GaussianBeamPropagation.ThinLens2D(2, 1; β=π / 2))
    @test GaussianBeamPropagation.RTM(GaussianBeamPropagation.ThinLens2D(1, 2; α=π / 12)) ≈
        GaussianBeamPropagation.RTM(GaussianBeamPropagation.ThinLens2D(2, 1; α=π / 12, β=π / 2))

    @test GaussianBeamPropagation.RTM(GaussianBeamPropagation.ThinLens2D(1, 2)) ≈ GaussianBeamPropagation.RTM(GaussianBeamPropagation.ThinLens2D(2, 1; β=π / 2))
    @test GaussianBeamPropagation.RTM(GaussianBeamPropagation.ThinLens2D(1, 2; α=π / 12)) ≈
        GaussianBeamPropagation.RTM(GaussianBeamPropagation.ThinLens2D(2, 1; α=π / 12, β=π / 2))
end
=#