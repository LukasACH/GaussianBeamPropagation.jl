using Test
using GaussianBeamPropagation
using Documenter
using StaticArrays

@testset verbose = true "GaussianBeamPropagation.jl" begin
    doctest(GaussianBeamPropagation)

    @testset "Descriptors" begin
        @testset "Curvature" begin
            @testset "Constructors" begin
                κ1 = Curvature(1.0)
                @test κ1.xx == 1.0
                @test κ1.xy == 0.0
                @test κ1.yy == 1.0

                κ2 = Curvature(1.0, 2.0, 3.0)
                @test κ2.xx == 1.0
                @test κ2.xy == 2.0
                @test κ2.yy == 3.0

                κ3 = Curvature(κ2)
                @test κ3.xx == 1.0
                @test κ3.xy == 2.0
                @test κ3.yy == 3.0

                κ4 = Curvature(1.0, 2.0; φ=deg2rad(45))
                @test κ4.xx ≈ 1.5
                @test κ4.xy ≈ -0.5
                @test κ4.yy ≈ 1.5
            end

            @testset "Type conversion and promotion" begin
                @test Curvature{Float64} == typeof(Curvature(1.0))
                @test Curvature{Float32} == typeof(Curvature(1.0f0))
                @test Curvature{Float64} == typeof(Curvature(1))
                @test Curvature{Float32} == typeof(Curvature(1, 2.0f0))
                @test Curvature{Float32} == typeof(Curvature(1, 2.0f0, Float16(3)))

                @test Curvature{Float32} == typeof(Curvature(1f0, 2f0; φ=Int(0)))
                @test Curvature{Float32} == typeof(Curvature(1f0, 2f0; φ=Float64(0)))
                @test Curvature{Float64} == typeof(Curvature(1, 2; φ=Float64(0)))

                @test Curvature{Float32} == typeof(Curvature{Float32}(Curvature(1.0)))
                @test Curvature{Float32} == typeof(convert(Curvature{Float32}, Curvature(1.0)))
            end

            @testset "Arithmetic" begin
                κ1 = Curvature(1.0, 2.0, 3.0)
                κ2 = Curvature(4.0, 5.0, 6.0)

                @test κ1 + κ2 == Curvature(5.0, 7.0, 9.0)
                @test κ1 - κ2 == Curvature(-3.0, -3.0, -3.0)
                @test 2 * κ1 == Curvature(2.0, 4.0, 6.0)
                @test κ1 * 2 == Curvature(2.0, 4.0, 6.0)
                @test κ1 / 2 == Curvature(0.5, 1.0, 1.5)
                @test -κ1 == Curvature(-1.0, -2.0, -3.0)
            end

            @testset "Show" begin
                κ = Curvature(1.0, 2.0; φ=deg2rad(30))
                io = IOBuffer()
                show(io, MIME"text/plain"(), κ)
                output = String(take!(io))
                @test occursin("Curvature{Float64}", output)
                @test occursin("κx = 1.0", output)
                @test occursin("κy = 2.0", output)
            end
        end

        @testset "Optical Power" begin
            @testset "Constructors" begin
                P1 = OpticalPower(1.0)
                @test P1.xx == 1.0
                @test P1.xy == 0.0
                @test P1.yy == 1.0

                P2 = OpticalPower(1.0, 2.0, 3.0)
                @test P2.xx == 1.0
                @test P2.xy == 2.0
                @test P2.yy == 3.0

                P3 = OpticalPower(P2)
                @test P3.xx == 1.0
                @test P3.xy == 2.0
                @test P3.yy == 3.0

                P4 = OpticalPower(1.0, 2.0; φ=deg2rad(45))
                @test P4.xx ≈ 1.5
                @test P4.xy ≈ -0.5
                @test P4.yy ≈ 1.5
            end

            @testset "Type conversion and promotion" begin
                @test OpticalPower{Float64} == typeof(OpticalPower(1.0))
                @test OpticalPower{Float32} == typeof(OpticalPower(1.0f0))
                @test OpticalPower{Float64} == typeof(OpticalPower(1))
                @test OpticalPower{Float32} == typeof(OpticalPower(1, 2.0f0))
                @test OpticalPower{Float32} == typeof(OpticalPower(1, 2.0f0, Float16(3)))

                @test OpticalPower{Float32} == typeof(OpticalPower(1f0, 2f0; φ=Int(0)))
                @test OpticalPower{Float32} == typeof(OpticalPower(1f0, 2f0; φ=Float64(0)))
                @test OpticalPower{Float64} == typeof(OpticalPower(1, 2; φ=Float64(0)))

                @test OpticalPower{Float32} == typeof(OpticalPower{Float32}(OpticalPower(1.0)))
                @test OpticalPower{Float32} == typeof(convert(OpticalPower{Float32}, OpticalPower(1.0)))
            end

            @testset "Arithmetic" begin
                P1 = OpticalPower(1.0, 2.0, 3.0)
                P2 = OpticalPower(4.0, 5.0, 6.0)

                @test P1 + P2 == OpticalPower(5.0, 7.0, 9.0)
                @test P1 - P2 == OpticalPower(-3.0, -3.0, -3.0)
                @test 2 * P1 == OpticalPower(2.0, 4.0, 6.0)
                @test P1 * 2 == OpticalPower(2.0, 4.0, 6.0)
                @test P1 / 2 == OpticalPower(0.5, 1.0, 1.5)
                @test -P1 == OpticalPower(-1.0, -2.0, -3.0)
            end

            @testset "Show" begin
                P = OpticalPower(1.0, 2.0; φ=deg2rad(30))
                io = IOBuffer()
                show(io, MIME"text/plain"(), P)
                output = String(take!(io))
                @test occursin("OpticalPower{Float64}", output)
                @test occursin("Px = 1.0", output)
                @test occursin("Py = 2.0", output)
            end
        end
    end

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
                        radius(Beam(λ; wx=2e-3, wy=4e-3, rx=Inf, ry=Inf, φ=deg2rad(30)),)
                    @test radius_ ≈ [2e-3, 4e-3, deg2rad(30)]
                end

                @testset "Normal 2" begin
                    radius_ =
                        radius(Beam(λ; wx=2e-3, wy=4e-3, rx=Inf, ry=Inf, φ=deg2rad(-15)),)
                    @test radius_ ≈ [2e-3, 4e-3, deg2rad(-15)]
                end
                @testset "Normal 3" begin
                    radius_ =
                        radius(Beam(λ; wx=4e-3, wy=2e-3, rx=Inf, ry=Inf, φ=deg2rad(15)),)
                    @test radius_ ≈ [4e-3, 2e-3, deg2rad(15)]
                end

                @testset "Normal 4" begin
                    radius_ =
                        radius(Beam(λ; wx=4e-3, wy=2e-3, rx=Inf, ry=Inf, φ=deg2rad(-30)),)
                    @test radius_ ≈ [4e-3, 2e-3, deg2rad(-30)]
                end

                @testset "45°, swapping x and y" begin
                    radius_ =
                        radius(Beam(λ; wx=2e-3, wy=4e-3, rx=Inf, ry=Inf, φ=deg2rad(45)),)
                    @test radius_ ≈ [4e-3, 2e-3, deg2rad(-45)]
                end

                @testset "-45°, swapping x and y" begin
                    radius_ =
                        radius(Beam(λ; wx=2e-3, wy=4e-3, rx=Inf, ry=Inf, φ=deg2rad(-45)),)
                    @test radius_ ≈ [4e-3, 2e-3, deg2rad(45)]
                end

                @testset "45°, not swapping x and y" begin
                    radius_ =
                        radius(Beam(λ; wx=4e-3, wy=2e-3, rx=Inf, ry=Inf, φ=deg2rad(45)),)
                    @test radius_ ≈ [4e-3, 2e-3, deg2rad(45)]
                end

                @testset "-45°, not swapping x and y" begin
                    radius_ =
                        radius(Beam(λ; wx=4e-3, wy=2e-3, rx=Inf, ry=Inf, φ=deg2rad(-45)),)
                    @test radius_ ≈ [4e-3, 2e-3, deg2rad(-45)]
                end
            end
        end

        @testset "Radius along a direction" begin
            λ = 632.8e-9
            beam = Beam(λ; wx=2e-3, wy=1e-3, rx=Inf, ry=Inf, φ=deg2rad(20))

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
            @test_throws ArgumentError Beam(632.8e-9; wx=1e-3, wy=1e-3)
        end

        @testset "Mirror" begin
            @testset "Simple mirror" begin
                λ = 632.8e-9
                power = OpticalPower(1.0, 1.0; φ=deg2rad(0))
                incidence = IncidenceAngle(ι=deg2rad(45), θ=deg2rad(0))
                beam_out =
                    FreeSpace(1) *
                        Mirror(power, incidence) *
                        Beam(λ; wx=1e-3, wy=1e-3, rx=Inf, ry=Inf, φ=deg2rad(0))

                @test radius(beam_out)[3] ≈ deg2rad(0)
            end

            @testset "Incidence angle" begin
                λ = 632.8e-9
                power = OpticalPower(1.0, 1.0; φ=deg2rad(10))
                incidence = IncidenceAngle(; ι=deg2rad(45), θ=deg2rad(-20))
                beam_out =
                    FreeSpace(1) *
                        Mirror(power, incidence) *
                        Beam(λ; wx=1e-3, wy=1e-3, rx=Inf, ry=Inf, φ=deg2rad(0))

                @test radius(beam_out)[3] ≈ deg2rad(-20)
            end
            @testset "Toric mirror" begin
                λ = 632.8e-9
                power = OpticalPower(2.0, 1.0; φ=deg2rad(10))
                incidence = IncidenceAngle(ι=0, θ=deg2rad(-20))

                beam_out =
                    FreeSpace(1) *
                        Mirror(power, incidence) *
                        Beam(λ; wx=1e-3, wy=1e-3, rx=Inf, ry=Inf, φ=deg2rad(0))

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

    @testset "Matrix Regression tests" begin
        @testset "Interface" begin
            n = 1.5
            incidence = IncidenceAngle(deg2rad(10), π/3)

            @test get_transfer_matrix(OpticalInterface(Curvature(), n)) ≈ SA[
                1.0 0.0 0.0 0.0;
                0.0 1.0 0.0 0.0;
                0.0 0.0 2/3 0.0;
                0.0 0.0 0.0 2/3;
            ]

            @test get_transfer_matrix(OpticalInterface(Curvature(1/2.0, 1/3.0), n)) ≈ SA[
                1.0 0.0 0.0 0.0;
                0.0 1.0 0.0 0.0;
                1/6 0.0 2/3 0.0;
                0.0 1/9 0.0 2/3;
            ]

            @test get_transfer_matrix(OpticalInterface(Curvature(1/2.0, 1/3.0; φ=π/6), n)) ≈
                SA[
                1.0 0.0 0.0 0.0;
                0.0 1.0 0.0 0.0;
                11/72 0.02405626121623441 2/3 0.0;
                0.02405626121623441 1/8 0.0 2/3;
            ]
            @test get_transfer_matrix(OpticalInterface(Curvature(), n, incidence)) ≈ SA[
                1.002149867536771 0.0037236798032304766 0.0 0.0;
                0.0037236798032305035 1.0064496026103134 0.0 0.0;
                0.0 0.0 0.6652456417038981 -0.002461287434338719;
                0.0 0.0 -0.002461287434338687 0.6624035917783612;
            ]

            @test get_transfer_matrix(OpticalInterface(
                Curvature(1/2.0, 1/3.0),
                n,
                incidence,
            )) ≈ SA[
                1.002149867536771 0.0037236798032304766 0.0 0.0;
                0.0037236798032305035 1.0064496026103134 0.0 0.0;
                0.16930659506108206 0.0014593981271599752 0.6652456417038981 -0.002461287434338719;
                0.0012489993990113287 0.1141244360464875 -0.002461287434338687 0.6624035917783612;
            ]

            @test get_transfer_matrix(OpticalInterface(
                Curvature(1/2.0, 1/3.0; φ=π/6),
                n,
                incidence,
            )) ≈ SA[
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
            @test radius(beam)[2] ≈ 7.235348837209305e-3 rtol=0.02 skip=false
            @test radius(beam)[1] ≈ 5.3600000000000065e-3 rtol=0.02 skip=false

            beam_p100 = FreeSpace(100e-3) * beam
            @test radius(beam_p100)[1] ≈ 10.973023255813956e-3 rtol=0.02 skip=false
            @test radius(beam_p100)[2] ≈ 5.5944186046511675e-3 rtol=0.02 skip=false

            beam_p200 = FreeSpace(200e-3) * beam
            @test radius(beam_p200)[1] ≈ 16.039069767441863e-3 rtol=0.02 skip=false
            @test radius(beam_p200)[2] ≈ 4.539534883720936e-3 rtol=0.02 skip=true

            beam_m100 = FreeSpace(-100e-3) * beam
            @test radius(beam_m100)[2] ≈ 8.081860465116284e-3 rtol=0.02 skip=false
            @test radius(beam_m100)[1] ≈ 0.5544186046511648e-3 rtol=0.02 skip=true

            beam_m200 = FreeSpace(-200e-3) * beam
            @test radius(beam_m200)[2] ≈ 9.293023255813958e-3 rtol=0.02 skip=false
            @test radius(beam_m200)[1] ≈ 4.4093023255814e-3 rtol=0.02 skip=false

            beam_waist = FreeSpace(-112e-3) * beam
            @test radius(beam_waist)[1] ≈ 0.0 atol=0.000001 skip=true

            # @info "" rms=sqrt(
            #     mapreduce(
            #         x->x^2,
            #         +,
            #         [
            #             radius(beam)[1] - 7.235348837209305e-3,
            #             radius(beam)[2] - 5.3600000000000065e-3,
            #             radius(beam_p100)[1] - 10.973023255813956e-3,
            #             radius(beam_p100)[2] - 5.5944186046511675e-3,
            #             radius(beam_p200)[1] - 16.039069767441863e-3,
            #             radius(beam_p200)[2] - 4.539534883720936e-3,
            #             radius(beam_m100)[1] - 8.081860465116284e-3,
            #             radius(beam_m100)[2] - 0.5544186046511648e-3,
            #             radius(beam_m200)[1] - 9.293023255813958e-3,
            #             radius(beam_m200)[2] - 4.4093023255814e-3,
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
                        OpticalPower(1/0.30012, 0; φ=deg2rad(0)),
                        lens_refractive_index,
                        IncidenceAngle(; ι=deg2rad(10), θ=deg2rad(0)),
                    ) *
                    FreeSpace(100e-3 - lens_thickness) *
                    FreeSpace(lens_thickness/1.5168) *
                    ThinLens(
                        OpticalPower(1/0.30012, 0; φ=deg2rad(45)),
                        lens_refractive_index,
                        IncidenceAngle(; ι=deg2rad(10), θ=deg2rad(45)),
                    ) *
                    FreeSpace(100e-3) *
                    Beam(800e-9; wx=10e-3, wy=10e-3, rx=Inf, ry=Inf, φ=0.0)
            @test radius(beam)[1]≈7.235348837209305e-3 rtol=0.01 skip=true
            @test radius(beam)[2]≈5.3600000000000065e-3 rtol=0.01 skip=true

            beam_p100 = FreeSpace(100e-3) * beam
            @test radius(beam_p100)[1]≈10.973023255813956e-3 rtol=0.01 skip=true
            @test radius(beam_p100)[2]≈5.5944186046511675e-3 rtol=0.01 skip=true

            beam_p200 = FreeSpace(200e-3) * beam
            @test radius(beam_p200)[1]≈16.039069767441863e-3 rtol=0.01 skip=true
            @test radius(beam_p200)[2]≈4.539534883720936e-3 rtol=0.01 skip=true

            beam_m100 = FreeSpace(-100e-3) * beam
            @test radius(beam_m100)[1]≈8.081860465116284e-3 rtol=0.01 skip=true
            @test radius(beam_m100)[2]≈0.5544186046511648e-3 rtol=0.01 skip=true

            beam_m200 = FreeSpace(-200e-3) * beam
            @test radius(beam_m200)[1]≈9.293023255813958e-3 rtol=0.01 skip=true
            @test radius(beam_m200)[2]≈4.4093023255814e-3 rtol=0.01 skip=true

            # beam_waist = FreeSpace(-109e-3) * beam
            # @test radius(beam_waist)[2] ≈ 0.0 atol=0.000001

            # @info "" rms=sqrt(
            #     mapreduce(
            #         x->x^2,
            #         +,
            #         [
            #             radius(beam)[1] - 7.235348837209305e-3,
            #             radius(beam)[2] - 5.3600000000000065e-3,
            #             radius(beam_p100)[1] - 10.973023255813956e-3,
            #             radius(beam_p100)[2] - 5.5944186046511675e-3,
            #             radius(beam_p200)[1] - 16.039069767441863e-3,
            #             radius(beam_p200)[2] - 4.539534883720936e-3,
            #             radius(beam_m100)[1] - 8.081860465116284e-3,
            #             radius(beam_m100)[2] - 0.5544186046511648e-3,
            #             radius(beam_m200)[1] - 9.293023255813958e-3,
            #             radius(beam_m200)[2] - 4.4093023255814e-3,
            #         ],
            #     ),
            # )
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

    @testset "Propagation fitting" begin
        # using Logging: with_logger, NullLogger
        λ = 1.03e-6
        ε2 = (λ / 4π)^2
        positions = collect(range(0.1, 2.0; length=8))

        # Beam sizes and orientation at the given positions, in the layout expected by `propagation_fitting`
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
                    rx=(wx / 2)^2,
                    ry=(wy / 2)^2,
                    φ_r,
                    θx=m2 * ε2 / (wx / 2)^2,
                    θy=m2 * ε2 / (wy / 2)^2,
                    φ_θ,
                    rθx,
                    rθy,
                    φ_rθ,
                    twist=0.0,
                ),
            )
        end

        fit(args...; kwargs...) =
        # with_logger(NullLogger()) do
            propagation_fitting(args...; kwargs...)
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

        @testset "observation" begin
            wx, wy, φ = [2.0e-3, 1.5e-3], [1.0e-3, 1.0e-3], [0.4, -0.2]
            xx, xy, yy = observation(wx, wy, φ)
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