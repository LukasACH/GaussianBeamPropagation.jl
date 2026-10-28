using Test
using GaussianBeamPropagation
using Documenter
using StaticArrays
using Unitful

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
        using GaussianBeamPropagation: tensor_from_principal, rotate_symmetric, rotate_general, principal_axes

        xx = 0.44793749401064553
        xy = 0.38922522287528594
        yx = 0.6485114409721432
        yy = 0.3645526281257251
        α = 1.6328107606377502

        @testset "tensor_from_principal" begin
            @test tensor_from_principal(1, 1, α)[1] ≈ 1.0
            @test tensor_from_principal(1, 1, α)[2] ≈ 0.0 atol=eps(Float64)
            @test tensor_from_principal(1, 1, α)[3] ≈ 1.0

            @test tensor_from_principal(xx, yy, 0)[1] ≈ xx
            @test tensor_from_principal(xx, yy, 0)[2] ≈ 0.0 atol=eps(Float64)
            @test tensor_from_principal(xx, yy, 0)[3] ≈ yy

            @test tensor_from_principal(xx, yy, π)[1] ≈ xx
            @test tensor_from_principal(xx, yy, π)[2] ≈ 0.0 atol=eps(Float64)
            @test tensor_from_principal(xx, yy, π)[3] ≈ yy

            @test tensor_from_principal(xx, yy, π/2)[1] ≈ yy
            @test tensor_from_principal(xx, yy, π/2)[2] ≈ 0.0 atol=eps(Float64)
            @test tensor_from_principal(xx, yy, π/2)[3] ≈ xx

            @test tensor_from_principal(xx, yy, π/4)[1] ≈ (xx+yy)/2
            @test tensor_from_principal(xx, yy, π/4)[2] ≈ (xx-yy)/2
            @test tensor_from_principal(xx, yy, π/4)[3] ≈ (xx+yy)/2
        end

        @testset "rotate_symmetric" begin
            @test rotate_symmetric(xx, 0.0, yy, α)[1] .≈ tensor_from_principal(xx, yy, α)[1]
            @test rotate_symmetric(xx, 0.0, yy, α)[2] .≈ tensor_from_principal(xx, yy, α)[2]
            @test rotate_symmetric(xx, 0.0, yy, α)[3] .≈ tensor_from_principal(xx, yy, α)[3]
        end

        @testset "rotate_general" begin
            @test rotate_general(xx, 0.0, 0.0, yy, α)[1] ≈ tensor_from_principal(xx, yy, α)[1]
            @test rotate_general(xx, 0.0, 0.0, yy, α)[2] ≈ tensor_from_principal(xx, yy, α)[2]
            @test rotate_general(xx, 0.0, 0.0, yy, α)[3] ≈ tensor_from_principal(xx, yy, α)[2]
            @test rotate_general(xx, 0.0, 0.0, yy, α)[4] ≈ tensor_from_principal(xx, yy, α)[3]

            @test rotate_general(xx, xy, xy, yy, α)[1] ≈ rotate_symmetric(xx, xy, yy, α)[1]
            @test rotate_general(xx, xy, xy, yy, α)[2] ≈ rotate_symmetric(xx, xy, yy, α)[2]
            @test rotate_general(xx, xy, xy, yy, α)[3] ≈ rotate_symmetric(xx, xy, yy, α)[2]
            @test rotate_general(xx, xy, xy, yy, α)[4] ≈ rotate_symmetric(xx, xy, yy, α)[3]

            @test rotate_general(xx, xy, yx, yy, α)[2] - rotate_general(xx, xy, yx, yy, α)[3] ≈ xy - yx
        end

        @testset "principal_axes" begin
            @test principal_axes(tensor_from_principal(1, 2, 0))[1] ≈ 1.0
            @test principal_axes(tensor_from_principal(1, 2, 0))[2] ≈ 2.0
            @test principal_axes(tensor_from_principal(1, 2, 0))[3] ≈ 0.0 atol=eps(Float64)

            @test principal_axes(tensor_from_principal(1, 2, π/2))[1] ≈ 2.0
            @test principal_axes(tensor_from_principal(1, 2, π/2))[2] ≈ 1.0
            @test principal_axes(tensor_from_principal(1, 2, π/2))[3] ≈ 0.0 atol=eps(Float64)

            @test principal_axes(tensor_from_principal(1, 2, π/4))[1] ≈ 2.0
            @test principal_axes(tensor_from_principal(1, 2, π/4))[2] ≈ 1.0
            @test principal_axes(tensor_from_principal(1, 2, π/4))[3] ≈ -π/4

            @test principal_axes(tensor_from_principal(1, 2, -π/4))[1] ≈ 2.0
            @test principal_axes(tensor_from_principal(1, 2, -π/4))[2] ≈ 1.0
            @test principal_axes(tensor_from_principal(1, 2, -π/4))[3] ≈ π/4

            @test principal_axes(tensor_from_principal(2, 1, π/4))[1] ≈ 2.0
            @test principal_axes(tensor_from_principal(2, 1, π/4))[2] ≈ 1.0
            @test principal_axes(tensor_from_principal(2, 1, π/4))[3] ≈ π/4

            @test principal_axes(tensor_from_principal(2, 1, -π/4))[1] ≈ 2.0
            @test principal_axes(tensor_from_principal(2, 1, -π/4))[2] ≈ 1.0
            @test principal_axes(tensor_from_principal(2, 1, -π/4))[3] ≈ -π/4
        end
    end

    @testset "Wigner Beam" begin
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
                    Beam(800e-9, BeamRadius(10e-3), BeamRadius(10e-3); φ=0.0)
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
                    Beam(800e-9, BeamRadius(10e-3), BeamRadius(10e-3); φ=0.0)
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

    @testset "Unitful" begin
        @test Beam(632.8u"nm", BeamRadius(1u"mm", 2u"m"), Waist(2u"mm"; z=10u"cm"); φ=0.3) ≈
              Beam(632.8e-9, BeamRadius(1e-3, 2.0), Waist(2e-3; z=0.1); φ=0.3)
        @test Beam(1u"μm", RayleighRange(1u"m"; z=1u"m")) ≈ Beam(1e-6, RayleighRange(1.0; z=1.0))
        @test Beam(1u"μm", ComplexBeamParameter((1 + 2im) * u"m")) ≈
              Beam(1e-6, ComplexBeamParameter(1.0 + 2.0im))
        @test OpticalPower(2u"m^-1") == OpticalPower(2)
        @test Curvature(1u"m^-1", 2u"m^-1"; φ=0.1) == Curvature(1, 2; φ=0.1)
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