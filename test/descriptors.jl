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
