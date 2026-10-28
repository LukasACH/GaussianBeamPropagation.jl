@testset "Tensors" begin
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
