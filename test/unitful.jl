@testset "Unitful" begin
    @test Beam(632.8u"nm", BeamRadius(1u"mm", 2u"m"), Waist(2u"mm"; z=10u"cm"); φ=0.3) ≈
          Beam(632.8e-9, BeamRadius(1e-3, 2.0), Waist(2e-3; z=0.1); φ=0.3)
    @test Beam(1u"μm", RayleighRange(1u"m"; z=1u"m")) ≈ Beam(1e-6, RayleighRange(1.0; z=1.0))
    @test Beam(1u"μm", ComplexBeamParameter((1 + 2im) * u"m")) ≈
          Beam(1e-6, ComplexBeamParameter(1.0 + 2.0im))
    @test OpticalPower(2u"m^-1") == OpticalPower(2)
    @test Curvature(1u"m^-1", 2u"m^-1"; φ=0.1) == Curvature(1, 2; φ=0.1)
end
