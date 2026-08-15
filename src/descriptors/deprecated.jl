@deprecate RadiusOfCurvature(Rx::Real, Ry::Real, φ::Real) Curvature(1/Rx, 1/Ry; φ)
@deprecate RadiusOfCurvature(Rx::Real, Ry::Real; φ::Real=0.0) Curvature(1/Rx, 1/Ry; φ)
@deprecate RadiusOfCurvature(R::Real) Curvature(1/R)
@deprecate RadiusOfCurvature(; Rx::Real=Inf, Ry::Real=Inf, φ::Real=0.0) Curvature(
    1/Rx, 1/Ry; φ,
)

@deprecate FocalPower(Dx::Real, Dy::Real, φ::Real) OpticalPower(Dx, Dy; φ)
@deprecate FocalPower(Dx::Real, Dy::Real; φ::Real=0.0) OpticalPower(Dx, Dy; φ)
@deprecate FocalPower(D::Real) OpticalPower(D)
@deprecate FocalPower(; Dx::Real=0.0, Dy::Real=0.0, φ::Real=0.0) OpticalPower(Dx, Dy; φ)

@deprecate FocalLength(fx::Real, fy::Real, φ::Real) OpticalPower(1/fx, 1/fy; φ)
@deprecate FocalLength(fx::Real, fy::Real; φ::Real=0.0) OpticalPower(1/fx, 1/fy; φ)
@deprecate FocalLength(f::Real) OpticalPower(1/f)
@deprecate FocalLength(; fx::Real=Inf, fy::Real=Inf, φ::Real=0.0) OpticalPower(
    1/fx, 1/fy; φ,
)