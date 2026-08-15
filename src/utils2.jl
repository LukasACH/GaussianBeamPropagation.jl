function lensmaker(front::Curvature, back::Curvature, n::Real)
    P = (n - 1) * (front - back)
    return OpticalPower(P.xx, P.xy, P.yy)
end
