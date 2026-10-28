module UnitfulExt

using GaussianBeamPropagation
using Unitful

const Length = Quantity{<:Real,Unitful.𝐋}
const InverseLength = Quantity{<:Real,Unitful.𝐋^-1}

# Lengths are given in metres inside the package
metres(x) = x isa Quantity ? ustrip(u"m", x) : x

GaussianBeamPropagation.OpticalPower(P_xx::InverseLength, P_xy::InverseLength, P_yy::InverseLength) =
    OpticalPower(ustrip(u"m^-1", P_xx), ustrip(u"m^-1", P_xy), ustrip(u"m^-1", P_yy))

GaussianBeamPropagation.OpticalPower(Px::InverseLength, Py::InverseLength; φ::Real=0.0f0) =
    OpticalPower(ustrip(u"m^-1", Px), ustrip(u"m^-1", Py); φ)

GaussianBeamPropagation.OpticalPower(P::InverseLength) = OpticalPower(ustrip(u"m^-1", P))

GaussianBeamPropagation.Curvature(κ_xx::InverseLength, κ_xy::InverseLength, κ_yy::InverseLength) =
    Curvature(ustrip(u"m^-1", κ_xx), ustrip(u"m^-1", κ_xy), ustrip(u"m^-1", κ_yy))

GaussianBeamPropagation.Curvature(κx::InverseLength, κy::InverseLength; φ::Real=0.0f0) =
    Curvature(ustrip(u"m^-1", κx), ustrip(u"m^-1", κy); φ)

GaussianBeamPropagation.Curvature(κ::InverseLength) = Curvature(ustrip(u"m^-1", κ))

GaussianBeamPropagation.BeamRadius(w::Length, R::Length=Inf * u"m") =
    BeamRadius(metres(w), metres(R))

GaussianBeamPropagation.Waist(w0::Length; z::Length=0.0u"m") = Waist(metres(w0); z=metres(z))

GaussianBeamPropagation.RayleighRange(zR::Length; z::Length=0.0u"m") =
    RayleighRange(metres(zR); z=metres(z))

GaussianBeamPropagation.ComplexBeamParameter(q::Quantity{<:Complex,Unitful.𝐋}) =
    ComplexBeamParameter(metres(q))

GaussianBeamPropagation.Beam(λ::Length, args...; kwargs...) =
    Beam(metres(λ), args...; (k => metres(v) for (k, v) in kwargs)...)

GaussianBeamPropagation.SecondOrderMoments(λ::Length; kwargs...) =
    SecondOrderMoments(metres(λ); (k => metres(v) for (k, v) in kwargs)...)

function GaussianBeamPropagation.SecondOrderMoments(
    rxrx::Quantity{<:Real,Unitful.𝐋^2},
    rxry::Quantity{<:Real,Unitful.𝐋^2},
    ryry::Quantity{<:Real,Unitful.𝐋^2},
    θxθx::Union{Real,Quantity{<:Real,NoDims}},
    θxθy::Union{Real,Quantity{<:Real,NoDims}},
    θyθy::Union{Real,Quantity{<:Real,NoDims}},
    rxθx::Length,
    rxθy::Length,
    ryθx::Length,
    ryθy::Length,
)
    return SecondOrderMoments(
        ustrip(u"m^2", rxrx),
        ustrip(u"m^2", rxry),
        ustrip(u"m^2", ryry),
        ustrip(NoUnits, θxθx),
        ustrip(NoUnits, θxθy),
        ustrip(NoUnits, θyθy),
        ustrip(u"m", rxθx),
        ustrip(u"m", rxθy),
        ustrip(u"m", ryθx),
        ustrip(u"m", ryθy),
    )
end

end
