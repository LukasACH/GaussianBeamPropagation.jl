"""
    SymmetricDescriptor{T}

Supertype of descriptors given by a symmetric 2×2 tensor `[xx xy; xy yy]`, such as
[`Curvature`](@ref) and [`OpticalPower`](@ref).

Subtypes have the fields `xx`, `xy` and `yy`, and get their constructors, arithmetic and
`show` method from `@symmetric_descriptor`.
"""
abstract type SymmetricDescriptor{T<:Real} end

"""
    @symmetric_descriptor D symbol inverse_symbol

Define the constructors, arithmetic and `show` method of the [`SymmetricDescriptor`](@ref) `D`.
`symbol` and `inverse_symbol` are the names of the quantity and its inverse used by `show`, for
example `"κ"` and `"R"` for a curvature and its radius of curvature.
"""
macro symmetric_descriptor(D, symbol, inverse_symbol)
    return esc(quote
        $D(d::$D) = $D(d.xx, d.xy, d.yy)
        $D{T}(d::$D) where {T<:Real} = $D{T}(T(d.xx), T(d.xy), T(d.yy))
        Base.convert(::Type{$D{T}}, d::$D) where {T<:Real} = $D{T}(d)

        $D(xx::Real, xy::Real, yy::Real) = $D(promote(xx, xy, yy)...)
        function $D(x::Real, y::Real; φ::Real=0.0f0)
            φ = convert(promote_type(typeof(x), typeof(y)), φ)
            return $D(tensor_from_principal(x, y, φ)...)
        end
        $D(d::Real) = $D(d, d)
        $D() = $D(0.0f0)

        Base.:-(d::$D) = $D(-d.xx, -d.xy, -d.yy)
        Base.:-(a::$D, b::$D) = $D(a.xx - b.xx, a.xy - b.xy, a.yy - b.yy)
        Base.:+(a::$D, b::$D) = $D(a.xx + b.xx, a.xy + b.xy, a.yy + b.yy)
        Base.:*(a::Real, d::$D) = $D(a * d.xx, a * d.xy, a * d.yy)
        Base.:*(d::$D, a::Real) = a * d
        Base.:/(d::$D, a::Real) = $D(d.xx / a, d.xy / a, d.yy / a)

        function Base.show(io::IO, ::MIME"text/plain", d::$D)
            x, y, φ = principal_axes(d.xx, d.xy, d.yy)
            return print(
                io, repr(d),
                " with ", $symbol, "x = ", x, " (", $inverse_symbol, "x = ", 1 / x, "),",
                " ", $symbol, "y = ", y, " (", $inverse_symbol, "y = ", 1 / y, "),",
                " φ = ", rad2deg(φ), "°",
            )
        end
    end)
end
