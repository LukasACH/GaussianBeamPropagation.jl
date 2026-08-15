import Base.:*;

function Rotation(α::Real)
    s, c = sincos(α)

    rtm = SA_F64[
        c; s; 0.0; 0.0;;
        -s; c; 0.0; 0.0;;
        0.0; 0.0; c; s;;
        0.0; 0.0; -s; c
    ]
    return rtm
end

# function rotate(x::Number, y::Number, α::Real)
#     s, c = sincos(α)
#     x_out = x * c^2 + y * s^2
#     xy_out = (x - y) * s * c
#     y_out = x * s^2 + y * c^2
#     return (x_out, xy_out, y_out)
# end
