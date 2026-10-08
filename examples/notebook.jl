### A Pluto.jl notebook ###
# v1.0.3

#> [frontmatter]
#> language = "en-GB"
#> 
#>     [[frontmatter.author]]
#>     name = "Lukas Affolter"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 26f1a4ae-ab0f-4230-95be-94a15062be45
begin
    import Pkg
    # activate a temporary environment
    Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(name="GaussianBeamPropagation", url="https://github.com/LukasACH/GaussianBeamPropagation.jl", rev="232beb0"),
        Pkg.PackageSpec(name="Unitful"),
        Pkg.PackageSpec(name="DSP"),
        Pkg.PackageSpec(name="Makie"),
        Pkg.PackageSpec(name="CairoMakie"),
        Pkg.PackageSpec(name="PlutoUI"),
    ])
    using Unitful, DSP, Makie, CairoMakie, PlutoUI
    using GaussianBeamPropagation
end

# ╔═╡ 9039be31-25cf-4192-adea-fdc1bc71c57c
λ = 1030u"nm"

# ╔═╡ 7f2fc98a-7a64-4c26-872e-0e8de626d4a6
D_disk = OpticalPower(1 / (f_disk = 5 * u"m"))

# ╔═╡ eaebc649-e56f-4fff-b271-942f6a251645
D_vex = OpticalPower(1 / (f_vex = -0.5u"m"))

# ╔═╡ 962ab9fe-fdc2-4adf-a802-beee37ab1874
D_4f = OpticalPower(1 / (f_4f = 0.75u"m"))

# ╔═╡ a9da072e-137e-494b-bc59-53b273269550
md"Number of round-trips: $(@bind n NumberField(1:20, default=2))"

# ╔═╡ 4705a3ad-2d69-43ce-b5e6-703f1a610ef3
δ_vex = 0.55u"m"

# ╔═╡ 18921edc-4607-49f9-bec6-a721c6786fd7
beam = Beam(λ; wx=5.8u"mm"/2, wy=5.6u"mm"/2, rx=Inf * u"m", ry=(Inf * u"m"))

# ╔═╡ ded8db82-0244-4df5-8d79-7f06d2ea820b
md"Fourier transform beam size: $(@bind w0_in_raw NumberField(3.00:0.05:7.00, default=5.7)) mm"

# ╔═╡ bfe681a4-b088-440a-b384-5d8c30d108d9
w0_in = w0_in_raw * u"mm" / 2

# ╔═╡ d53a337a-68b5-4366-905d-3f479b2c4bef
zR_in = u"m"(π * w0_in^2 / λ)

# ╔═╡ 3c302bf8-8fb3-4de0-84c1-73118ba91d31
δ_4f = f_4f^2 * (1 / zR_in - 1 / f_disk + 1 / (δ_vex - f_vex))

# ╔═╡ 4adb1395-344a-40a9-8eb5-8f12131419b7
δ_fourier = zR_in * f_vex^2 / (2(δ_vex - f_vex)^2) + δ_vex * f_vex / (δ_vex - f_vex)

# ╔═╡ 2aa89eae-aa50-4217-85a4-4996dce10dbb
md"""
### Layout:
| Distance | Length |
|:-------- | ------- |
| 4f-imaging | $(round(u"cm", δ_4f; sigdigits=3)) |
| disk to vex | $(round(u"cm", δ_vex; sigdigits=3)) |
| vex to end | $(round(u"cm", δ_fourier; sigdigits=3)) |
"""

# ╔═╡ b96e83d3-6917-4924-b530-c40f8de43964
md"Astigmatism compensated: $(@bind compensated Switch())"

# ╔═╡ 83c75062-091a-4180-9f4b-f7799356eadd
md"""
Astigmatic surface of thin disk:
- Variance: $(@bind disk_variance_change PlutoUI.Slider(0.0:0.001:0.02, show_value=true))
- Angle: $(@bind disk_variance_angle PlutoUI.Slider(0:180, show_value=true))
"""

# ╔═╡ 0a107c1e-3c84-40a4-8b7b-cd2d75ad5d4f
begin
    local fig = Figure()
    local axis = Axis(fig[1,1])
    
        local xvals = @. (4:2:22)
        scatter!(axis, xvals, [6791.54, 5779.78, 5808.05, 5876.84, 5854.01, 5701.75, 5904.95, 5806.29, 6110.66, 8597.03] ./ 2000)
    
        scatter!(axis, xvals, [11947.04, 9984.71, 8936.01, 7238.42, 6241.51, 5934.08, 5730.87,5716.14, 5696.1, 6521.51] ./ 2000)

    fig
end

# ╔═╡ adbc4a99-02f4-4dac-abad-3b32b6150cde
IncidenceAngle(ι=atan(sqrt(64^2 + 4^2) / 2, 330), θ=0)

# ╔═╡ 49e5a4f1-6d31-4fc2-b055-4c9f89a4904a
incidence_vex = incidence_vex = IncidenceAngle(ι=atan(50 / 2, 330), θ=0) # IncidenceAngle(ι=atan(sqrt(64^2 + 4^2) / 2, 330), θ=0)

# ╔═╡ 352e37a2-b912-40cc-94c4-a945bd86019e
begin
    if compensated
        # (dx1, dy1) = (24.5, 24.5) # (29.2190664122, 29.2190664122) # (49, 49)
        # (dx2, dy2) = (27.5, 27.5) # (30, 30)
        (dx1, dy1) = (30, 125-95) # (29.2190664122, 29.2190664122) # (49, 49)
        (dx2, dy2) = (38, 125-95) # (30, 30)

        # Limiting is d2:
        d1 = sqrt(dx1^2 + dy1^2)
        d2 = sqrt(dx2^2 + dy2^2) # 39 # from calculation uncompensated

        incidence_4f_1 = IncidenceAngle(ι=atan(d1 / 2, 480), θ=atan(dy1, dx1))
        incidence_4f_2 = IncidenceAngle(ι=atan(d2 / 2, 530), θ=atan(-dy2, dx2))

        # incidence_4f_1 = IncidenceAngle(ι=atan(50 / 2, 480), θ=deg2rad(30))
        # incidence_4f_2 = IncidenceAngle(ι=atan(39 / 2, 530), θ=-deg2rad(30))

        # incidence_4f_1 = IncidenceAngle(ι=atan(40 / 2, 480), θ=deg2rad(00))
        # incidence_4f_2 = IncidenceAngle(ι=atan(40 / 2, 530), θ=-deg2rad(00))

        # incidence_4f_1 = IncidenceAngle(ι=incidence_4f_2.ι * (1 + δ_4f / f_4f), θ=deg2rad(45)) # 

        d1p = 2tan(incidence_4f_1.ι)
        d2p = 2tan(incidence_4f_2.ι)

        @info "" d1 d2 d1p d2p rad2deg(incidence_4f_1.ι) rad2deg(incidence_4f_2.ι)

        @info "" sincosd(30) .* 50 sincosd(30) .* 39
    else
        # One is sized 1.5 inches 480 away, the other 1 inch 750 away
        # the main is 2 inch
        # the futher away one at the position of the closer one is
        dx1 = begin
            r_m = 2 * 25.4 / 2
            d_i = 750
            r_i = 1 * 25.4 / 2
            d_o = 480
            r_o = 1.5 * 25.4 / 2

            r_1 = (min(d_i, d_o) * r_i + (d_i - min(d_i, d_o)) * r_m) / d_i
            r_2 = (min(d_i, d_o) * r_o + (d_o - min(d_i, d_o)) * r_m) / d_o

            r_1 + r_2 + 5 # 5mm for mount on closer mirror
        end
        dx2 = begin
            r_m = 2 * 25.4 / 2
            d_i = 530
            r_i = 1.5 * 25.4 / 2
            d_o = 750
            r_o = 21 / 2

            r_1 = (min(d_i, d_o) * r_i + (d_i - min(d_i, d_o)) * r_m) / d_i
            r_2 = (min(d_i, d_o) * r_o + (d_o - min(d_i, d_o)) * r_m) / d_o

            r_1 + r_2 + 5 # 5mm for mount on closer mirror
        end
        # @info "" dx1 dx2

        incidence_4f_1 = IncidenceAngle(ι=atan(dx1 / 2, 480) * 2, θ=0)
        incidence_4f_2 = IncidenceAngle(ι=atan(dx2 / 2, 530) * 2, θ=0)
    end
end

# ╔═╡ 2a15961b-e028-4692-9097-b7bb75558e08


# ╔═╡ fff03b9f-474b-42fd-9f8d-2914c8ed641a
d1_4f = FreeSpace(ustrip(u"m", f_4f))

# ╔═╡ 2e403cf7-a811-45e0-b77c-3cace0d4a28e
dc_4f = FreeSpace(ustrip(u"m", 2f_4f + δ_4f))

# ╔═╡ b1576b17-f8e7-45c5-b8d5-ab91f88b65f3
d2_4f = FreeSpace(ustrip(u"m", f_4f))

# ╔═╡ c0ab9bb5-6dd8-4f0a-afc8-b1063c961d0a
d_vex = FreeSpace(ustrip(u"m", δ_vex)) #FreeSpace2D(0.582103)

# ╔═╡ 1da9a051-6af1-4eda-957d-6fb1dca94a5f
d_fourier = FreeSpace(ustrip(u"m", δ_fourier)) #FreeSpace2D(2.66144)

# ╔═╡ 10bf3584-2608-4dbe-bffd-c93b205404dc
D_detuning = OpticalPower(δ_4f / f_4f^2)

# ╔═╡ 3fa53dff-4631-4fa4-b9a4-83029e4b135b
D_fourier = D_detuning + D_disk

# ╔═╡ 5b1a8ab7-1844-449e-892b-cd2996051420
italic(x) = rich(x; font=:italic)

# ╔═╡ 2156bbcd-1fd2-434a-a9dc-be9166d174a7
function lims_pad(; low, high, pad=0.05)
    return (; low=low - (high - low) * pad, high=high + (high - low) * pad)
end

# ╔═╡ e345e713-830d-4820-80b3-5348e3ca69ee
incidence_disk = IncidenceAngle(ι=atan(48 / 2, 580), θ=0)

# ╔═╡ eea0920f-b457-465b-b263-478f3bbc7418
m4f1 = Mirror(D_4f, incidence=incidence_4f_1)

# ╔═╡ ab65d764-0e97-4891-aee3-c0d4480fe3c3
m4f2 = Mirror(D_4f, incidence=incidence_4f_2)

# ╔═╡ 926545c7-40fc-470f-85e1-5754978af885
mvex = Mirror(D_vex, incidence=incidence_vex)

# ╔═╡ d9094462-c1bb-41b2-954d-4f0995d3f946
part5 = d1_4f * m4f1 * dc_4f

# ╔═╡ ab6aaa86-91b7-48fb-b887-5749a9942020
part3 = d_vex * mvex * d_fourier * d_fourier * mvex * d_vex

# ╔═╡ d3576f04-9c0c-4763-b4b5-1958b31b4ffe
part1 = dc_4f * m4f1 * d1_4f

# ╔═╡ 0097f3bd-989d-42b5-aace-4c328a3c8638
# ╠═╡ disabled = true
#=╠═╡
beam = Beam(λ; wx=w0_in, wy=w0_in, rx=Inf * u"m", ry=(-Inf * u"m"))
  ╠═╡ =#

# ╔═╡ 61661f96-a646-4383-8810-3567461a3991
Ds = range(-0.06, 0.06, 1_00)

# ╔═╡ 878b2db6-1703-47e2-b5c6-5fe75a39827b
sample_points = Ds .* ustrip(u"m", f_4f)^2

# ╔═╡ 6201acae-04a9-43b1-9940-11894e6589d3
CairoMakie.activate!(type = "svg")

# ╔═╡ 99580059-46f8-4eca-afb1-f0bd863bd0b9
md"## Plotting"

# ╔═╡ 351909be-23ba-4f99-b19d-38626aeb459b
disk_variance = OpticalPower(disk_variance_change, -disk_variance_change; φ=deg2rad(disk_variance_angle))

# ╔═╡ 7fb71e8b-2085-478b-be03-9ee11fde9ac2
disk = Mirror(D_disk + disk_variance, incidence=incidence_disk)

# ╔═╡ 8e43b31c-9dde-4d97-8e25-d46e99db3191
part4 = m4f2 * d2_4f * disk

# ╔═╡ dd1f81cb-fa26-4b60-a1e9-deea16823669
part2 = disk * d2_4f * m4f2

# ╔═╡ a39bce74-278b-4f76-a457-6594961f643b
function amplifier(; ΔD=0.0, Δ4f=0.0)
    Δ4f = FreeSpace(Δ4f)
    Δdisk = Mirror(FocalPower(ΔD))
    return (part5 * Δ4f * part4 * Δdisk * part3 * Δdisk * part2 * Δ4f * part1)^n
end

# ╔═╡ c5a3348c-1534-4292-aa90-aa6ea3bc0b6a
md"## Support functions"

# ╔═╡ c933d42d-8ab4-49ae-b505-b6ed4fbcc26e
sample_beam(beams::AbstractVector{<:Beam}) = sample_beam(iso_radius, beams)

# ╔═╡ 0d5352e9-ee3f-47d7-af14-1047aa6266c4
function sample_beam_continous(beams::AbstractVector{<:Beam})
    return sample_beam_continous(iso_radius, beams)
end

# ╔═╡ b5592863-a581-4a75-8959-e0d0278b6ada
mod_pi(x) = mod2pi(2x) / 2

# ╔═╡ b26ed39c-c7fd-4c72-b7ca-ec2737069cdc
function sample_beam(func, beams::AbstractVector{<:Beam})
    r = map(func, beams)
    x = @. getfield(r, :x)
    y = @. getfield(r, :y)
    φ = @. getfield(r, :φ)

    return (; rx=x, ry=y, φx=mod_pi.(φ), φy=mod_pi.(φ .+ π / 2))
end

# ╔═╡ 659cf583-5d3d-43b3-9cf2-8878abbace85
mod_pi2(x) = mod_pi(2x) / 2

# ╔═╡ 4f757570-60c0-46fb-b152-a6f793fd077a
mod_mpi2_pi2(x) = mod_pi(x + π / 2) - π / 2

# ╔═╡ e8b355e3-19c0-4006-8687-37bc761fa157
mod_mpi4_pi4(x) = mod_mpi2_pi2(2x) / 2

# ╔═╡ f285b189-a0aa-45dd-8bfc-3e309c0cfcb4
function find_continous_angle(ang::AbstractVector{T}) where {T<:AbstractFloat}
    angle_90 = @. mod2pi(4ang + π) / 4 - π / 4
    angle_90_unwrapped = unwrap(angle_90; range=π / 2)
    ints = @. Int(round((ang - angle_90_unwrapped) / (π / 2)))
    # return ints

    angle_clamped = @. mod(2ang / π + 0.5, 1) - 0.5
    angle_clamped_2 = @. mod(2ang / π + 1, 2) - 1
    angle_unwrapped = unwrap(angle_clamped; range=1)
    ints = @. Int(round(angle_clamped_2 - angle_unwrapped))
    return ints

    # The angle is between 0 and 2π. We only care about 0 and π. We want to detect if the angle changed by π/2.
    # When we take the difference between two consecutive angles, if they are offset by more than π/4, we consider this a change.
    # For that, we want to first calculate the differences, and then check if they are between π/4 and 3π/4 or 5π/4 and 7π/4.
    angle_differece = @views ang[begin:(end - 1)] - ang[(begin + 1):end]
    # Now we map everything onto 0 and π.
    angle_differece_mod = @. mod2pi(2angle_differece) / 2

    bang = @. (1π / 4 < angle_differece_mod < 3π / 4)

    flip = Vector{Int}(undef, length(ang))
    fill!(flip, 0)
    # flip[begin] = 0
    # flip_view = @view flip[begin+1:end]
    # cumsum!(flip_view, bang)
    return flip
end

# ╔═╡ 7c348e71-22e6-4f2b-98e5-73ab5b29c935
function sample_beam_continous(func, beams::AbstractVector{<:Beam})
    (; rx, ry, φx, φy) = sample_beam(func, beams)

    sw_φ = find_continous_angle(φx)

    rx_out = @. ifelse(iseven(sw_φ), rx, ry)
    ry_out = @. ifelse(iseven(sw_φ), ry, rx)

    φx_in = @. mod_pi(ifelse(iseven(sw_φ), φx, φy))
    φy_in = @. mod_pi(ifelse(iseven(sw_φ), φy, φx))

    φx_out = unwrap(φx_in; range=pi) # @. mod2pi(ifelse(iseven(nang), 2ang, 2ang + π) + π/2) / 2 - π/4
    φy_out = unwrap(φy_in; range=pi) # @. mod2pi(ifelse(isodd(nang), 2ang, 2ang + π) + π/2) / 2 - π/4

    return (; rx_out, ry_out, φx_out, φy_out)
end

# ╔═╡ b0edc5a7-45cb-4dfd-942d-c0f88d2de2ce
begin
    rx, ry, φ_rx, φ_ry = sample_beam_continous(
        radius,
        # map(Δ4f -> amplifier(; Δ4f) * beam, sample_points),
        map(ΔD -> amplifier(; ΔD) * beam, Ds),
    )
    
    clamp!(rx, 0, 1000)
    clamp!(ry, 0, 1000)
end

# ╔═╡ 68b9bc70-9470-4c47-a7a6-cb6ccb1db4e2
Rx, Ry, φ_Rx, φ_Ry = sample_beam_continous(
        phase_curvature,
        # map(Δ4f -> amplifier(; Δ4f) * beam, sample_points),
        map(ΔD -> amplifier(; ΔD) * beam, Ds),
    )

# ╔═╡ c2155af3-2479-462c-8348-f82975ba38d0
macro px_str(value)
    return parse(Float64, value) * 0.1
end

# ╔═╡ 9748ed29-2425-47f4-bc3a-10797fb89b92
m_px() = 1 * 3 # units_per_pixel

# ╔═╡ 9b76b722-65ce-4755-9cf8-1bd46f3051c0
m_in() = 96 * m_px() # pt / in

# ╔═╡ 109645f8-2132-4cd7-beab-30ae0b0e46c0
m_cm() = m_in() / 2.54

# ╔═╡ 335d8235-55b7-4ce8-a92d-7326ff8bfb6c
m_mm() = m_cm() / 10

# ╔═╡ b826e551-f94d-41ef-bc6c-7c86c67ff973
begin
        degrees = false
        style_wx = (; color=Makie.wong_colors()[1], linestyle=(:solid))
        style_wy = (; color=Makie.wong_colors()[1], linestyle=(:dot, :dense))
        style_Rx = (; color=Makie.wong_colors()[2], linestyle=(:solid))
        style_Ry = (; color=Makie.wong_colors()[2], linestyle=(:dot, :dense))

        fig = Figure()

        common_options = (;
            xticks=collect(-0.05:0.025:0.05),
            xlabel=rich(
                "variation of disk focal power  Δ",
                italic("D"),
                subscript("disk"),
                " ∕ m",
                superscript("−1"),
            ),
            limits=((-0.06, 0.06), nothing),
            xminorgridvisible=true,
            yminorgridvisible=true,
            xminorticksvisible=true,
            yminorticksvisible=true,
            xminorticks=IntervalsBetween(5),
        )

        axis_waist = Axis(
            fig[1, 1];
            ylabel=rich("radius  ", italic("w"), " ∕ mm"),
            yticks=collect(0:0.5:5),
            yminorticks=IntervalsBetween(5),
            common_options...,
        )
        ylims!(axis_waist; lims_pad(; low=2.5, high=4.0)...)

    

        # local xvals = @. (-9:2:9) * 2 / 1000 / 0.75^2 - 0.01
        # scatter!(axis_waist, xvals, [6791.54, 5779.78, 5808.05, 5876.84, 5854.01, 5701.75, 5904.95, 5806.29, 6110.66, 8597.03] ./ 2000)
        # scatter!(axis_waist, xvals, [11947.04, 9984.71, 8936.01, 7238.42, 6241.51, 5934.08, 5730.87,5716.14, 5696.1, 6521.51] ./ 2000)
    
        axis_phase = Axis(
            fig[2, 1];
            ylabel=rich(
                "phase-front curvature  ",
                italic("R"),
                superscript("−1"),
                " ∕ km",
                superscript("−1"),
                "",
            ),
            # yticks=collect(-20:20:20),
            yminorticks=IntervalsBetween(4),
            common_options...,
        )
        ylims!(axis_phase; lims_pad(; low=-20, high=20)...)

        yticks_angle = (
            (collect((-(π/2)):(π/4):(π/2)) .- 0.000001),
            ["−π∕2", "−π∕4", "0", "π∕4", "π∕2"],
        )
        #yticks_angle = (0:(π / 4):(π / 2), ["0", "π/4", "π/2"])

        axis_angle = Axis(
            fig[3, 1];
            ylabel="principal axes orientation",
            yticks=yticks_angle,
            yminorticks=IntervalsBetween(3),
            common_options...,
        )
        ylims!(axis_angle; lims_pad(; low=0, high=π / 2, pad=0.2)...)

        lines!(axis_waist, Ds, 1000 .* rx; label="X", style_wx...)
        lines!(axis_waist, Ds, 1000 .* ry; label="Y", style_wy...)
        lines!(axis_phase, Ds, 1000 ./ Rx; label="X", style_Rx...)
        lines!(axis_phase, Ds, 1000 ./ Ry; label="Y", style_Ry...)


        lines!(axis_waist, Ds, 1000 .* (rx .* cos.(φ_rx .- π/4) .^ 2 .+ ry .* sin.(φ_rx .- π/4) .^ 2); label="X", style_Rx...)
        lines!(axis_waist, Ds, 1000 .* (ry .* cos.(φ_rx .- π/4) .^ 2 .+ rx .* sin.(φ_rx .- π/4) .^ 2); label="Y", style_Ry...)
    
        begin
            mod_φ_Rx = @. mod_pi(φ_Rx + (π / 4)) - (π / 4)
            mod_φ_Ry = @. mod_pi(φ_Ry + (π / 4)) - (π / 4)
            int_φ_Rx = @. Int(round((φ_Rx - mod_φ_Rx) / π))
            int_φ_Ry = @. Int(round((φ_Ry - mod_φ_Ry) / π))

            mod_φ_rx = @. mod_pi(φ_rx + (π / 4)) - (π / 4)
            mod_φ_ry = @. mod_pi(φ_ry + (π / 4)) - (π / 4)
            int_φ_rx = @. Int(round((φ_rx - mod_φ_rx) / (π)))
            int_φ_ry = @. Int(round((φ_ry - mod_φ_ry) / (π)))
            for i in unique(int_φ_Rx)
                lines!(axis_angle, Ds, φ_Rx .- i * π; label="X", style_Rx...)
            end
            for i in unique(int_φ_Ry)
                lines!(axis_angle, Ds, φ_Ry .- i * π; label="X", style_Ry...)
            end
            for i in unique(int_φ_rx)
                lines!(axis_angle, Ds, φ_rx .- i * π; label="X", style_wx...)
            end
            for i in unique(int_φ_ry)
                lines!(axis_angle, Ds, φ_ry .- i * π; label="Y", style_wy...)
            end
        end

        linkxaxes!(axis_waist, axis_phase, axis_angle)

        yticklabelspace = maximum(
            tight_yticklabel_spacing!, [axis_waist, axis_phase, axis_angle]
        )
        axis_waist.yticklabelspace = yticklabelspace
        axis_phase.yticklabelspace = yticklabelspace
        axis_angle.yticklabelspace = yticklabelspace

        hidexdecorations!(axis_waist; grid=false, minorgrid=false)
        hidexdecorations!(axis_phase; grid=false, minorgrid=false)
        # hidexdecorations!(axis_divergence; grid=false, minorgrid=false)
        linkxaxes!(axis_waist, axis_phase, axis_angle)

        rowgap!(fig.layout, 3m_mm())
        colsize!(fig.layout, 1, (6m_cm()))
        rowsize!(fig.layout, 1, (3m_cm()))
        rowsize!(fig.layout, 2, (3m_cm()))
        rowsize!(fig.layout, 3, (3m_cm()))
        resize_to_layout!(fig)

    fig
    end

# ╔═╡ 111a3c9f-fcff-4b7a-9ff1-cf33303e1ec6
m_pt() = m_in() / 72

# ╔═╡ b89fca52-0fe5-4f8a-956c-22dd6d0e3415
const fontsize_12pt = (;
    normal=12m_pt(), small=11m_pt(), footnote=10m_pt(), script=8m_pt(), tiny=6m_pt()
)

# ╔═╡ 405bf3ee-b122-4aee-9578-b151d542041a
const fontsize_11pt = (;
    normal=11m_pt(), small=10m_pt(), footnote=9m_pt(), script=8m_pt(), tiny=6m_pt()
)

# ╔═╡ bdcdf051-5d35-4408-8406-3104ef52d52c
const fontsize_10pt = (;
    normal=10m_pt(), small=9m_pt(), footnote=8m_pt(), script=7m_pt(), tiny=5m_pt()
)

# ╔═╡ d44fc394-9607-4a56-be18-05ecb01686a5
get_plot_theme() = get_plot_theme(fontsize_12pt);

# ╔═╡ 7f2124ae-bc43-45a8-b8b5-3fe776261a71
#let small = (9pt, 10pt, 11pt)
#let footnotesize = (8pt, 9pt, 10pt)
#let scriptsize = (7pt, 8pt, 8pt)
#let tiny = (5pt, 6pt, 6pt)

function get_plot_theme(fontsize)
    return Theme(;
        fonts=(;
            regular="EB Garamond Regular",
            italic="EB Garamond Italic",
            bold="EB Garamond Bold",
            bold_italic="EB Garamond Bold Italic",
        ),
        fontsize=fontsize.normal,
        figure_padding=5m_mm(),#2m_mm(),
        backgroundcolor=:white,
        rowgap=3m_mm(),
        colgap=3m_mm(),
        markersize=2m_mm(),
        markerstrokewidth=0m_px(),
        linewidth=1m_pt(),
        Label=Attributes(; fontsize=fontsize.footnote),
        Axis=Attributes(;
            alignmode=Inside(),
            aspect=nothing,
            autolimitaspect=nothing,
            backgroundcolor=:white,
            dim1_conversion=nothing,
            dim2_conversion=nothing,
            flip_ylabel=false,
            halign=:center,
            height=nothing,
            tellheight=true,
            tellwidth=true,
            valign=:center,
            width=nothing,

            # Axis settings
            limits=(nothing, nothing),
            xautolimitmargin=(0.05f0, 0.05f0),
            xaxisposition=:bottom,
            xreversed=false,
            xscale=identity,
            yautolimitmargin=(0.05f0, 0.05f0),
            yaxisposition=:left,
            yreversed=false,
            yscale=identity,

            # Panning & Zooming
            xrectzoom=true,
            yrectzoom=true,
            zoombutton=true,
            xzoomlock=false,
            yzoomlock=false,
            xzoomkey=Makie.Keyboard.x,
            yzoomkey=Makie.Keyboard.y,
            panbutton=Makie.Mouse.right,
            xpanlock=false,
            ypanlock=false,
            xpankey=Makie.Keyboard.x,
            ypankey=Makie.Keyboard.y,

            # Title
            title="",
            titlealign=:center,
            # titlecolor=(Makie.@inherit :textcolor :black),
            titlefont=:bold,
            titlegap=0.25 * fontsize.normal,
            titlelineheight=1,
            titlesize=fontsize.normal,
            titlevisible=true,

            # Subtitle
            subtitle="",
            # subtitlecolor=(Makie.@inherit :textcolor :black),
            subtitlefont=:regular,
            subtitlegap=0m_px(),
            subtitlelineheight=1,
            subtitlesize=fontsize.small,

            # Spine
            spinewidth=0.8m_pt(),
            topspinevisible=true,
            topspinecolor=:black,
            rightspinevisible=true,
            rightspinecolor=:black,
            bottomspinevisible=true,
            bottomspinecolor=:black,
            leftspinevisible=true,
            leftspinecolor=:black,
            xtrimspine=false,
            ytrimspine=false,

            # Axislabels
            xlabel="",
            ylabel="",
            # xlabelcolor=(@inherit :textcolor :black),
            # ylabelcolor=(@inherit :textcolor :black),
            xlabelfont=:regular,
            ylabelfont=:regular,
            xlabelpadding=3 / 16 * fontsize.footnote, # 3m_px(),
            ylabelpadding=5 / 16 * fontsize.footnote, # 5m_px(),
            xlabelrotation=Makie.automatic,
            ylabelrotation=Makie.automatic,
            xlabelsize=fontsize.footnote,
            ylabelsize=fontsize.footnote, # 12pt
            xlabelvisible=true,
            ylabelvisible=true,

            # Ticklabels
            xticklabelalign=Makie.automatic,
            yticklabelalign=Makie.automatic,
            # xticklabelcolor=(@inherit :textcolor :black),
            # yticklabelcolor=(@inherit :textcolor :black),
            xticklabelfont=:regular,
            yticklabelfont=:regular,
            xticklabelpad=2 / 16 * fontsize.footnote, # 0m_px(), # 2m_px(),
            yticklabelpad=4 / 16 * fontsize.footnote, # 0m_px(), # 4m_px(),
            xticklabelrotation=0,
            yticklabelrotation=π / 2,
            xticklabelsize=fontsize.footnote,
            yticklabelsize=fontsize.footnote, # 12pt
            xticklabelspace=Makie.automatic,
            yticklabelspace=Makie.automatic,
            xticklabelsvisible=true,
            yticklabelsvisible=true,

            # Grid
            xgridcolor=RGBAf(0, 0, 0, 0.12),
            ygridcolor=RGBAf(0, 0, 0, 0.12),
            xgridstyle=nothing,
            ygridstyle=nothing,
            xgridvisible=true,
            ygridvisible=true,
            xgridwidth=0.6m_pt(),
            ygridwidth=0.6m_pt(),
            xminorgridcolor=RGBAf(0, 0, 0, 0.05),
            yminorgridcolor=RGBAf(0, 0, 0, 0.05),
            xminorgridstyle=nothing,
            yminorgridstyle=nothing,
            xminorgridvisible=false,
            yminorgridvisible=false,
            xminorgridwidth=0.4m_pt(),
            yminorgridwidth=0.4m_pt(),

            # Ticks
            xtickalign=0,
            ytickalign=0,
            xtickcolor=RGBf(0, 0, 0),
            ytickcolor=RGBf(0, 0, 0),
            xtickformat=Makie.automatic,
            ytickformat=Makie.automatic,
            xticks=Makie.automatic,
            yticks=Makie.automatic,
            xticksize=1m_mm(), # 5m_px(),
            yticksize=1m_mm(), # 5m_px(),
            xticksmirrored=false,
            yticksmirrored=false,
            xticksvisible=true,
            yticksvisible=true,
            xtickwidth=0.6m_pt(), # 1m_px(),
            ytickwidth=0.6m_pt(), # 1m_px(),

            # Minorticks
            xminortickalign=0,
            yminortickalign=0,
            xminortickcolor=:black,
            yminortickcolor=:black,
            xminorticks=IntervalsBetween(2),
            yminorticks=IntervalsBetween(2),
            xminorticksize=0.5m_mm(),
            yminorticksize=0.5m_mm(),
            xminorticksvisible=false,
            yminorticksvisible=false,
            xminortickwidth=0.4m_pt(),
            yminortickwidth=0.4m_pt(),

            # yticklabelspace=:max_auto,
            #xtickformat=values -> [L"%$(value)" for value in Makie.showoff_minus(values)],
            #ytickformat=values -> [L"%$(value)" for value in Makie.showoff_minus(values)],
        ),
        Legend=Attributes(;
            # alignmode=Inside(),
            # alpha=1,
            backgroundcolor=:white,
            # bgcolor=nothing,
            colgap=1 * fontsize.script,
            framecolor=:transparent,
            # framevisible=false,
            framewidth=1m_px(),
            # gridshalign=:centre,
            # gridsvalign=:centre,
            groupgap=1 * fontsize.script,
            # halign=:centre,
            # heatmapcolorrange=automatic,
            # heatmaplimits=(0 .. 1, 0 .. 1),
            # heatmapvalues=[0 0.3;0.6 1],
            # height=Auto(),
            # imagecolorrange=automatic,
            # imagelimits=(0 .. 1, 0 .. 1),
            # imagevalues=[0 0.3;0.6 1],
            # label="undefined",
            # labelcolor=@inherit :textcolor :black,
            # labelfont=:regular,
            # labelhalign=:left,
            # labeljustification=automatic,
            labelsize=1 * fontsize.script, # (Makie.@inherit :fontsize 12m_pt()),
            # labelvalign=:centre,
            # linecolor=theme(scene, :linecolor),
            # linecolormap=theme(scene, :colormap),
            # linecolorrange=automatic,
            # linepoints=[Point2f(0, 0.5), Point2f(1, 0.5)],
            # linestyle=:solid,
            # linewidth=theme(scene, :linewidth),
            margin=(0.0f0, 0.0f0, 0.0f0, 0.0f0),
            # marker=theme(scene, :marker),
            # markercolor=theme(scene, :markercolor),
            # markercolormap=theme(scene, :colormap),
            # markercolorrange=automatic,
            # markerpoints=[Point2f(0.5, 0.5)],
            # markersize=theme(scene, :markersize),
            # markerstrokecolor=theme(scene, :markerstrokecolor),
            # markerstrokewidth=theme(scene, :markerstrokewidth),
            # mesh=Rect3f(Point3f(-0.7), Vec3f(1.4)),
            # meshcolor=(wong_colors())[1],
            # meshcolormap=theme(scene, :colormap),
            # meshcolorrange=automatic,
            # meshscattercolor=theme(scene, :markercolor),
            # meshscattercolormap=theme(scene, :colormap),
            # meshscattercolorrange=automatic,
            meshscattermarker=Sphere(Point3f(0), 1.0m_px()),
            # meshscatterpoints=[Point3f(0)],
            # meshscatterrotation=Quaternionf(0, 0, 0, 1),
            meshscattersize=0.8m_px(),
            # nbanks=1,
            # orientation=:vertical,
            padding=(0m_pt(), 0m_pt(), 0m_pt(), 0m_pt()), # (6.0f0, 6.0f0, 6.0f0, 6.0f0)
            patchcolor=:transparent,
            patchlabelgap=1 * 5 / 16 * fontsize.script,
            patchsize=(20 / 16 * fontsize.script, 20 / 16 * fontsize.script),
            # patchstrokecolor=:transparent,
            patchstrokewidth=0.8m_mm(),
            # polycolor=theme(scene, :patchcolor),
            # polycolormap=theme(scene, :colormap),
            # polycolorrange=automatic,
            # polypoints=[Point2f(0, 0), Point2f(1, 0), Point2f(1, 1), Point2f(0, 1)],
            # polystrokecolor=theme(scene, :patchstrokecolor),
            # polystrokewidth=theme(scene, :patchstrokewidth),
            rowgap=0, # -0.5 * fontsize.script, #-8m_px(),
            # surfacecolormap=theme(scene, :colormap),
            # surfacecolorrange=automatic,
            # surfacedata=(-0.7 .. 0.7, -0.7 .. 0.7, [-0.007 * x ^ 3 * (1 - 0.05 * y ^ 2) for x = -5:5, y = -5:5]),
            # surfacevalues=automatic,
            # tellheight=automatic,
            # tellwidth=automatic,
            # titlecolor=@inherit :textcolor :black,
            titlefont=:bold,
            titlegap=0m_px(),
            # titlehalign=:centre,
            # titleposition=:top,
            titlesize=fontsize.script, # (Makie.@inherit :fontsize 12m_pt()),
            # titlevalign=:centre,
            # titlevisible=true,
            # valign=:centre,
            # width=Auto(),
        ),
        CairoMakie=(px_per_unit=2.0 / m_px(), pt_per_unit=1 / m_pt()),
        px_per_unit=2.0 / m_px(),
        pt_per_unit=1.0 / m_pt(),
    )
end

# ╔═╡ 766bfe5b-75fb-4afa-aef8-93f6dc1ab385
set_theme!(get_plot_theme())

# ╔═╡ Cell order:
# ╠═26f1a4ae-ab0f-4230-95be-94a15062be45
# ╟─9039be31-25cf-4192-adea-fdc1bc71c57c
# ╟─d53a337a-68b5-4366-905d-3f479b2c4bef
# ╠═7f2fc98a-7a64-4c26-872e-0e8de626d4a6
# ╠═eaebc649-e56f-4fff-b271-942f6a251645
# ╠═962ab9fe-fdc2-4adf-a802-beee37ab1874
# ╟─2aa89eae-aa50-4217-85a4-4996dce10dbb
# ╟─a9da072e-137e-494b-bc59-53b273269550
# ╟─3c302bf8-8fb3-4de0-84c1-73118ba91d31
# ╟─4705a3ad-2d69-43ce-b5e6-703f1a610ef3
# ╠═4adb1395-344a-40a9-8eb5-8f12131419b7
# ╠═18921edc-4607-49f9-bec6-a721c6786fd7
# ╟─bfe681a4-b088-440a-b384-5d8c30d108d9
# ╟─ded8db82-0244-4df5-8d79-7f06d2ea820b
# ╟─b96e83d3-6917-4924-b530-c40f8de43964
# ╟─83c75062-091a-4180-9f4b-f7799356eadd
# ╟─0a107c1e-3c84-40a4-8b7b-cd2d75ad5d4f
# ╠═b826e551-f94d-41ef-bc6c-7c86c67ff973
# ╠═adbc4a99-02f4-4dac-abad-3b32b6150cde
# ╠═49e5a4f1-6d31-4fc2-b055-4c9f89a4904a
# ╠═352e37a2-b912-40cc-94c4-a945bd86019e
# ╠═2a15961b-e028-4692-9097-b7bb75558e08
# ╠═fff03b9f-474b-42fd-9f8d-2914c8ed641a
# ╠═2e403cf7-a811-45e0-b77c-3cace0d4a28e
# ╠═b1576b17-f8e7-45c5-b8d5-ab91f88b65f3
# ╠═c0ab9bb5-6dd8-4f0a-afc8-b1063c961d0a
# ╠═1da9a051-6af1-4eda-957d-6fb1dca94a5f
# ╠═10bf3584-2608-4dbe-bffd-c93b205404dc
# ╠═3fa53dff-4631-4fa4-b9a4-83029e4b135b
# ╟─5b1a8ab7-1844-449e-892b-cd2996051420
# ╟─2156bbcd-1fd2-434a-a9dc-be9166d174a7
# ╠═e345e713-830d-4820-80b3-5348e3ca69ee
# ╠═eea0920f-b457-465b-b263-478f3bbc7418
# ╠═ab65d764-0e97-4891-aee3-c0d4480fe3c3
# ╠═7fb71e8b-2085-478b-be03-9ee11fde9ac2
# ╠═926545c7-40fc-470f-85e1-5754978af885
# ╟─d9094462-c1bb-41b2-954d-4f0995d3f946
# ╟─8e43b31c-9dde-4d97-8e25-d46e99db3191
# ╟─ab6aaa86-91b7-48fb-b887-5749a9942020
# ╠═dd1f81cb-fa26-4b60-a1e9-deea16823669
# ╠═d3576f04-9c0c-4763-b4b5-1958b31b4ffe
# ╠═a39bce74-278b-4f76-a457-6594961f643b
# ╠═0097f3bd-989d-42b5-aace-4c328a3c8638
# ╠═61661f96-a646-4383-8810-3567461a3991
# ╠═878b2db6-1703-47e2-b5c6-5fe75a39827b
# ╟─b0edc5a7-45cb-4dfd-942d-c0f88d2de2ce
# ╟─68b9bc70-9470-4c47-a7a6-cb6ccb1db4e2
# ╠═766bfe5b-75fb-4afa-aef8-93f6dc1ab385
# ╠═6201acae-04a9-43b1-9940-11894e6589d3
# ╟─99580059-46f8-4eca-afb1-f0bd863bd0b9
# ╠═351909be-23ba-4f99-b19d-38626aeb459b
# ╟─c5a3348c-1534-4292-aa90-aa6ea3bc0b6a
# ╠═c933d42d-8ab4-49ae-b505-b6ed4fbcc26e
# ╠═b26ed39c-c7fd-4c72-b7ca-ec2737069cdc
# ╠═0d5352e9-ee3f-47d7-af14-1047aa6266c4
# ╠═7c348e71-22e6-4f2b-98e5-73ab5b29c935
# ╠═b5592863-a581-4a75-8959-e0d0278b6ada
# ╠═659cf583-5d3d-43b3-9cf2-8878abbace85
# ╠═4f757570-60c0-46fb-b152-a6f793fd077a
# ╠═e8b355e3-19c0-4006-8687-37bc761fa157
# ╠═f285b189-a0aa-45dd-8bfc-3e309c0cfcb4
# ╠═c2155af3-2479-462c-8348-f82975ba38d0
# ╠═9748ed29-2425-47f4-bc3a-10797fb89b92
# ╠═9b76b722-65ce-4755-9cf8-1bd46f3051c0
# ╠═109645f8-2132-4cd7-beab-30ae0b0e46c0
# ╠═335d8235-55b7-4ce8-a92d-7326ff8bfb6c
# ╠═111a3c9f-fcff-4b7a-9ff1-cf33303e1ec6
# ╠═b89fca52-0fe5-4f8a-956c-22dd6d0e3415
# ╠═405bf3ee-b122-4aee-9578-b151d542041a
# ╠═bdcdf051-5d35-4408-8406-3104ef52d52c
# ╠═d44fc394-9607-4a56-be18-05ecb01686a5
# ╟─7f2124ae-bc43-45a8-b8b5-3fe776261a71
