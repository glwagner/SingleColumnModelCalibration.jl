using Oceananigans
using Oceananigans.Operators: Δzᶜᶜᶜ
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity

using ParameterEstimocean
using ParameterEstimocean.InverseProblems: BatchedInverseProblem, inverting_forward_map
using ParameterEstimocean.PseudoSteppingSchemes: Kovachki2018InitialConvergenceRatio
using ParameterEstimocean.Parameters: build_parameters_named_tuple

using Printf
using JLD2
using LinearAlgebra

using CairoMakie
using ElectronDisplay

set_theme!(Theme(fontsize=16))

architecture = CPU()

Nhours = 72
resolution = "1m"
suite = "$(Nhours)_hour_suite"
stop_time = Nhours * hours
times = [2hours, stop_time]
Nt = length(times)
case_path(case) = joinpath("..", "data", suite, resolution, case * "_instantaneous_statistics.jld2")

observation_library = []

cases = [
    "free_convection",
    "weak_wind_strong_cooling",
    "med_wind_med_cooling",
    "strong_wind_weak_cooling",
    "strong_wind",
    "strong_wind_no_rotation",
]

titles = [
    "Free convection",
    "Weak wind \n strong cooling",
    "Medium wind \n medium cooling",
    "Strong wind \n weak cooling",
    "Strong wind \n no cooling",
    "Strong wind \n no rotation",
]

for case in cases
    if case == "free_convection"
        field_names = (:b, :e)
    elseif case == "strong_wind_no_rotation"
        field_names = (:b, :u, :e)
    else
        field_names = (:b, :u, :v, :e)
    end

    obs = SyntheticObservations(case_path(case); times, field_names)
    push!(observation_library, obs)
end

#####
##### Figure
#####

grid = observation_library[1].grid
z = znodes(Center, grid)

sim_color = (:seagreen, 0.6)
zlim      = -256
les_str   = string("LES at t = ", prettytime(times[end]))

fig = Figure(resolution=(1200, 600))

ax_b = []
ax_u = []

for c = 1:6
    if c < 5
        yaxisposition=:left
    else
        yaxisposition=:right
    end

    Label(fig[1, c], titles[c], tellwidth=false)

    ax_bc = Axis(fig[2, c]; ylabel="z (m)", xlabel="Buoyancy (m s⁻²)", yaxisposition, xticks=[0.0386, 0.039, 0.0394])
    push!(ax_b, ax_bc)

    ax_uc = c == 1 ? nothing : Axis(fig[3, c], ylabel="z (m)", xlabel="Velocities (m s⁻¹)"; yaxisposition, xticks=-0.1:0.1:0.3)
    push!(ax_u, ax_uc)

    b_init   = interior(observation_library[c].field_time_serieses.b[1], 1, 1, :)
    b_obs    = interior(observation_library[c].field_time_serieses.b[Nt], 1, 1, :)

    lines!(ax_b[c], b_init, z, linewidth=2, label="Initial condition at t = 2 hours", color=sim_color, linestyle=:dot)
    lines!(ax_b[c], b_obs,  z, linewidth=8, label=les_str, color=sim_color)

    xlims!(ax_b[c], 0.0386, maximum(b_init) + 2e-5)
    ylims!(ax_b[c], zlim, 0)

    if c > 1
        ylims!(ax_u[c], zlim, 0)

        u_obs    = interior(observation_library[c].field_time_serieses.u[Nt], 1, 1, :)
        lines!(ax_u[c], u_obs, z, linewidth=8, label="u, LES",                          color=sim_color)

        if :v ∈ keys(observation_library[c].field_time_serieses)
            v_obs    = interior(observation_library[c].field_time_serieses.v[Nt], 1, 1, :)
            lines!(ax_u[c], v_obs, z, color=sim_color, linewidth=4, linestyle=:dash) #, label="v, LES")
        end
    end

    hidespines!(ax_b[c], :t)
    c != 1 && hidespines!(ax_u[c], :t)

    c != 1 && hidespines!(ax_b[c], :l)
    c != 1 && c != 6 && hideydecorations!(ax_b[c], grid=false)
    c != 6 && hidespines!(ax_b[c], :r)

    c > 2 && hidespines!(ax_u[c], :l)
    c != 1 && c != 6 && hidespines!(ax_u[c], :r)
    c > 2 && c != 6 && hideydecorations!(ax_u[c], grid=false)
end

Legend(fig[3, 1], ax_b[2])
text!(ax_u[2], +0.1, -50.0, text="u")
text!(ax_u[2], -0.14, -110.0, text="v")
#Legend(fig[5, 1], ax_u[2])

display(fig)

