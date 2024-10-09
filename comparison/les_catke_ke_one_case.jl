using Oceananigans
using Oceananigans.Operators: Δzᶜᶜᶜ
using Oceananigans.Units
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities:
    CATKEVerticalDiffusivity,
    TKEDissipationVerticalDiffusivity,
    TKEDissipationEquations

using ParameterEstimocean
using ParameterEstimocean.InverseProblems: BatchedInverseProblem, inverting_forward_map
using ParameterEstimocean.PseudoSteppingSchemes: Kovachki2018InitialConvergenceRatio
using ParameterEstimocean.Parameters: build_parameters_named_tuple

using Printf
using JLD2
using NCDatasets
using LinearAlgebra
using GLMakie

using SingleColumnModelCalibration:
    dependent_parameter_sets,
    build_batched_inverse_problem,
    prior_library

set_theme!(Theme(fontsize=22))

catke = CATKEVerticalDiffusivity()
tke_dissipation_equations = TKEDissipationEquations(Cᵂu★ = 3.2)
k_epsilon = TKEDissipationVerticalDiffusivity(; tke_dissipation_equations)
resolution = "1m"
kecolor = (:royalblue1, 0.8)
catkecolor = (:black, 0.9)
LEScolor = (:seagreen, 0.6)

cases = ["free_convection",
         "weak_wind_strong_cooling",
         "med_wind_med_cooling",
         "strong_wind_weak_cooling",
         "strong_wind",
         "strong_wind_no_rotation"]
         #"strong_wind_and_sunny"]

titles = [
          "6 hour simulation (extreme forcing)",
          "12 hour simulation (strong forcing)",
          "24 hour simulation (medium forcing)",
          "48 hour simulation (weak forcing)",
          "72 hour simulation (very weak forcing)",
]

labels = Dict(
    "SMCLT" => "SMC-LT (Harcourt 2015)",
    "EPBL-LT" => "ePBL (Reichl and Li 2019)",
    "KPP-CVMix" => "KPP (Large et al. 1994)"
)

grid_parameters = [(size=128, z=(-256, 0))]
k₀ = 28
cc = 5
case = cases[cc]

#####
##### Figure
#####

fig = Figure(size=(800, 500))

d = 4
dash = Linestyle([0.0, d, 1.6d, 2.6d])
ax_b = []
ax_e = []

suites = [6, 12, 24, 48, 72]
resolutions = ["1m" for c = 1:7]

zlims = [-200 for c = 1:7]
zlims[1] = -150
zlims[2] = -160
zlims[3] = -165
zlims[4] = -180
zlims[5] = -180
zlims[6] = -150
zlims[7] = -150

xticks = [([1e-4, 2e-4, 3e-4, 4e-4, 5e-4], ["1", "2", "3", "4", "5"]) for c = 1:7]

zlim = zlims[cc]

c = 4
resolution = resolutions[c]
suitehrs = suites[c]

# Batch the inverse problems
suite = "$(suitehrs)_hour_suite"
suite_parameters = (; name=suite, resolution, stop_time=suitehrs*hours)

batched_ip = build_batched_inverse_problem(catke,
                                           "extended_stability_conv_adj";
                                           #"variable_stabilities";
                                           Nensemble = 1,
                                           Δt = 1minutes,
                                           start_time = 10minutes,
                                           grid_parameters,
                                           suite_parameters)

parameters = (; minimum_turbulent_kinetic_energy=1e-9)
forward_run!(batched_ip, [parameters])

ip = batched_ip[1]
times = observation_times(ip.observations)
Nt = length(times)
LES_str = string("Large eddy simulation")

grid = ip.simulation.model.grid
z = znodes(grid, Center())
Δz = Δzᶜᶜᶜ(1, 1, grid.Nz, grid)
Δz_str = @sprintf("%d", Δz)

Label(fig[2, 1:2], titles[c], tellwidth=false)
ax_b = Axis(fig[3, 1]; ylabel="z (m)", xlabel="Buoyancy \n (10⁻⁴ × m s⁻²)", yaxisposition=:left, xticks=xticks[cc])
ax_e = Axis(fig[3, 2]; ylabel="z (m)", xlabel="Turbulent kinetic \n energy (cm² s⁻²)", yaxisposition=:right)
σe = 1e4 # scaling for tke

ylims!(ax_b, zlim, 5)
ylims!(ax_e, zlim, 5)

b_init = interior(ip.observations[cc].field_time_serieses.b[1], 1, 1, :)
b_obs  = interior(ip.observations[cc].field_time_serieses.b[Nt], 1, 1, :)
e_obs  = interior(ip.observations[cc].field_time_serieses.e[Nt], 1, 1, :)
b₀ = b_init[k₀]

lines!(ax_b, b_obs .- b₀, z, linewidth=8, label=LES_str, color=LEScolor)
lines!(ax_e, σe .* e_obs, z, linewidth=8, label=LES_str, color=LEScolor)

b = interior(ip.time_series_collector.field_time_serieses.b[Nt], 1, cc, :)
e = interior(ip.simulation.model.tracers.e, 1, cc, :)
lines!(ax_b, b .- b₀, z, linewidth=2; label="CATKE", color=catkecolor)
lines!(ax_e, σe .* e, z, linewidth=2; label="CATKE", color=catkecolor)

batched_ip = build_batched_inverse_problem(k_epsilon,
                                           "variable_stabilities";
                                           Nensemble = 1,
                                           Δt = 1minutes,
                                           start_time = 10minutes,
                                           grid_parameters,
                                           suite_parameters)

parameters = NamedTuple()
forward_run!(batched_ip, [parameters])
ip = batched_ip[1]
b = interior(ip.time_series_collector.field_time_serieses.b[Nt], 1, cc, :)
e = interior(ip.simulation.model.tracers.e, 1, cc, :)
lines!(ax_b, b .- b₀, z, linewidth=2; label="k-ϵ (Umlauf and Buchard 2003)", color=kecolor)
lines!(ax_e, σe .* e, z, linewidth=2; label="k-ϵ (Umlauf and Buchard 2003)", color=kecolor)

hidespines!(ax_b, :t)
hidespines!(ax_e, :t)
hidespines!(ax_e, :l)
hidespines!(ax_b, :r)

Legend(fig[1, 1:2], ax_b, nbanks=4, framevisible=false)

xlims!(ax_b, 3e-5, 5e-4)

rowsize!(fig.layout, 1, Relative(0.1))
save("gotm_catke_k_epsilon_tke_comparison_one_forcing_$case.png", fig)

display(fig)

