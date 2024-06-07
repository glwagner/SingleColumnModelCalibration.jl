using Oceananigans
using Oceananigans.Operators: Δzᶜᶜᶜ
using Oceananigans.Units
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    CATKEVerticalDiffusivity,
    MixingLength,
    TurbulentKineticEnergyEquation

using ParameterEstimocean
using ParameterEstimocean.InverseProblems: BatchedInverseProblem, inverting_forward_map
using ParameterEstimocean.PseudoSteppingSchemes: Kovachki2018InitialConvergenceRatio
using ParameterEstimocean.Parameters: build_parameters_named_tuple

using Printf
using JLD2
using NCDatasets
using LinearAlgebra
using CairoMakie

using SingleColumnModelCalibration:
    dependent_parameter_sets,
    build_batched_inverse_problem,
    prior_library

set_theme!(Theme(fontsize=22))

closure_label = "CATKE"

@load "optimal_catke.jld2" optimal_catke
closure = optimal_catke

cases = ["free_convection",
         "weak_wind_strong_cooling",
         "med_wind_med_cooling",
         "strong_wind_weak_cooling",
         "strong_wind",
         "strong_wind_no_rotation",
         "strong_wind_and_sunny"]

titles = [
          "6 hour simulation \n (extreme forcing)",
          "12 hour simulation \n (strong forcing)",
          "24 hour simulation \n (medium forcing)",
          "48 hour simulation \n (weak forcing)",
          "72 hour simulation \n (very weak forcing)",
]

labels = Dict(
    "SMCLT" => "SMC-LT (Harcourt 2015)",
    "EPBL-LT" => "ePBL (Reichl and Li 2019)",
    "KPP-CVMix" => "KPP (Large et al. 1994)"
)

grid_parameters = [(size=128, z=(-256, 0))]
k₀ = 28

#for cc = 1:7
cc = 7

#####
##### Figure
#####

fig = Figure(size=(1200, 500))

d = 4
dash = Linestyle([0.0, d, 1.6d, 2.6d])
ax_c = []

suites = [6, 12, 24, 48, 72]
resolutions = ["1m" for c = 1:5]

zlims = [-200 for c = 1:7]
zlims[1] = -150
zlims[2] = -160
zlims[3] = -165
zlims[4] = -180
zlims[5] = -180
zlims[6] = -150
zlims[7] = -150

zlim = zlims[cc]

for c = 1:5
    resolution = resolutions[c]
    suitehrs = suites[c]

    # Batch the inverse problems
    suite = "$(suitehrs)_hour_suite"
    resolution = "1m"
    #catkecolor = (:royalblue1, 0.8)
    catkecolor = (:black, 0.9)
    LEScolor = (:seagreen, 0.6)
    suite_parameters = (; name=suite, resolution, stop_time=suitehrs*hours)

    batched_ip = build_batched_inverse_problem(closure,
                                               #"variable_Pr_conv_adj";
                                               "extended_stability_conv_adj";
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

    case = cases[cc]
    yaxisposition = c < 5 ? :left : :right
    Label(fig[2, c], titles[c], tellwidth=false)

    ax_cc = Axis(fig[3, c]; ylabel="z (m)", xlabel="Passive \n tracer", yaxisposition, xticks=[-2, 0, 2, 4, 6, 8])
    push!(ax_c, ax_cc)
    ylims!(ax_c[c], zlim, 5)

    c_obs = interior(ip.observations[cc].field_time_serieses.c[Nt], 1, 1, :)
    lines!(ax_c[c], c_obs, z, linewidth=8, label=LES_str, color=LEScolor)

    c_catke = interior(ip.time_series_collector.field_time_serieses.c[Nt], 1, cc, :)
    lines!(ax_c[c], c_catke, z, linewidth=2; label="CATKE", color=catkecolor)

    hidespines!(ax_c[c], :t)
    c != 1 && hidespines!(ax_c[c], :l)
    c != 1 && c != 5 && hideydecorations!(ax_c[c], grid=false)
    c != 5 && hidespines!(ax_c[c], :r)

    Legend(fig[1, 1:5], ax_c[1], nbanks=5, framevisible=false)
end

if cc == 1 # free_convection
    [xlims!(ax, -2.5, 5.5) for ax in ax_c]
elseif cc == 2 # weak wind, strong cooling
    [xlims!(ax, -2.5, 4.5) for ax in ax_c]
elseif cc == 3 # weak wind, strong cooling
    [xlims!(ax, -2.5, 3.8) for ax in ax_c]
elseif cc == 4 # weak wind, strong cooling
    [xlims!(ax, -2.5, 3.8) for ax in ax_c]
elseif cc == 5 # weak wind, strong cooling
    [xlims!(ax, -2.5, 4.8) for ax in ax_c]
elseif cc == 6 # weak wind, strong cooling
    [xlims!(ax, -2.5, 6.2) for ax in ax_c]
elseif cc == 7 # weak wind, strong cooling
    [xlims!(ax, -2.5, 9.0) for ax in ax_c]
end

xtxt = -2
ztxt = -3
txts = ["(a)", "(b)", "(c)", "(d)", "(e)"]
for a = 1:5
    ax = ax_c[a]
    txt = txts[a]
    text!(ax, xtxt, ztxt, text=txt, align=(:left, :top))
end

rowsize!(fig.layout, 1, Relative(0.1))

display(fig)
case = cases[cc]
save("les_catke_forced_tracer_$case.pdf", fig)

#end

