using Oceananigans
using Oceananigans.Operators: Δzᶜᶜᶜ
using Oceananigans.Units
#using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities:
    CATKEVerticalDiffusivity,
    TKEDissipationVerticalDiffusivity

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

@load "optimal_catke.jld2" optimal_catke
closure = optimal_catke
#closure = TKEDissipationVerticalDiffusivity()

cases = ["free_convection",
         "weak_wind_strong_cooling",
         "med_wind_med_cooling",
         "strong_wind_weak_cooling",
         "strong_wind",
         "strong_wind_no_rotation"]

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
#for cc = 1:6
cc = 2

#####
##### Figure
#####

fig = Figure(size=(1200, 500))

d = 4
dash = Linestyle([0.0, d, 1.6d, 2.6d])
ax_b = []

suites = [6, 12, 24, 48, 72]
resolutions = ["1m" for c = 1:5]

zlims = [-200 for c = 1:6]
zlims[1] = -150
zlims[2] = -160
zlims[3] = -165
zlims[4] = -180
zlims[5] = -180
zlims[6] = -150

xticks = [([1e-4, 2e-4, 3e-4, 4e-4, 5e-4], ["1", "2", "3", "4", "5"]) for c = 1:6]

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
                                               #"extended_stability_conv_adj";
                                               "variable_stabilities";
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

    ax_bc = Axis(fig[3, c]; ylabel="z (m)", xlabel="Buoyancy \n (10⁻⁴ × m s⁻²)", yaxisposition, xticks=xticks[cc])
    push!(ax_b, ax_bc)
    ylims!(ax_b[c], zlim, 5)

    b_init = interior(ip.observations[cc].field_time_serieses.b[1], 1, 1, :)
    b_obs  = interior(ip.observations[cc].field_time_serieses.b[Nt], 1, 1, :)
    b₀ = b_init[k₀]

    #lines!(ax_b[c], b_init .- b₀, z, linewidth=2, label="Initial condition", color=LEScolor)
    lines!(ax_b[c], b_obs  .- b₀, z, linewidth=8, label=LES_str, color=LEScolor)

    #for parameterization_str in ["SMCLT", "EPBL-LT", "KPP-CVMix"]
    gotm_closures = ["SMCLT", "KPP-CVMix"]
    gotm_tsteps = [1, 120]
    for (parameterization_str, tstep) in zip(gotm_closures, gotm_tsteps)
        basedir = "/Users/gregorywagner/Projects/SingleColumnModelCalibration/comparison/gotm-new"
        casedir = string(case, "_", suitehrs, "hours_model_", parameterization_str, "_tstep", tstep)
        filename = "gotm_output_file.nc"
        filepath = joinpath(basedir, casedir, filename)
        ds = Dataset(filepath)
        α = 2e-4
        g = 9.81
        T_gotm = ds["temp"][1, 1, :, :]
        b_gotm = ds["buoy"][1, 1, :, :]
        z_gotm = ds["z"][1, 1, :, 1]
        label = labels[parameterization_str]
        lines!(ax_b[c], b_gotm[:, end] .- b_gotm[k₀, 1], z_gotm; linestyle=dash, linewidth=2, label)
    end

    b = interior(ip.time_series_collector.field_time_serieses.b[Nt], 1, cc, :)
    #lines!(ax_b[c], b .- b₀, z, linewidth=2; label="CATKE", color=catkecolor)
    lines!(ax_b[c], b .- b₀, z, linewidth=2; label="k-ϵ", color=catkecolor)

    hidespines!(ax_b[c], :t)
    c != 1 && hidespines!(ax_b[c], :l)
    c != 1 && c != 5 && hideydecorations!(ax_b[c], grid=false)
    c != 5 && hidespines!(ax_b[c], :r)

    Legend(fig[1, 1:5], ax_b[1], nbanks=5, framevisible=false)
end

if cc == 1 # free_convection
    xlims!(ax_b[1], 9e-5, 2.3e-4)
    xlims!(ax_b[2], 9e-5, 2.3e-4)
    xlims!(ax_b[3], 9e-5, 2.3e-4)
    xlims!(ax_b[4], 9e-5, 2.2e-4)
    xlims!(ax_b[5], 9e-5, 2.1e-4)
    xtxt, ztxt = 1e-4, -10
elseif cc == 2 # weak wind, strong cooling
    xlims!(ax_b[1], 9e-5, 2.7e-4)
    xlims!(ax_b[2], 9e-5, 2.6e-4)
    xlims!(ax_b[3], 9e-5, 2.6e-4)
    xlims!(ax_b[4], 9e-5, 2.5e-4)
    xlims!(ax_b[5], 9e-5, 2.3e-4)
    xtxt, ztxt = 1e-4, -5
elseif cc == 3 # medium wind, medium cooling
    xlims!(ax_b[1], 8e-5, 3.1e-4)
    xlims!(ax_b[2], 8e-5, 3.0e-4)
    xlims!(ax_b[3], 8e-5, 3.0e-4)
    xlims!(ax_b[4], 8e-5, 3.0e-4)
    xlims!(ax_b[5], 8e-5, 2.8e-4)
    xtxt, ztxt = 9e-5, -5
elseif cc == 4 # strong wind, weak cooling
    xlims!(ax_b[1], 3e-5, 3.6e-4)
    xlims!(ax_b[2], 3e-5, 3.6e-4)
    xlims!(ax_b[3], 3e-5, 3.6e-4)
    xlims!(ax_b[4], 3e-5, 3.6e-4)
    xlims!(ax_b[5], 3e-5, 3.6e-4)
    xtxt, ztxt = 4e-5, -3
elseif cc == 5 # strong wind
    xlims!(ax_b[1], 3e-5, 5.0e-4)
    xlims!(ax_b[2], 3e-5, 4.8e-4)
    xlims!(ax_b[3], 3e-5, 4.6e-4)
    xlims!(ax_b[4], 3e-5, 4.8e-4)
    xlims!(ax_b[5], 3e-5, 4.8e-4)
    xtxt, ztxt = 4e-5, -3
elseif cc == 6 # strong wind, no rotation
    xlims!(ax_b[1], 9e-5, 5.1e-4)
    xlims!(ax_b[2], 9e-5, 4.8e-4)
    xlims!(ax_b[3], 9e-5, 4.4e-4)
    xlims!(ax_b[4], 9e-5, 4.4e-4)
    xlims!(ax_b[5], 9e-5, 4.3e-4)
    xtxt, ztxt = 1e-4, -3
end

txts = ["(a)", "(b)", "(c)", "(d)", "(e)"]
for a = 1:5
    ax = ax_b[a]
    txt = txts[a]
    text!(ax, xtxt, ztxt, text=txt, align=(:left, :top))
end

rowsize!(fig.layout, 1, Relative(0.1))

display(fig)
case = cases[cc]
save("gotm_catke_forcing_comparison_$case.pdf", fig)
#save("gotm_k_epsilon_forcing_comparison_$case.pdf", fig)

#end

