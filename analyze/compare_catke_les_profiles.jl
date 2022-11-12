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

using SingleColumnModelCalibration:
    lesbrary_inverse_problem,
    get_best_parameters,
    dependent_parameter_sets,
    prior_library,
    ecco_vertical_grid,
    load_summaries

set_theme!(Theme(fontsize=16))

# dir = "." #calibration_results/round1"
# name = :ri_based
# closure = RiBasedVerticalDiffusivity()

dir = "../results"
include("results.jl")
name = :complex_dissipation_conv_adj
closure = CATKEVerticalDiffusivity()

filename = joinpath(dir, calibration_filenames[name])
summaries = load_summaries(filename)
dependent_parameters = dependent_parameter_sets[string(name)]
#dependent_parameters = NamedTuple()

mean_θ₀ = summaries[0].ensemble_mean
mean_θₙ = summaries[end].ensemble_mean
best_θₙ = get_best_parameters(summaries)

@show summaries[end]

parameter_names = keys(best_θₙ)
free_parameters = FreeParameters(prior_library; names=parameter_names, dependent_parameters)

mean_θ₀ = build_parameters_named_tuple(free_parameters, mean_θ₀)
mean_θₙ = build_parameters_named_tuple(free_parameters, mean_θₙ)
best_θₙ = build_parameters_named_tuple(free_parameters, best_θₙ)

savename = string("catke_", name, "_parameters.jld2")
@save savename mean=mean_θₙ best=best_θₙ

# Two grids: "coarse" with ECCO vertical resolution to z=-256 m, and a fine grid with 4m resolution
Nz_ecco = length(ecco_vertical_grid) - 1
#coarse_regrid = RectilinearGrid(size=Nz_ecco, z=ecco_vertical_grid, topology=(Flat, Flat, Bounded))
coarse_regrid = RectilinearGrid(size=32, z=(-256, 0), topology=(Flat, Flat, Bounded))
med_regrid    = RectilinearGrid(size=48, z=(-256, 0), topology=(Flat, Flat, Bounded))
fine_regrid   = RectilinearGrid(size=64; z=(-256, 0), topology=(Flat, Flat, Bounded))

# Batch the inverse problems
Nensemble = 3
architecture = CPU()

suite = "48_hour_suite"
times = [2hours, 48hours]

#times = [2hours, 24hours]
#suite = "24_hour_suite"

#times = [2hours, 12hours]
#suite = "12_hour_suite"

inverse_problem_kwargs = (; suite, free_parameters, Nensemble, architecture, closure)
coarse_ip = lesbrary_inverse_problem(coarse_regrid; times, Δt=20minutes, inverse_problem_kwargs...)
med_ip    = lesbrary_inverse_problem(med_regrid; times,    Δt=20minutes, inverse_problem_kwargs...)
fine_ip   = lesbrary_inverse_problem(fine_regrid; times,   Δt=20minutes, inverse_problem_kwargs...)

forward_run!(fine_ip, [mean_θ₀, mean_θₙ, best_θₙ])
forward_run!(med_ip, [mean_θ₀, mean_θₙ, best_θₙ])
forward_run!(coarse_ip, [mean_θ₀, mean_θₙ, best_θₙ])

@show coarse_ip.simulation.model.closure[3]

#####
##### Figure
#####

fig = Figure(resolution=(1200, 600))

z_fine = znodes(Center, fine_regrid)
z_coarse = znodes(Center, coarse_regrid)
z_med = znodes(Center, med_regrid)
times = observation_times(fine_ip.observations)
Nt = length(times)

ax_b = []
ax_u = []

titles = [
    "Free convection",
    "Weak wind \n strong cooling",
    "Medium wind \n medium cooling",
    "Strong wind \n weak cooling",
    "Strong wind \n no cooling",
    "Strong wind \n no rotation",
]

sim_color    = (:seagreen, 0.6)
fine_color   = (:darkred, 0.8)
med_color   = (:orange, 0.8)
coarse_color = (:royalblue1, 0.8)
k_plot       = 3
zlim         = -192

lesstr = string("LES at t = ", prettytime(times[end]))

Δz_coarse = Δzᶜᶜᶜ(1, 1, coarse_regrid.Nz, coarse_regrid)
Δz_coarse_str = @sprintf("%.1f", Δz_coarse)
Δz_fine = Δzᶜᶜᶜ(1, 1, fine_regrid.Nz, fine_regrid)
Δz_fine_str = @sprintf("%.1f", Δz_fine)
Δz_med = Δzᶜᶜᶜ(1, 1, med_regrid.Nz, med_regrid)
Δz_med_str = @sprintf("%.1f", Δz_med)

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

    b_init   = interior(fine_ip.observations[c].field_time_serieses.b[1], 1, 1, :)
    b_obs    = interior(fine_ip.observations[c].field_time_serieses.b[Nt], 1, 1, :)
    b_fine   = interior(  fine_ip.time_series_collector.field_time_serieses.b[Nt], k_plot, c, :)
    b_med    = interior(   med_ip.time_series_collector.field_time_serieses.b[Nt], k_plot, c, :)
    b_coarse = interior(coarse_ip.time_series_collector.field_time_serieses.b[Nt], k_plot, c, :)

    lines!(ax_b[c], b_init,   z_fine,   linewidth=2, label="Initial condition at t = 2 hours", color=sim_color, linestyle=:dot)
    lines!(ax_b[c], b_obs,    z_fine,   linewidth=8, label=lesstr, color=sim_color)
    lines!(ax_b[c], b_fine,   z_fine,   linewidth=3, label="CATKE, Δz ≈ $Δz_fine_str m", color=fine_color)
    lines!(ax_b[c], b_med,    z_med,    linewidth=3, label="CATKE, Δz = $Δz_med_str m", color=med_color)
    lines!(ax_b[c], b_coarse, z_coarse, linewidth=3, label="CATKE, Δz = $Δz_coarse_str m", color=coarse_color)

    xlims!(ax_b[c], 0.0386, maximum(b_init) + 2e-5)
    ylims!(ax_b[c], zlim, 0)

    if c > 1
        #xlims!(ax_u[c], -0.15, 0.35)
        ylims!(ax_u[c], zlim, 0)

        u_obs    = interior(fine_ip.observations[c].field_time_serieses.u[Nt], 1, 1, :)
        u_fine   = interior(  fine_ip.time_series_collector.field_time_serieses.u[Nt], k_plot, c, :)
        u_med    = interior(  med_ip.time_series_collector.field_time_serieses.u[Nt], k_plot, c, :)
        u_coarse = interior(coarse_ip.time_series_collector.field_time_serieses.u[Nt], k_plot, c, :)

        lines!(ax_u[c], u_obs,    z_fine,   linewidth=8, label="u, LES",                          color=sim_color)
        lines!(ax_u[c], u_fine,   z_fine,   linewidth=3, label="u, CATKE, Δz ≈ $Δz_fine_str m",   color=fine_color)
        lines!(ax_u[c], u_med,    z_med,    linewidth=3, label="u, CATKE, Δz ≈ $Δz_med_str m",    color=med_color)
        lines!(ax_u[c], u_coarse, z_coarse, linewidth=3, label="u, CATKE, Δz = $Δz_coarse_str m", color=coarse_color)

        if :v ∈ keys(fine_ip.observations[c].field_time_serieses)
            v_obs    = interior(fine_ip.observations[c].field_time_serieses.v[Nt], 1, 1, :)
            v_fine   = interior(  fine_ip.time_series_collector.field_time_serieses.v[Nt], k_plot, c, :)
            v_med    = interior(   med_ip.time_series_collector.field_time_serieses.v[Nt], k_plot, c, :)
            v_coarse = interior(coarse_ip.time_series_collector.field_time_serieses.v[Nt], k_plot, c, :)

            lines!(ax_u[c], v_obs,    z_fine,   color=sim_color,    linewidth=4, linestyle=:dash) #, label="v, LES")
            lines!(ax_u[c], v_fine,   z_fine,   color=fine_color,   linewidth=2, linestyle=:dash, label="v") #, fine resolution CATKE")
            lines!(ax_u[c], v_med,    z_med,    color=med_color,    linewidth=2, linestyle=:dash) #, med resolution CATKE")
            lines!(ax_u[c], v_coarse, z_coarse, color=coarse_color, linewidth=2, linestyle=:dash) #, label="v, coarse resolution CATKE")
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

save("$(name)_catke_LES_comparison.png", fig)
