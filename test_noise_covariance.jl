using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using ParameterEstimocean.InverseProblems: BatchedInverseProblem
using ParameterEstimocean.PseudoSteppingSchemes: Kovachki2018InitialConvergenceRatio
using Printf
using JLD2
using Statistics: cov

# Import "lesbrary_inverse_problem" + plotting utility "calibration_progress_figure"
include("multi_resolution_calibration_utilities.jl")
include("parameter_sets.jl")

dir = "/home/greg/Projects/LocalOceanClosureCalibration/calibration_summaries"
name = "complex_conv_adj"

non_ensemble_closure = nothing # convective_adjustment

free_parameters = FreeParameters(prior_library,
                                 names = parameter_sets[name],
                                 dependent_parameters = dependent_parameter_sets[name])
@show free_parameters

# Two grids: "coarse" with ECCO vertical resolution to z=-256 m, and a fine grid with 4m resolution
Nz_ecco = length(ecco_vertical_grid) - 1
coarse_regrid = RectilinearGrid(size=Nz_ecco, z=ecco_vertical_grid, topology=(Flat, Flat, Bounded))
fine_regrid   = RectilinearGrid(size=48; z=(-256, 0), topology=(Flat, Flat, Bounded))

fine_Δt = 5minutes
coarse_Δt = 10minutes

# Batch the inverse problems
times = [6hours, 12hours, 18hours, 24hours]
tke_weight = 1.0
#obs_1m = batched_lesbrary_observations(coarse_regrid; times, tke_weight, resolution="1m")
#obs_2m = batched_lesbrary_observations(coarse_regrid; times, tke_weight, resolution="2m")
#obs_4m = batched_lesbrary_observations(coarse_regrid; times, tke_weight, resolution="4m")

obs_1m = batched_lesbrary_observations(coarse_regrid; times, tke_weight, resolution="1m")
obs_2m = batched_lesbrary_observations(coarse_regrid; times, tke_weight, resolution="2m")
obs_4m = batched_lesbrary_observations(coarse_regrid; times, tke_weight, resolution="4m")

fig = Figure(resolution=(1200, 1200))

for (c, case) in enumerate(cases)
    # Make axes
    label = replace(case, "_" => "\n")
    axs = make_axes(fig, c, label)

    # Plot observed data for each field
    case_obs = obs_4m[c]
    case_dataset = case_obs.field_time_serieses
    grid = case_obs.grid
    Nt = length(times)
    case_names = keys(case_dataset)
    initial_case_field_data = NamedTuple(n => interior(getproperty(case_dataset, n)[1])[1, 1, :] for n in case_names)
    final_case_field_data = NamedTuple(n => interior(getproperty(case_dataset, n)[Nt])[1, 1, :] for n in case_names)
    plot_fields!(axs, "Observed at t = " * prettytime(times[1]), (:black, 0.6), grid, initial_case_field_data...;
                 linewidth=2, linestyle=:dash)
    plot_fields!(axs, "Observed at t = " * prettytime(times[Nt]), (:gray23, 0.4), grid, final_case_field_data...; linewidth=6)
end

display(fig)

Γ = cov([obs_1m, obs_2m, obs_4m])
d = [Γ[n, n] for n=1:size(Γ, 1)]
ϵ = 1e-2 * mean(abs, d)
Γ .+= ϵ * Diagonal(randn(size(Γ, 1)))
fig = Figure(resolution=(1200, 1200))
ax = Axis(fig[1, 1])
heatmap!(ax, Γ, colormap=:balance, colorrange=(-0.01, 0.01))
#heatmap!(ax, inv(Γ), colormap=:balance, colorrange=(-0.01, 0.01))
display(fig)

inv_Γ = inv(Γ)
@show maximum(inv_Γ) minimum(inv_Γ)

lines([Γ[n, n] for n=1:size(Γ, 1)])
