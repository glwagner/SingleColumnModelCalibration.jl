using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using JLD2
using Printf
using Statistics
using GLMakie

using ParameterEstimocean.Parameters: transform_to_constrained
using SingleColumnModelCalibration: calibrate_parameter_set, build_batched_inverse_problem
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity

set_theme!(Theme(fontsize=24))

closure = CATKEVerticalDiffusivity()
name = "constant_Pr_no_shear"
suffix = "Nens400_Δt1200_τ10000_Nz32_Nz64_12_hour_suite_24_hour_suite_48_hour_suite.jld2"
dataset_filename = "calibration_summary_" * suffix
@load dataset_filename dataset

# Find best run
data = dataset[name]

Φ★, best_run = findmin(data[:final_minimum_objectives])
Φ_best              = data[:mean_parameter_objective_serieses][best_run]
optimal_parameters  = data[:final_best_parameters][best_run]
iteration_summaries = data[:iteration_summaries][best_run]

mean_parameters = parent([summ.ensemble_mean for summ in iteration_summaries])
Niters = length(mean_parameters)
Nsample = 10
Nens = max(Nsample, Niters)

#####
##### Generate ensemble mean data
#####

grid_parameters = [
    #(size=32, z=(-256, 0)),
    (size=64, z=(-256, 0))
]

suite_parameters = [
    #(name = "12_hour_suite", stop_time=12hours),
    #(name = "24_hour_suite", stop_time=24hours),
    #(name = "48_hour_suite", stop_time=48hours),
    (name = "72_hour_suite", stop_time=48hours),
]

ip = build_batched_inverse_problem(closure, name;
                                   Nensemble = Nens,
                                   grid_parameters,
                                   suite_parameters)

forward_run!(ip, mean_parameters)

pseudotimes = parent([summ.pseudotime for summ in iteration_summaries])

fig = Figure()
axh = Axis(fig[1, 2], xlabel="Pseudotime", ylabel="z (m)", xscale=log10)
axb = Axis(fig[1, 1], xlabel="Buoyancy (m s⁻²)", ylabel="z (m)")

c = 5 # case
m = 1 # batch member

b_obs = interior(ip[m].observations[c].field_time_serieses.b[end], 1:1, 1, :) 
b_catke = interior(ip[m].simulation.model.tracers.b, 1:Niters, c, :)
b′ = b_catke .- b_obs

z = znodes(Center, ip[m].simulation.model.grid)
blim = 4 * maximum(abs, b′[Niters, :])

hm = heatmap!(axh, pseudotimes[2:end], z, b′[2:end, :], colormap=:balance, colorrange=(-blim, blim))
Colorbar(fig[1, 3], hm, label="b_CATKE - b_obs")

lines!(axb, b_obs[1, :], z, linewidth=10, color=(:black, 0.4))

selected_iters = [10, 20, 40, Niters-1]

colors = [
    (:red,        0.05),
    (:seagreen,   0.05),
    (:black,      0.05),
    (:royalblue1, 0.05),
]

for (n, iter) in enumerate(selected_iters)
    summary = iteration_summaries[iter]
    μ = mean(summary.unconstrained_parameters, dims=2)[:]
    Σ = cov(summary.unconstrained_parameters, dims=2)
    X = rand(MvNormal(μ, Σ), Nsample)
    C = transform_to_constrained(priors, X)
    forward_run!(ip, C)
    b_catke = collect(interior(ip[m].simulation.model.tracers.b, 1:Nsample, c, :))
    b_sorted = sort(b_catke, dims=1) # (x, z)

    for α = 1:Nsample
        lines!(axb, b_catke[α, :],  z, linewidth=2,  color=colors[n])
    end
end

vlines!(axh, pseudotimes[selected_iters])

display(fig)

