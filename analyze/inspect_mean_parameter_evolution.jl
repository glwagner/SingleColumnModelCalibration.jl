using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using JLD2
using Printf
using Statistics
using GLMakie
using Distributions

using ParameterEstimocean.Parameters: transform_to_constrained
using SingleColumnModelCalibration: calibrate_parameter_set, build_batched_inverse_problem
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity

set_theme!(Theme(fontsize=36))

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
Nsample = 100
Nens = max(Nsample, Niters)

#####
##### Generate ensemble mean data
#####

grid_parameters = [
    #(size=32, z=(-256, 0)),
    (size=64, z=(-256, 0))
    #(size=128, z=(-256, 0))
]

suite_parameters = [
    (name = "6_hour_suite", resolution="0.75m", stop_time=6hours),
    #(name = "12_hour_suite", stop_time=12hours),
    #(name = "24_hour_suite", stop_time=24hours),
    #(name = "48_hour_suite", stop_time=48hours),
    (name = "72_hour_suite", stop_time=72hours),
]

ip = build_batched_inverse_problem(closure, name;
                                   Nensemble = Nens,
                                   Δt = 20minutes, #Δt = 1minutes,
                                   grid_parameters,
                                   suite_parameters)

forward_run!(ip, mean_parameters)
pseudotimes = parent([summ.pseudotime for summ in iteration_summaries])

T = [0.1, 1, 100]
selected_iters = [findfirst(τ -> τ >= τn, pseudotimes) for τn in T]

cases = [1, 5]
batch_members = [1, 2]
rows = 2

colors = [
    (:yellow,      0.3),
    (:red,         0.4),
    (:royalblue1,  0.8),
    (:black,       0.8),
]

fig = Figure(resolution=(1800, 1200))

for row = 1:rows
    c = cases[row]
    m = batch_members[row]

    axh = Axis(fig[row, 2], xlabel=L"Pseudotime, $\mathcal{T}$", ylabel=L"$z$ (m)", xscale=log10, yaxisposition=:right)
    hideydecorations!(axh)

    b_obs = interior(ip[m].observations[c].field_time_serieses.b[end], 1:1, 1, :) 
    b_catke = interior(ip[m].simulation.model.tracers.b, 1:Niters, c, :)
    Δb = b_catke .- b_obs

    z = znodes(Center, ip[m].simulation.model.grid)
    max_Δb = maximum(abs, Δb[Niters, :])

    hm = heatmap!(axh, pseudotimes[2:end], z, Δb[2:end, :] ./ max_Δb, colormap=:balance, colorrange=(-2, 2))
    Colorbar(fig[row, 3], hm, label=L"$\Delta b / \max[\Delta b(\mathcal{T}_f]$")

    for (n, iter) in enumerate(selected_iters)
        vlines!(axh, pseudotimes[iter], color=(first(colors[n]), 0.8), linewidth=6)
    end

    if row == 1
        hidexdecorations!(axh)
    end
end

for row = 1:rows
    c = cases[row]
    m = batch_members[row]

    axb = Axis(fig[row, 1], xlabel="Buoyancy (m s⁻²)", ylabel="z (m)")
    xlims!(axb, 0.0385, 0.0391)

    Nz = ip[m].simulation.model.grid.Nz
    z = znodes(Center, ip[m].simulation.model.grid)
    b_obs = interior(ip[m].observations[c].field_time_serieses.b[end], 1:1, 1, :) 

    lines!(axb, b_obs[1, :], z, linewidth=10, color=:black, label=L"\text{LES}")

    priors = ip.free_parameters.priors

    for (n, iter) in enumerate(selected_iters)
        summary = iteration_summaries[iter-1]
        pseudotime = summary.pseudotime
        μ = mean(summary.unconstrained_parameters, dims=2)[:]
        Σ = cov(summary.unconstrained_parameters, dims=2)
        X = rand(MvNormal(μ, Σ), Nsample)
        C = transform_to_constrained(priors, X)
        forward_run!(ip, C)
        b_catke = collect(interior(ip[m].simulation.model.tracers.b, 1:Nsample, c, :))

        min_b = minimum(b_catke, dims=1)
        max_b = maximum(b_catke, dims=1)

        left  = [Point2f(min_b[k], z[k]) for k = 1:Nz]
        right = [Point2f(max_b[k], z[k]) for k = 1:Nz]

        exponent = round(Int, log10(pseudotime) * 10) / 10
        tstr = exponent % 1 == 0 ? @sprintf("%d", 10^exponent) : @sprintf("%.1f", pseudotime)
        label = L"CATKE, $\mathcal{T} \approx %$tstr"
        band!(axb, left, right, color=colors[n]; label)
    end

    axislegend(axb, position=:rb)

    if row == 1
        hidexdecorations!(axb, grid=false)
    end
end

colsize!(fig.layout, 1, Relative(0.35))

display(fig)

