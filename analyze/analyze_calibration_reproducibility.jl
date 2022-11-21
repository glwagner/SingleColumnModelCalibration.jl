using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using SingleColumnModelCalibration: calibrate_parameter_set, finitefindmin
using CairoMakie
using ElectronDisplay
using JLD2
using Printf
using Statistics
using LinearAlgebra

import ParameterEstimocean.EnsembleKalmanInversions: eki_objective, UninitializedForwardMapOutput

using Oceananigans.TurbulenceClosures:
    RiBasedVerticalDiffusivity,
    CATKEVerticalDiffusivity

function eki_objective(y, G, Γ⁻¹²)
    Nens = size(G, 2)
    Φ = [norm(Γ⁻¹² * (y .- view(G, :, k)))^2 / 2 for k = 1:Nens]
    return Φ
end

function collect_calibration_data(name, suffix;
                                  closure = CATKEVerticalDiffusivity(),
                                  Nrepeats = 10,
                                  dir = joinpath("..", "results"))
    iteration_summaries = []
    final_minimum_objectives = Float64[]
    minimum_objective_serieses = []
    final_best_parameters = []
    final_mean_parameters = []
    mean_parameter_serieses = []

    for i = 1:Nrepeats
        filename = string(name, "_", i, "_", suffix)
        filepath = joinpath(dir, filename)
        file = jldopen(filepath)
        push!(iteration_summaries, file["iteration_summaries"])
        close(file)

        # Extract the last summary from each run
        last_summary = iteration_summaries[end][end]

        Φ★, k★ = finitefindmin(map(obj -> first(obj), last_summary.objective_values))
        push!(final_minimum_objectives, Φ★)
        push!(final_best_parameters, last_summary.parameters[k★])
        push!(final_mean_parameters, last_summary.ensemble_mean)

        push!(mean_parameter_serieses, [deepcopy(summ.ensemble_mean) for summ in iteration_summaries[end]])

        @printf "%.3e " Φ★
        for (k, p) in zip(keys(last_summary.parameters[k★]),
                          values(last_summary.parameters[k★]))

            @printf "%s: %.2e " k p
        end
        println(" ")

        run_objectives = zeros(length(iteration_summaries[end]))
        for (i, summary) in enumerate(iteration_summaries[end])
            min_iter_objective, k = finitefindmin(map(obj -> first(obj), summary.objective_values))
            run_objectives[i] = min_iter_objective
        end

        push!(minimum_objective_serieses, run_objectives)
    end

    Nens = 0
    for parameter_series in mean_parameter_serieses
        Nens = max(Nens, length(parameter_series))
    end

    #####
    ##### Generate ensemble mean data
    #####

    grid_parameters = [
        (size=32, z=(-256, 0)),
        (size=64, z=(-256, 0))
    ]

    suite_parameters = [
        (name = "12_hour_suite", stop_time=12hours),
        (name = "24_hour_suite", stop_time=24hours),
        (name = "48_hour_suite", stop_time=48hours),
    ]

    eki = calibrate_parameter_set(closure, name;
                                  Nensemble = Nens,
                                  stop_pseudotime = 0.0,
                                  stop_iteration = 0,
                                  resampler = nothing,
                                  plot_progress = false,
                                  forward_map_output = UninitializedForwardMapOutput(),
                                  grid_parameters,
                                  suite_parameters)

    Γ⁻¹² = eki.precomputed_arrays[:inv_sqrt_Γy]
    y = observation_map(eki.inverse_problem)
    Nobs = length(y)
    final_mean_parameter_objectives = Float64[]
    mean_parameter_objective_serieses = []

    for parameter_series in mean_parameter_serieses
        G = forward_map(eki.inverse_problem, parent(parameter_series))
        Nt = length(parameter_series)
        series_objectives = eki_objective(y, view(G, :, 1:Nt), Γ⁻¹²)
        push!(mean_parameter_objective_serieses, map(obj -> first(obj), series_objectives))
        push!(final_mean_parameter_objectives, first(series_objectives[end]))
    end

    #####
    ##### Human-readable summary
    #####

    names = keys(first(final_best_parameters))

    mean_best_parameters = []
    std_best_parameters = []
    mean_mean_parameters = []
    std_mean_parameters = []

    for name in names
        best_param = [final_best_parameters[i][name] for i = 1:Nrepeats]
        mean_param = [final_mean_parameters[i][name] for i = 1:Nrepeats]
        push!(mean_best_parameters, mean(best_param))
        push!(std_best_parameters, std(best_param))
        push!(mean_mean_parameters, mean(mean_param))
        push!(std_mean_parameters, std(mean_param))
    end

    bests_table = [collect(names) mean_best_parameters std_best_parameters std_best_parameters ./ mean_best_parameters]
    means_table = [collect(names) mean_mean_parameters std_mean_parameters std_mean_parameters ./ mean_mean_parameters]

    data = Dict()

    # Some metadata
    data[:name] = name
    data[:suffix] = suffix
    data[:bests_table] = bests_table
    data[:mean_table] = means_table

    # All the iteration summaries from all EKI iterations from every run
    data[:iteration_summaries] = iteration_summaries

    # Data from the end of each run
    data[:final_best_parameters] = final_best_parameters
    data[:final_minimum_objectives] = final_minimum_objectives
    data[:final_mean_parameters] = final_mean_parameters
    data[:final_mean_parameter_objectives] = final_mean_parameter_objectives
    data[:final_minimum_objectives] = final_minimum_objectives

    # Series for EKI evolution for every run
    data[:mean_parameter_serieses] = mean_parameter_serieses
    data[:minimum_objective_serieses] = minimum_objective_serieses
    data[:mean_parameter_objective_serieses] = mean_parameter_objective_serieses

    return data
end

name = "constant_Pr"
#name = "constant_Pr_conv_adj"
#name = "variable_Pr_conv_adj"
#suffix = "Nens1000_Δt1200_Nz32_Nz64_12_hour_suite_24_hour_suite_48_hour_suite.jld2"
suffix = "Nens100_Δt1200_Nz32_Nz64_12_hour_suite_24_hour_suite_48_hour_suite.jld2"
Nrepeats = 10

names = [
    "constant_Pr",
    #"variable_Pr",
    #"constant_Pr_conv_adj",
    #"variable_Pr_conv_adj"
]

suffixes = [
    "Nens100_Δt1200_τ100_Nz32_Nz64_12_hour_suite_24_hour_suite_48_hour_suite.jld2",
    "Nens100_Δt1200_Nz32_Nz64_12_hour_suite_24_hour_suite_48_hour_suite.jld2",
    "Nens100_Δt1200_Nz32_Nz64_12_hour_suite_24_hour_suite_48_hour_suite.jld2",
    "Nens100_Δt1200_τ10_Nz32_Nz64_12_hour_suite_24_hour_suite_48_hour_suite.jld2",
]

labels = [
    "Constant \n Pr",
    #"Variable \n Pr",
    #"Conv. adj. \n Constant \n Pr",
    #"Conv. Adj. \n Variable \n Pr",
]

vals = collect(1:length(names))

#=
dataset = []
for i = 1:length(names)
    run_data = collect_calibration_data(names[i], suffixes[i]; Nrepeats, dir="../results")
    push!(dataset, run_data)
end
=#

for (name, run_data) in zip(names, dataset)
    Φ★, best_run = findmin(run_data[:final_minimum_objectives])
    optimal_parameters = run_data[:final_best_parameters][best_run]

    savedir = joinpath("..", "parameters")
    savename = string(name, "_best_parameters.jld2")
    savepath = joinpath(savedir, savename)

    rm(savepath; force=true)
    file = jldopen(savepath, "a+")
    file["optimal_parameters"] = optimal_parameters
    file["final_best_parameters"] = run_data[:final_best_parameters]
    file["final_mean_parameters"] = run_data[:final_mean_parameters]
    close(file)
end

#####
##### Plot
#####

fig = Figure(resolution=(1200, 400))

ax1 = Axis(fig[1, 1], ylabel="EKI objective", xticks=(vals, labels))#, yscale=log10)
ax2 = Axis(fig[1, 2], xlabel="EKI iteration", ylabel="EKI objective", yscale=log10)

colors = [:orange1, :black, :seagreen, :royalblue]

for (n, d) in enumerate(dataset)
    name = names[n]
    color = (colors[n], 0.6)
    scatter!(ax1, n * ones(Nrepeats), d[:final_minimum_objectives]; color, marker=:star5, label=string(name, ", best"))
    scatter!(ax1, n * ones(Nrepeats), d[:final_mean_parameter_objectives]; color, marker=:circle, label=string(name, ", mean"))
    
    objective_scale = 1.0
    
    for i = 1:Nrepeats
        #@show d[:minimum_objective_serieses][i])
        @show d[:mean_parameter_objective_serieses][i]

        τ = [summary.pseudotime for summary in parent(d[:iteration_summaries][i])]
        color = (colors[n], 0.3)
        #ls = lines!(ax2, τ, d[:minimum_objective_serieses][i]; color)
        #ls = lines!(ax2, τ, d[:mean_parameter_objective_serieses][i][1:length(τ)]; color, linewidth=3)
        ls = lines!(ax2, d[:minimum_objective_serieses][i]; color)
        ls = lines!(ax2, d[:mean_parameter_objective_serieses][i]; color, linewidth=3)
    end
end

#Legend(fig[1, 3], ax1)

display(fig)

