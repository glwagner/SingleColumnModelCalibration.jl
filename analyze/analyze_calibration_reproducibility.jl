using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using SingleColumnModelCalibration: calibrate_parameter_set, finitefindmin
using GLMakie
using JLD2
using Printf
using Statistics
using LinearAlgebra

import ParameterEstimocean.EnsembleKalmanInversions: eki_objective, UninitializedForwardMapOutput

using Oceananigans.TurbulenceClosures:
    RiBasedVerticalDiffusivity,
    CATKEVerticalDiffusivity

#suffix = "Nens1000_Δt300_τ10000_Nz32_Nz64_Nz128_12_hour_suite_24_hour_suite_48_hour_suite_nonsplit_tke_stepping_tight_priors.jld2"
#suffix = "Nens1000_Δt300_τ10000_Nz32_Nz64_Nz128_12_hour_suite_24_hour_suite_48_hour_suite_nonsplit_tke_stepping.jld2"
#suffix = "Nens1000_Δt300_τ10000_Nz32_Nz64_Nz128_12_hour_suite_24_hour_suite_48_hour_suite_nonsplit_tke_stepping_surface.jld2"
#suffix = "Nens100_Δt60_τ10000_Nz32_Nz64_12_hour_suite_24_hour_suite_72_hour_suite_new_sp_plus_scale_nonsplit_tke_stepping.jld2"
suffix = "Nens1000_Δt60_τ10000_Nz32_Nz64_Nz128_12_hour_suite_24_hour_suite_48_hour_suite_with_tracer.jld2"
Nrepeats = 1

dataset_filename = "calibration_summary_" * suffix

names = [
    "extended_stability_conv_adj",
    #"variable_Pr_conv_adj",
    #"ri_based",
]

closure = CATKEVerticalDiffusivity(turbulent_kinetic_energy_time_step=nothing)

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

        @printf "%.3e \n" Φ★
        for (k, p) in zip(keys(last_summary.parameters[k★]),
                          values(last_summary.parameters[k★]))

            @printf "%s: %.2e \n" k p
        end

        run_objectives = zeros(length(iteration_summaries[end]))
        for (i, summary) in enumerate(iteration_summaries[end])
            min_iter_objective, k = finitefindmin(map(obj -> first(obj), summary.objective_values))
            run_objectives[i] = min_iter_objective
        end

        push!(minimum_objective_serieses, run_objectives)
    end

    Nens = 0
    for parameter_series in mean_parameter_serieses
        Nens = min(max(Nens, length(parameter_series)), 100)
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

    @show eki

    Γ⁻¹² = eki.precomputed_arrays[:inv_sqrt_Γy]
    @show size(Γ⁻¹²)
    y = observation_map(eki.inverse_problem)
    Nobs = length(y)
    final_mean_parameter_objectives = Float64[]
    mean_parameter_objective_serieses = []

    for parameter_series in mean_parameter_serieses
        parameter_series = parent(parameter_series)[end-Nens+1:end]
        G = forward_map(eki.inverse_problem, parameter_series)
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

dataset = Dict()

for n = 1:length(names)
    name = names[n]
    data = collect_calibration_data(name, suffix; closure, Nrepeats, dir="../results")
    dataset[name] = data

    # Find best run
    Φ★, best_run = findmin(data[:final_minimum_objectives])
    Φþ, worst_run = findmax(data[:final_minimum_objectives])
    Φ_best = data[:mean_parameter_objective_serieses][best_run]
    iteration_summaries = data[:iteration_summaries][best_run]
    optimal_parameters = data[:final_best_parameters][best_run]

    @info name
    for (k, v) in zip(keys(optimal_parameters), values(optimal_parameters))
        @printf "% 6s %.5f \n" k v
    end

    best_first_particles = []
    best_first_objectives = []
    for r = 1:Nrepeats
        first_iteration_summary = first(data[:iteration_summaries][r])
        Φ₁, k₁ = findmin(o -> isnan(first(o)) ? Inf : first(o), first_iteration_summary.objective_values)
        θ = first_iteration_summary.parameters[k₁]
        push!(best_first_particles, θ)
        push!(best_first_objectives, Φ₁)
    end

    Φ₁★, r₁★ = findmin(best_first_objectives)
    θ₁★ = best_first_particles[r₁★]

    savedir = joinpath("..", "parameters")
    savename = string(name, "_best_parameters.jld2")
    savepath = joinpath(savedir, savename)

    rm(savepath; force=true)
    file = jldopen(savepath, "a+")
    file["optimal_parameters"] = optimal_parameters
    file["example_first_mean_parameters"] = first(data[:iteration_summaries][1]).ensemble_mean
    file["first_best_parameters"] = θ₁★
    file["final_best_parameters"] = data[:final_best_parameters]
    file["final_mean_parameters"] = data[:final_mean_parameters]
    file["first_mean_parameters"] = data[:final_mean_parameters]
    close(file)
end

# Save summary data
@save dataset_filename dataset

include("plot_calibration_reproducibility.jl")
include("compare_les_parameterization.jl")

