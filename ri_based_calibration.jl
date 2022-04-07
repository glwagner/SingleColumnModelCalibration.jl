using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using LinearAlgebra
using GLMakie
using DataDeps
using Distributions

using Oceananigans.TurbulenceClosures: RiBasedVerticalDiffusivity

include("utils.jl")

#####
##### Setup
#####

cases = ["free_convection",
         "strong_wind_weak_cooling",
         "weak_wind_strong_cooling",
         "strong_wind",
         "strong_wind_no_rotation"]

datapaths = [@datadep_str("two_day_suite_1m/$(case)_instantaneous_statistics.jld2") for case in cases]

times = [2hours, 12hours, 36hours]
field_names = (:b, :u, :v)
transformation = ZScore()

z = [-256.0,
     -238.3,
     -207.2,
     -182.3,
     -162.5,
     -146.5,
     -133.0,
     -121.3,
     -110.5,
     -100.2,
     - 90.0,
     - 80.0,
     - 70.0,
     - 60.0,
     - 50.0,
     - 40.0,
     - 30.0,
     - 20.0,
     - 10.0,
        0.0]

regrid = RectilinearGrid(size=length(z)-1; z, topology=(Flat, Flat, Bounded))

observations = [SyntheticObservations(path; field_names, times, transformation, regrid)
                for path in datapaths]

colorcycle = [:black, :red, :darkblue, :orange, :pink1, :seagreen, :magenta2]
markercycle = [:rect, :utriangle, :star5, :circle, :cross, :+, :pentagon]

function make_figure_axes(n=1)
    fig = Figure(resolution=(1200, n*400))
    axs = []
    for i = 1:n
        ax_b = Axis(fig[i, 1], xlabel = "Buoyancy \n[cm s⁻²]", ylabel = "z [m]")
        ax_u = Axis(fig[i, 2], xlabel = "x-velocity, u \n[cm s⁻¹]")
        ax_v = Axis(fig[i, 3], xlabel = "y-velocity, v \n[cm s⁻¹]")
        push!(axs, (ax_b, ax_u, ax_v))
    end
    n == 1 && (axs = first(axs))
    return fig, axs
end

function plot_fields!(axs, b, u, v, label, color, grid=first(observations).grid; kwargs...)
    z = znodes(Center, grid)
    ## Note unit conversions below, eg m s⁻² -> cm s⁻²:
    lines!(axs[1], 1e2 * b, z; color, label, kwargs...)
    lines!(axs[2], 1e2 * u, z; color, label, kwargs...)
    lines!(axs[3], 1e2 * v, z; color, label, kwargs...)
    return nothing
end

#####
##### Calibration
#####

ri_based_closure = RiBasedVerticalDiffusivity()

simulation = ensemble_column_model_simulation(observations;
                                              Nensemble = 100,
                                              architecture = CPU(),
                                              tracers = (:b, :e),
                                              closure = ri_based_closure)

Qᵘ = simulation.model.velocities.u.boundary_conditions.top.condition
Qᵇ = simulation.model.tracers.b.boundary_conditions.top.condition
N² = simulation.model.tracers.b.boundary_conditions.bottom.condition

simulation.Δt = 10minutes

for (i, obs) in enumerate(observations)
    view(Qᵘ, :, i) .= obs.metadata.parameters.momentum_flux
    view(Qᵇ, :, i) .= obs.metadata.parameters.buoyancy_flux
    view(N², :, i) .= obs.metadata.parameters.N²_deep
end

priors = (ν₀   = lognormal(mean=0.1, std=0.1),
          κ₀   = lognormal(mean=0.1, std=0.1),
          Ri₀ν = Normal(-0.5, 1.0),
          Ri₀κ = Normal(-0.5, 1.0),
          Riᵟν = lognormal(mean=1.0,  std=1.0),
          Riᵟκ = lognormal(mean=1.0,  std=1.0))

free_parameters = FreeParameters(priors)
calibration = InverseProblem(observations, simulation, free_parameters)
eki = EnsembleKalmanInversion(calibration; convergence_rate=0.95)

#####
##### Visualization
#####

modeled_time_serieses = calibration.time_series_collector.field_time_serieses 
observed, mean_modeled, best_modeled, worst_modeled = [], [], [], []
Nt = length(first(observations).times)
for (c, obs) in enumerate(observations)
    push!(observed,      map(name -> interior(obs.field_time_serieses[name][Nt], 1, 1, :), field_names))
    push!(mean_modeled,  map(name -> interior(  modeled_time_serieses[name][Nt], 1, c, :), field_names))
    push!(best_modeled,  map(name -> interior(  modeled_time_serieses[name][Nt], 2, c, :), field_names))
    push!(worst_modeled, map(name -> interior(  modeled_time_serieses[name][Nt], 3, c, :), field_names))
end

function compare_model_observations(model_label="modeled")
    fig, axs = make_figure_axes(length(observations))
    mean_model_label = "mean " * model_label
    best_model_label = "best " * model_label
    worst_model_label = "worst " * model_label
    for (c, obs) in enumerate(observations)
        plot_fields!(axs[c], observed[c]..., "observed at t = " * prettytime(times[end]), (:black, 0.6), linewidth=5)
        plot_fields!(axs[c], mean_modeled[c]..., mean_model_label, :darkblue, linewidth=2)
        plot_fields!(axs[c], best_modeled[c]..., best_model_label, :orange, linewidth=3)
        plot_fields!(axs[c], worst_modeled[c]..., worst_model_label, :seagreen, linewidth=1)
        [axislegend(ax, position=:rb, merge=true, fontsize=10) for ax in axs[c]]
    end
    return fig
end

#####
##### Calibration
#####

while true
    try
        iterate!(eki, iterations=4)
        latest_summary = eki.iteration_summaries[end]
        @show latest_summary

        # Plot best current parameters
        Niter = eki.iteration

        min_error, k_min = finitefindmin(latest_summary.mean_square_errors)
        max_error, k_max = finitefindmax(latest_summary.mean_square_errors)

        mean_parameters = latest_summary.ensemble_mean
        best_parameters = latest_summary.parameters[k_min]
        worst_parameters = latest_summary.parameters[k_max]

        forward_run!(calibration, [mean_parameters, best_parameters, worst_parameters])
        fig = compare_model_observations("($Niter iters)")
        save("multi_case_model_observation_comparison_iteration_$Niter.png", fig)
        display(fig)

    catch err
        err isa InterruptException && break
        throw(err)
    end
end

#=
ensemble_means = NamedTuple(n => map(summary -> summary.ensemble_mean[n], eki.iteration_summaries)
                            for n in calibration.free_parameters.names)

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Ensemble Kalman iteration", ylabel = "Parameter value")

for (i, name) in enumerate(calibration.free_parameters.names)
    label = string(name)
    marker = markercycle[i]
    color = colorcycle[i]
    scatterlines!(ax, 0:Niter, parent(ensemble_means[name]); marker, color, label)
end
=#

