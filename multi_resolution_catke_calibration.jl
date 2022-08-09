using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using ParameterEstimocean: Transformation
using ParameterEstimocean.InverseProblems: BatchedInverseProblem
using LinearAlgebra
using DataDeps
using Distributions
using ParameterEstimocean.Observations: forward_map_names
using GLMakie

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity

#####
##### Compile LESbrary
#####

times = [44hours, 48hours] #6hours, 12hours, 24hours, 36hours, 48hours]
field_names = (:b, :e, :u, :v)
Nensemble = 100
Δt = 10minutes
closure = CATKEVerticalDiffusivity()

include("prior_library.jl")

parameter_names = (:CᵂwΔ,  :Cᵂu★, :Cᴰ⁻, :Cᴰʳ,
                   :Cˢc,   :Cˢu,  :Cˢe,
                   :Cᵇc,   :Cᵇu,  :Cᵇe,
                   :Cᴷc⁻,  :Cᴷu⁻, :Cᴷe⁻,
                   :Cᴷcʳ,  :Cᴷuʳ, :Cᴷeʳ,
                   :CᴰRiᶜ, :CᴰRiʷ,
                   :CᴷRiᶜ, :CᴷRiʷ)

#:Cᴬc,   :Cᴬe,
#:Cᴬc,  :Cᴬu, :Cᴬe, :Cᴬˢc, :Cᴬˢu, :Cᴬˢe,
#:Cᴬc,  :Cᴬˢc,

free_parameters = FreeParameters(prior_library, names=parameter_names)

cases = ["free_convection",
         "strong_wind_weak_cooling",
         "med_wind_med_cooling",
         "weak_wind_strong_cooling",
         "strong_wind",
         "strong_wind_no_rotation"]

suite = "two_day_suite"
case_path(case) = joinpath("data", suite, case * "_instantaneous_statistics.jld2")

transformation = (b = ZScore(),
                  u = ZScore(),
                  v = ZScore(),
                  e = RescaledZScore(0.05))

function build_inverse_problem(regrid)

    observation_library = Dict()

    # Don't optimize u, v for free_convection
    observation_library["free_convection"] =
        SyntheticObservations(case_path("free_convection"); transformation, times, regrid,
                              field_names = (:b, :e))
                                                                    
    # Don't optimize v for non-rotating cases
    observation_library["strong_wind_no_rotation"] =
        SyntheticObservations(case_path("strong_wind_no_rotation"); transformation, times, regrid,
                              field_names = (:b, :e, :u))

    # The rest are standard
    for case in ["strong_wind", "med_wind_med_cooling", "strong_wind_weak_cooling", "weak_wind_strong_cooling"]
        observation_library[case] = SyntheticObservations(case_path(case); field_names,
                                                          transformation, times, regrid)
    end

    observations = [observation_library[case] for case in cases]

    simulation = ensemble_column_model_simulation(observations; closure, Nensemble,
                                                  architecture = CPU(),
                                                  tracers = (:b, :e))

    simulation.Δt = Δt    

    Qᵘ = simulation.model.velocities.u.boundary_conditions.top.condition
    Qᵇ = simulation.model.tracers.b.boundary_conditions.top.condition
    N² = simulation.model.tracers.b.boundary_conditions.bottom.condition
    
    for (case, obs) in enumerate(observations)
        f = obs.metadata.parameters.coriolis_parameter
        view(Qᵘ, :, case) .= obs.metadata.parameters.momentum_flux
        view(Qᵇ, :, case) .= obs.metadata.parameters.buoyancy_flux
        view(N², :, case) .= obs.metadata.parameters.N²_deep
        view(simulation.model.coriolis, :, case) .= Ref(FPlane(f=f))
    end

    batched_observations = BatchedSyntheticObservations(observations)
    ip = InverseProblem(batched_observations, simulation, free_parameters)

    return ip
end

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

coarse_regrid = RectilinearGrid(size=length(z)-1; z, topology=(Flat, Flat, Bounded))
fine_regrid = RectilinearGrid(size=128; z=(-256, 0), topology=(Flat, Flat, Bounded))

coarse_ip = build_inverse_problem(coarse_regrid)
fine_ip = build_inverse_problem(fine_regrid)

calibration = BatchedInverseProblem(coarse_ip, fine_ip)

resampler = Resampler(resample_failure_fraction=0.0, acceptable_failure_fraction=1.0)
#eki = EnsembleKalmanInversion(coarse_ip; resampler, pseudo_stepping = ConstantConvergence(0.9))
eki = EnsembleKalmanInversion(calibration; resampler, pseudo_stepping = ConstantConvergence(0.9))

#####
##### Plot utils
#####

function finitefind(a, val, find)
    b = deepcopy(a)
    b[.!isfinite.(a)] .= val
    return find(b)
end

finitefindmin(a) = finitefind(a, Inf, findmin)
finitefindmax(a) = finitefind(a, -Inf, findmax)

function make_axes(fig, row=1, label=nothing)
    ax_b = Axis(fig[row, 1], xlabel = "Buoyancy \n[cm s⁻²]", ylabel = "z [m]")
    ax_u = Axis(fig[row, 2], xlabel = "x-velocity \n[cm s⁻¹]")
    ax_v = Axis(fig[row, 3], xlabel = "y-velocity \n[cm s⁻¹]")
    ax_e = Axis(fig[row, 4], xlabel = "Turbulent kinetic energy \n[cm² s⁻²]")
    if !isnothing(label)
        ax_t = Axis(fig[row, 5])
        xlims!(0, 1)
        ylims!(0, 1)
        hidespines!(ax_t)
        hidedecorations!(ax_t)
        text!(ax_t, label, justification=:left, align=(:left, :center), position=(0, 0.5))
    end
    return (ax_b, ax_u, ax_v, ax_e)
end

function get_modeled_case(ip, c, name, k=1)
    model_time_serieses = ip.time_series_collector.field_time_serieses 
    field = getproperty(model_time_serieses, name)[Nt]
    return interior(field, k, c, :)
end

function plot_fields!(axs, label, color, grid, b, e, u=zeros(size(b)), v=zeros(size(b)); linewidth=2, linestyle=:solid)
    z = znodes(Center, grid)
    b, u, v, e = Tuple(Array(f) for f in (b, u, v, e))

    for (q, name) in zip((b, u, v, e), ("b", "u", "v", "e"))
        any(isnan.(q)) && @warn("NaNs found in $label $(name)!")
    end

    ## Note unit conversions below, eg m s⁻¹ -> cm s⁻¹:cyan
    lines!(axs[1], 1e2 * b, z; color, linestyle, label, linewidth) 
    lines!(axs[2], 1e2 * u, z; color, linestyle, label, linewidth)
    lines!(axs[3], 1e2 * v, z; color, linestyle, label, linewidth)
    lines!(axs[4], 1e4 * e, z; color, linestyle, label, linewidth)

    return nothing
end


function make_fig(eki)
    high_res_ip = eki.inverse_problem[2]
    times = first(high_res_ip.observations).times
    field_names = forward_map_names(high_res_ip.observations)
    Nt = length(times)

    latest_summary = eki.iteration_summaries[end]
    min_error, k_min = finitefindmin(latest_summary.mean_square_errors)
    # max_error, k_max = finitefindmax(latest_summary.mean_square_errors)

    fig = Figure(resolution=(1200, 1200))

    # Plot case by case
    for (c, case) in enumerate(cases)
        # Make axes
        label = replace(case, "_" => "\n")
        axs = make_axes(fig, c, label)

        # Plot observed data for each field
        case_obs = high_res_ip.observations[c]
        case_dataset = case_obs.field_time_serieses
        case_names = keys(case_dataset)
        case_field_data = NamedTuple(n => interior(getproperty(case_dataset, n)[Nt])[1, 1, :] for n in case_names)
        plot_fields!(axs, "Observed at t = " * prettytime(times[Nt]), (:gray23, 0.6), case_field_data...; linewidth=4)

        # Plot model case with minimum error
        ip = eki.inverse_problem[1] # low res 
        obs = ip.observations
        grid = ip.observations.grid
        min_error_data = NamedTuple(n => get_modeled_case(ip, c, n, k_min) for n in keys(obs.field_time_serieses))
        plot_fields!(axs, "min (low res)", :navy, grid, min_error_data...)

        ip = eki.inverse_problem[2] # high res 
        obs = ip.observations
        grid = ip.observations.grid
        min_error_data = NamedTuple(n => get_modeled_case(ip, c, n, k_min) for n in keys(obs.field_time_serieses))
        plot_fields!(axs, "min (hi res)", :orange, grid, min_error_data...)

        fig[1, 6] = Legend(fig, axs[1]) 
    end

    return fig
end

#####
##### Calibrate
#####

fig = make_fig(eki)
display(fig)

iterate!(eki, iterations=10)

fig = make_fig(eki)
display(fig)

