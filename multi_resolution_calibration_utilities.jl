using DataDeps
using Oceananigans
using CairoMakie
using ElectronDisplay

using ParameterEstimocean
using ParameterEstimocean: Transformation
using ParameterEstimocean.Observations: forward_map_names

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    CATKEVerticalDiffusivity,
    MixingLength,
    TurbulentKineticEnergyEquation

include("prior_library.jl")

# Coarse grid used by ECCO
ecco_vertical_grid = [-256.0,
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

cases = ["free_convection",
         "strong_wind_weak_cooling",
         "med_wind_med_cooling",
         "weak_wind_strong_cooling",
         "strong_wind",
         "strong_wind_no_rotation"]
    
function lesbrary_inverse_problem(regrid;
                                  free_parameters,
                                  times = [48hours - 10minutes, 48hours],
                                  Nensemble = 500,
                                  Δt = 10minutes,
                                  field_names = (:b, :e, :u, :v),
                                  closure = CATKEVerticalDiffusivity(),
                                  non_ensemble_closure = nothing,
                                  suite = "one_day_suite",
                                  tke_weight = 0.1,
                                  architecture = GPU())

    normalizations = (b = ZScore(),
                      u = ZScore(),
                      v = ZScore(),
                      e = RescaledZScore(tke_weight))

    # Retain only the upper half
    Nz = regrid.Nz
    Nh = floor(Int, Nz/3)
    space = SpaceIndices(z=Nh:Nz)

    transformation = NamedTuple(n => Transformation(; space, normalization=normalizations[n])
                                for n in keys(normalizations))

    case_path(case) = joinpath("data", suite, case * "_instantaneous_statistics.jld2")

    observation_library = Dict()

    # Don't optimize u, v for free_convection
    free_convection_names = filter(n -> n == :u || n == :v, field_names)
    observation_library["free_convection"] =
        SyntheticObservations(case_path("free_convection"); transformation, times, regrid,
                              field_names = free_convection_names)
                                                                    
    # Don't optimize v for non-rotating cases
    strong_wind_no_rotation_names = filter(n -> n == :v, field_names)
    observation_library["strong_wind_no_rotation"] =
        SyntheticObservations(case_path("strong_wind_no_rotation"); transformation, times, regrid,
                              field_names = strong_wind_no_rotation_names)

    # The rest are standard
    for case in ["strong_wind", "med_wind_med_cooling", "strong_wind_weak_cooling", "weak_wind_strong_cooling"]
        observation_library[case] = SyntheticObservations(case_path(case); transformation, times, regrid, field_names)
    end

    observations = [observation_library[case] for case in cases]

    simulation = ensemble_column_model_simulation(observations; closure, Nensemble, architecture,  non_ensemble_closure,
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

    weights = [1.0 for _ in observations]
    #weights[1] = 0.1 # downweight free convection case
    batched_observations = BatchedSyntheticObservations(observations; weights)

    ip = InverseProblem(batched_observations, simulation, free_parameters)

    return ip
end

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
    times = ip.time_series_collector.times
    Nt = length(times)
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

linestyles = [nothing,
              :dash,
              :dot,
              :dashdot,
              :dashdotdot]

function calibration_progress_figure(eki; Nparticles=2)
    high_res_ip = eki.inverse_problem isa BatchedInverseProblem ?
        eki.inverse_problem[2] : eki.inverse_problem
    times = first(high_res_ip.observations).times
    field_names = forward_map_names(high_res_ip.observations)
    Nt = length(times)

    latest_summary = eki.iteration_summaries[end]
    min_error, k_min = finitefindmin(latest_summary.mean_square_errors)
    # max_error, k_max = finitefindmax(latest_summary.mean_square_errors)
    
    errors = deepcopy(latest_summary.mean_square_errors)
    notnans = isfinite.(errors)
    errors[.!notnans] .= +Inf
    kk = sortperm(errors)

    fig = Figure(resolution=(1200, 1200))

    # Plot case by case
    for (c, case) in enumerate(cases)
        # Make axes
        label = replace(case, "_" => "\n")
        axs = make_axes(fig, c, label)

        # Plot observed data for each field
        case_obs = high_res_ip.observations[c]
        case_dataset = case_obs.field_time_serieses
        grid = case_obs.grid
        case_names = keys(case_dataset)
        initial_case_field_data = NamedTuple(n => interior(getproperty(case_dataset, n)[1])[1, 1, :] for n in case_names)
        final_case_field_data = NamedTuple(n => interior(getproperty(case_dataset, n)[Nt])[1, 1, :] for n in case_names)
        plot_fields!(axs, "Observed at t = " * prettytime(times[1]), (:black, 0.6), grid, initial_case_field_data...;
                     linewidth=2, linestyle=:dash)
        plot_fields!(axs, "Observed at t = " * prettytime(times[Nt]), (:gray23, 0.4), grid, final_case_field_data...; linewidth=6)

        if eki.inverse_problem isa BatchedInverseProblem
            # Plot model case with minimum error
            ip = eki.inverse_problem[1] # low res 
            obs = ip.observations[c]
            grid = obs.grid
            iter = eki.iteration

            for p in 1:Nparticles
                data = NamedTuple(n => get_modeled_case(ip, c, n, kk[p]) for n in keys(obs.field_time_serieses))
                plot_fields!(axs, "rank $p, iter $iter (low res)", :navy, grid, data...; linestyle=linestyles[p])
            end

            ip = eki.inverse_problem[2] # high res 
            obs = ip.observations[c]
            grid = obs.grid

            for p in 1:Nparticles
                data = NamedTuple(n => get_modeled_case(ip, c, n, kk[p]) for n in keys(obs.field_time_serieses))
                plot_fields!(axs, "rank $p, iter $iter (hi res)", :orange, grid, data...; linestyle=linestyles[p])
            end

        else
            ip = eki.inverse_problem
            obs = ip.observations[c]
            grid = obs.grid
            iter = eki.iteration

            for p in 1:Nparticles
                data = NamedTuple(n => get_modeled_case(ip, c, n, kk[p]) for n in keys(obs.field_time_serieses))
                plot_fields!(axs, "rank $p, iter $iter", :navy, grid, data...; linestyle=linestyles[p])
            end
        end


        fig[1, 6] = Legend(fig, axs[1]) 
    end

    return fig
end

