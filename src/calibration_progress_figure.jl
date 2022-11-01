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
        plot_fields!(axs, "Obs at t = " * prettytime(times[1]), (:black, 0.6), grid, initial_case_field_data...;
                     linewidth=2, linestyle=:dash)
        plot_fields!(axs, "Obs at t = " * prettytime(times[Nt]), (:gray23, 0.4), grid, final_case_field_data...; linewidth=6)

        if eki.inverse_problem isa BatchedInverseProblem
            # Plot model case with minimum error
            ip = eki.inverse_problem[1] # low res 
            Nz = size(ip.simulation.model.grid, 3)
            obs = ip.observations[c]
            grid = obs.grid
            iter = eki.iteration

            for p in 1:Nparticles
                data = NamedTuple(n => get_modeled_case(ip, c, n, kk[p]) for n in keys(obs.field_time_serieses))
                plot_fields!(axs, "rank $p, iter $iter (Nz = $Nz)", :navy, grid, data...; linestyle=linestyles[p])
            end

            ip = eki.inverse_problem[2] # high res 
            Nz = size(ip.simulation.model.grid, 3)
            obs = ip.observations[c]
            grid = obs.grid

            for p in 1:Nparticles
                data = NamedTuple(n => get_modeled_case(ip, c, n, kk[p]) for n in keys(obs.field_time_serieses))
                plot_fields!(axs, "rank $p, iter $iter (Nz = $Nz)", :orange, grid, data...; linestyle=linestyles[p])
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

