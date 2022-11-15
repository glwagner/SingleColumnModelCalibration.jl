vectorize(a) = [a]
vectorize(a::Vector) = a

function calibrate_parameter_set(name, closure;
                                 prefixname = name,
                                 savedir = @__DIR__,
                                 Nensemble = 100,
                                 tke_weight = 0.0,
                                 pseudotime_limit = 1.0,
                                 max_iterations = 1000,
                                 resample_failure_fraction = 0.2,
                                 acceptable_failure_fraction = 1.0,
                                 initial_convergence_ratio = 0.7,
                                 plot_progress = true,
                                 mark_failed_particles = ObjectiveLossThreshold(2.0),
                                 free_parameters = get_free_parameters(name),
                                 architecture = CPU(),
                                 Ntimes = 3,
                                 Δt = 20minutes,
                                 start_time = 2hours,
                                 grid_parameters = [(size=19, z=ecco_vertical_grid), (size=64, z=(-256, 0))],
                                 suite_parameters = [(name="24_hour_suite", stop_time=24hours)])

    suite_parameters = vectorize(suite_parameters)
    grid_parameters = vectorize(grid_parameters)

    suitestr = ""
    for (n, s) in enumerate(suite_parameters)
        suitestr *= s.name
        n < length(suite_parameters) && (suitestr *= "_")
    end

    gridstr = ""
    for (n, p) in enumerate(grid_parameters)
        gridstr *= string("Nz", p.size)
        n < length(grid_parameters) && (gridstr *= "_")
    end

    prefix = @sprintf("%s_Nens%d_Δt%d_%s", prefixname, Nensemble, Δt, gridstr)

    #####
    ##### Build grids and noise covariances for each grid
    #####

    inverse_problems = []
    weights = []
    noise_covariances = []

    grids = [RectilinearGrid(; p..., topology=(Flat, Flat, Bounded)) for p in grid_parameters]
    inverse_problem_kwargs = (; free_parameters, Nensemble, architecture, closure, tke_weight)

    for p in suite_parameters
        suite = p.name
        stop_time = p.stop_time
        times = [start_time, stop_time] #collect(range(start_time, stop=stop_time, length=Ntimes))

        grid_inverse_problems = []
        grid_weights = []
        grid_Γ = []

        for grid in grids
            # Estimate noise covariance based on discrepency between LES with different resolution
            Γ = estimate_noise_covariance(grid; times, tke_weight, suite)
            ip = lesbrary_inverse_problem(grid; times, Δt, suite, inverse_problem_kwargs...)
            weight = 1 / grid.Nz 

            push!(grid_inverse_problems, ip)
            push!(grid_weights, weight)
            push!(grid_Γ, Γ)
        end

        push!(inverse_problems, grid_inverse_problems)
        push!(noise_covariances, grid_Γ)
        push!(weights, grid_weights)
    end

    Ngrids = length(grids)
    Nsuites = length(suite_parameters)

    ip_list      = tuple([inverse_problems[s][g]  for g=1:Ngrids, s=1:Nsuites]...)
    weights_list = tuple([weights[s][g]           for g=1:Ngrids, s=1:Nsuites]...)
    Γ_list       = tuple([noise_covariances[s][g] for g=1:Ngrids, s=1:Nsuites]...)

    batched_ip = BatchedInverseProblem(ip_list...; weights=weights_list)
    batched_Γ = Matrix(BlockDiagonal([wi * Γi for (wi, Γi) in zip(weights_list, Γ_list)]))

    resampler = Resampler(; resample_failure_fraction, acceptable_failure_fraction)
    pseudo_stepping = Kovachki2018InitialConvergenceRatio(; initial_convergence_ratio)

    eki = EnsembleKalmanInversion(batched_ip;
                                  noise_covariance = batched_Γ,
                                  pseudo_stepping,
                                  mark_failed_particles,
                                  resampler)

    if plot_progress
        figs = calibration_progress_figure(eki; Nsuites, Ngrids)
        for fig in figs
            display(fig)
        end
    end

    @info "Now for the main event..."

    while eki.pseudotime < pseudotime_limit && eki.iteration < max_iterations
        iterate!(eki)

        if (eki.iteration % 10) == 0.0

            if plot_progress
                figs = calibration_progress_figure(eki; Nsuites, Ngrids)
                for fig in figs
                    display(fig)
                end
            end

            latest_summary = eki.iteration_summaries[end]
            @show latest_summary
        end
    end

    @show eki.iteration_summaries[end]

    datapath = joinpath(savedir, string(prefix, "_", suitestr, ".jld2"))
    @save datapath iteration_summaries=eki.iteration_summaries

    if plot_progress
        figs = calibration_progress_figure(eki; Nsuites, Ngrids)

        for (suite, fig) in zip(suite_parameters, figs)
            display(fig)
            suitename = suite.name
            figpath = joinpath(savedir, string(prefix, "_", suitename, ".png"))
            save(figpath, fig)
        end
    end

    return eki
end

