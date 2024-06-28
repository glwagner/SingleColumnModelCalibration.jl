vectorize(a) = [a]
vectorize(a::Vector) = a

rectilinear_grids_from_parameters(grid_parameters) =
    [RectilinearGrid(; p..., topology=(Flat, Flat, Bounded)) for p in grid_parameters]

default_start_time = 10minutes
default_Ntimes = 2

function generate_filepath(; suite_parameters,
                             grid_parameters,
                             Nensemble,
                             stop_pseudotime,
                             filename,
                             dir,
                             Δt)

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

    prefix = @sprintf("%s_Nens%d_Δt%d_τ%d_%s", filename, Nensemble,
                      Δt, stop_pseudotime, gridstr)

    filepath = joinpath(dir, string(prefix, "_", suitestr, ".jld2"))

    return filepath
end

function build_batched_inverse_problem(closure, name="";
                                       free_parameters = get_free_parameters(name),
                                       grid_parameters,
                                       suite_parameters,
                                       Nensemble = 100,
                                       start_time = default_start_time,
                                       Ntimes = default_Ntimes,
                                       Δt = 20minutes,
                                       architecture = CPU(),
                                       overwrite_existing = false,
                                       tke_weight = 0.0)

    suite_parameters = vectorize(suite_parameters)
    grid_parameters = vectorize(grid_parameters)

    #####
    ##### Build grids and noise covariances for each grid
    #####

    inverse_problems = []
    weights = []

    grids = rectilinear_grids_from_parameters(grid_parameters)
    inverse_problem_kwargs = (; free_parameters, Nensemble, architecture, closure, tke_weight)

    for p in suite_parameters
        suite = p.name
        stop_time = p.stop_time
        times = Ntimes == 2 ? [start_time, stop_time] :
                Ntimes == 3 ? [start_time, 2stop_time/3, stop_time] :
                Ntimes == 4 ? [start_time, stop_time/3, 2stop_time/3, stop_time] :
                              collect(range(start_time, stop=stop_time, length=Ntimes))
        observations_resolution = :resolution ∈ keys(p) ? p.resolution : "1m" 

        grid_inverse_problems = []
        grid_weights = []

        for grid in grids
            ip = lesbrary_inverse_problem(grid; times, Δt, suite, observations_resolution,
                                          inverse_problem_kwargs...)

            # Assume the grid is not "too stretched" for this to be useful
            weight = 1 / grid.Nz

            push!(grid_inverse_problems, ip)
            push!(grid_weights, weight)
        end

        push!(inverse_problems, grid_inverse_problems)
        push!(weights, grid_weights)
    end

    Ngrids = length(grids)
    Nsuites = length(suite_parameters)

    ip_list      = tuple([inverse_problems[s][g]  for g=1:Ngrids, s=1:Nsuites]...)
    weights_list = tuple([weights[s][g]           for g=1:Ngrids, s=1:Nsuites]...)

    batched_ip = BatchedInverseProblem(ip_list...; weights=weights_list)

    return batched_ip
end

function build_ensemble_kalman_inversion(closure, name="";
                                         grid_parameters,
                                         suite_parameters,
                                         start_time = default_start_time,
                                         Ntimes = default_Ntimes,
                                         tke_weight = 0.0,
                                         # EnsembleKalmanInverion parameters
                                         noise_covariance = nothing,
                                         resample_failure_fraction = 0.1,
                                         acceptable_failure_fraction = 1.0,
                                         forward_map_output = nothing,
                                         resampler = Resampler(; resample_failure_fraction, acceptable_failure_fraction),
                                         initial_convergence_ratio = 0.7,
                                         mark_failed_particles = ObjectiveLossThreshold(3.0),
                                         other_kw...)

    batched_ip = build_batched_inverse_problem(closure, name; Ntimes, grid_parameters, suite_parameters, other_kw...)
    grids = rectilinear_grids_from_parameters(grid_parameters)
    noise_covariances = []

    for p in suite_parameters
        suite = p.name
        stop_time = p.stop_time
        times = Ntimes == 2 ? [start_time, stop_time] : collect(range(start_time, stop=stop_time, length=Ntimes))
        grid_Γ = []

        for grid in grids
            # Estimate noise covariance based on discrepency between LES with different resolution
            Γ = estimate_noise_covariance(grid; times, tke_weight, suite)
            push!(grid_Γ, Γ)
        end

        push!(noise_covariances, grid_Γ)
    end

    Ngrids = length(grids)
    Nsuites = length(suite_parameters)
    Γ_list = tuple([noise_covariances[s][g] for g=1:Ngrids, s=1:Nsuites]...)

    if isnothing(noise_covariance)
        batched_Γ = Matrix(BlockDiagonal([wi * Γi for (wi, Γi) in zip(batched_ip.weights, Γ_list)]))
    else
        batched_Γ = noise_covariance
    end

    pseudo_stepping = Kovachki2018InitialConvergenceRatio(; initial_convergence_ratio)

    eki = EnsembleKalmanInversion(batched_ip;
                                  noise_covariance = batched_Γ,
                                  pseudo_stepping,
                                  forward_map_output,
                                  mark_failed_particles,
                                  resampler)

    return eki
end

function calibrate_parameter_set(closure, name="";
                                 prefixname = name,
                                 savedir = @__DIR__,
                                 stop_pseudotime = 1.0,
                                 stop_iteration = 1000,
                                 plot_progress = true,
                                 eki_kwargs...)

    eki = build_ensemble_kalman_inversion(closure, name; eki_kwargs...)

    if plot_progress
        figs = calibration_progress_figure(eki; Nsuites, Ngrids)
        for fig in figs
            display(fig)
        end
    end

    @info "Now for the main event..."

    while eki.pseudotime < stop_pseudotime && eki.iteration < stop_iteration
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

