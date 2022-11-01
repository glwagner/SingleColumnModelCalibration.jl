function iterate_until!(eki, pseudotime_limit, max_iterations=1000)
    while eki.pseudotime < pseudotime_limit && eki.iteration < max_iterations
        iterate!(eki)

        if eki.iteration % 10 == 0
            @show eki.iteration_summaries[end]
            display(calibration_progress_figure(eki))
        end
    end

    @show eki.iteration_summaries[end]
    display(calibration_progress_figure(eki))

    return nothing
end

function calibrate_parameter_set(name, closure;
                                 savedir = @__DIR__,
                                 Nensemble = 200,
                                 tke_weight = 1e-2,
                                 suite = "one_day_suite",
                                 pseudotime_limit = 4.0,
                                 max_iterations = 1000,
                                 resample_failure_fraction = 0.2,
                                 acceptable_failure_fraction = 1.0,
                                 initial_convergence_ratio = 0.7,
                                 coarse_calibration = true,
                                 fine_calibration = true,
                                 combined_calibration = true,
                                 fine_Nz = 32,
                                 mark_failed_particles = ObjectiveLossThreshold(3.0),
                                 free_parameters = get_free_parameters(name),
                                 fine_Δt = 5minutes,
                                 coarse_Δt = 5minutes,
                                 architecture = CPU())

    prefix = @sprintf("%s_Nens%d_%s", name, Nensemble, suite)

    @info "Calibrating free parameters $free_parameters..."

    #####
    ##### Build grids and noise covariances for each grid
    #####

    times = collect(range(2hours, stop=24hours, length=4))

    # Two grids: "coarse" and "fine"
    fine_regrid = RectilinearGrid(size=fine_Nz; z=(-256, 0), topology=(Flat, Flat, Bounded))

    # Coarse grid has ECCO vertical resolution to z=-256 m
    Nz_ecco = length(ecco_vertical_grid) - 1
    coarse_regrid = RectilinearGrid(size=Nz_ecco, z=ecco_vertical_grid, topology=(Flat, Flat, Bounded))

    # Estimate noise covariance based on discrepency between LES with different resolution
    coarse_Γ = estimate_noise_covariance(coarse_regrid; times, tke_weight, suite)
    fine_Γ = estimate_noise_covariance(fine_regrid; times, tke_weight, suite)

    #####
    ##### Single-resolution calibrations
    #####

    coarse_eki = nothing
    fine_eki = nothing
    eki = nothing

    inverse_problem_kwargs = (; suite, free_parameters, Nensemble, architecture, closure)
    resampler = Resampler(; resample_failure_fraction, acceptable_failure_fraction)
    pseudo_stepping = Kovachki2018InitialConvergenceRatio(; initial_convergence_ratio)
    eki_kwargs = (; pseudo_stepping, mark_failed_particles, resampler)

    if fine_calibration
        @info "Calibrating $name parameters to the fine grid..."
        fine_ip = lesbrary_inverse_problem(fine_regrid; times, Δt=fine_Δt, inverse_problem_kwargs...)
        fine_eki = EnsembleKalmanInversion(fine_ip; noise_covariance=fine_Γ, eki_kwargs...)
        iterate_until!(fine_eki, pseudotime_limit, max_iterations)

        name = string(prefix, "_Nz", fine_regrid.Nz, "_calibration")
        datapath = joinpath(savedir, name * ".jld2")
        @save datapath iteration_summaries=fine_eki.iteration_summaries

        progress_fig = calibration_progress_figure(fine_eki)
        display(progress_fig)
        figpath = joinpath(savedir, name * ".png")
        save(figpath, progress_fig)
    end

    if coarse_calibration
        @info "Calibrating $name parameters to the coarse grid..."
        coarse_ip = lesbrary_inverse_problem(coarse_regrid; times, Δt=coarse_Δt, inverse_problem_kwargs...)
        coarse_eki = EnsembleKalmanInversion(coarse_ip; noise_covariance=coarse_Γ, eki_kwargs...)
        iterate_until!(coarse_eki, pseudotime_limit, max_iterations)

        name = string(prefix, "_Nz", coarse_regrid.Nz, "_calibration")
        datapath = joinpath(savedir, name * ".jld2")
        @save datapath iteration_summaries=coarse_eki.iteration_summaries

        progress_fig = calibration_progress_figure(coarse_eki)
        display(progress_fig)
        figpath = joinpath(savedir, name * ".png")
        save(figpath, progress_fig)
    end

    #####
    ##### "Combined" calibration to multiple resolutions
    #####

    if combined_calibration
        pseudo_stepping = Kovachki2018InitialConvergenceRatio(; initial_convergence_ratio)
        coarse_ip = lesbrary_inverse_problem(coarse_regrid; times, Δt=coarse_Δt, inverse_problem_kwargs...)
        fine_ip = lesbrary_inverse_problem(fine_regrid; times, Δt=fine_Δt, inverse_problem_kwargs...)

        weights = (1.0, coarse_regrid.Nz / fine_regrid.Nz)
        batched_ip = BatchedInverseProblem(coarse_ip, fine_ip; weights)

        combined_Γ = Matrix(BlockDiagonal([weights[1] * coarse_Γ, weights[2] * fine_Γ]))
        pseudo_stepping = Kovachki2018InitialConvergenceRatio(; initial_convergence_ratio)

        fine_ip = lesbrary_inverse_problem(fine_regrid; times, Δt=fine_Δt, inverse_problem_kwargs...)
        coarse_ip = lesbrary_inverse_problem(coarse_regrid; times, Δt=coarse_Δt, inverse_problem_kwargs...)

        eki = EnsembleKalmanInversion(batched_ip;
                                      noise_covariance = combined_Γ,
                                      eki_kwargs...)

        progress = calibration_progress_figure(eki)
        display(progress)
        @show eki.iteration_summaries[end]

        @info "Now for the main event..."

        while eki.pseudotime < pseudotime_limit && eki.iteration < max_iterations
            iterate!(eki)
            progress = calibration_progress_figure(eki)
            latest_summary = eki.iteration_summaries[end]
            display(progress)
            @show latest_summary
        end

        @show eki.iteration_summaries[end]

        name = string(prefix, "_combined_Nz", coarse_regrid.Nz, "_Nz", fine_regrid.Nz, "_calibration")
        datapath = joinpath(savedir, name * ".jld2")
        @save datapath iteration_summaries=eki.iteration_summaries

        progress_fig = calibration_progress_figure(eki)
        display(progress_fig)
        figpath = joinpath(savedir, name * ".png")
        save(figpath, progress_fig)
    end

    return coarse_eki, fine_eki, eki
end

