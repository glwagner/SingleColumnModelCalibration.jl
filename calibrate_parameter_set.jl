using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using ParameterEstimocean.InverseProblems: BatchedInverseProblem
using ParameterEstimocean.PseudoSteppingSchemes: Kovachki2018InitialConvergenceRatio
using Printf
using JLD2

# Import "lesbrary_inverse_problem" + plotting utility "calibration_progress_figure"
include("multi_resolution_calibration_utilities.jl")
include("parameter_sets.jl")

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

function calibrate_parameter_set(name;
                                 dir = "/home/greg/Projects/LocalOceanClosureCalibration/calibration_summaries",
                                 Nensemble = 200,
                                 tke_weight = 1e-2,
                                 suite = "one_day_suite",
                                 pseudotime_limit = 1.0,
                                 max_iterations = 1000,
                                 resample_failure_fraction = 0.2,
                                 acceptable_failure_fraction = 1.0,
                                 initial_convergence_ratio = 0.7,
                                 fine_Nz = 32,
                                 fine_Δt = 10minutes,
                                 architecture = GPU())

    savename = @sprintf("%s_Nens%d_%s_fine_Nz%d", name, Nensemble, suite, fine_Nz)

    #####
    ##### Build closure and free parameters
    #####

    mixing_length = MixingLength(; neutral_default_mixing_length_parameters...)
    turbulent_kinetic_energy_equation = TurbulentKineticEnergyEquation(; neutral_default_tke_parameters...)
    closure = CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation)

    free_parameters = FreeParameters(prior_library,
                                     names = parameter_sets[name],
                                     dependent_parameters = dependent_parameter_sets[name])

    @info "Calibrating free parameters $free_parameters..."

    #####
    ##### Build grids and noise covariances for each grid
    #####

    times = range(2hours, stop=24hours, length=4)

    # Two grids: "coarse" and "fine"
    fine_regrid = RectilinearGrid(size=fine_Nz; z=(-256, 0), topology=(Flat, Flat, Bounded))

    # Coarse grid has ECCO vertical resolution to z=-256 m
    Nz_ecco = length(ecco_vertical_grid) - 1
    coarse_Δt = 10minutes
    coarse_regrid = RectilinearGrid(size=Nz_ecco, z=ecco_vertical_grid, topology=(Flat, Flat, Bounded))

    # Estimate noise covariance based on discrepency between LES with different resolution
    coarse_Γ = estimate_noise_covariance(coarse_regrid; times, tke_weight, suite)
    fine_Γ = estimate_noise_covariance(fine_regrid; times, tke_weight, suite)

    #####
    ##### Single-resolution calibrations
    #####

    times = [2hours, 12hours, 18hours, 24hours]
    inverse_problem_kwargs = (; suite, free_parameters, Nensemble, architecture, closure)
    resampler = Resampler(; resample_failure_fraction, acceptable_failure_fraction)

    fine_ip = lesbrary_inverse_problem(fine_regrid; times, Δt=fine_Δt, inverse_problem_kwargs...)
    coarse_ip = lesbrary_inverse_problem(coarse_regrid; times, Δt=coarse_Δt, inverse_problem_kwargs...)

    pseudo_stepping = Kovachki2018InitialConvergenceRatio(; initial_convergence_ratio)
    #pseudo_stepping = ConstantConvergence(initial_convergence_ratio)
    fine_eki = EnsembleKalmanInversion(fine_ip; resampler, pseudo_stepping, noise_covariance=fine_Γ)
    coarse_eki = EnsembleKalmanInversion(coarse_ip; resampler, pseudo_stepping, noise_covariance=coarse_Γ)

    @info "Calibrating $name parameters to the fine grid..."
    iterate_until!(fine_eki, pseudotime_limit, max_iterations)

    datapath = joinpath(dir, "$(savename)_fine_calibration.jld2")
    @save datapath iteration_summaries=fine_eki.iteration_summaries

    @info "Calibrating $name parameters to the coarse grid..."
    iterate_until!(coarse_eki, pseudotime_limit, max_iterations)

    datapath = joinpath(dir, "$(savename)_coarse_calibration.jld2")
    @save datapath iteration_summaries=coarse_eki.iteration_summaries

    #####
    ##### "Combined" calibration to multiple resolutions
    #####

    pseudo_stepping = Kovachki2018InitialConvergenceRatio(; initial_convergence_ratio)
    fine_ip = lesbrary_inverse_problem(fine_regrid; times, Δt=fine_Δt, inverse_problem_kwargs...)
    preliminary_eki = EnsembleKalmanInversion(fine_ip; resampler, pseudo_stepping, noise_covariance=fine_Γ)

    @info "Performing a preliminary calibration of $name to the fine grid..."
    iterate_until!(preliminary_eki, pseudotime_limit / 4, max_iterations)

    final_pseudotime = preliminary_eki.pseudotime + pseudotime_limit
    times = collect(range(2hours, stop=24hours, length=4))
    weights = (1.0, coarse_regrid.Nz / fine_regrid.Nz)
    batched_ip = BatchedInverseProblem(coarse_ip, fine_ip; weights)

    combined_Γ = Matrix(BlockDiagonal([weights[1] * coarse_Γ, weights[2] * fine_Γ]))
    pseudo_stepping = Kovachki2018InitialConvergenceRatio(; initial_convergence_ratio)

    fine_ip = lesbrary_inverse_problem(fine_regrid; times, Δt=fine_Δt, inverse_problem_kwargs...)
    coarse_ip = lesbrary_inverse_problem(coarse_regrid; times, Δt=coarse_Δt, inverse_problem_kwargs...)

    eki = EnsembleKalmanInversion(batched_ip; resampler, pseudo_stepping,
                                  noise_covariance = combined_Γ,
                                  unconstrained_parameters = preliminary_eki.unconstrained_parameters)

    eki.pseudotime = preliminary_eki.pseudotime

    progress = calibration_progress_figure(eki)
    display(progress)
    @show eki.iteration_summaries[end]

    @info "Now for the main event..."

    while eki.pseudotime < final_pseudotime && eki.iteration < max_iterations
        iterate!(eki)
        progress = calibration_progress_figure(eki)
        latest_summary = eki.iteration_summaries[end]
        display(progress)
        @show latest_summary
    end

    display(calibration_progress_figure(eki))
    @show eki.iteration_summaries[end]

    datapath = joinpath(dir, "$(savename)_combined_calibration.jld2")
    @save datapath iteration_summaries=eki.iteration_summaries

    return nothing
end

