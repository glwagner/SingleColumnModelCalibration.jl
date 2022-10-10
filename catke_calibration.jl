using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using ParameterEstimocean.InverseProblems: BatchedInverseProblem
using ParameterEstimocean.PseudoSteppingSchemes: Kovachki2018InitialConvergenceRatio
using Printf
using JLD2

# using GLMakie

# Import "lesbrary_inverse_problem" + plotting utility "calibration_progress_figure"
include("multi_resolution_calibration_utilities.jl")
include("parameter_sets.jl")

dir = "/home/greg/Projects/LocalOceanClosureCalibration/calibration_summaries"
name = "complex_conv_adj"

non_ensemble_closure = nothing # convective_adjustment

free_parameters = FreeParameters(prior_library,
                                 names = parameter_sets[name],
                                 dependent_parameters = dependent_parameter_sets[name])
@show free_parameters

# Two grids: "coarse" with ECCO vertical resolution to z=-256 m, and a fine grid with 4m resolution
Nz_ecco = length(ecco_vertical_grid) - 1
coarse_regrid = RectilinearGrid(size=Nz_ecco, z=ecco_vertical_grid, topology=(Flat, Flat, Bounded))
fine_regrid   = RectilinearGrid(size=48; z=(-256, 0), topology=(Flat, Flat, Bounded))

fine_Δt = 5minutes
coarse_Δt = 10minutes

# Batch the inverse problems
times = [6hours, 12hours, 18hours, 24hours]
Nensemble = 400
tke_weight = 0.01
architecture = GPU()
inverse_problem_kwargs = (; free_parameters, Nensemble, architecture, closure, non_ensemble_closure)
ip = lesbrary_inverse_problem(fine_regrid; times, Δt=fine_Δt, inverse_problem_kwargs...)
resampler = Resampler(resample_failure_fraction=0.0, acceptable_failure_fraction=1.0)

@info "Performing some preliminary iterations with the coarse model..."
pseudo_stepping = ConstantConvergence(0.9)
preliminary_eki = EnsembleKalmanInversion(ip; resampler, pseudo_stepping)
iterate!(preliminary_eki, iterations=40)

display(calibration_progress_figure(preliminary_eki))
@show preliminary_eki.iteration_summaries[end]

@info "Now for the main event..."
times = [6hours, 12hours, 18hours, 24hours]
weights = (1.0, coarse_regrid.Nz / fine_regrid.Nz)
coarse_ip = lesbrary_inverse_problem(coarse_regrid; times, inverse_problem_kwargs...)
fine_ip   = lesbrary_inverse_problem(fine_regrid;   times, Δt=fine_Δt, inverse_problem_kwargs...)
batched_ip = BatchedInverseProblem(coarse_ip, fine_ip; weights)

pseudo_stepping = Kovachki2018InitialConvergenceRatio(initial_convergence_ratio=0.8)
eki = EnsembleKalmanInversion(batched_ip; resampler, pseudo_stepping,
                              unconstrained_parameters = preliminary_eki.unconstrained_parameters)

eki.pseudotime = preliminary_eki.pseudotime

display(calibration_progress_figure(eki))
@show eki.iteration_summaries[end]

# Hierarchical calibration, moving the start time back
pseudotimes = [1e-3, 1e-2, 1e1]
start_times = [23hours, 12hours, 2hours]

for (n, start_time) in enumerate(start_times)
    times .= range(start_time, stop=24hours, length=4)
    #times[1] = start_time
    coarse_ip = lesbrary_inverse_problem(coarse_regrid; times, inverse_problem_kwargs...)
    fine_ip   = lesbrary_inverse_problem(fine_regrid;   times, Δt=fine_Δt, inverse_problem_kwargs...)
    batched_ip = BatchedInverseProblem(coarse_ip, fine_ip; weights)
    eki.inverse_problem = batched_ip

    while eki.pseudotime < pseudotimes[n] || eki.iteration > 1000
        @time iterate!(eki)
        progress = calibration_progress_figure(eki)
        display(progress)

        # Figure
        saveprefix = @sprintf("%s_progress_%02d", name, eki.iteration)
        figname = saveprefix * ".png"
        figpath = joinpath(dir, figname)
        save(figpath, progress)

        # Iteration data
        latest_summary = eki.iteration_summaries[end]
        dataname = saveprefix * ".jld2"
        datapath = joinpath(dir, dataname)
        @save datapath latest_summary
    
        @show eki.iteration_summaries[end]
    end
end

datapath = joinpath(dir, "$(name)_calibration.jld2")
@save datapath iteration_summaries=eki.iteration_summaries

