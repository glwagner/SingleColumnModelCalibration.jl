using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using ParameterEstimocean.InverseProblems: BatchedInverseProblem
using Printf
using JLD2

# using GLMakie

# Import "lesbrary_inverse_problem" + plotting utility "calibration_progress_figure"
include("multi_resolution_calibration_utilities.jl")

# Coarse grid used by ECCO
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

# "Fine" grid
fine_regrid = RectilinearGrid(size=48; z=(-256, 0), topology=(Flat, Flat, Bounded))

# Batch the inverse problems
times = [6hours, 24hours]
Nensemble = 400
architecture = GPU()
ip = lesbrary_inverse_problem(fine_regrid; times, Nensemble, architecture)
resampler = Resampler(resample_failure_fraction=0.0, acceptable_failure_fraction=1.0)

@info "Performing some preliminary iterations with the coarse model..."
preliminary_eki = EnsembleKalmanInversion(ip; resampler, pseudo_stepping = ConstantConvergence(0.8))
iterate!(preliminary_eki, iterations=10)

display(calibration_progress_figure(preliminary_eki))
@show preliminary_eki.iteration_summaries[end]

@info "Now for the main event..."
times = [12hours, 24hours]
weights = (2, 1)
coarse_ip = lesbrary_inverse_problem(coarse_regrid; times, Nensemble, architecture)
fine_ip   = lesbrary_inverse_problem(fine_regrid;   times, Nensemble, architecture)
batched_ip = BatchedInverseProblem(coarse_ip, fine_ip; weights)
eki = EnsembleKalmanInversion(batched_ip; resampler, pseudo_stepping = ConstantConvergence(0.8),
                              unconstrained_parameters = preliminary_eki.unconstrained_parameters)

eki.pseudotime = preliminary_eki.pseudotime

display(calibration_progress_figure(eki))
@show eki.iteration_summaries[end]

# Hierarchical calibration, moving the start time back
for start_time in [22hours, 18hours, 12hours, 6hours, 2hours]
    times[1] = start_time
    coarse_ip = lesbrary_inverse_problem(coarse_regrid; times, Nensemble, architecture)
    fine_ip   = lesbrary_inverse_problem(fine_regrid;   times, Nensemble, architecture)
    batched_ip = BatchedInverseProblem(coarse_ip, fine_ip; weights)
    eki.inverse_problem = batched_ip

    for i = 1:10
        @time iterate!(eki)
        progress = calibration_progress_figure(eki)
        display(progress)

        dir = "/home/greg/Projects/LocalOceanClosureCalibration/calibration_summaries"

        # Figure
        saveprefix = @sprintf("progress_%02d", eki.iteration)
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

