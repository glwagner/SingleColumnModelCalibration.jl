using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using ParameterEstimocean.InverseProblems: BatchedInverseProblem
using ParameterEstimocean.PseudoSteppingSchemes: Kovachki2018InitialConvergenceRatio
using Printf
using JLD2

# Import "lesbrary_inverse_problem" + plotting utility "calibration_progress_figure"
include("multi_resolution_calibration_utilities.jl")

#####
##### Setup
#####

dir = "/home/greg/Projects/LocalOceanClosureCalibration/calibration_summaries"
name = "goldilocks_conv_adj"
Nensemble = 200
tke_weight = 0.01
suite = "one_day_suite"
architecture = GPU()

#####
##### Build closure and free parameters
#####

include("parameter_sets.jl")

mixing_length = MixingLength(; neutral_default_mixing_length_parameters...)
turbulent_kinetic_energy_equation = TurbulentKineticEnergyEquation(; neutral_default_tke_parameters...)
closure = CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation)

free_parameters = FreeParameters(prior_library,
                                 names = parameter_sets[name],
                                 dependent_parameters = dependent_parameter_sets[name])

@show free_parameters

#####
##### Build grids and noise covariances for each grid
#####

times = range(2hours, stop=24hours, length=4)

# Two grids: "coarse" and "fine"
fine_Δt = 10minutes
fine_regrid = RectilinearGrid(size=24; z=(-256, 0), topology=(Flat, Flat, Bounded))

# Coarse grid has ECCO vertical resolution to z=-256 m
Nz_ecco = length(ecco_vertical_grid) - 1
coarse_Δt = 10minutes
coarse_regrid = RectilinearGrid(size=Nz_ecco, z=ecco_vertical_grid, topology=(Flat, Flat, Bounded))

# Estimate noise covariance based on discrepency between LES with different resolution
coarse_Γ = estimate_noise_covariance(coarse_regrid; times, tke_weight, suite)
fine_Γ = estimate_noise_covariance(fine_regrid; times, tke_weight, suite)

#####
##### Preliminary single-resolution calibration
#####

times = [2hours, 12hours, 18hours, 24hours]
inverse_problem_kwargs = (; suite, free_parameters, Nensemble, architecture, closure)
resampler = Resampler(resample_failure_fraction=0.2, acceptable_failure_fraction=1.0)

fine_ip = lesbrary_inverse_problem(fine_regrid; times, Δt=fine_Δt, inverse_problem_kwargs...)
coarse_ip = lesbrary_inverse_problem(coarse_regrid; times, Δt=coarse_Δt, inverse_problem_kwargs...)

pseudo_stepping = Kovachki2018InitialConvergenceRatio(initial_convergence_ratio=0.2)
fine_eki = EnsembleKalmanInversion(fine_ip; resampler, pseudo_stepping, noise_covariance=fine_Γ)

@info "Performing some preliminary iterations..."
while fine_eki.pseudotime < 0.2
    iterate!(fine_eki)
    fine_eki.iteration % 10 == 0 && @show(fine_eki.iteration_summaries[end])
end

display(calibration_progress_figure(fine_eki))
@show fine_eki.iteration_summaries[end]

#####
##### Second "transfer" calibration to multiple resolutions
#####

final_pseudotime = fine_eki.pseudotime + 1.0
times = collect(range(2hours, stop=24hours, length=4))
weights = (1.0, coarse_regrid.Nz / fine_regrid.Nz)
batched_ip = BatchedInverseProblem(coarse_ip, fine_ip; weights)

combined_Γ = Matrix(BlockDiagonal([weights[1] * coarse_Γ, weights[2] * fine_Γ]))
pseudo_stepping = Kovachki2018InitialConvergenceRatio(initial_convergence_ratio=0.2)

fine_ip = lesbrary_inverse_problem(fine_regrid; times, Δt=fine_Δt, inverse_problem_kwargs...)
coarse_ip = lesbrary_inverse_problem(coarse_regrid; times, Δt=coarse_Δt, inverse_problem_kwargs...)

eki = EnsembleKalmanInversion(batched_ip; resampler, pseudo_stepping,
                              noise_covariance = combined_Γ,
                              unconstrained_parameters = fine_eki.unconstrained_parameters)

eki.pseudotime = fine_eki.pseudotime

progress = calibration_progress_figure(eki)
display(progress)
@show eki.iteration_summaries[end]

@info "Now for the main event..."

while eki.pseudotime < final_pseudotime || eki.iteration < 1000
    iterate!(eki)
    progress = calibration_progress_figure(eki)
    latest_summary = eki.iteration_summaries[end]
    display(progress)
    @show latest_summary

    #=
    # Save figure
    saveprefix = @sprintf("%s_progress_%02d", name, eki.iteration)
    figname = saveprefix * ".png"
    figpath = joinpath(dir, figname)
    save(figpath, progress)

    # Iteration data
    dataname = saveprefix * ".jld2"
    datapath = joinpath(dir, dataname)
    @save datapath latest_summary
    =#
end


#=
# Hierarchical calibration, moving the start time back
pseudodurations = [1e-1, 1e-1, 1e0]
start_times = [23hours, 12hours, 2hours]

for (n, start_time) in enumerate(start_times)
    times .= range(start_time, stop=24hours, length=4)

    coarse_ip = lesbrary_inverse_problem(coarse_regrid; times, Δt=coarse_Δt, inverse_problem_kwargs...)
    fine_ip   = lesbrary_inverse_problem(fine_regrid;   times, Δt=fine_Δt, inverse_problem_kwargs...)

    batched_ip = BatchedInverseProblem(coarse_ip, fine_ip; weights)
    eki.inverse_problem = batched_ip

    next_pseudotime = pseudodurations[n] + eki.pseudotime
    while eki.pseudotime < next_pseudotime || eki.iteration > 1000
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
=#

display(calibration_progress_figure(eki))
@show eki.iteration_summaries[end]

datapath = joinpath(dir, "$(name)_calibration.jld2")
@save datapath iteration_summaries=eki.iteration_summaries
