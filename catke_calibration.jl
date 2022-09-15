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

parameter_names = (
    :CᵂwΔ,  :Cᵂu★,
    :Cᴰ⁻,
    :Cˢc,   :Cˢu,  :Cˢe,
    :Cᵇc,   :Cᵇu,  :Cᵇe,
    :Cᴷc⁻,  :Cᴷu⁻, :Cᴷe⁻,
    # :Ku_adjustment,
    # :Kc_adjustment,
    # :Ke_adjustment,
    :Cᴰʳ, :CᴰRiᶜ, :CᴰRiʷ,
    #:Cᵟc,   :Cᵟu,  :Cᵟe,
    :Cᴷcʳ,  :Cᴷuʳ, :Cᴷeʳ,
    :CᴷRiᶜ, :CᴷRiʷ,
)

neutral_default_mixing_length_parameters =
    Dict(
         :Cᵟc => 0.5,
         :Cᵟu => 0.5,
         :Cᵟe => 0.5,
         :Cᴷcʳ => 0.0,
         :Cᴷuʳ => 0.0,
         :Cᴷeʳ => 0.0,
        )

neutral_default_tke_parameters =
    Dict(
         :Cᴰʳ => 0.0,
         :CᴰRiᶜ => 0.0,
         :CᴰRiʷ => 0.0,
        )

mixing_length = MixingLength(; neutral_default_mixing_length_parameters...)
turbulent_kinetic_energy_equation = TurbulentKineticEnergyEquation(; neutral_default_tke_parameters...)
closure = CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation)

convective_adjustment = ConvectiveAdjustmentVerticalDiffusivity(convective_κz=0.0)
non_ensemble_closure = nothing # convective_adjustment

free_parameters = FreeParameters(prior_library, names=parameter_names)

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
Nensemble = 800
architecture = GPU()
ip = lesbrary_inverse_problem(fine_regrid; free_parameters, times, Nensemble, architecture, closure, non_ensemble_closure)
resampler = Resampler(resample_failure_fraction=0.0, acceptable_failure_fraction=1.0)

@info "Performing some preliminary iterations with the coarse model..."
pseudo_stepping = ConstantConvergence(0.8)
preliminary_eki = EnsembleKalmanInversion(ip; resampler, pseudo_stepping)
iterate!(preliminary_eki, iterations=20)

display(calibration_progress_figure(preliminary_eki))
@show preliminary_eki.iteration_summaries[end]

@info "Now for the main event..."
times = [6hours, 24hours]
weights = (2, 1)
coarse_ip = lesbrary_inverse_problem(coarse_regrid; free_parameters, times, Nensemble, architecture, closure, non_ensemble_closure)
fine_ip   = lesbrary_inverse_problem(fine_regrid;   free_parameters, times, Nensemble, architecture, closure, non_ensemble_closure)
batched_ip = BatchedInverseProblem(coarse_ip, fine_ip; weights)

pseudo_stepping = Kovachki2018InitialConvergenceRatio(initial_convergence_ratio=0.8)
eki = EnsembleKalmanInversion(batched_ip; resampler, pseudo_stepping,
                              unconstrained_parameters = preliminary_eki.unconstrained_parameters)

eki.pseudotime = preliminary_eki.pseudotime

display(calibration_progress_figure(eki))
@show eki.iteration_summaries[end]

# Hierarchical calibration, moving the start time back
for start_time in [22hours, 18hours, 12hours, 6hours, 2hours]
    times[1] = start_time
    coarse_ip = lesbrary_inverse_problem(coarse_regrid; free_parameters, times, Nensemble, architecture, closure, non_ensemble_closure)
    fine_ip   = lesbrary_inverse_problem(fine_regrid;   free_parameters, times, Nensemble, architecture, closure, non_ensemble_closure)
    batched_ip = BatchedInverseProblem(coarse_ip, fine_ip; weights)
    eki.inverse_problem = batched_ip

    for i = 1:100
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

