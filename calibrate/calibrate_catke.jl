using Oceananigans
using Oceananigans.Units
using JLD2
using Printf

using Oceananigans.TurbulenceClosures:
    RiBasedVerticalDiffusivity,
    TKEDissipationVerticalDiffusivity,
    CATKEVerticalDiffusivity

using ParameterEstimocean: iterate!

using SingleColumnModelCalibration:
    build_ensemble_kalman_inversion,
    generate_filepath,
    parameter_sets

grid_parameters = [
    (size=32, z=(-256, 0)),
    (size=64, z=(-256, 0)),
    (size=128, z=(-256, 0)),
]

suite_parameters = [
    (name = "12_hour_suite", stop_time=12hours),
    (name = "24_hour_suite", stop_time=24hours),
    (name = "48_hour_suite", stop_time=48hours),
]

closure = CATKEVerticalDiffusivity(tke_time_step=nothing)
name = "extended_stability_conv_adj"
label = "uq_run"

architecture = GPU()
resample_failure_fraction = 0.1
stop_pseudotime = 10
max_iterations = 1000
Nensemble = 1000
Δt = 1minutes
start_time = time_ns()

eki = build_ensemble_kalman_inversion(closure, name;
                                      start_time = 10minutes,
                                      architecture,
                                      Nensemble,
                                      tke_weight = 0.0,
                                      Ntimes = 2,
                                      Δt,
                                      grid_parameters,
                                      suite_parameters,
                                      resample_failure_fraction)

filepath = generate_filepath(; Δt,
                             suite_parameters,
                             grid_parameters,
                             stop_pseudotime,
                             Nensemble,
                             dir = ".",
                             filename=name)
    
prefix = filepath[1:end-5]
filepath = prefix * "_$label.jld2"
rm(filepath, force=true)

file = jldopen(filepath, "a+")
file["stop_pseudotime"] = stop_pseudotime
file["noise_covariance"] = eki.noise_covariance
close(file)

while (eki.pseudotime < stop_pseudotime) && (eki.iteration < max_iterations)
    latest_summary = eki.iteration_summaries[end]
    @show latest_summary

    @info "Saving data to $filepath..."
    i = eki.iteration
    file = jldopen(filepath, "a+")
    file["iteration_summaries/$i"] = latest_summary
    file["forward_map_outputs/$i"] = eki.forward_map_output
    close(file)

    @time iterate!(eki)
end

elapsed = 1e-9 * (time_ns() - start_time)
@info "Calibrating $name parameters took " * prettytime(elapsed)

using Statistics
using Distributions
using LinearAlgebra
using ParameterEstimocean.EnsembleKalmanInversions: resampling_forward_map!, IterationSummary
using OffsetArrays

for update = 1:3
    start_time = time_ns()
    # Generate new ensemble with non-zero mean and unit variance in unconstrained space
    Xⁿ = eki.unconstrained_parameters
    μ = [mean(Xⁿ, dims=2)...]
    mn = MvNormal(μ, I)
    X₀ = rand(mn, Nensemble)
    eki.unconstrained_parameters .= X₀

    # Re-initialize EKI
    G = resampling_forward_map!(eki, X₀)
    summary = IterationSummary(eki, X₀, G)
    eki.iteration_summaries = OffsetArray([summary], -1)
    eki.iteration = 0
    eki.pseudotime = 0

    # Iterate again
    filepath = prefix * "_update$update" * "_$label.jld2"
    rm(filepath, force=true)

    file = jldopen(filepath, "a+")
    file["stop_pseudotime"] = stop_pseudotime
    file["unconstrained_prior_mean"] = μ
    file["noise_covariance"] = eki.noise_covariance
    close(file)

    while (eki.pseudotime < stop_pseudotime) && (eki.iteration < max_iterations)
        latest_summary = eki.iteration_summaries[end]
        @show latest_summary

        @info "Saving data to $filepath..."
        i = eki.iteration
        file = jldopen(filepath, "a+")
        file["iteration_summaries/$i"] = latest_summary
        file["forward_map_outputs/$i"] = eki.forward_map_output
        close(file)

        @time iterate!(eki)
    end

    elapsed = 1e-9 * (time_ns() - start_time)
    @info "Calibrating $name parameters took " * prettytime(elapsed)
end

