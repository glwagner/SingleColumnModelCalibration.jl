module SingleColumnModelCalibration

using LinearAlgebra
using BlockDiagonals
using Distributions
using Printf
using DataDeps
using JLD2
using OffsetArrays

#using CairoMakie
#using ElectronDisplay

using Oceananigans
using Oceananigans.Units

using ParameterEstimocean
using ParameterEstimocean: Transformation
using ParameterEstimocean.Observations: forward_map_names
using ParameterEstimocean.InverseProblems: BatchedInverseProblem
using ParameterEstimocean.PseudoSteppingSchemes: Kovachki2018InitialConvergenceRatio
using ParameterEstimocean.EnsembleKalmanInversions: best_next_best, nanminimum, IterationSummary

# Coarse grid used by ECCO
ecco_vertical_grid = [-256.0,
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

cases = ["free_convection",
         "strong_wind_weak_cooling",
         "med_wind_med_cooling",
         "weak_wind_strong_cooling",
         "strong_wind",
         "strong_wind_no_rotation"]

# Weight by number of fields, excluding TKE
observation_weights = Dict(
    "free_convection"          => 1.0,
    "strong_wind_weak_cooling" => 1/3,
    "med_wind_med_cooling"     => 1/3,
    "weak_wind_strong_cooling" => 1/3,
    "strong_wind"              => 1/3,
    "strong_wind_no_rotation"  => 1/2,
)

include("parameter_sets.jl")
include("lesbrary_inverse_problem.jl")
include("calibration_progress_figure.jl")
include("calibrate_parameter_set.jl")

#####
##### Utils
#####

function load_summaries(filepath)
    file = jldopen(filepath)
    summaries = file["iteration_summaries"]
    close(file)
    return summaries
end

function get_best_parameters(summary::IterationSummary)
    _, k = findmin(summary.mean_square_errors)
    return summary.parameters[k]
end

# Fine best *global* parameters
function get_best_parameters(summaries::AbstractVector{<:IterationSummary})
    best_parameters = summaries[0].mean_square_errors[1]
    best_mse = Inf

    for summary in summaries
        mse, k = findmin(summary.mean_square_errors)
        parameters = summary.parameters[k]
        best_parameters = ifelse(mse < best_mse, parameters, best_parameters)
    end

    return best_parameters
end

end # module
