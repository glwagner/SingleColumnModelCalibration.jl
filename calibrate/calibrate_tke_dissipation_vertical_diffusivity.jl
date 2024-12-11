using Oceananigans
using Oceananigans.Units
using JLD2
using Printf

using Oceananigans.TurbulenceClosures: TKEDissipationVerticalDiffusivity

using ParameterEstimocean: iterate!

using SingleColumnModelCalibration:
    build_ensemble_kalman_inversion,
    generate_filepath,
    parameter_sets

#####
##### Set up the problem
#####

# Parameters of the calibration problem
closure         = TKEDissipationVerticalDiffusivity()
architecture    = GPU()
Nensemble       = 1000
simulation_Δt   = 30.0
stop_pseudotime = 1e3

# Choose the parameter set to calibrate.
# There must be corresponding entries in the following dictionaries, which
# are defined in src/parameter_sets.jl:
#
#   * SingleColumnModelCalibration.parameter_sets
#   * SingleColumnModelCalibration.dependent_parameter_sets
#   * SingleColumnModelCalibration.boundary_library
#
parameter_set = "pretty_simple_stabilities"
label = ""

# Data to calibrate against
suite_parameters = [
    (name = "12_hour_suite", stop_time=12hours),
    (name = "24_hour_suite", stop_time=24hours),
    (name = "48_hour_suite", stop_time=48hours),
]

# Grids to include
grid_parameters = [
    (size=32, z=(-256, 0)),
    (size=64, z=(-256, 0)),
    (size=128, z=(-256, 0)),
]

filename = "TKEDissipationVerticalDiffusivity_$(parameter_set)$(label)"
filepath = generate_filepath(; filename,
                             Nensemble,
                             suite_parameters,
                             grid_parameters,
                             stop_pseudotime,
                             Δt = simulation_Δt,
                             dir = ".")
 
@info "Saving data to $filepath"

#####
##### Build and run Ensemble Kalman Inversion
#####

eki = build_ensemble_kalman_inversion(closure, parameter_set;
                                      architecture,
                                      Nensemble,
                                      grid_parameters,
                                      suite_parameters,
                                      Δt = simulation_Δt)

   
# Run the calibration
max_iterations = 1000

while (eki.pseudotime < stop_pseudotime) && (eki.iteration < max_iterations)
    @time iterate!(eki)

    if eki.iteration % 5 == 0
        @show eki.iteration_summaries[end]

        rm(filepath; force=true)
        
        @info "Saving data to $filepath..."
        file = jldopen(filepath, "a+")
        file["stop_pseudotime"] = stop_pseudotime
        file["iteration_summaries"] = eki.iteration_summaries
        close(file)
    end
end

