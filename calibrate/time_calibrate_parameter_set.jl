using Oceananigans
using Oceananigans.Units
using JLD2

using Oceananigans.TurbulenceClosures:
    RiBasedVerticalDiffusivity,
    CATKEVerticalDiffusivity

using ParameterEstimocean: iterate!

using SingleColumnModelCalibration:
    build_ensemble_kalman_inversion,
    generate_filepath,
    parameter_sets

#name = "ri_based"
#closure = RiBasedVerticalDiffusivity()

# Other names:
# "constant_Pr"
# "constant_Pr_conv_adj"

name =  "constant_Pr_no_shear"
   
grid_parameters = [
    (size=32, z=(-256, 0)),
    (size=64, z=(-256, 0)),
]

suite_parameters = [
    (name = "12_hour_suite", stop_time=12hours),
    (name = "24_hour_suite", stop_time=24hours),
    (name = "48_hour_suite", stop_time=48hours),
]

resultsdir = "../results"

closure = CATKEVerticalDiffusivity()
architecture = CPU()
Nensemble = 400
resample_failure_fraction = 0.05
stop_pseudotime = 1.0
Δt = 20minutes

start_time = time_ns()
    
eki = build_ensemble_kalman_inversion(closure, name;
                                      architecture,
                                      Nensemble,
                                      Δt,
                                      grid_parameters,
                                      suite_parameters,
                                      resample_failure_fraction)

iterate!(eki)
@time iterate!(eki)

