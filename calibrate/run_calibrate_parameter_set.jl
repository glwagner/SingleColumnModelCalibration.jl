using Oceananigans
using Oceananigans.Units

using Oceananigans.TurbulenceClosures:
    RiBasedVerticalDiffusivity,
    CATKEVerticalDiffusivity

using SingleColumnModelCalibration: calibrate_parameter_set, parameter_sets

#name = "ri_based"
#closure = RiBasedVerticalDiffusivity()

name = "shear_constant_Pr_conv_adj"
closure = CATKEVerticalDiffusivity()

grid_parameters = [
    (size=32, z=(-256, 0)),
    (size=64, z=(-256, 0))
]

suite_parameters = [
    (name = "12_hour_suite", stop_time=10hours),
    (name = "24_hour_suite", stop_time=24hours),
    (name = "48_hour_suite", stop_time=48hours),
]

resultsdir = "../results"

start_time = time_ns()

for i = 1:10
    eki = calibrate_parameter_set(name, closure;
                                  prefixname = string(name, "_", i),
                                  Nensemble = 100,
                                  pseudotime_limit = 1.0,
                                  plot_progress = false,
                                  resample_failure_fraction = 0.0,
                                  savedir = resultsdir,
                                  grid_parameters,
                                  suite_parameters)
end

elapsed = 1e-9 * (time_ns() - start_time)

@info "Calibrating $name parameters took " * prettytime(elapsed)

