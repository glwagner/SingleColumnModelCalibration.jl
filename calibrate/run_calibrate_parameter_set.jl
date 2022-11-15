using Oceananigans
using Oceananigans.Units

using Oceananigans.TurbulenceClosures:
    RiBasedVerticalDiffusivity,
    CATKEVerticalDiffusivity

using SingleColumnModelCalibration: calibrate_parameter_set, parameter_sets

#name = "ri_based"
#closure = RiBasedVerticalDiffusivity()

name = "constant_Pr_conv_adj"
#name = "variable_Pr"
closure = CATKEVerticalDiffusivity()
architecture = CPU()

grid_parameters = [
    (size=32, z=(-256, 0)),
    (size=64, z=(-256, 0))
]

suite_parameters = [
    (name = "12_hour_suite", stop_time=12hours),
    (name = "24_hour_suite", stop_time=24hours),
    (name = "48_hour_suite", stop_time=48hours),
]

resultsdir = "../results"

start_time = time_ns()

@sync begin
    for i = 1:10
        @async begin
            eki = calibrate_parameter_set(name, closure; architecture,
                                          prefixname = string(name, "_", i),
                                          Nensemble = 1000,
                                          pseudotime_limit = 4.0,
                                          plot_progress = false,
                                          resample_failure_fraction = 0.2,
                                          savedir = resultsdir,
                                          grid_parameters,
                                          suite_parameters)
        end
    end
end

elapsed = 1e-9 * (time_ns() - start_time)

@info "Calibrating $name parameters took " * prettytime(elapsed)

