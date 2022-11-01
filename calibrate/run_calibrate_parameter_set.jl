using Oceananigans
using Oceananigans.Units

using Oceananigans.TurbulenceClosures:
    RiBasedVerticalDiffusivity,
    CATKEVerticalDiffusivity

using SingleColumnModelCalibration: calibrate_parameter_set

#name = "ri_based"
#closure = RiBasedVerticalDiffusivity()

name = "nemo_like_conv_adj"
closure = CATKEVerticalDiffusivity()

start_time = time_ns()

coarse_eki, fine_eki, eki = calibrate_parameter_set(name, closure;
                                                    Nensemble = 100,
                                                    architecture = CPU(),
                                                    tke_weight = 0.0,
                                                    pseudotime_limit = 1.0,
                                                    coarse_Δt = 20minutes,
                                                    fine_Δt = 20minutes,
                                                    fine_Nz = 64,
                                                    fine_calibration = false,
                                                    coarse_calibration = false,
                                                    combined_calibration = true)

elapsed = 1e-9 * (time_ns() - start_time)

@info "Calibration took "* prettytime(elapsed)

