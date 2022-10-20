include("calibrate_parameter_set.jl")

#=
start_time = time_ns()
calibrate_parameter_set("shear_nemo_like_conv_adj")
elapsed = 1e-9 * (time_ns() - start_time)
@info "Calibration took "* prettytime(elapsed)
=#

sets_to_calibrate = [
    "shear_nemo_like",
    "shear_nemo_like_conv_adj",
    "goldilocks",
    "goldilocks_conv_adj",
    "complex",
    "complex_conv_adj",
]

for set in sets_to_calibrate
    calibrate_parameter_set(set)
end

