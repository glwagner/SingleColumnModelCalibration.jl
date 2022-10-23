include("calibrate_parameter_set.jl")

using Oceananigans.TurbulenceClosures: RiBasedVerticalDiffusivity

start_time = time_ns()

calibrate_parameter_set("basic_conv_adj",
                        Nensemble = 200,
                        tke_weight = 1e-2,
                        coarse_Δt = 10minutes,
                        fine_Nz = 64,
                        fine_calibration = false,
                        preliminary_combined_calibration = false,
                        coarse_calibration = false)

elapsed = 1e-9 * (time_ns() - start_time)
@info "Calibration took "* prettytime(elapsed)

#=
sets_to_calibrate = [
    "shear_nemo_like",
    "shear_nemo_like_conv_adj",
    "basic",
    "basic_conv_adj",
    "goldilocks",
    "goldilocks_conv_adj",
    "complex",
    "complex_conv_adj",
]

for set in sets_to_calibrate
    calibrate_parameter_set(set)
end
=#

#=
# Ri-based calibration
Ri_based_bounds = Dict()
Ri_based_bounds[:ν₀] = (0.0, 0.1)
Ri_based_bounds[:κ₀] = (0.0, 1.0)
Ri_based_bounds[:νᶜ] = (0.0, 2.0)
Ri_based_bounds[:κᶜ] = (0.0, 2.0)
Ri_based_bounds[:Ri₀] = (0.0, 1.0)
Ri_based_bounds[:Riᵟ] = (0.0, 1.0)

Ri_based_priors = Dict()
for p in keys(Ri_based_bounds)
    Ri_based_priors[p] = ScaledLogitNormal(; bounds=Ri_based_bounds[p])
end

name = "Ri_based_calibration"
#free_parameters = FreeParameters(Ri_based_priors, names=(:ν₀, :κ₀, :νᶜ, :κᶜ, :Ri₀ν, :Ri₀κ, :Riᵟν, :Riᵟκ))
free_parameters = FreeParameters(Ri_based_priors, names=(:ν₀, :κ₀, :κᶜ, :Ri₀, :Riᵟ))
closure = RiBasedVerticalDiffusivity()
calibrate_parameter_set(name; free_parameters, closure, preliminary_combined_calibration=false)
=#
