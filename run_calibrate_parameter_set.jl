include("calibrate_parameter_set.jl")

using Oceananigans.TurbulenceClosures: RiBasedVerticalDiffusivity, viscosity_location

#sets_to_calibrate = ["shear_nemo_like"]

sets_to_calibrate = [
    #"shear_nemo_like",
    #"shear_nemo_like_conv_adj",
    #"basic",
    #"basic_conv_adj",
    #"goldilocks",
    #"goldilocks_conv_adj",
    "complex",
    "complex_conv_adj",
    "complex_dissipation_length",
    "complex_dissipation_length_conv_adj",
]

for set in sets_to_calibrate
    start_time = time_ns()
    
    calibrate_parameter_set(set,
                            Nensemble = 200,
                            tke_weight = 1e-2,
                            pseudotime_limit = 1.0,
                            coarse_Δt = 20minutes,
                            fine_Nz = 64,
                            fine_calibration = true,
                            coarse_calibration = true,
                            preliminary_combined_calibration = false)
    
    elapsed = 1e-9 * (time_ns() - start_time)

    @info "Calibration took "* prettytime(elapsed)
end

#=
# Ri-based calibration
Ri_based_bounds = Dict()
Ri_based_bounds[:ν₀]  = (0.0,  1.0)
Ri_based_bounds[:κ₀]  = (0.0,  1.0)
Ri_based_bounds[:κᶜ]  = (0.0, 10.0)
Ri_based_bounds[:Cᵉ]  = (0.0, 10.0)
Ri_based_bounds[:Ri₀] = (0.0,  1.0)
Ri_based_bounds[:Riᵟ] = (0.0,  1.0)

Ri_based_priors = Dict()
for p in keys(Ri_based_bounds)
    Ri_based_priors[p] = ScaledLogitNormal(; bounds=Ri_based_bounds[p])
end

name = "Ri_based_calibration"
free_parameters = FreeParameters(Ri_based_priors, names=(:ν₀, :κ₀, :κᶜ, :Cᵉ, :Ri₀, :Riᵟ))
#free_parameters = FreeParameters(Ri_based_priors, names=(:ν₀, :κ₀, :κᶜ, :Ri₀, :Riᵟ))
closure = RiBasedVerticalDiffusivity()

#@show viscosity_location(closure)

calibrate_parameter_set(name; free_parameters, closure,
                        Nensemble = 200,
                        tke_weight = 0.0,
                        pseudotime_limit = 1.0,
                        coarse_Δt = 20minutes,
                        fine_Nz = 64,
                        fine_calibration = false,
                        coarse_calibration = true,
                        combined_calibration = true,
                        preliminary_combined_calibration = false)

=#
