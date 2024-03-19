using Oceananigans
using Oceananigans.Units
using JLD2

using Oceananigans.TurbulenceClosures:
    RiBasedVerticalDiffusivity,
    CATKEVerticalDiffusivity

using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities:
    CATKEMixingLength,
    CATKEEquation

using ParameterEstimocean: iterate!

using SingleColumnModelCalibration:
    build_ensemble_kalman_inversion,
    generate_filepath,
    parameter_sets

grid_parameters = [
    # (size=24, z=(-256, 0)),
    (size=32, z=(-256, 0)),
    (size=64, z=(-256, 0)),
    (size=128, z=(-256, 0)),
]

suite_parameters = [
    (name = "12_hour_suite", stop_time=12hours),
    (name = "24_hour_suite", stop_time=24hours),
    (name = "48_hour_suite", stop_time=48hours),
]

resultsdir = "../results"

turbulent_kinetic_energy_equation = CATKEEquation(
    CˡᵒD  = 1.0,
    CʰⁱD  = 1.0,
    CᶜD   = 0.0,
    CᵉD   = 0.0,
    Cᵂu★  = 1.0,
    CᵂwΔ  = 1.0,
    Cᵂϵ   = 0.0,
)

mixing_length = CATKEMixingLength(
    Cˢ   = Inf, 
    Cᵇ   = Inf, 
    Cᶜc  = 0.0,
    Cᶜe  = 0.0,
    Cᵉc  = 0.0,
    Cᵉe  = 0.0,
    Cˢᵖ  = 0.0,
    Cˡᵒu = 1.0,
    Cʰⁱu = 1.0,
    Cˡᵒc = 1.0,
    Cʰⁱc = 1.0,
    Cˡᵒe = 1.0,
    Cʰⁱe = 1.0,
    CRi⁰ = 1.0,
    CRiᵟ = 0.0,
)

minimum_turbulent_kinetic_energy = 1e-9
minimum_convective_buoyancy_flux = 1e-15
closure = CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation)

#                                   minimum_turbulent_kinetic_energy,
#                                   minimum_convective_buoyancy_flux)

name = "variable_Pr_conv_adj"
#name = "fixed_Ric"

# Other names:
# "variable_Pr"
# "constant_Pr"
# "constant_Pr_conv_adj"
#
# Another closure:
# closure = RiBasedVerticalDiffusivity()
# name = "ri_based"

architecture = CPU()
resample_failure_fraction = 0.1
stop_pseudotime = 1e4
max_iterations = Inf
Nensemble = 400
Δt = 5minutes
irepeat = try ARGS[1]; catch; 1; end
start_time = time_ns()

eki = build_ensemble_kalman_inversion(closure, name;
                                      start_time = 1hours,
                                      architecture,
                                      Nensemble,
                                      tke_weight = 0.0,
                                      Δt,
                                      grid_parameters,
                                      suite_parameters,
                                      resample_failure_fraction)

label = "convective_depth_default_dimensional_tight_priors"
logname = string(name, "_Nens", Nensemble, "_", irepeat, "_", label, ".txt")

filename = string(name, "_", irepeat)
filepath = generate_filepath(; Δt,
                             dir = resultsdir,
                             suite_parameters,
                             grid_parameters,
                             stop_pseudotime,
                             Nensemble,
                             filename)
    
filepath = filepath[1:end-5] * "_$label.jld2"

while (eki.pseudotime < stop_pseudotime) && (eki.iteration < max_iterations)
    iterate!(eki)

    if eki.iteration % 10 == 0
        open(logname, "a") do io
             show(io, "text/plain", eki.iteration_summaries[end])
             write(io, '\n')
             write(io, '\n')
        end

        @show eki.iteration_summaries[end]
    end

    if eki.iteration % 100 == 0
        rm(filepath; force=true)
        
        @info "Saving data to $filepath..."
        file = jldopen(filepath, "a+")
        file["resample_failure_fraction"] = resample_failure_fraction
        file["stop_pseudotime"] = stop_pseudotime
        file["iteration_summaries"] = eki.iteration_summaries
        close(file)
    end
end

elapsed = 1e-9 * (time_ns() - start_time)

@info "Calibrating $name parameters took " * prettytime(elapsed)
=#
