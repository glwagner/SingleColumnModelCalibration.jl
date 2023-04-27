using Oceananigans
using Oceananigans.Units
using JLD2

using Oceananigans.TurbulenceClosures:
    RiBasedVerticalDiffusivity,
    CATKEVerticalDiffusivity

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    MixingLength,
    TurbulentKineticEnergyEquation

using ParameterEstimocean: iterate!

using SingleColumnModelCalibration:
    build_ensemble_kalman_inversion,
    generate_filepath,
    parameter_sets

grid_parameters = [
    (size=32, z=(-256, 0)),
    (size=128, z=(-256, 0)),
]

suite_parameters = [
    (name = "12_hour_suite", stop_time=12hours),
    (name = "24_hour_suite", stop_time=24hours),
    (name = "48_hour_suite", stop_time=48hours),
]

resultsdir = "../results"

# Other names:
# "constant_Pr"
# "constant_Pr_conv_adj"

#closure = RiBasedVerticalDiffusivity()
#name = "ri_based"
#name = "constant_Pr_no_shear"
#name = "variable_Pr"
name = "variable_Pr_conv_adj"

turbulent_kinetic_energy_equation = TurbulentKineticEnergyEquation(
    C⁻D   = 1.0,
    C⁺D   = 1.0,
    CᶜD   = 0.0,
    CᵉD   = 0.0,
    Cᵂu★  = 1.0,
    CᵂwΔ  = 1.0,
)

mixing_length = MixingLength(
    Cᵇ   = Inf, 
    Cᶜc  = 0.0,
    Cᶜe  = 0.0,
    Cᵉc  = 0.0,
    Cᵉe  = 0.0,
    Cˢᶜ  = 0.0,
    C⁻u  = 1.0,
    C⁺u  = 1.0,
    C⁻c  = 1.0,
    C⁺c  = 1.0,
    C⁻e  = 1.0,
    C⁺e  = 1.0,
    CRiʷ = 1.0,
    CRiᶜ = 0.0,
)

minimum_turbulent_kinetic_energy = 1e-6
minimum_convective_buoyancy_flux = 1e-11
closure = CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation,
                                   minimum_turbulent_kinetic_energy,
                                   minimum_convective_buoyancy_flux)

architecture = CPU()
resample_failure_fraction = 0.05
stop_pseudotime = 1e3
Nensemble = 400
Δt = 10minutes
Nrepeats = 10

repeat_ekis = []
for i = 1:Nrepeats
    eki = build_ensemble_kalman_inversion(closure, name;
                                          architecture,
                                          Nensemble,
                                          tke_weight = 0.0,
                                          Δt,
                                          grid_parameters,
                                          suite_parameters,
                                          resample_failure_fraction)
    push!(repeat_ekis, eki)
end

start_time = time_ns()

for i = 1:Nrepeats
    eki = repeat_ekis[i]

    logname = string(name, "_", i, ".txt")

    while eki.pseudotime < stop_pseudotime
        iterate!(eki)

        if eki.iteration % 10 == 0
            open(logname, "a") do io
                 show(io, "text/plain", eki.iteration_summaries[end])
                 write(io, '\n')
                 write(io, '\n')
            end

            @show eki.iteration_summaries[end]
        end
    end

    filename = string(name, "_", i)
    filepath = generate_filepath(; Δt,
                                 dir = resultsdir,
                                 suite_parameters,
                                 grid_parameters,
                                 stop_pseudotime,
                                 Nensemble,
                                 filename)


    rm(filepath; force=true)

    @info "Saving data to $filepath..."
    file = jldopen(filepath, "a+")
    file["resample_failure_fraction"] = resample_failure_fraction
    file["stop_pseudotime"] = stop_pseudotime
    file["iteration_summaries"] = eki.iteration_summaries
    close(file)
end

elapsed = 1e-9 * (time_ns() - start_time)

@info "Calibrating $name parameters took " * prettytime(elapsed)

