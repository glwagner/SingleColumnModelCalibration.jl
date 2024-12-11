using Oceananigans
using Oceananigans.Units
using JLD2
using ParameterEstimocean
using SingleColumnModelCalibration
using Printf

using Oceananigans.TurbulenceClosures:
    RiBasedVerticalDiffusivity,
    TKEDissipationVerticalDiffusivity,
    CATKEVerticalDiffusivity

using ParameterEstimocean: iterate!

using SingleColumnModelCalibration:
    build_ensemble_kalman_inversion,
    generate_filepath,
    parameter_sets

grid_parameters = [
    (size=32, z=(-256, 0)),
    #(size=64, z=(-256, 0)),
    #(size=128, z=(-256, 0)),
]

suite_parameters = [
    (name = "12_hour_suite", stop_time=12hours),
    (name = "24_hour_suite", stop_time=24hours),
    (name = "48_hour_suite", stop_time=48hours),
]

resultsdir = "../results"

closure = CATKEVerticalDiffusivity(tke_time_step=nothing)
name = "extended_stability_conv_adj"

architecture = CPU()
resample_failure_fraction = 0.1
Nensemble = 1
Δt = 30.0

eki = build_ensemble_kalman_inversion(closure, name;
                                      start_time = 10minutes,
                                      architecture,
                                      Nensemble,
                                      tke_weight = 0.0,
                                      Ntimes = 2,
                                      Δt,
                                      grid_parameters,
                                      suite_parameters,
                                      resample_failure_fraction)


#=
dir = joinpath("..", "results")
filename = "extended_stability_conv_adj_1_Nens1000_Δt60_τ10000_Nz32_Nz64_Nz128_12_hour_suite_24_hour_suite_48_hour_suite_with_tracer.jld2"
filepath = joinpath(dir, filename)

file = jldopen(filepath)
summaries = file["iteration_summaries"]
close(file)

s1 = first(summaries)
Ne = length(s1.parameters)
Ni = length(summaries)
Np = length(first(s1.parameters))
parameters = zeros(Ni, Ne, Np)
objectives = zeros(Ni, Ne)

for i = 1:Ni
    s = summaries[i-1]
    objectives[i, :] .= map(o -> first(o), s.objective_values)
    
    for e = 1:Ne
        parameters[i, e, :] .= values(s.parameters[e])
    end
end

parametersfilename = "catke_parameters.jld2"
rm(parametersfilename, force=true)
file = jldopen(parametersfilename, "a+")
file["parameters"] = parameters
file["objectives"] = objectives
file["observations"] = eki.mapped_observations
file["noise_covariance"] = eki.noise_covariance
close(file)
=#
