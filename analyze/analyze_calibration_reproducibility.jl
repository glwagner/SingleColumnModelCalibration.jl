using Oceananigans
using ParameterEstimocean
using CairoMakie
using ElectronDisplay
using JLD2
using Printf
using Statistics

dir = joinpath("..", "results") #, "archive")
prefix = "shear_constant_Pr_conv_adj_"
suffix = "_Nens100_Δt1200_Nz32_Nz64_12_hour_suite_24_hour_suite_48_hour_suite.jld2"
Nruns = 10

summaries = []
min_errors = Float64[]
error_series = []
best_parameters = []

for i = 1:Nruns
    filename = string(prefix, i, suffix)
    filepath = joinpath(dir, filename)
    file = jldopen(filepath)
    push!(summaries, file["iteration_summaries"])
    close(file)

    last_summary = summaries[end][end]
    min_err, k_best = findmin(last_summary.mean_square_errors)
    push!(min_errors, min_err)

    push!(best_parameters, last_summary.parameters[k_best])

    @printf "%.3e " min_err
    for (k, p) in zip(keys(last_summary.parameters[k_best]),
                      values(last_summary.parameters[k_best]))

        @printf "%s: %.2e " k p
    end
    println(" ")

    errors = zeros(length(summaries[end]))
    for (i, summary) in enumerate(summaries[end])
        min_err, k = findmin(filter(!isnan, summary.mean_square_errors))
        errors[i] = min_err
    end

    push!(error_series, errors)
end

names = keys(first(best_parameters))
@show names

mean_parameters = []
std_parameters = []
for name in names
    best_param = [best_parameters[i][name] for i = 1:Nruns]
    push!(mean_parameters, mean(best_param))
    push!(std_parameters, std(best_param))
end

@show [collect(names) mean_parameters std_parameters std_parameters ./ mean_parameters]

fig = Figure()
ax = Axis(fig[1, 1])
scatter!(ax, min_errors)

ax2 = Axis(fig[1, 2])
for i = 1:Nruns
    lines!(ax2, error_series[i])
    scatter!(ax2, error_series[i])
end

display(fig)
