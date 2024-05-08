using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using SingleColumnModelCalibration: calibrate_parameter_set, finitefindmin
using JLD2
using Printf
using Statistics
using LinearAlgebra
using GLMakie

set_theme!(Theme(fontsize=24))

name = "constant_Pr_no_shear"
suffix = "Nens400_Δt1200_τ10000_Nz32_Nz64_12_hour_suite_24_hour_suite_48_hour_suite.jld2"
dataset_filename = "calibration_summary_" * suffix

@load dataset_filename dataset

# Find best run
data = dataset[name]
Nrepeats = length(data[:final_minimum_objectives])
Φ★, best_run = findmin(data[:final_minimum_objectives])
Φþ, worst_run = findmax(data[:final_minimum_objectives])
Φ_best = data[:mean_parameter_objective_serieses][best_run]
optimal_parameters = data[:final_best_parameters][best_run]

best_first_particles = []
best_first_objectives = []
for r = 1:Nrepeats
    first_iteration_summary = first(data[:iteration_summaries][r])
    Φ₁, k₁ = findmin(o -> isnan(first(o)) ? Inf : first(o), first_iteration_summary.objective_values)
    θ = first_iteration_summary.parameters[k₁]
    push!(best_first_particles, θ)
    push!(best_first_objectives, Φ₁)
end

Φ₁★, r₁★ = findmin(best_first_objectives)
θ₁★ = best_first_particles[r₁★]

savedir = joinpath("..", "parameters")
savename = string(name, "_best_parameters.jld2")
savepath = joinpath(savedir, savename)

rm(savepath; force=true)
file = jldopen(savepath, "a+")
file["optimal_parameters"] = optimal_parameters
file["example_first_mean_parameters"] = first(data[:iteration_summaries][1]).ensemble_mean
file["first_best_parameters"] = θ₁★
file["final_best_parameters"] = data[:final_best_parameters]
file["final_mean_parameters"] = data[:final_mean_parameters]
file["first_mean_parameters"] = data[:final_mean_parameters]
close(file)

rr = sortperm(data[:final_minimum_objectives])

iteration_summaries = data[:iteration_summaries][best_run]
Ni = length(iteration_summaries)
Nc = length(iteration_summaries[0].ensemble_mean)
Ne = length(iteration_summaries[0].parameters)

iterations = 0:(Ni-1)

detC = []

for r = 1:Nrepeats
    summaries = data[:iteration_summaries][r]
    detCr = parent([det(cov(s.unconstrained_parameters, dims=2))^(1/Nc) for s in summaries])
    push!(detC, detCr)
end

X = [s.unconstrained_parameters for s in iteration_summaries]
X = parent(X)

K = [cov(x, dims=2) for x in X]
E = [eigen(k) for k in K]

Ni = length(E)
Nλ = length(first(E).values)
λ = [E[j].values for i=1:Nλ, j=1:Ni]
λ = [E[j].values[i] for i=1:Nλ, j=1:Ni]

fig = Figure(resolution=(1200, 700))

xticks = LogTicks(-3:4)
axV = Axis(fig[1, 1]; xticks, xscale=log10, xlabel="Pseudotime", ylabel="Ensemble size")
axE = Axis(fig[1, 2]; xticks, xscale=log10, xlabel="Pseudotime", ylabel=L"\lambda / \max(\lambda)")

xlims!(axV, 4e-3, 1e4)
xlims!(axE, 4e-3, 1e4)

for r in rr
    pseudotimes = parent([summ.pseudotime for summ in data[:iteration_summaries][r]])
    ln = lines!(axV, pseudotimes[2:end], detC[r][2:end], linewidth=3)
    ln.color[] = (ln.color.val, 0.6)
end

pseudotimes = parent([summ.pseudotime for summ in data[:iteration_summaries][best_run]])
lines!(axV, pseudotimes[2:end], detC[best_run][2:end], linewidth=8, color=(:black, 0.4))

for i = 1:Nλ-1
    ln = lines!(axE, pseudotimes[2:end], λ[end-i, 2:end] ./ λ[end, 2:end], linewidth=8) 
    ln.color[] = (ln.color.val, 0.8)
end

display(fig)
