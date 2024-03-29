using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using SingleColumnModelCalibration: calibrate_parameter_set, finitefindmin
using JLD2
using Printf
using Statistics
using LinearAlgebra

using CairoMakie
using ElectronDisplay

# using GLMakie

set_theme!(Theme(fontsize=24))

name = "constant_Pr_no_shear"
suffix = "Nens400_Δt1200_τ1000_Nz32_Nz64_12_hour_suite_24_hour_suite_48_hour_suite.jld2"
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

max_Ni = maximum(length(summaries) for summaries in data[:iteration_summaries])

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

fig = Figure(resolution=(1200, 700))

xticks = LogTicks(-3:4)
axΦ = Axis(fig[1, 1:2]; xticks, xscale=log10, xlabel="Pseudotime", ylabel="Normalized \n EKI Objective")
axV = Axis(fig[1, 3:4]; xticks, xscale=log10, xlabel="Pseudotime", ylabel="Ensemble size")

axΦ2 = Axis(fig[1, 1:2], alignmode=Outside(10), tellheight=false, height=150, width=280, valign=:top, halign=:right,
            xticks=LogTicks(0:4), xscale=log10, backgroundcolor=(:black, 0.1),
            xticklabelsize = 18,  yticklabelsize = 18,
            xlabelsize = 18,  ylabelsize = 18,
            xlabel = "Pseudotime", ylabel="Normalized \n EKI Objective")
            
ylims!(axΦ2, 0.9, 1.5)
xlims!(axΦ2, 1, 1e3)

τlims = (4e-3, 2e3)
xticks = LogTicks([-2, 0, 2, 4])
axCc = Axis(fig[2, 1]; xticks, xscale=log10, ylabel=L"\mathbb{C}^{\ell}_c")
axCu = Axis(fig[3, 1]; xticks, xscale=log10, xlabel="Pseudotime", ylabel=L"\mathbb{C}^{\ell}_u")
axCe = Axis(fig[2, 2]; xticks, xscale=log10, ylabel=L"\mathbb{C}^{\ell}_e")
axCD = Axis(fig[3, 2]; xticks, xscale=log10, xlabel="Pseudotime", ylabel=L"\mathbb{C}^{\ell}_D")

axCb = Axis(fig[2, 3]; xticks, xscale=log10, ylabel=L"\mathbb{C}^{b}")
axCw = Axis(fig[3, 4]; xticks, xscale=log10, xlabel="Pseudotime", ylabel=L"\mathbb{C}^{||}_Q")
axC★ = Axis(fig[3, 3]; xticks, xscale=log10, xlabel="Pseudotime", ylabel=L"\mathbb{C}^{=}_Q")

hidexdecorations!(axCc; grid=false)
hidexdecorations!(axCb; grid=false)
hidexdecorations!(axCe; grid=false)

hidespines!(axΦ, :t, :r)
hidespines!(axV, :t, :r)

hidespines!(axCc, :t, :r, :b)
hidespines!(axCu, :t, :r)
hidespines!(axCb, :t, :r, :b)

hidespines!(axCe, :t, :r, :b)
hidespines!(axCD, :t, :r)
hidespines!(axCw, :t, :r)
hidespines!(axC★, :t, :r)

max_Ni = 80
xlims!(axΦ,  τlims...)
xlims!(axV,  τlims...)

xlims!(axCc, τlims...)
xlims!(axCu, τlims...)
xlims!(axCb, τlims...)
xlims!(axCe, τlims...)
xlims!(axCD, τlims...)
xlims!(axCw, τlims...)
xlims!(axC★, τlims...)

ylims!(axCc, 0.0, 1.0)
ylims!(axCu, 0.0, 1.0)
ylims!(axCb, 0.0, 1.0)

ylims!(axCe, 0.0, 10.0)
ylims!(axCD, 0.0, 10.0)
ylims!(axCw, 0.0, 20.0)
ylims!(axC★, 0.0, 2.0)

for r in rr #[[1, 2, 3, 4, 5]]
    Φ = data[:mean_parameter_objective_serieses][r]
    run_iterations = 0:(length(Φ)-1)
    pseudotimes = parent([summ.pseudotime for summ in data[:iteration_summaries][r]])

    ln = lines!(axΦ, pseudotimes[2:end], Φ[2:end] ./ Φ★, linewidth=3)
    ln.color[] = (ln.color.val, 0.6)

    ln = lines!(axΦ2, pseudotimes[2:end], Φ[2:end] ./ Φ★, linewidth=3)
    ln.color[] = (ln.color.val, 0.8)

    ln = lines!(axV, pseudotimes[2:end], detC[r][2:end], linewidth=3)
    ln.color[] = (ln.color.val, 0.6)
end

pseudotimes = parent([summ.pseudotime for summ in data[:iteration_summaries][best_run]])
lines!(axΦ, pseudotimes[2:end], Φ_best[2:end] ./ Φ★, linewidth=8, color=(:black, 0.4))
lines!(axΦ2, pseudotimes[2:end], Φ_best[2:end] ./ Φ★, linewidth=8, color=(:black, 0.6))
lines!(axV, pseudotimes[2:end], detC[best_run][2:end], linewidth=8, color=(:black, 0.4))

markersize = 5
scatter_α = 0.05
scatter_colors = [(:royalblue1, scatter_α),
                  (:red,        scatter_α),
                  (:seagreen,   scatter_α),
                  (:purple,   scatter_α),
                  (:purple,   scatter_α),
                  (:purple,   scatter_α),
                  (:purple,   scatter_α),
                  (:purple,   scatter_α),
                  (:purple,   scatter_α),
                  (:purple,   scatter_α)]

bcolor1 = (:royalblue1, 0.5)
bcolor2 = (:red,        0.5)

for (ir, r) in enumerate(rr[[1]])
    color = scatter_colors[ir]

    iteration_summaries = data[:iteration_summaries][r]
    Ni = length(iteration_summaries)
    
    C⁺u  = [iteration_summaries[i-1].parameters[p][:C⁺u]  for p = 1:Ne, i = 1:Ni]
    C⁺c  = [iteration_summaries[i-1].parameters[p][:C⁺c]  for p = 1:Ne, i = 1:Ni]
    C⁺D  = [iteration_summaries[i-1].parameters[p][:C⁺D]  for p = 1:Ne, i = 1:Ni]
    C⁺e  = [iteration_summaries[i-1].parameters[p][:C⁺e]  for p = 1:Ne, i = 1:Ni]
    CᵂwΔ = [iteration_summaries[i-1].parameters[p][:CᵂwΔ] for p = 1:Ne, i = 1:Ni]
    Cᵂu★ = [iteration_summaries[i-1].parameters[p][:Cᵂu★] for p = 1:Ne, i = 1:Ni]
    Cᵇ   = [iteration_summaries[i-1].parameters[p][:Cᵇ]   for p = 1:Ne, i = 1:Ni]
    
    for i = 2:Ni
        τ = iteration_summaries[i-1].pseudotime
        scatter!(axCu, τ * ones(Ne), C⁺u[:, i];  markersize, marker=:circle, color)
        scatter!(axCc, τ * ones(Ne), C⁺c[:, i];  markersize, marker=:circle, color)
        scatter!(axCD, τ * ones(Ne), C⁺D[:, i];  markersize, marker=:circle, color)
        scatter!(axCe, τ * ones(Ne), C⁺e[:, i];  markersize, marker=:circle, color)
        scatter!(axCw, τ * ones(Ne), CᵂwΔ[:, i]; markersize, marker=:circle, color)
        scatter!(axC★, τ * ones(Ne), Cᵂu★[:, i]; markersize, marker=:circle, color)
        scatter!(axCb, τ * ones(Ne), Cᵇ[:, i];   markersize, marker=:circle, color)
    end
end

display(fig)

save("summarize_minimalist_calibration.png", fig)

for (k, v) in zip(keys(optimal_parameters), values(optimal_parameters))
    @printf "% 6s %.5f \n" k v
end

