using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using SingleColumnModelCalibration: calibrate_parameter_set, finitefindmin
using GLMakie
#using CairoMakie
#using ElectronDisplay
using JLD2
using Printf
using Statistics
using LinearAlgebra

names = [
    #"constant_Pr_no_shear",
    #"variable_Pr",
    #"variable_Pr_conv_adj",
    "extended_stability_conv_adj",
    #"ri_based"
]

labels_dict = Dict(
    "constant_Pr_no_shear" => "Constant \n Pr \n (no shear)",
    "constant_Pr" => "Constant \n Pr",
    "variable_Pr" => "Variable \n Pr",
    "constant_Pr_no_shear_simple_conv_adj" => "Conv. adj. \n Constant \n Pr",
    "variable_Pr_conv_adj" => "Conv. Adj. \n Variable \n Pr",
    "extended_stability_conv_adj" => "Conv. Adj. \n Extended \n Stab",
    "ri_based" => "Ri based",
)

labels = [labels_dict[n] for n in names]

suffixes = [suffix]
# Nrepeats = 1

# Extract best parameters
for n = 1:length(names)
    suffix = suffixes[n]
    name = names[n]

    dataset_filename = "calibration_summary_" * suffix
    @load dataset_filename dataset
    run_data = dataset[name]

    Φ★, best_run = findmin(run_data[:final_minimum_objectives])
    optimal_parameters = run_data[:final_best_parameters][best_run]

    best_particles = []
    best_objectives = []
    for r = 1:Nrepeats
        first_iteration_summary = first(run_data[:iteration_summaries][r])
        Φ₁, k₁ = findmin(o -> isnan(first(o)) ? Inf : first(o), first_iteration_summary.objective_values)
        θ = first_iteration_summary.parameters[k₁]
        push!(best_particles, θ)
        push!(best_objectives, Φ₁)
    end

    _, r₁★ = findmin(best_objectives)
    θ₁★ = best_particles[r₁★]
    @show θ₁★ 

    savedir = joinpath("..", "parameters")
    savename = string(name, "_best_parameters.jld2")
    savepath = joinpath(savedir, savename)

    rm(savepath; force=true)
    file = jldopen(savepath, "a+")
    file["optimal_parameters"] = optimal_parameters
    file["example_first_mean_parameters"] = first(run_data[:iteration_summaries][1]).ensemble_mean
    file["first_best_parameters"] = θ₁★
    file["final_best_parameters"] = run_data[:final_best_parameters]
    file["final_mean_parameters"] = run_data[:final_mean_parameters]
    file["first_mean_parameters"] = run_data[:final_mean_parameters]
    close(file)
end


#####
##### Plot
#####

fig = Figure(resolution=(1200, 400))

tickindices = collect(1:length(names))

ax1 = Axis(fig[1, 1], ylabel="EKI objective", xticks=(tickindices, labels))#, yscale=log10)
ax2 = Axis(fig[1, 2], xlabel="EKI iteration", ylabel="EKI objective", yscale=log10)

colors = [:orange1, :black, :seagreen, :royalblue, :red]

for n = 1:length(names)
    suffix = suffixes[n]
    name = names[n]

    dataset_filename = "calibration_summary_" * suffix
    @load dataset_filename dataset
    d = dataset[name]

    color = (colors[n], 0.6)
    scatter!(ax1, n * ones(Nrepeats), d[:final_minimum_objectives]; color, marker=:star5, label=string(name, ", best"))
    scatter!(ax1, n * ones(Nrepeats), d[:final_mean_parameter_objectives]; color, marker=:circle, label=string(name, ", mean"))
    
    objective_scale = 1.0
    
    for i = 1:Nrepeats
        #@show d[:minimum_objective_serieses][i])
        @show d[:mean_parameter_objective_serieses][i]

        τ = [summary.pseudotime for summary in parent(d[:iteration_summaries][i])]
        Nτ = length(τ) 
        color = (colors[n], 0.3)
        ls = scatter!(ax2, [1, Nτ], d[:minimum_objective_serieses][i][[1, Nτ]]; marker=:star5, color)
        ls = lines!(ax2, d[:mean_parameter_objective_serieses][i]; color, linewidth=3)
    end
end

display(fig)

#####
##### Plot
#####

fig2 = Figure(resolution=(1200, 800))

ax1 = Axis(fig2[1, 1], ylabel="Cˡᵒu, Cʰⁱu", xticks=(tickindices, labels))
ax2 = Axis(fig2[1, 2], ylabel="Cˡᵒe, Cʰⁱe", xticks=(tickindices, labels))
ax3 = Axis(fig2[2, 1], ylabel="Cˡᵒc, Cʰⁱc", xticks=(tickindices, labels))
ax4 = Axis(fig2[2, 2], ylabel="CˡᵒD, CʰⁱD", xticks=(tickindices, labels))
ax5 = Axis(fig2[3, 1], ylabel="Cˢ", xticks=(tickindices, labels))
ax6 = Axis(fig2[3, 2], ylabel="Cˢ", xticks=(tickindices, labels))
ax7 = Axis(fig2[4, 1], ylabel="Cᶜ", xticks=(tickindices, labels))
ax8 = Axis(fig2[4, 2], ylabel="Cᵉ", xticks=(tickindices, labels))

#ylims!(ax1, 0.3, 1.3)
#ylims!(ax2, 0.3, 1.3)

#ylims!(ax3, 0.2, 0.8)
#ylims!(ax4, 0.2, 0.8)

for n = 1:length(names)
    suffix = suffixes[n]
    name = names[n]

    dataset_filename = "calibration_summary_" * suffix
    @load dataset_filename dataset
    d = dataset[name]

    color = (colors[n], 0.6)

    scatter!(ax1, n * ones(Nrepeats), map(C -> C.Cʰⁱu, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
    scatter!(ax1, n * ones(Nrepeats), map(C -> C.Cʰⁱu, d[:final_best_parameters]); color, marker=:star5, label=string(name, ", best"))

    scatter!(ax2, n * ones(Nrepeats), map(C -> C.Cʰⁱe, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
    scatter!(ax2, n * ones(Nrepeats), map(C -> C.Cʰⁱe, d[:final_best_parameters]); color, marker=:star5, label=string(name, ", best"))

    scatter!(ax3, n * ones(Nrepeats), map(C -> C.Cʰⁱc, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
    scatter!(ax3, n * ones(Nrepeats), map(C -> C.Cʰⁱc, d[:final_best_parameters]); color, marker=:star5, label=string(name, ", best"))

    scatter!(ax4, n * ones(Nrepeats), map(C -> C.CʰⁱD, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
    scatter!(ax4, n * ones(Nrepeats), map(C -> C.CʰⁱD, d[:final_best_parameters]); color, marker=:star5, label=string(name, ", best"))

    scatter!(ax5, n * ones(Nrepeats), map(C -> C.Cˢ, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
    scatter!(ax5, n * ones(Nrepeats), map(C -> C.Cˢ, d[:final_best_parameters]); color, marker=:star5, label=string(name, ", best"))

    try
        scatter!(ax6, n * ones(Nrepeats), map(C -> C.Cˢ, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax6, n * ones(Nrepeats), map(C -> C.Cˢ, d[:final_best_parameters]); color, marker=:star5, label=string(name, ", best"))
    catch
    end

    try
        scatter!(ax1, n * ones(Nrepeats) .- 0.25, map(C -> C.Cˡᵒu, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax1, n * ones(Nrepeats) .- 0.25, map(C -> C.Cˡᵒu, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    catch 
        scatter!(ax1, n * ones(Nrepeats) .- 0.25, map(C -> C.Cʰⁱu, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax1, n * ones(Nrepeats) .- 0.25, map(C -> C.Cʰⁱu, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    end

    try
        scatter!(ax2, n * ones(Nrepeats) .- 0.25, map(C -> C.Cˡᵒe, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax2, n * ones(Nrepeats) .- 0.25, map(C -> C.Cˡᵒe, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    catch 
        scatter!(ax2, n * ones(Nrepeats) .- 0.25, map(C -> C.Cʰⁱe, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax2, n * ones(Nrepeats) .- 0.25, map(C -> C.Cʰⁱe, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    end

    try
        scatter!(ax3, n * ones(Nrepeats) .- 0.25, map(C -> C.Cˡᵒc, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax3, n * ones(Nrepeats) .- 0.25, map(C -> C.Cˡᵒc, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    catch 
        scatter!(ax3, n * ones(Nrepeats) .- 0.25, map(C -> C.Cʰⁱc, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax3, n * ones(Nrepeats) .- 0.25, map(C -> C.Cʰⁱc, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    end

    try
        scatter!(ax4, n * ones(Nrepeats) .- 0.25, map(C -> C.CˡᵒD, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax4, n * ones(Nrepeats) .- 0.25, map(C -> C.CˡᵒD, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    catch 
        scatter!(ax4, n * ones(Nrepeats) .- 0.25, map(C -> C.CʰⁱD, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax4, n * ones(Nrepeats) .- 0.25, map(C -> C.CʰⁱD, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    end

    try
        scatter!(ax7, n * ones(Nrepeats), map(C -> C.Cᶜc, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax7, n * ones(Nrepeats), map(C -> C.Cᶜc, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))

        scatter!(ax7, n * ones(Nrepeats) .- 0.25, map(C -> C.Cᶜe, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax7, n * ones(Nrepeats) .- 0.25, map(C -> C.Cᶜe, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))

        scatter!(ax7, n * ones(Nrepeats) .- 0.5, map(C -> C.CᶜD, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax7, n * ones(Nrepeats) .- 0.5, map(C -> C.CᶜD, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    catch 
    end

    try
        scatter!(ax8, n * ones(Nrepeats), map(C -> C.Cᵉc, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax8, n * ones(Nrepeats), map(C -> C.Cᵉc, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
        scatter!(ax8, n * ones(Nrepeats) .- 0.25, map(C -> C.Cᵉe, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax8, n * ones(Nrepeats) .- 0.25, map(C -> C.Cᵉe, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
        scatter!(ax8, n * ones(Nrepeats) .- 0.5, map(C -> C.CᵉD, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax8, n * ones(Nrepeats) .- 0.5, map(C -> C.CᵉD, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    catch 
    end

end
  
display(fig2)

