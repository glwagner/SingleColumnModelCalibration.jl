using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using SingleColumnModelCalibration: calibrate_parameter_set, finitefindmin
using CairoMakie
using ElectronDisplay
using JLD2
using Printf
using Statistics
using LinearAlgebra

names = [
    "constant_Pr",
    "constant_Pr_conv_adj",
    "variable_Pr",
    "variable_Pr_conv_adj"
]

labels_dict = Dict(
    "constant_Pr" => "Constant \n Pr",
    "variable_Pr" => "Variable \n Pr",
    "constant_Pr_conv_adj" => "Conv. adj. \n Constant \n Pr",
    "variable_Pr_conv_adj" => "Conv. Adj. \n Variable \n Pr",
)

labels = [labels_dict[n] for n in names]

#suffix = "Nens100_Δt1200_τ1000_Nz32_Nz64_Nz128_12_hour_suite_24_hour_suite_48_hour_suite.jld2"
suffix = "Nens400_Δt1200_τ1000_Nz32_Nz64_12_hour_suite_24_hour_suite_48_hour_suite.jld2"
dataset_filename = "calibration_summary_" * suffix

@load dataset_filename dataset

# Save best parameters

# Extract best parameters
for n = 1:length(names)
    name = names[n]
    run_data = dataset[name]

    Φ★, best_run = findmin(run_data[:final_minimum_objectives])
    optimal_parameters = run_data[:final_best_parameters][best_run]

    savedir = joinpath("..", "parameters")
    savename = string(name, "_best_parameters.jld2")
    savepath = joinpath(savedir, savename)

    rm(savepath; force=true)
    file = jldopen(savepath, "a+")
    file["optimal_parameters"] = optimal_parameters
    file["final_best_parameters"] = run_data[:final_best_parameters]
    file["final_mean_parameters"] = run_data[:final_mean_parameters]
    close(file)
end


#####
##### Plot
#####

fig = Figure(resolution=(1200, 400))

tickindices = collect(1:length(names))

ax1 = Axis(fig[1, 1], ylabel="EKI objective", xticks=(tickindices, labels))#, yscale=log10)
ax2 = Axis(fig[1, 2], xlabel="EKI iteration", ylabel="EKI objective", yscale=log10)

colors = [:orange1, :black, :seagreen, :royalblue]

for n = 1:length(names)
    name = names[n]
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
        #ls = lines!(ax2, τ, d[:minimum_objective_serieses][i]; color)
        #ls = lines!(ax2, τ, d[:mean_parameter_objective_serieses][i][1:length(τ)]; color, linewidth=3)
        ls = scatter!(ax2, [1, Nτ], d[:minimum_objective_serieses][i][[1, Nτ]]; marker=:star5, color)
        ls = lines!(ax2, d[:mean_parameter_objective_serieses][i]; color, linewidth=3)
    end
end

display(fig)

#####
##### Plot
#####

fig2 = Figure(resolution=(1200, 800))

ax1 = Axis(fig2[1, 1], ylabel="C⁻u, C⁺u", xticks=(tickindices, labels))
ax2 = Axis(fig2[1, 2], ylabel="C⁻e, C⁺e", xticks=(tickindices, labels))
ax3 = Axis(fig2[2, 1], ylabel="C⁻c, C⁺c", xticks=(tickindices, labels))
ax4 = Axis(fig2[2, 2], ylabel="C⁻D, C⁺D", xticks=(tickindices, labels))
ax5 = Axis(fig2[3, 1], ylabel="Cᵇ", xticks=(tickindices, labels))
ax6 = Axis(fig2[3, 2], ylabel="Cˢ", xticks=(tickindices, labels))
ax7 = Axis(fig2[4, 1], ylabel="Cᶜ", xticks=(tickindices, labels))
ax8 = Axis(fig2[4, 2], ylabel="Cᵉ", xticks=(tickindices, labels))

#ylims!(ax1, 0.3, 1.3)
#ylims!(ax2, 0.3, 1.3)

#ylims!(ax3, 0.2, 0.8)
#ylims!(ax4, 0.2, 0.8)

for n = 1:length(names)
    name = names[n]
    d = dataset[name]
    color = (colors[n], 0.6)

    scatter!(ax1, n * ones(Nrepeats), map(C -> C.C⁺u, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
    scatter!(ax1, n * ones(Nrepeats), map(C -> C.C⁺u, d[:final_best_parameters]); color, marker=:star5, label=string(name, ", best"))

    scatter!(ax2, n * ones(Nrepeats), map(C -> C.C⁺e, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
    scatter!(ax2, n * ones(Nrepeats), map(C -> C.C⁺e, d[:final_best_parameters]); color, marker=:star5, label=string(name, ", best"))

    scatter!(ax3, n * ones(Nrepeats), map(C -> C.C⁺c, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
    scatter!(ax3, n * ones(Nrepeats), map(C -> C.C⁺c, d[:final_best_parameters]); color, marker=:star5, label=string(name, ", best"))

    scatter!(ax4, n * ones(Nrepeats), map(C -> C.C⁺D, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
    scatter!(ax4, n * ones(Nrepeats), map(C -> C.C⁺D, d[:final_best_parameters]); color, marker=:star5, label=string(name, ", best"))

    scatter!(ax5, n * ones(Nrepeats), map(C -> C.Cᵇ, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
    scatter!(ax5, n * ones(Nrepeats), map(C -> C.Cᵇ, d[:final_best_parameters]); color, marker=:star5, label=string(name, ", best"))

    try
        scatter!(ax6, n * ones(Nrepeats), map(C -> C.Cˢ, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax6, n * ones(Nrepeats), map(C -> C.Cˢ, d[:final_best_parameters]); color, marker=:star5, label=string(name, ", best"))
    catch
    end

    try
        scatter!(ax1, n * ones(Nrepeats) .- 0.25, map(C -> C.C⁻u, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax1, n * ones(Nrepeats) .- 0.25, map(C -> C.C⁻u, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    catch 
        scatter!(ax1, n * ones(Nrepeats) .- 0.25, map(C -> C.C⁺u, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax1, n * ones(Nrepeats) .- 0.25, map(C -> C.C⁺u, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    end

    try
        scatter!(ax2, n * ones(Nrepeats) .- 0.25, map(C -> C.C⁻e, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax2, n * ones(Nrepeats) .- 0.25, map(C -> C.C⁻e, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    catch 
        scatter!(ax2, n * ones(Nrepeats) .- 0.25, map(C -> C.C⁺e, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax2, n * ones(Nrepeats) .- 0.25, map(C -> C.C⁺e, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    end

    try
        scatter!(ax3, n * ones(Nrepeats) .- 0.25, map(C -> C.C⁻c, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax3, n * ones(Nrepeats) .- 0.25, map(C -> C.C⁻c, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    catch 
        scatter!(ax3, n * ones(Nrepeats) .- 0.25, map(C -> C.C⁺c, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax3, n * ones(Nrepeats) .- 0.25, map(C -> C.C⁺c, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    end

    try
        scatter!(ax4, n * ones(Nrepeats) .- 0.25, map(C -> C.C⁻D, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax4, n * ones(Nrepeats) .- 0.25, map(C -> C.C⁻D, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
    catch 
        scatter!(ax4, n * ones(Nrepeats) .- 0.25, map(C -> C.C⁺D, d[:final_mean_parameters]); color, marker=:circle, label=string(name, ", mean"))
        scatter!(ax4, n * ones(Nrepeats) .- 0.25, map(C -> C.C⁺D, d[:final_best_parameters]); color, marker=:star5,  label=string(name, ", best"))
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

