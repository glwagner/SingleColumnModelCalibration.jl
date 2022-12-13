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

set_theme!(Theme(fontsize=20))

name = "variable_Pr_conv_adj"
suffix = "Nens400_Δt1200_τ1000_Nz32_Nz64_12_hour_suite_24_hour_suite_48_hour_suite.jld2"
dataset_filename = "calibration_summary_" * suffix

@load dataset_filename dataset

# Find best run
data = dataset[name]
Nrepeats = length(data[:final_minimum_objectives])
Φ★, best_run = findmin(data[:final_minimum_objectives])
iteration_summaries = data[:iteration_summaries][best_run]
Φ_best = data[:mean_parameter_objective_serieses][best_run]

max_Ni = maximum(length(summaries) for summaries in data[:iteration_summaries])

Ni = length(iteration_summaries)
Nc = length(iteration_summaries[0].ensemble_mean)

iterations = 0:(Ni-1)

detC = []

for r = 1:Nrepeats
    summaries = data[:iteration_summaries][r]
    detCr = parent([det(cov(s.unconstrained_parameters, dims=2))^(1/Nc) for s in summaries])
    push!(detC, detCr)
end

C⁺u = [iteration_summaries[i-1].parameters[p][:C⁺u] for p = 1:Nc, i = 1:Ni]
C⁻u = [iteration_summaries[i-1].parameters[p][:C⁻u] for p = 1:Nc, i = 1:Ni]

C⁺c = [iteration_summaries[i-1].parameters[p][:C⁺c] for p = 1:Nc, i = 1:Ni]
C⁻c = [iteration_summaries[i-1].parameters[p][:C⁻c] for p = 1:Nc, i = 1:Ni]

C⁺D = [iteration_summaries[i-1].parameters[p][:C⁺D] for p = 1:Nc, i = 1:Ni]
C⁻D = [iteration_summaries[i-1].parameters[p][:C⁻D] for p = 1:Nc, i = 1:Ni]

C⁺e = [iteration_summaries[i-1].parameters[p][:C⁺e] for p = 1:Nc, i = 1:Ni]
C⁻e = [iteration_summaries[i-1].parameters[p][:C⁻e] for p = 1:Nc, i = 1:Ni]

Cᶜc = [iteration_summaries[i-1].parameters[p][:Cᶜc] for p = 1:Nc, i = 1:Ni]
Cᵉc = [iteration_summaries[i-1].parameters[p][:Cᵉc] for p = 1:Nc, i = 1:Ni]

CᵂwΔ = [iteration_summaries[i-1].parameters[p][:CᵂwΔ] for p = 1:Nc, i = 1:Ni]
Cᵂu★ = [iteration_summaries[i-1].parameters[p][:Cᵂu★] for p = 1:Nc, i = 1:Ni]

Cᵇ = [iteration_summaries[i-1].parameters[p][:Cᵇ] for p = 1:Nc, i = 1:Ni]

fig = Figure(resolution=(1200, 700))

axΦ = Axis(fig[1, 1:2], xlabel="EKI iteration", ylabel=L"\Phi \left ( \langle \mathcal{C} \rangle \right )")
axV = Axis(fig[1, 3:4], xlabel="EKI iteration", ylabel=L"\mathrm{det}|\mathcal{K}(\mathbf{C}, \mathbf{C})|^{1/P}")

axCc = Axis(fig[2, 1], ylabel=L"\mathbb{C}^{+}_c \, , \, \, \mathbb{C}^{-}_c")
axCu = Axis(fig[3, 1], xlabel="EKI iteration", ylabel=L"\mathbb{C}^{+}_u \, , \, \, \mathbb{C}^{-}_u")
axCb = Axis(fig[2, 2], ylabel=L"\mathbb{C}^{b}")
axCC = Axis(fig[3, 2], xlabel="EKI iteration", ylabel=L"\mathbb{C}^{c}_c \, , \, \, \mathbb{C}^{e}_c")

axCe = Axis(fig[2, 3], ylabel=L"\mathbb{C}^{+}_e \, , \, \, \mathbb{C}^{-}_e")
axCD = Axis(fig[3, 3], xlabel="EKI iteration", ylabel=L"\mathbb{C}^{+}_D \, , \, \, \mathbb{C}^{-}_D")
axCw = Axis(fig[2, 4], ylabel=L"\mathbb{C}^{W}_{wΔ}")
axC★ = Axis(fig[3, 4], xlabel="EKI iteration", ylabel=L"\mathbb{C}^{W}_{u★}")

hidexdecorations!(axCc; grid=false)
#hidexdecorations!(axCu; grid=false)
hidexdecorations!(axCb; grid=false)
#hidexdecorations!(axCC; grid=false)

hidexdecorations!(axCe; grid=false)
#hidexdecorations!(axCD; grid=false)
hidexdecorations!(axCw; grid=false)
#hidexdecorations!(axC★; grid=false)

hidespines!(axΦ, :t, :r)
hidespines!(axV, :t, :r)

hidespines!(axCc, :t, :r, :b)
hidespines!(axCu, :t, :r)
hidespines!(axCb, :t, :r, :b)
hidespines!(axCC, :t, :r)

hidespines!(axCe, :t, :r, :b)
hidespines!(axCD, :t, :r)
hidespines!(axCw, :t, :r, :b)
hidespines!(axC★, :t, :r)

max_Ni = 50
xlims!(axΦ, -1, max_Ni)
xlims!(axV, -1, max_Ni)
xlims!(axCc, -1, max_Ni)
xlims!(axCu, -1, max_Ni)
xlims!(axCC, -1, max_Ni)
xlims!(axCb, -1, max_Ni)

xlims!(axCe, -1, max_Ni)
xlims!(axCD, -1, max_Ni)
xlims!(axCw, -1, max_Ni)
xlims!(axC★, -1, max_Ni)

ylims!(axCc, 0.0, 1.0)
ylims!(axCu, 0.0, 1.0)
ylims!(axCC, 0.0, 10.0)
ylims!(axCb, 0.0, 2.0)

ylims!(axCe, 0.0, 10.0)
ylims!(axCD, 0.0, 10.0)
ylims!(axCw, 0.0, 10.0)
ylims!(axC★, 0.0, 10.0)

for r = 1:Nrepeats
    Φ = data[:mean_parameter_objective_serieses][r]
    run_iterations = 0:(length(Φ)-1)

    ln = lines!(axΦ, run_iterations, Φ, linewidth=3)
    ln.color[] = (ln.color.val, 0.6)

    ln = lines!(axV, run_iterations, detC[r], linewidth=3)
    ln.color[] = (ln.color.val, 0.6)
end

lines!(axΦ, iterations, Φ_best, linewidth=8, color=(:black, 0.4))
lines!(axV, iterations, detC[best_run], linewidth=8, color=(:black, 0.4))

function parameter_band!(ax, C, color; kw...)
    C_mean = mean(C, dims=1)[:]
    C_std = std(C, dims=1)[:]

    # C_hi = maximum(C, dims=1)[:] 
    # C_lo = minimum(C, dims=1)[:]
    #
    C_hi = @. C_mean + C_std
    C_lo = @. C_mean - C_std

    return band!(ax, iterations, C_hi, C_lo; color, overdraw=true, shading=true, transparency=true, kw...)
end

markersize = 5
scolor1 = (:royalblue1, 0.3)
scolor2 = (:red,        0.3)

bcolor1 = (:royalblue1, 0.5)
bcolor2 = (:red,        0.5)

for i = 1:Ni
    scatter!(axCu, (i-1) * ones(Nc), C⁺u[:, i]; markersize, marker=:circle,    color=scolor1)
    scatter!(axCu, (i-1) * ones(Nc), C⁻u[:, i]; markersize, marker=:utriangle, color=scolor2)
    
    scatter!(axCc, (i-1) * ones(Nc), C⁺c[:, i]; markersize, marker=:circle,    color=scolor1)
    scatter!(axCc, (i-1) * ones(Nc), C⁻c[:, i]; markersize, marker=:utriangle, color=scolor2)

    scatter!(axCD, (i-1) * ones(Nc), C⁺D[:, i]; markersize, marker=:circle,    color=scolor1)
    scatter!(axCD, (i-1) * ones(Nc), C⁻D[:, i]; markersize, marker=:utriangle, color=scolor2)

    scatter!(axCe, (i-1) * ones(Nc), C⁺e[:, i]; markersize, marker=:circle,    color=scolor1)
    scatter!(axCe, (i-1) * ones(Nc), C⁻e[:, i]; markersize, marker=:utriangle, color=scolor2)

    scatter!(axCC, (i-1) * ones(Nc), Cᶜc[:, i]; markersize, marker=:circle,    color=scolor1)
    scatter!(axCC, (i-1) * ones(Nc), Cᵉc[:, i]; markersize, marker=:utriangle, color=scolor2)

    scatter!(axCw, (i-1) * ones(Nc), CᵂwΔ[:, i]; markersize, marker=:circle, color=scolor1)

    scatter!(axC★, (i-1) * ones(Nc), Cᵂu★[:, i]; markersize, marker=:circle, color=scolor1)

    scatter!(axCb, (i-1) * ones(Nc), Cᵇ[:, i]; markersize, marker=:circle, color=scolor1)
end

bn1 = parameter_band!(axCc, C⁺c, bcolor1, label=L"\mathbb{C}^+_c")
bn2 = parameter_band!(axCc, C⁻c, bcolor2, label=L"\mathbb{C}^-_c")

Legend(fig[2, 1],
       [bn1, bn2],
       [L"\mathbb{C}^+_c", L"\mathbb{C}^-_c"],
       tellheight=false, tellwidth=false,
       halign=:right, valign=:top,
       orientation=:horizontal, framevisible=false)

bn1 = parameter_band!(axCu, C⁺u, bcolor1, label=L"\mathbb{C}^+_u")
bn2 = parameter_band!(axCu, C⁻u, bcolor2, label=L"\mathbb{C}^-_u")

Legend(fig[3, 1],
       [bn1, bn2],
       [L"\mathbb{C}^+_u", L"\mathbb{C}^-_u"],
       tellheight=false, tellwidth=false,
       halign=:right, valign=:top,
       orientation=:horizontal, framevisible=false)

bn1 = parameter_band!(axCD, C⁺D, bcolor1, label=L"\mathbb{C}^+_D")
bn2 = parameter_band!(axCD, C⁻D, bcolor2, label=L"\mathbb{C}^-_D")

Legend(fig[3, 3],
       [bn1, bn2],
       [L"\mathbb{C}^+_D", L"\mathbb{C}^-_D"],
       tellheight=false, tellwidth=false,
       halign=:right, valign=:top,
       orientation=:horizontal, framevisible=false)

bn1 = parameter_band!(axCe, C⁺e, bcolor1, label=L"\mathbb{C}^+_e")
bn2 = parameter_band!(axCe, C⁻e, bcolor2, label=L"\mathbb{C}^-_e")

Legend(fig[2, 3],
       [bn1, bn2],
       [L"\mathbb{C}^+_e", L"\mathbb{C}^-_e"],
       tellheight=false, tellwidth=false,
       halign=:right, valign=:top,
       orientation=:horizontal, framevisible=false)

bn1 = parameter_band!(axCC, Cᶜc, bcolor1, label=L"\mathbb{C}^c_c")
bn2 = parameter_band!(axCC, Cᵉc, bcolor2, label=L"\mathbb{C}^e_c")

Legend(fig[3, 2],
       [bn1, bn2],
       [L"\mathbb{C}^c_c", L"\mathbb{C}^e_c"],
       tellheight=false, tellwidth=false,
       halign=:right, valign=:top,
       orientation=:horizontal, framevisible=false)

parameter_band!(axCw, CᵂwΔ, bcolor1, label=L"\mathbb{C}^W_{wΔ}")
parameter_band!(axC★, Cᵂu★, bcolor1, label=L"\mathbb{C}^W_{u★}")
parameter_band!(axCb, Cᵇ,   bcolor1, label=L"\mathbb{C}^b")

display(fig)

save("summarize_best_calibration.png", fig)

