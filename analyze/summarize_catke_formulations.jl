using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using ParameterEstimocean
using ParameterEstimocean.Parameters: build_parameters_named_tuple
using LinearAlgebra
using Statistics
using SingleColumnModelCalibration:
    ecco_vertical_grid,
    cases,
    lesbrary_inverse_problem,
    prior_library,
    dependent_parameter_sets,
    parameter_sets,
    get_best_parameters,
    load_summaries

using GLMakie

set_theme!(Theme(fontsize=24))

dir = "../results"
include("results.jl")
closure = CATKEVerticalDiffusivity()

name = :ri_dependent_dissipation
filename = calibration_filenames[name]
filepath = joinpath(dir, filename)
summaries = load_summaries(filepath)
mean_θₙ = summaries[end].ensemble_mean
parameter_names = keys(mean_θₙ)
free_parameters = FreeParameters(prior_library, names=parameter_names)

θ_dict = Dict()

for name in keys(calibration_filenames)
    filename = calibration_filenames[name]
    filepath = joinpath(dir, filename)
    summaries = load_summaries(filepath)

    # mean or best?
    θ = summaries[end].ensemble_mean
    #θ = get_best_parameters(summaries)

    # Don't forget the dependent parameters
    free_parameters = FreeParameters(prior_library; names=keys(θ), dependent_parameters=dependent_parameter_sets[string(name)])
    θ_dict[name] = build_parameters_named_tuple(free_parameters, θ)
end

θ_vec = [θ_dict[name] for name in keys(calibration_filenames)]

# Two grids: "coarse" with ECCO vertical resolution to z=-256 m, and a fine grid with 4m resolution
Nz_ecco = length(ecco_vertical_grid) - 1
coarse_regrid = RectilinearGrid(size=Nz_ecco, z=ecco_vertical_grid, topology=(Flat, Flat, Bounded))
med_regrid   = RectilinearGrid(size=32; z=(-256, 0), topology=(Flat, Flat, Bounded))
fine_regrid   = RectilinearGrid(size=64; z=(-256, 0), topology=(Flat, Flat, Bounded))

# Batch the inverse problems
times = [2hours, 48hours]
Nensemble = length(θ_vec)
architecture = CPU()
inverse_problem_kwargs = (; free_parameters, Nensemble, architecture, closure, suite="two_day_suite", tke_weight=0.0)

Δt1 = 20minutes
Δt2 = 1minutes

coarse_ip_dt1 = lesbrary_inverse_problem(coarse_regrid; times, Δt=Δt1, inverse_problem_kwargs...)
coarse_ip_dt2 = lesbrary_inverse_problem(coarse_regrid; times, Δt=Δt2, inverse_problem_kwargs...)
med_ip_dt1    = lesbrary_inverse_problem(med_regrid; times,    Δt=Δt1, inverse_problem_kwargs...)
med_ip_dt2    = lesbrary_inverse_problem(med_regrid; times,    Δt=Δt2, inverse_problem_kwargs...)
fine_ip_dt1   = lesbrary_inverse_problem(fine_regrid; times,   Δt=Δt1, inverse_problem_kwargs...)
fine_ip_dt2   = lesbrary_inverse_problem(fine_regrid; times,   Δt=Δt2, inverse_problem_kwargs...)

G_fine_dt1   = forward_map(fine_ip_dt1, θ_vec)
G_fine_dt2   = forward_map(fine_ip_dt2, θ_vec)
G_med_dt1    = forward_map(med_ip_dt1, θ_vec)
G_med_dt2    = forward_map(med_ip_dt2, θ_vec)
G_coarse_dt1 = forward_map(coarse_ip_dt1, θ_vec)
G_coarse_dt2 = forward_map(coarse_ip_dt2, θ_vec)

y_fine   = observation_map(fine_ip_dt1)
y_med = observation_map(med_ip_dt1)
y_coarse = observation_map(coarse_ip_dt1)

# Number of data scenarios
Ncases = length(coarse_ip_dt1.observations.observations) # Number of observational cases

# Number of parameter sets
Nsets = length(calibration_filenames) # Number of parameter sets considered

"""
    analyze_errors(ip, y, G)

Return a "global" estimate of errors for each G[:, k].
"""
function analyze_errors(ip, y, G)

    Ny, Nθ = size(G)
    Nunits = sum(length(obs.forward_map_names) for obs in ip.observations.observations)

    representative_obs = ip.observations[1]
    Nz = size(representative_obs.grid, 3)

    e_cases = [zeros(Ncases, 1) for _ in 1:Nsets]
    iy₁ = 1 # First y-index pertaining to case

    for case in 1:Ncases
        obs = ip.observations[case]
        Nfields = length(obs.forward_map_names)
        Ny_fields = Ny / Nunits
        Ncase = Int(Nfields * Ny_fields)

        y_case = y[iy₁:iy₁ + Ncase - 1, :]
        G_case = G[iy₁:iy₁ + Ncase - 1, :]
        e_case = [norm(G_case[:, i] .- y_case) for i = 1:size(G, 2)]

        for s in 1:Nsets
            # Normalize by "amount" of data?
            # Or, we should actually do this normalization during training...
            e_cases[s][case] = e_case[s] / Nz
            #e_cases[s][case] = e_case[s]
        end
        
        iy₁ += Ncase
    end

    # Normalize by "amount" of data?
    # Or, we should actually do this normalization during training...
    e_global = [norm(G[:, i] .- y) for i = 1:Nsets] ./ Nz
    #e_global = [norm(G[:, i] .- y) for i = 1:Nsets]

    return e_global, e_cases
end

coarse_global_errors_dt1, coarse_case_errors_dt1 = analyze_errors(coarse_ip_dt1, y_coarse, G_coarse_dt1)
coarse_global_errors_dt2, coarse_case_errors_dt2 = analyze_errors(coarse_ip_dt2, y_coarse, G_coarse_dt2)
med_global_errors_dt1,    med_case_errors_dt1    = analyze_errors(med_ip_dt1, y_med, G_med_dt1)
med_global_errors_dt2,    med_case_errors_dt2    = analyze_errors(med_ip_dt2, y_med, G_med_dt2)
fine_global_errors_dt1, fine_case_errors_dt1     = analyze_errors(fine_ip_dt1, y_fine, G_fine_dt1)
fine_global_errors_dt2, fine_case_errors_dt2     = analyze_errors(fine_ip_dt2, y_fine, G_fine_dt2)

case_errors = [zeros(Ncases, 6) for _ in 1:Nsets]

for s = 1:Nsets
    case_errors[s][:, 1] .= coarse_case_errors_dt1[s]
    case_errors[s][:, 2] .= coarse_case_errors_dt2[s]
    case_errors[s][:, 3] .= med_case_errors_dt1[s]
    case_errors[s][:, 4] .= med_case_errors_dt2[s]
    case_errors[s][:, 5] .= fine_case_errors_dt1[s]
    case_errors[s][:, 6] .= fine_case_errors_dt2[s]
end

#####
##### Figure
#####

case_names = [replace(name, "_" => " ") for name in cases]
elims = (0.9, 1.1) .* mean(coarse_global_errors_dt1)

symbol_names = keys(calibration_filenames)
string_names = [replace(string(name), "_" => " ") for name in symbol_names]

scale = identity
cmin = minimum(minimum(filter(!isnan, scale.(e))) for e in case_errors)
cmax = maximum(maximum(filter(!isnan, scale.(e))) for e in case_errors)

@show cmin cmax

xlabels = [
    #"Δz ≈ 10 m \n Δt = 20 min",
    #"Δz ≈ 10 m \n Δt = 1 min",
    #"Δz = 8 m \n Δt = 20 min",
    #"Δz = 8 m \n Δt = 1 min",
    #"Δz = 4 m \n Δt = 20 min",
    #"Δz = 4 m \n Δt = 1 min",
    "≈10 m \n 20 min",
    "≈10 m \n 1 min",
    "8 m \n 20 min",
    "8 m \n 1 min",
    "4 m \n 20 min",
    "4 m \n 1 min",
]

fig = Figure(resolution=(2400, 1200))

for s = 1:Nsets
    yaxisposition = s < Nsets ? :left : :right
    ylabel = "Scaled error"
    ax = Axis(fig[1, s]; xticks=(3:3, string_names[s:s]), yaxisposition, ylabel)

    barplot!(ax, 0.5, coarse_global_errors_dt1[s])
    barplot!(ax, 1.5, coarse_global_errors_dt2[s])
    barplot!(ax, 2.5,    med_global_errors_dt1[s])
    barplot!(ax, 3.5,    med_global_errors_dt2[s])
    barplot!(ax, 4.5,   fine_global_errors_dt1[s])
    barplot!(ax, 5.5,   fine_global_errors_dt2[s])
    s > 1 && s < Nsets && hideydecorations!(ax)
    s > 1 && hidespines!(ax, :l)
    s < Nsets && hidespines!(ax, :r)
    hidespines!(ax, :t)

    hidedecorations!(ax, label=false, ticklabels=false, ticks=false, minorticks=false,
                     minorgrid=true, grid=true)

    Δ = cmax - cmin
    colorrange = (cmin + 0.1 * Δ, cmax - 0.1 * Δ)

    axc = Axis(fig[2, s], yticks=(1:Ncases, case_names), xticks=(1:6, xlabels))
    heatmap!(axc, scale.(case_errors[s]'); colormap=:reds, colorrange)

    s > 1 && s < Ncases && hideydecorations!(axc)
    #s == 1 && hidexdecorations!(axc)
end

display(fig)

