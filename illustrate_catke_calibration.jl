using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using ParameterEstimocean.InverseProblems: BatchedInverseProblem, inverting_forward_map
using ParameterEstimocean.PseudoSteppingSchemes: Kovachki2018InitialConvergenceRatio
using Printf
using JLD2
using LinearAlgebra

# using GLMakie

# Import "lesbrary_inverse_problem" + plotting utility "calibration_progress_figure"
include("multi_resolution_calibration_utilities.jl")

calibration_filenames = (
    full = "full_catke_calibration.jld2",
    simple = "simple_catke_calibration.jld2",
    simple_with_conv_adj = "simple_catke_conv_adj_calibration.jld2",
)

function load_summaries(filename)
    file = jldopen(filename)
    summaries = file["iteration_summaries"]
    close(file)
    return summaries
end

function get_best_parameters(summary)
    _, k = findmin(summary.mean_square_errors)
    return summary.parameters[k]
end

parameter_names = (
    :CᵂwΔ,  :Cᵂu★,
    :Cᴰ⁻, :Cᴰʳ, :CᴰRiᶜ, :CᴰRiʷ,
    :Cˢc,   :Cˢu,  :Cˢe,
    :Cᵇc,   :Cᵇu,  :Cᵇe,
    :Cᴷc⁻,  :Cᴷu⁻, :Cᴷe⁻,
    :Cᴷcʳ,  :Cᴷuʳ, :Cᴷeʳ,
    :CᴷRiᶜ, :CᴷRiʷ,
    :Cᵟc,   :Cᵟu,  :Cᵟe,
    #:Ku_adjustment, :Kc_adjustment, :Ke_adjustment,
)

neutral_parameter_values = (
    Cᴰʳ   = 0.0,
    CᴰRiᶜ = 0.0,
    CᴰRiʷ = 0.0,
    Cᴷcʳ  = 0.0,
    Cᴷuʳ  = 0.0,
    Cᴷeʳ  = 0.0,
    CᴷRiᶜ = 0.0,
    CᴷRiʷ = 0.0,
    Cᵟc   = 0.5,
    Cᵟu   = 0.5,
    Cᵟe   = 0.5,
    Ku_adjustment = 0.0,
    Kc_adjustment = 0.0,
    Ke_adjustment = 0.0)

closure = CATKEVerticalDiffusivity()
free_parameters = FreeParameters(prior_library, names=parameter_names)

# Two grids: "coarse" with ECCO vertical resolution to z=-256 m, and a fine grid with 4m resolution
Nz_ecco = length(ecco_vertical_grid) - 1
coarse_regrid = RectilinearGrid(size=Nz_ecco, z=ecco_vertical_grid, topology=(Flat, Flat, Bounded))
fine_regrid   = RectilinearGrid(size=48; z=(-256, 0), topology=(Flat, Flat, Bounded))

# Batch the inverse problems
times = [2hours, 24hours]
Nensemble = 3
architecture = CPU()
inverse_problem_kwargs = (; free_parameters, Nensemble, architecture, closure)
coarse_ip = lesbrary_inverse_problem(coarse_regrid; times, inverse_problem_kwargs...)
fine_ip   = lesbrary_inverse_problem(fine_regrid; times, inverse_problem_kwargs...)

y = observation_map(fine_ip)
summaries = load_summaries(calibration_filenames.full)
#summaries = load_summaries(calibration_filenames.simple)
mean_θ₀ = summaries[0].ensemble_mean
mean_θₙ = summaries[end].ensemble_mean
best_θₙ = get_best_parameters(summaries[end])

# Add neutral defaults (second in `merge` takes precedence)
# mean_θ₀ = merge(neutral_parameter_values, mean_θ₀)
# mean_θₙ = merge(neutral_parameter_values, mean_θₙ)
# best_θₙ = merge(neutral_parameter_values, best_θₙ)

G = forward_map(fine_ip, [mean_θ₀, mean_θₙ, best_θₙ])

Niterations = length(summaries)
Nensemble = length(summaries[0].parameters)
errors = zeros(Niterations, Nensemble)
det_Cᶿᶿ = zeros(Niterations)
pseudotimes = zeros(Niterations)

for (n, summary) in enumerate(summaries)
    errors[n, :] .= summary.mean_square_errors
    X = summary.unconstrained_parameters
    det_Cᶿᶿ[n] = det(cov(X, dims=2))
    pseudotimes[n] = summary.pseudotime
end

mean_iteration_ip = lesbrary_inverse_problem(fine_regrid; Nensemble=Niterations, times, free_parameters, architecture, closure)
θᵢ = [summaries[i-1].ensemble_mean for i = 1:Niterations]

@info "Evaluating mean parameters from every iteration..."
Gᵢ = @time forward_map(mean_iteration_ip, θᵢ)
mean_θ_errors = [norm(y .- Gᵢ[:, i]) for i in 1:Niterations]

#####
##### Figure
#####

fig = Figure(resolution=(1200, 600))
ax = Axis(fig[1:2, 1], ylabel="y", xlabel="Observation index")
lines!(ax, y[:], linewidth=6, color=(:indigo, 0.6))

lines!(ax, G[:, 1], linewidth=3, linestyle=:dash, color=(:red, 0.6), label="Initial mean")
lines!(ax, G[:, 2], linewidth=3, color=(:black, 0.6), label="Final mean")
lines!(ax, G[:, 3], linewidth=3, color=(:seagreen, 0.6), label="Final best")

axislegend(ax)

ax_e = Axis(fig[1, 2], xlabel="Pseudotime", ylabel="Mean square error", yscale=log10)
ylims!(ax_e, 1e-1, 1e3)

for k = 1:Nensemble
    #lines!(ax_e, pseudotimes, errors[:, k], color=(:black, 0.1))
    scatter!(ax_e, pseudotimes, errors[:, k], markersize=5, color=(:black, 0.1))
    scatter!(ax_e, pseudotimes, mean_θ_errors, linewidth=4, color=(:black, 0.6))
end

ax_C = Axis(fig[2, 2], xlabel="Pseudotime", ylabel="det(cov(X))", yscale=log10)
lines!(ax_C, pseudotimes, det_Cᶿᶿ)
scatter!(ax_C, pseudotimes, det_Cᶿᶿ)

display(fig)


