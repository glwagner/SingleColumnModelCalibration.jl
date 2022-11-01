#=
using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using ParameterEstimocean.InverseProblems: BatchedInverseProblem, inverting_forward_map
using ParameterEstimocean.PseudoSteppingSchemes: Kovachki2018InitialConvergenceRatio
using Printf
using JLD2
using LinearAlgebra
using GLMakie

set_theme!(Theme(fontsize=16))

# Import "lesbrary_inverse_problem" + plotting utility "calibration_progress_figure"
include("multi_resolution_calibration_utilities.jl")

calibration_filenames = (
    complex                  = "complex_catke_calibration.jld2",
    simple                   = "simple_catke_calibration.jld2",
    simple_with_conv_adj     = "simple_catke_conv_adj_calibration.jld2",
    goldilocks               = "goldilocks_catke_calibration.jld2",
    goldilocks_with_conv_adj = "goldilocks_catke_conv_adj_calibration.jld2",
    #nemo                     = "nemo_catke_calibration.jld2",
    #nemo_with_conv_adj       = "nemo_catke_conv_adj_calibration.jld2",
    #sanity_check             = "sanity_check_calibration.jld2",
)

function load_summaries(filename)
    dir = "calibration_results"
    file = jldopen(joinpath(dir, filename))
    summaries = file["iteration_summaries"]
    close(file)
    return summaries
end

function get_best_parameters(summary)
    _, k = findmin(summary.mean_square_errors)
    return summary.parameters[k]
end

neutral_default_mixing_length_parameters = Dict(
    :Cᵟc => 0.5,
    :Cᵟu => 0.5,
    :Cᵟe => 0.5,
    :Cᵇu => Inf,
    :Cᵇc => Inf,
    :Cᵇe => Inf,
    :Cˢu => Inf,
    :Cˢc => Inf,
    :Cˢe => Inf,
    :Cᴷcʳ => 0.0,
    :Cᴷuʳ => 0.0,
    :Cᴷeʳ => 0.0,
)

neutral_default_tke_parameters = Dict(
    :CᵂwΔ => 0.0,
    :Cᵂu★ => 0.0,
    :Cᴰ⁻ => 0.0,
    :Cᴰʳ => 0.0,
    :CᴰRiᶜ => 0.0,
    :CᴰRiʷ => 0.0,
)

other_default_parameters = Dict(
    :Ku_adjustment => 0.0,
    :Kc_adjustment => 0.0,
    :Ke_adjustment => 0.0,
)

mixing_length = MixingLength(; neutral_default_mixing_length_parameters...)
turbulent_kinetic_energy_equation = TurbulentKineticEnergyEquation(; neutral_default_tke_parameters...)
closure = CATKEVerticalDiffusivity(; other_default_parameters...,
                                   mixing_length, turbulent_kinetic_energy_equation)

name = :complex
filename = calibration_filenames[name]
summaries = load_summaries(filename)
mean_θₙ = summaries[end].ensemble_mean
parameter_names = keys(mean_θₙ)
free_parameters = FreeParameters(prior_library, names=parameter_names)

θ_dict = Dict()

for name in keys(calibration_filenames)
    filename = calibration_filenames[name]
    last_summary = load_summaries(filename)[end]
    θ = last_summary.ensemble_mean
    #θ = get_best_parameters(last_summary)
    θ_dict[name] = θ
end

θ_vec = [θ_dict[name] for name in keys(calibration_filenames)]

# Two grids: "coarse" with ECCO vertical resolution to z=-256 m, and a fine grid with 4m resolution
Nz_ecco = length(ecco_vertical_grid) - 1
coarse_regrid = RectilinearGrid(size=Nz_ecco, z=ecco_vertical_grid, topology=(Flat, Flat, Bounded))
fine_regrid   = RectilinearGrid(size=48; z=(-256, 0), topology=(Flat, Flat, Bounded))

# Batch the inverse problems
times = collect(range(2hours, stop=24hours, length=4))
Nensemble = length(θ_vec)
architecture = CPU()
inverse_problem_kwargs = (; free_parameters, Nensemble, architecture, closure, Δt=1minute)
coarse_ip = lesbrary_inverse_problem(coarse_regrid; times, inverse_problem_kwargs...)
fine_ip   = lesbrary_inverse_problem(fine_regrid; times, inverse_problem_kwargs...)

G_fine   = forward_map(fine_ip, θ_vec)
G_coarse = forward_map(coarse_ip, θ_vec)

y_fine   = observation_map(fine_ip)
y_coarse = observation_map(coarse_ip)
=#

# Number of data scenarios
Ncases = length(coarse_ip.observations.observations) # Number of observational cases

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
            #e_cases[s][case] = e_case[s] / (Nz * (Nfields - 1))
            e_cases[s][case] = e_case[s]
        end
        
        iy₁ += Ncase
    end

    # Normalize by "amount" of data?
    # Or, we should actually do this normalization during training...
    #global_errors = [norm(G[:, i] .- y) for i = 1:Nsets] ./ (Nunits * Nz)
    e_global = [norm(G[:, i] .- y) for i = 1:Nsets]

    return e_global, e_cases
end

coarse_global_errors, coarse_case_errors = analyze_errors(coarse_ip, y_coarse, G_coarse)
fine_global_errors, fine_case_errors = analyze_errors(fine_ip, y_fine, G_fine)

case_errors = [zeros(Ncases, 2) for _ in 1:Nsets]

for s = 1:Nsets
    case_errors[s][:, 1] .= coarse_case_errors[s]
    case_errors[s][:, 2] .= fine_case_errors[s]
end

#####
##### Figure
#####

case_names = [replace(name, "_" => " ") for name in cases]
elims = (0.9, 1.1) .* mean(coarse_global_errors)

symbol_names = keys(calibration_filenames)
string_names = [replace(string(name), "_" => " ") for name in symbol_names]

cmin = minimum(minimum(e) for e in case_errors)
cmax = maximum(maximum(e) for e in case_errors)

@show cmin cmax

fig = Figure(resolution=(1200, 600))

for s = 1:Nsets
    yaxisposition = s < Nsets ? :left : :right
    ylabel = "Scaled error"
    ax = Axis(fig[1, s]; xticks=(1:1, string_names[s:s]), yaxisposition, ylabel)

    #ylims!(ax, elims)
    barplot!(ax, 0.5, coarse_global_errors[s])
    barplot!(ax, 1.5, fine_global_errors[s])
    s > 1 && s < Nsets && hideydecorations!(ax)
    s > 1 && hidespines!(ax, :l)
    s < Nsets && hidespines!(ax, :r)
    hidespines!(ax, :t)

    hidedecorations!(ax, label=false, ticklabels=false, ticks=false, minorticks=false,
                     minorgrid=true, grid=true)

    axc = Axis(fig[2, s], yticks=(1:Ncases, case_names))
    heatmap!(axc, case_errors[s]', colormap=:reds) #, clims=(0, 1))

    s > 1 && hidedecorations!(axc)
    s == 1 && hidexdecorations!(axc)
end


display(fig)

#=
z_fine = znodes(Center, fine_regrid)
z_coarse = znodes(Center, coarse_regrid)
times = observation_times(fine_ip.observations)
Nt = length(times)

ax_b = []
ax_u = []

titles = [
    "Free convection",
    "Weak wind \n strong cooling",
    "Medium wind \n medium cooling",
    "Strong wind \n weak cooling",
    "Strong wind \n no cooling",
    "Strong wind \n no rotation",
]

sim_color    = (:seagreen, 0.6)
fine_color   = (:darkred, 0.8)
coarse_color = (:royalblue1, 0.8)
k_plot       = 2

for c = 1:6
    if c < 5
        yaxisposition=:left
    else
        yaxisposition=:right
    end

    Label(fig[1, c], titles[c], tellwidth=false)

    ax_bc = Axis(fig[2, c]; ylabel="z (m)", xlabel="Buoyancy (m s⁻²)", yaxisposition, xticks=[0.0386, 0.039, 0.0394])
    push!(ax_b, ax_bc)

    ax_uc = c == 1 ? nothing : Axis(fig[3, c], ylabel="z (m)", xlabel="Velocities (m s⁻¹)"; yaxisposition, xticks=-0.1:0.1:0.3)
    push!(ax_u, ax_uc)

    b_init   = interior(fine_ip.observations[c].field_time_serieses.b[1], 1, 1, :)
    b_obs    = interior(fine_ip.observations[c].field_time_serieses.b[Nt], 1, 1, :)
    b_fine   = interior(  fine_ip.time_series_collector.field_time_serieses.b[Nt], k_plot, c, :)
    b_coarse = interior(coarse_ip.time_series_collector.field_time_serieses.b[Nt], k_plot, c, :)

    lines!(ax_b[c], b_init,   z_fine,   linewidth=2, label="Initial condition at t = 2 hours", color=sim_color, linestyle=:dot)
    lines!(ax_b[c], b_obs,    z_fine,   linewidth=8, label="LES at t = 1 day", color=sim_color)
    lines!(ax_b[c], b_fine,   z_fine,   linewidth=3, label="CATKE, Δz ≈ 10 m", color=fine_color)
    lines!(ax_b[c], b_coarse, z_coarse, linewidth=3, label="CATKE, Δz = 4 m", color=coarse_color)

    xlims!(ax_b[c], 0.0386, maximum(b_init) + 2e-5)
    ylims!(ax_b[c], -192, 0)

    if c > 1
        #xlims!(ax_u[c], -0.15, 0.35)
        ylims!(ax_u[c], -192, 0)

        u_obs    = interior(fine_ip.observations[c].field_time_serieses.u[Nt], 1, 1, :)
        u_fine   = interior(  fine_ip.time_series_collector.field_time_serieses.u[Nt], k_plot, c, :)
        u_coarse = interior(coarse_ip.time_series_collector.field_time_serieses.u[Nt], k_plot, c, :)

        lines!(ax_u[c], u_obs,    z_fine,   linewidth=8, label="u, LES",              color=sim_color)
        lines!(ax_u[c], u_fine,   z_fine,   linewidth=3, label="u, CATKE, Δz ≈ 10 m", color=fine_color)
        lines!(ax_u[c], u_coarse, z_coarse, linewidth=3, label="u, CATKE, Δz = 4 m",  color=coarse_color)

        if :v ∈ keys(fine_ip.observations[c].field_time_serieses)
            v_obs    = interior(fine_ip.observations[c].field_time_serieses.v[Nt], 1, 1, :)
            v_fine   = interior(  fine_ip.time_series_collector.field_time_serieses.v[Nt], k_plot, c, :)
            v_coarse = interior(coarse_ip.time_series_collector.field_time_serieses.v[Nt], k_plot, c, :)

            lines!(ax_u[c], v_obs,    z_fine,   color=sim_color,    linewidth=4, linestyle=:dash) #, label="v, LES")
            lines!(ax_u[c], v_fine,   z_fine,   color=fine_color,   linewidth=2, linestyle=:dash, label="v") #, fine resolution CATKE")
            lines!(ax_u[c], v_coarse, z_coarse, color=coarse_color, linewidth=2, linestyle=:dash) #, label="v, coarse resolution CATKE")
        end
    end

    hidespines!(ax_b[c], :t)
    c != 1 && hidespines!(ax_u[c], :t)

    c != 1 && hidespines!(ax_b[c], :l)
    c != 1 && c != 6 && hideydecorations!(ax_b[c], grid=false)
    c != 6 && hidespines!(ax_b[c], :r)

    c > 2 && hidespines!(ax_u[c], :l)
    c != 1 && c != 6 && hidespines!(ax_u[c], :r)
    c > 2 && c != 6 && hideydecorations!(ax_u[c], grid=false)
end

Legend(fig[3, 1], ax_b[2])
text!(ax_u[2], +0.1, -50.0, text="u")
text!(ax_u[2], -0.14, -110.0, text="v")
#Legend(fig[5, 1], ax_u[2])

display(fig)

save("$(name)_catke_LES_comparison.png", fig)
=#
