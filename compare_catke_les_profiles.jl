using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using ParameterEstimocean.InverseProblems: BatchedInverseProblem, inverting_forward_map
using ParameterEstimocean.PseudoSteppingSchemes: Kovachki2018InitialConvergenceRatio
using Printf
using JLD2
using LinearAlgebra

# using GLMakie

set_theme!(Theme(fontsize=16))

# Import "lesbrary_inverse_problem" + plotting utility "calibration_progress_figure"
include("multi_resolution_calibration_utilities.jl")

calibration_filenames = (
    complex = "full_catke_calibration.jld2",
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

name = :complex
filename = calibration_filenames[name]
summaries = load_summaries(filename)

mean_θ₀ = summaries[0].ensemble_mean
mean_θₙ = summaries[end].ensemble_mean
best_θₙ = get_best_parameters(summaries[end])

closure = CATKEVerticalDiffusivity()

parameter_names = keys(best_θₙ)
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


forward_run!(fine_ip, [mean_θ₀, mean_θₙ, best_θₙ])
forward_run!(coarse_ip, [mean_θ₀, mean_θₙ, best_θₙ])

#####
##### Figure
#####

fig = Figure(resolution=(1200, 600))

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
