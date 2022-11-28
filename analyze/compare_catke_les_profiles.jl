using Oceananigans
using Oceananigans.Operators: Δzᶜᶜᶜ
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity

using ParameterEstimocean
using ParameterEstimocean.InverseProblems: BatchedInverseProblem, inverting_forward_map
using ParameterEstimocean.PseudoSteppingSchemes: Kovachki2018InitialConvergenceRatio
using ParameterEstimocean.Parameters: build_parameters_named_tuple

using Printf
using JLD2
using LinearAlgebra

using CairoMakie
using ElectronDisplay

using SingleColumnModelCalibration:
    dependent_parameter_sets,
    build_batched_inverse_problem,
    prior_library

set_theme!(Theme(fontsize=16))

# dir = "." #calibration_results/round1"
# name = :ri_based
# closure = RiBasedVerticalDiffusivity()

dir = "../parameters"
name = "constant_Pr"
#name = "variable_Pr"
#name = "constant_Pr_conv_adj"
#name = "variable_Pr_conv_adj"
closure = CATKEVerticalDiffusivity()

filepath = joinpath(dir, string(name) * "_best_parameters.jld2")
file = jldopen(filepath)
optimal_parameters = file["optimal_parameters"]
close(file)

dependent_parameters = dependent_parameter_sets[string(name)]
parameter_names = keys(optimal_parameters)
free_parameters = FreeParameters(prior_library; names=parameter_names, dependent_parameters)
optimal_parameters = build_parameters_named_tuple(free_parameters, optimal_parameters)
@show optimal_parameters

# Batch the inverse problems
grid_parameters = [
    (size=128, z=(-256, 0)),
    (size=64,  z=(-256, 0)),
    (size=32,  z=(-256, 0)),
]

grid_colors = [
    (:black, 0.8),
    (:darkred, 0.8),
    (:royalblue1, 0.8)
]

suite_parameters = [
    (name = "6_hour_suite",  resolution="0.75m", stop_time=6hours),
    (name = "12_hour_suite", resolution="1m", stop_time=12hours),
    (name = "18_hour_suite", resolution="1m", stop_time=18hours),
    (name = "24_hour_suite", resolution="1m", stop_time=24hours),
    (name = "36_hour_suite", resolution="1m", stop_time=36hours),
    (name = "48_hour_suite", resolution="1m", stop_time=48hours),
    (name = "72_hour_suite", resolution="1m", stop_time=72hours),
]

batched_ip = build_batched_inverse_problem(closure, name;
                                           Nensemble = 1,
                                           Δt = 1minute,
                                           grid_parameters,
                                           suite_parameters)

forward_run!(batched_ip, [optimal_parameters])

titles = [
    "Free convection",
    "Weak wind \n strong cooling",
    "Medium wind \n medium cooling",
    "Strong wind \n weak cooling",
    "Strong wind \n no cooling",
    "Strong wind \n no rotation",
]

#####
##### Figure
#####

LES_color = (:seagreen, 0.6)

k★      = 1
zlim    = -192
Ngrids  = length(grid_parameters)
Nsuites = length(suite_parameters)

suite_names = [suite_parameters[s].name for s = 1:length(suite_parameters)]

for (s, suite) in enumerate(suite_names)

    fig = Figure(resolution=(1200, 600))

    i1 = Ngrids*(s-1) + 1
    ip1 = batched_ip[i1]
    times = observation_times(ip1.observations)
    Nt = length(times)

    ax_b = []
    ax_u = []

    LES_str = string("LES at t = ", prettytime(times[end]))

    grid1 = ip1.simulation.model.grid
    z1 = znodes(Center, grid1)
    Δz1 = Δzᶜᶜᶜ(1, 1, grid1.Nz, grid1)
    Δz1_str = @sprintf("%d", Δz1)
    
    for c = 1:6
        yaxisposition = c < 6 ? :left : :right
        Label(fig[1, c], titles[c], tellwidth=false)

        ax_bc = Axis(fig[2, c]; ylabel="z (m)", xlabel="Buoyancy (m s⁻²)", yaxisposition, xticks=[0.0386, 0.039, 0.0394])
        push!(ax_b, ax_bc)

        ax_uc = c == 1 ? nothing : Axis(fig[3, c], ylabel="z (m)", xlabel="Velocities (m s⁻¹)"; yaxisposition, xticks=-0.1:0.1:0.3)
        push!(ax_u, ax_uc)

        b_init = interior(ip1.observations[c].field_time_serieses.b[1], 1, 1, :)
        b_obs  = interior(ip1.observations[c].field_time_serieses.b[Nt], 1, 1, :)
        start_time = ip1.observations[c].times[1]

        lines!(ax_b[c], b_init, z1, linewidth=2, label="Initial condition at t = " * prettytime(start_time), color=LES_color, linestyle=:dot)
        lines!(ax_b[c], b_obs,  z1, linewidth=8, label=LES_str, color=LES_color)
        
        xlims!(ax_b[c], 0.0386, maximum(b_init) + 2e-5)
        ylims!(ax_b[c], zlim, 0)

        if c > 1
            ylims!(ax_u[c], zlim, 0)
            u_obs = interior(ip1.observations[c].field_time_serieses.u[Nt], 1, 1, :)
            lines!(ax_u[c], u_obs, z1, linewidth=8, label="u, LES",  color=LES_color)

            if :v ∈ keys(ip1.observations[c].field_time_serieses)
                v_obs = interior(ip1.observations[c].field_time_serieses.v[Nt], 1, 1, :)
                lines!(ax_u[c], v_obs, z1, color=LES_color, linewidth=4, linestyle=:dash) #, label="v, LES")
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

    iᵇ = Ngrids*(s-1)

    for c = 1:6
        for iᵍ = 1:Ngrids
            ip = batched_ip[iᵇ + iᵍ]
            grid = ip.simulation.model.grid
            z = znodes(Center, grid)
            Δz = Δzᶜᶜᶜ(1, 1, grid.Nz, grid)
            Δz_str = @sprintf("%d", Δz)
            color = grid_colors[iᵍ]

            b = interior(ip.time_series_collector.field_time_serieses.b[Nt], k★, c, :)
            lines!(ax_b[c], b, z, linewidth=3; label="CATKE, Δz = $Δz_str m", color)

            if c > 1
                u = interior(ip.time_series_collector.field_time_serieses.u[Nt], k★, c, :)
                lines!(ax_u[c], u, z, linewidth=3; label="u, CATKE, Δz = $Δz_str m", color)

                if :v ∈ keys(ip1.observations[c].field_time_serieses)
                    v = interior(ip.time_series_collector.field_time_serieses.v[Nt], k★, c, :)
                    lines!(ax_u[c], v, z; color, linewidth=2, linestyle=:dash, label="v")
                end
            end
        end
    end

    Legend(fig[3, 1], ax_b[2])
    text!(ax_u[2], +0.1, -50.0, text="u")
    text!(ax_u[2], -0.14, -110.0, text="v")

    display(fig)

    #save("$(name)_catke_LES_comparison.png", fig)
end
