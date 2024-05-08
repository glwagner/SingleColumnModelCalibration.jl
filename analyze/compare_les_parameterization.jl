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

using SingleColumnModelCalibration:
    dependent_parameter_sets,
    build_batched_inverse_problem,
    prior_library

set_theme!(Theme(fontsize=19))

# dir = "." #calibration_results/round1"
# name = "ri_based"
# closure = RiBasedVerticalDiffusivity()
# closure_label = "RiBased"

dir = "../parameters"
#name = "constant_Pr_no_shear"
#name = "variable_Pr"
#name = "constant_Pr_conv_adj"
name = "variable_Pr_conv_adj"
#name = "fixed_Ric"

# suffix = "Nens2000_Δt20_τ10000_Nz32_Nz64_Nz128_12_hour_suite_24_hour_suite_48_hour_suite_convective_momentum_mixing_min_tke.jld2"
closure = CATKEVerticalDiffusivity(turbulent_kinetic_energy_time_step=5minute)
#closure = CATKEVerticalDiffusivity(turbulent_kinetic_energy_time_step=nothing)
closure_label = "CATKE"

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
    #(size=256, z=(-256, 0)),
    (size=128, z=(-256, 0)),
    (size=64,  z=(-256, 0)),
    #(size=32,  z=(-256, 0)),
    (size=16,  z=(-256, 0)),
]

grid_colors = [
    (:royalblue1, 0.8),
    (:darkred, 0.8),
    (:black, 0.8),
]

suite_parameters = [
    (name = "6_hour_suite",  resolution="0.75m", stop_time=6hours),
    (name = "12_hour_suite", resolution="1m", stop_time=12hours),
#    (name = "18_hour_suite", resolution="1m", stop_time=18hours),
    (name = "24_hour_suite", resolution="1m", stop_time=24hours),
#    (name = "36_hour_suite", resolution="1m", stop_time=36hours),
    (name = "48_hour_suite", resolution="1m", stop_time=48hours),
    (name = "72_hour_suite", resolution="1m", stop_time=72hours),
]

batched_ip = build_batched_inverse_problem(closure, name;
                                           Nensemble = 1,
                                           Δt = 20minutes,
                                           grid_parameters,
                                           suite_parameters)

forward_run!(batched_ip, [optimal_parameters])

optimal_catke = batched_ip[1].simulation.model.closure[1]
@save "optimal_catke.jld2" optimal_catke

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

    fig = Figure(size=(1200, 600))

    i1 = Ngrids*(s-1) + 1
    ip1 = batched_ip[i1]
    times = observation_times(ip1.observations)
    Nt = length(times)

    ax_b = []
    ax_u = []

    LES_str = string("LES at t = ", prettytime(times[end]))

    grid1 = ip1.simulation.model.grid
    z1 = znodes(grid1, Center())
    Δz1 = Δzᶜᶜᶜ(1, 1, grid1.Nz, grid1)
    Δz1_str = @sprintf("%d", Δz1)
    
    for c = 1:6
        yaxisposition = c < 6 ? :left : :right
        Label(fig[1, c], titles[c], tellwidth=false)

        if c == 1
            ax_bc = Axis(fig[2, c]; ylabel="z (m)", xlabel="Buoyancy \n (m s⁻²)", yaxisposition, xticks=[0.0386, 0.0389, 0.0392])
        elseif c < 5
            ax_bc = Axis(fig[2, c]; ylabel="z (m)", xlabel="Buoyancy \n (m s⁻²)", yaxisposition, xticks=[0.0386, 0.0389])
        else
            ax_bc = Axis(fig[2, c]; ylabel="z (m)", xlabel="Buoyancy \n (m s⁻²)", yaxisposition, xticks=[0.0386, 0.0390])
        end
        push!(ax_b, ax_bc)

        xticks = [-0.2, 0.0, 0.2, 0.4]
        ax_uc = c == 1 ? nothing : Axis(fig[3, c], ylabel="z (m)", xlabel="Velocities \n (m s⁻¹)"; yaxisposition, xticks)
        push!(ax_u, ax_uc)

        b_init = interior(ip1.observations[c].field_time_serieses.b[1], 1, 1, :)
        b_obs  = interior(ip1.observations[c].field_time_serieses.b[Nt], 1, 1, :)
        start_time = ip1.observations[c].times[1]

        lines!(ax_b[c], b_init, z1, linewidth=2, label="Initial condition", color=LES_color, linestyle=:dot)
        lines!(ax_b[c], b_obs,  z1, linewidth=8, label=LES_str, color=LES_color)
        
        if c == 6
            if s == 1
                !isnothing(ax_u[c]) && xlims!(ax_u[c], -0.1, 0.6)
            else
                !isnothing(ax_u[c]) && xlims!(ax_u[c], -0.1, 0.4)
            end
        else
            if s == 1
                !isnothing(ax_u[c]) && xlims!(ax_u[c], -0.25, 0.4)
            else
                !isnothing(ax_u[c]) && xlims!(ax_u[c], -0.07, 0.2)
            end
        end

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
            z = znodes(grid, Center())
            Δz = Δzᶜᶜᶜ(1, 1, grid.Nz, grid)
            Δz_str = @sprintf("%d", Δz)
            color = grid_colors[iᵍ]

            b = interior(ip.time_series_collector.field_time_serieses.b[Nt], k★, c, :)
            lines!(ax_b[c], b, z, linewidth=3; label="$closure_label, Δz = $Δz_str m", color)

            if c > 1
                u = interior(ip.time_series_collector.field_time_serieses.u[Nt], k★, c, :)
                lines!(ax_u[c], u, z, linewidth=3; label="u, $closure_label, Δz = $Δz_str m", color)

                if :v ∈ keys(ip1.observations[c].field_time_serieses)
                    v = interior(ip.time_series_collector.field_time_serieses.v[Nt], k★, c, :)
                    lines!(ax_u[c], v, z; color, linewidth=2, linestyle=:dash, label="v")
                end
            end
        end
    end

    Legend(fig[3, 1], ax_b[2])
    if s == 1
        text!(ax_u[2], +0.18, -50.0, text="u")
        text!(ax_u[2], -0.2, -130.0, text="v")
    else
        text!(ax_u[2], +0.09, -35.0, text="u")
        text!(ax_u[2], -0.06, -30.0, text="v")
    end

    text!(ax_b[1], 0.03905, -185, text="(a)")
    text!(ax_b[2], 0.03905, -185, text="(b)")
    text!(ax_b[3], 0.03905, -185, text="(c)")
    text!(ax_b[4], 0.03905, -185, text="(d)")
    text!(ax_b[5], 0.03905, -185, text="(e)")
    text!(ax_b[6], 0.03905, -185, text="(f)")

    if suite == "6_hour_suite"
        text!(ax_u[2], 0.27, -185, text="(g)")
        text!(ax_u[3], 0.27, -185, text="(h)")
        text!(ax_u[4], 0.27, -185, text="(i)")
        text!(ax_u[5], 0.27, -185, text="(j)")
        text!(ax_u[6], 0.37, -185, text="(k)")
    else
        text!(ax_u[2], 0.14, -185, text="(g)")
        text!(ax_u[3], 0.14, -185, text="(h)")
        text!(ax_u[4], 0.14, -185, text="(i)")
        text!(ax_u[5], 0.14, -185, text="(j)")
        text!(ax_u[6], 0.27, -185, text="(k)")
    end

    display(fig)

    save("$(name)_$(suite)_assessment_$(suffix)_conservative.pdf", fig)
end

