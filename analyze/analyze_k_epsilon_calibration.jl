using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity, TKEDissipationVerticalDiffusivity
using Oceananigans.Operators: Δzᶜᶜᶜ
using ParameterEstimocean
using JLD2
using Printf
using GLMakie

using SingleColumnModelCalibration:
    dependent_parameter_sets,
    build_batched_inverse_problem,
    prior_library

dir = joinpath("..", "results")
closure = TKEDissipationVerticalDiffusivity()
closure_label = "k-ϵ"
name = "variable_stabilities"
filename = "variable_stabilities_1_Nens2000_Δt60_τ10000_Nz32_Nz64_Nz128_12_hour_suite_24_hour_suite_48_hour_suite_k_epsilon.jld2"
filepath = joinpath(dir, filename)

file = jldopen(filepath)
summaries = file["iteration_summaries"]
close(file)

# Reshape data
Niters = length(summaries)
Nens = length(summaries[1].parameters)

parameters = [s.parameters[ω] for s in summaries, ω = 1:Nens]
objectives = [s.objective_values[ω][1] for s in summaries, ω = 1:Nens]

best_objective, i★ = findmin(objectives)
optimal_parameters = parameters[i★]

default_parameters = (
     Cσe = closure.stability_functions.Cσe, 
     Cσϵ = closure.stability_functions.Cσϵ,
     Cu₀ = closure.stability_functions.Cu₀,
     Cu₁ = closure.stability_functions.Cu₁,
     Cu₂ = closure.stability_functions.Cu₂,
     Cc₀ = closure.stability_functions.Cc₀,
     Cc₁ = closure.stability_functions.Cc₁,
     Cc₂ = closure.stability_functions.Cc₂,
     Cd₁ = closure.stability_functions.Cd₁,
     Cd₂ = closure.stability_functions.Cd₂,
     Cd₃ = closure.stability_functions.Cd₃,
     Cd₄ = closure.stability_functions.Cd₄,
     Cd₅ = closure.stability_functions.Cd₅,
     Cᵋϵ = closure.tke_dissipation_equations.Cᵋϵ,
     Cᴾϵ = closure.tke_dissipation_equations.Cᴾϵ,
     Cᵇϵ = closure.tke_dissipation_equations.Cᵇϵ,
    Cᵂu★ = closure.tke_dissipation_equations.Cᵂu★,
    CᵂwΔ = closure.tke_dissipation_equations.CᵂwΔ,
)

for name in keys(optimal_parameters)
    C★ = optimal_parameters[name]
    C₀ = default_parameters[name]
    print(@sprintf("% 4s: %6.3f (%.6f)", name, C★, C₀), '\n')
end

grid_parameters = [(size=128, z=(-256, 0))]

suite_parameters = [
    (name = "6_hour_suite",  resolution="1m", stop_time=6hours),
    (name = "12_hour_suite", resolution="1m", stop_time=12hours),
    (name = "24_hour_suite", resolution="1m", stop_time=24hours),
    (name = "48_hour_suite", resolution="1m", stop_time=48hours),
    (name = "72_hour_suite", resolution="1m", stop_time=72hours),
]

batched_ip = build_batched_inverse_problem(closure, name;
                                           Nensemble = 2,
                                           Δt = 10.0,
                                           grid_parameters,
                                           suite_parameters)

forward_run!(batched_ip, [optimal_parameters, NamedTuple()])

#####
##### Figure
#####

LES_color = (:seagreen, 0.6)
grid_colors = (:black, 0.8)

zlim    = -192
Ngrids  = length(grid_parameters)
Nsuites = length(suite_parameters)

suite_names = [suite_parameters[s].name for s = 1:length(suite_parameters)]

for (s, suite) in enumerate(suite_names)

    fig = Figure(size=(1800, 1200))

    i1 = Ngrids*(s-1) + 1
    ip1 = batched_ip[i1]
    grid1 = ip1.simulation.model.grid
    z1 = znodes(grid1, Center())
    times = observation_times(ip1.observations)
    Nt = length(times)

    ax_b = []
    ax_u = []
    ax_c = []

    LES_str = string("LES at t = ", prettytime(times[end]))
    iᵇ = Ngrids*(s-1)

    for c = 1:7
        yaxisposition = c < 7 ? :left : :right
        # Label(fig[1, c], titles[c], tellwidth=false)

        ax_bc = Axis(fig[2, c]; ylabel="z (m)", xlabel="Buoyancy \n (m s⁻²)", yaxisposition)
        push!(ax_b, ax_bc)

        ax_uc = c == 1 ? nothing : Axis(fig[3, c], ylabel="z (m)", xlabel="Velocities \n (m s⁻¹)"; yaxisposition)
        push!(ax_u, ax_uc)

        ax_cc = Axis(fig[4, c]; ylabel="z (m)", xlabel="Passive \n tracer", yaxisposition)
        push!(ax_c, ax_cc)

        b_obs = interior(ip1.observations[c].field_time_serieses.b[Nt], 1, 1, :)
        c_obs = interior(ip1.observations[c].field_time_serieses.c[Nt], 1, 1, :)
        start_time = ip1.observations[c].times[1]

        lines!(ax_b[c], b_obs, z1, linewidth=8, label=LES_str, color=LES_color)
        lines!(ax_c[c], c_obs, z1, linewidth=8, label=LES_str, color=LES_color)

        if c > 1
            u_obs = interior(ip1.observations[c].field_time_serieses.u[Nt], 1, 1, :)
            lines!(ax_u[c], u_obs, z1, linewidth=8, label="u, LES",  color=LES_color)

            if :v ∈ keys(ip1.observations[c].field_time_serieses)
                v_obs = interior(ip1.observations[c].field_time_serieses.v[Nt], 1, 1, :)
                lines!(ax_u[c], v_obs, z1, color=LES_color, linewidth=4, linestyle=:dash)
            end

            ylims!(ax_u[c], zlim, 0)
        end

        for iᵍ = 1:Ngrids
            ip = batched_ip[iᵇ + iᵍ]
            grid = ip.simulation.model.grid
            z = znodes(grid, Center())
            Δz = Δzᶜᶜᶜ(1, 1, grid.Nz, grid)
            Δz_str = @sprintf("%d", Δz)
            color = grid_colors[iᵍ]

            b = interior(ip.time_series_collector.field_time_serieses.b[Nt], 1, c, :)
            lines!(ax_b[c], b, z, linewidth=3; label="$closure_label, Δz = $Δz_str m", color)

            C = interior(ip.time_series_collector.field_time_serieses.c[Nt], 1, c, :)
            lines!(ax_c[c], C, z, linewidth=3; label="$closure_label, Δz = $Δz_str m", color)

            if c > 1
                u = interior(ip.time_series_collector.field_time_serieses.u[Nt], 1, c, :)
                lines!(ax_u[c], u, z, linewidth=3; label="u, $closure_label, Δz = $Δz_str m", color)

                if :v ∈ keys(ip1.observations[c].field_time_serieses)
                    v = interior(ip.time_series_collector.field_time_serieses.v[Nt], 1, c, :)
                    lines!(ax_u[c], v, z; color, linewidth=2, linestyle=:dash, label="v")
                end
            end
        end

        ylims!(ax_b[c], zlim, 0)
        hidespines!(ax_b[c], :t)
    end

    Legend(fig[3, 1], ax_b[2])

    display(fig)
end


