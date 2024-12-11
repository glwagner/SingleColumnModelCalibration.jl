using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using ParameterEstimocean
using JLD2
using Printf

using SingleColumnModelCalibration:
    dependent_parameter_sets,
    build_batched_inverse_problem,
    prior_library

dir = joinpath("..", "results")
closure = CATKEVerticalDiffusivity()
closure_label = "CATKE"
name = "extended_stability"
figdir = joinpath("figs", "CATKE_" * name)
#filename = "extended_stability_conv_adj_1_Nens1000_Δt60_τ10000_Nz32_Nz64_Nz128_12_hour_suite_24_hour_suite_48_hour_suite_with_tracer.jld2"
filename = "extended_stability_1_Nens1000_Δt30_τ10000_Nz32_Nz64_Nz128_12_hour_suite_24_hour_suite_48_hour_suite_no_convection.jld2"
filepath = joinpath(dir, filename)

try
    mkdir(figdir)
catch; end


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

Cˡᵒc = optimal_parameters.Cˡᵒc
Cˡᵒu = optimal_parameters.Cˡᵒu
Cˡᵒe = optimal_parameters.Cˡᵒe
CˡᵒD = optimal_parameters.CˡᵒD
Cʰⁱc = optimal_parameters.Cʰⁱc
Cʰⁱu = optimal_parameters.Cʰⁱu
Cʰⁱe = optimal_parameters.Cʰⁱe
# Cᶜc = optimal_parameters.Cᶜc
# Cᶜu = optimal_parameters.Cᶜu
# Cᶜe = optimal_parameters.Cᶜe
Cᵘⁿc = optimal_parameters.Cᵘⁿc
Cᵘⁿu = optimal_parameters.Cᵘⁿu
Cᵘⁿe = optimal_parameters.Cᵘⁿe
Cˢ = optimal_parameters.Cˢ

derived_parameters = (
    Pr₀ = Cˡᵒu / Cˡᵒc,
    Pr∞ = Cʰⁱu / Cʰⁱc,
    Pr⁻ = Cᵘⁿu / Cᵘⁿc,
    # Prᶜ = Cᶜu / Cᶜc,
    Sc₀ = Cˡᵒu / Cˡᵒe,
    Sc∞ = Cʰⁱu / Cʰⁱe,
    Sc⁻ = Cᵘⁿu / Cᵘⁿe,
    # Scᶜ = Cᶜu / Cᶜe,
    Ri★ = Cˡᵒu / (Cˡᵒc + CˡᵒD),
    ϰvk = Cˡᵒu * Cˢ,
)

all_parameters = merge(optimal_parameters, derived_parameters)

for name in keys(all_parameters)
    C = all_parameters[name]
    print(@sprintf("% 4s: %6.3f", name, C), '\n')
end

grid_parameters = [(size=128, z=(-256, 0))]

grid_colors = [
    (:royalblue1, 0.8),
    (:darkred, 0.8),
    (:black, 0.8),
]

suite_parameters = [
    (name = "6_hour_suite",  resolution="1m", stop_time=6hours),
    (name = "12_hour_suite", resolution="1m", stop_time=12hours),
    (name = "24_hour_suite", resolution="1m", stop_time=24hours),
    (name = "48_hour_suite", resolution="1m", stop_time=48hours),
    (name = "72_hour_suite", resolution="1m", stop_time=72hours),
]

batched_ip = build_batched_inverse_problem(closure, name;
                                           #Nensemble = 1,
                                           Δt = 1minutes,
                                           grid_parameters,
                                           suite_parameters)

forward_run!(batched_ip, [optimal_parameters])

#####
##### Figure
#####

LES_color = (:seagreen, 0.6)

k★      = 1
zlim    = -192
Ngrids  = length(grid_parameters)
Nsuites = length(suite_parameters)

suite_names = [suite_parameters[s].name for s = 1:length(suite_parameters)]

titles = [
    "Free convection",
    "Weak wind, strong cooling",
    "Mid wind, mid cooling",
    "Strong wind, weak cooling",
    "Strong wind",
    "Strong wind, no rotating",
    "Strong wind and sunny",
]

for (s, suite) in enumerate(suite_names)

    fig = Figure(size=(1800, 1200))

    i1 = Ngrids*(s-1) + 1
    ip1 = batched_ip[i1]
    times = observation_times(ip1.observations)
    Nt = length(times)

    ax_b = []
    ax_u = []
    ax_c = []

    LES_str = string("LES at t = ", prettytime(times[end]))

    grid1 = ip1.simulation.model.grid
    z1 = znodes(grid1, Center())
    Δz1 = Δzᶜᶜᶜ(1, 1, grid1.Nz, grid1)
    Δz1_str = @sprintf("%d", Δz1)
    
    for c = 1:6
        yaxisposition = c < 6 ? :left : :right
        Label(fig[1, c], titles[c], tellwidth=false)

        ax_bc = Axis(fig[2, c]; ylabel="z (m)", xlabel="Buoyancy \n (m s⁻²)", yaxisposition)
        push!(ax_b, ax_bc)

        xticks = [-0.2, 0.0, 0.2, 0.4]
        ax_uc = c == 1 ? nothing : Axis(fig[3, c], ylabel="z (m)", xlabel="Velocities \n (m s⁻¹)"; yaxisposition, xticks)
        push!(ax_u, ax_uc)

        ax_cc = Axis(fig[4, c]; ylabel="z (m)", xlabel="Passive \n tracer", yaxisposition)
        push!(ax_c, ax_cc)

        b_init = interior(ip1.observations[c].field_time_serieses.b[1], 1, 1, :)
        b_obs  = interior(ip1.observations[c].field_time_serieses.b[Nt], 1, 1, :)
        c_obs  = interior(ip1.observations[c].field_time_serieses.c[Nt], 1, 1, :)
        start_time = ip1.observations[c].times[1]

        lines!(ax_b[c], b_obs,  z1, linewidth=8, label=LES_str, color=LES_color)
        lines!(ax_c[c], c_obs,  z1, linewidth=8, label=LES_str, color=LES_color)
        
        xlims!(ax_b[c], -6e-4, -1e-4)
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

            C = interior(ip.time_series_collector.field_time_serieses.c[Nt], k★, c, :)
            lines!(ax_c[c], C, z, linewidth=3; label="$closure_label, Δz = $Δz_str m", color)

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

    # text!(ax_b[1], 0.03905, -185, text="(a)")
    # text!(ax_b[2], 0.03905, -185, text="(b)")
    # text!(ax_b[3], 0.03905, -185, text="(c)")
    # text!(ax_b[4], 0.03905, -185, text="(d)")
    # text!(ax_b[5], 0.03905, -185, text="(e)")
    # text!(ax_b[6], 0.03905, -185, text="(f)")

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

    figpath = joinpath(figdir, suite * ".png")
    save(figpath, fig)

    # save("$(name)_$(suite)_assessment_$(suffix)_conservative.pdf", fig)
end

