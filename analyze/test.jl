using Oceananigans.Operators: Δzᶜᶜᶜ

# Fine best *global* parameters
function get_best_parameters(summaries)
    best_parameters = summaries[0].mean_square_errors[1]
    best_mse = Inf

    for summary in summaries
        mse, k = findmin(summary.mean_square_errors)
        parameters = summary.parameters[k]
        best_parameters = ifelse(mse < best_mse, parameters, best_parameters)
    end

    return best_parameters
end

summaries = coarse_eki.iteration_summaries

mean_θ₀ = summaries[0].ensemble_mean
mean_θₙ = summaries[end].ensemble_mean
best_θₙ = get_best_parameters(summaries)

coarse_ip = coarse_eki.inverse_problem
grid = coarse_ip.simulation.model.grid

forward_run!(coarse_ip, [mean_θ₀, mean_θₙ, best_θₙ])

@show coarse_ip.simulation.model.closure[3]

#####
##### Figure
#####

fig = Figure(resolution=(1200, 600))

z_coarse = znodes(Center, grid)
times = observation_times(coarse_ip.observations)
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
coarse_color = (:royalblue1, 0.8)
k_plot       = 3

Δz_coarse = Δzᶜᶜᶜ(1, 1, grid.Nz, grid)
Δz_coarse_str = @sprintf("%.1f", Δz_coarse)
Δz_coarse_str = @sprintf("%.1f", Δz_coarse)

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

    b_init   = interior(coarse_ip.observations[c].field_time_serieses.b[1], 1, 1, :)
    b_obs    = interior(coarse_ip.observations[c].field_time_serieses.b[Nt], 1, 1, :)
    b_coarse = interior(coarse_ip.time_series_collector.field_time_serieses.b[Nt], k_plot, c, :)

    lines!(ax_b[c], b_init,   z_coarse,   linewidth=2, label="Initial condition at t = 2 hours", color=sim_color, linestyle=:dot)
    lines!(ax_b[c], b_obs,    z_coarse,   linewidth=8, label="LES at t = 1 day", color=sim_color)
    lines!(ax_b[c], b_coarse, z_coarse, linewidth=3, label="CATKE, Δz = $Δz_coarse_str m", color=coarse_color)

    xlims!(ax_b[c], 0.0386, maximum(b_init) + 2e-5)
    ylims!(ax_b[c], -192, 0)

    if c > 1
        #xlims!(ax_u[c], -0.15, 0.35)
        ylims!(ax_u[c], -192, 0)

        u_obs    = interior(coarse_ip.observations[c].field_time_serieses.u[Nt], 1, 1, :)
        u_coarse = interior(coarse_ip.time_series_collector.field_time_serieses.u[Nt], k_plot, c, :)

        lines!(ax_u[c], u_obs,    z_coarse,   linewidth=8, label="u, LES",              color=sim_color)
        lines!(ax_u[c], u_coarse, z_coarse, linewidth=3, label="u, CATKE, Δz = $Δz_coarse_str m",  color=coarse_color)

        if :v ∈ keys(coarse_ip.observations[c].field_time_serieses)
            v_obs    = interior(coarse_ip.observations[c].field_time_serieses.v[Nt], 1, 1, :)
            v_coarse = interior(coarse_ip.time_series_collector.field_time_serieses.v[Nt], k_plot, c, :)

            lines!(ax_u[c], v_obs,    z_coarse,   color=sim_color,    linewidth=4, linestyle=:dash) #, label="v, LES")
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

