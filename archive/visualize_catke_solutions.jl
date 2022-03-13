using Oceananigans
using OceanTurbulenceParameterEstimation
using Oceananigans.Units
using DataDeps
using GLMakie
using Printf

suite = "six_day_suite"
case_path(case) = @datadep_str("$(suite)_1m/$(case)_instantaneous_statistics.jld2")

cases = [
         "free_convection",
         "weak_wind_strong_cooling",
         "strong_wind_weak_cooling",
         "strong_wind",
         "strong_wind_no_rotation",
        ]

filename = string("catke_simulation_", suite, "_", Δz, "m.jld2")
u_catke_time_series = FieldTimeSeries(filename, "u", boundary_conditions=nothing)
v_catke_time_series = FieldTimeSeries(filename, "v", boundary_conditions=nothing)
T_catke_time_series = FieldTimeSeries(filename, "T", boundary_conditions=nothing)
e_catke_time_series = FieldTimeSeries(filename, "e", boundary_conditions=nothing)

Δz = 4
regrid_size = (1, 1, Int(256/Δz))
observations = [SyntheticObservations(case_path(case); field_names=(:u, :v, :T, :e), regrid_size) for case in cases]

#####
##### Animate!
#####

fig = Figure(resolution=(1800, 800))

T_axs = []
u_axs = []
#e_axs = []

for (i, case) in enumerate(cases)
    ax_l = Label(fig[1, i+1], replace(case, "_" => "\n"), tellwidth=false)

    if i == 1
        ax_T = Axis(fig[2, i+1], ylabel = "z (m)")
        ax_u = Axis(fig[3, i+1], ylabel = "z (m)")
    else
        ax_T = Axis(fig[2, i+1])
        ax_u = Axis(fig[3, i+1])
    end

    #ax_e = Axis(fig[4, i+1])

    #xlims!(ax_e, 0, 0.005) 
    xlims!(ax_u, -0.3, 0.5)

    #hidexdecorations!(ax_T)
    #hidexdecorations!(ax_u)
    #hidexdecorations!(ax_e)

    if i > 1
        hideydecorations!(ax_T)
        hideydecorations!(ax_u)
    end

    #hidexdecorations!(ax_e)

    if i > 1
        hidespines!(ax_T, :r, :t, :l)
        hidespines!(ax_u, :r, :t, :l)
    else
        hidespines!(ax_T, :t, :r)
        hidespines!(ax_u, :t, :r)
    end

    #hidespines!(ax_e)

    push!(T_axs, ax_T)
    push!(u_axs, ax_u)
    #push!(e_axs, ax_e)
end

Nt = length(first(observations).times)
slider = Slider(fig[5, :], range=1:Nt, startvalue=1)
n = slider.value

###
### Vertical velocity
###

z = znodes(T_catke_time_series)

labels = ["constant Pr", "variable Pr", "convective adjustment"]

for (j, case) in enumerate(cases)
    obs = observations[j]

    nobs = @lift begin
        model_time = u_catke_time_series.times[$n]
        findfirst(t -> t ≈ model_time, obs.times)
    end

    u_les = @lift interior(obs.field_time_serieses.u[$nobs])[1, 1, :]
    v_les = @lift interior(obs.field_time_serieses.v[$nobs])[1, 1, :]
    T_les = @lift interior(obs.field_time_serieses.T[$nobs])[1, 1, :]
    e_les = @lift interior(obs.field_time_serieses.e[$nobs])[1, 1, :]

    ax_T = T_axs[j]
    ax_u = u_axs[j]
    #ax_e = e_axs[j]

    linewidth = 5
    color = :indigo
    lines!(ax_T, T_les, z; linewidth, color=(color, 0.8))
    lines!(ax_u, u_les, z; linewidth, color=(color, 0.8), label = "u, LES")
    lines!(ax_u, v_les, z; linewidth, color=(color, 0.5), label = "v, LES")
    #lines!(ax_e, e_les, z; linewidth, color=(color, 0.8))

    linewidth = 2
    colors = [:darkorange3, :black, :lightseagreen]
    for i = 1:3
        u_catke = @lift interior(u_catke_time_series[$n])[i, j, :]
        v_catke = @lift interior(v_catke_time_series[$n])[i, j, :]
        T_catke = @lift interior(T_catke_time_series[$n])[i, j, :]
        e_catke = @lift interior(e_catke_time_series[$n])[i, j, :]

        color = colors[i]
        label = labels[i]
        lines!(ax_T, T_catke, z; linewidth, color)
        lines!(ax_u, u_catke, z; linewidth, color, label = "u, $label")
        lines!(ax_u, v_catke, z; linewidth, color=(color, 0.5), label = "v, $label")
        #lines!(ax_e, e_catke, z; linewidth, color)
    end

end

times = u_catke_time_series.times
txt = @lift string("Five boundary layers after ", prettytime(times[$n]))
Label(fig[0, 2:5], txt)

Label(fig[3, 1], "Temperature \n (ᵒC)", tellheight=false)
Label(fig[4, 1], "Velocity \n (m s⁻¹)", tellheight=false)
#Label(fig[4, 1], "Turbulent \n kinetic \n energy", tellheight=false)
Legend(fig[4, 7], u_axs[1], tellheight=false, framevisible=false)

display(fig)

#=
record(fig, string("catke_", suite, ".mp4"), 10:800; framerate=16) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end
=#
