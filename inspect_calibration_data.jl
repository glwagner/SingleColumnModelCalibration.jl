using ParameterEstimocean
using Oceananigans
using GLMakie

cases = ["free_convection",
         "strong_wind_weak_cooling",
         "med_wind_med_cooling",
         "weak_wind_strong_cooling",
         "strong_wind",
         "strong_wind_no_rotation"]

suite = "two_day_suite"
case_path(case) = joinpath("data", suite, case * "_instantaneous_statistics.jld2")

fontsize_theme = Theme(fontsize = 28, linewidth=4)
set_theme!(fontsize_theme)

field_names = (:u, :v, :b, :e)
observations = [SyntheticObservations(case_path(case); field_names) for case in cases]

obs1 = first(observations)
z = znodes(Center, obs1.grid)
t = obs1.times
Nt = length(t)

@show obs1

fig = Figure(resolution=(1400, 600))
axb = Axis(fig[1, 1], xlabel="Buoyancy (m s⁻²)", ylabel="z (m)")
axe = Axis(fig[1, 2], xlabel="Turbulent kinetic energy (m² s⁻²)", ylabel="z (m)", yaxisposition=:right)
# axu = Axis(fig[1, 3], xlabel="Velocities (m s⁻¹)", ylabel="z (m)", yaxisposition=:right, xscale=log10)

hidespines!(axb, :r, :t)
hidespines!(axe, :l, :t, :l)
#hidespines!(axu, :r, :t)

blines = []

for (case, obs) in zip(cases, observations)
    bt = obs.field_time_serieses.b
    ut = obs.field_time_serieses.u
    vt = obs.field_time_serieses.v
    et = obs.field_time_serieses.e
    b = bt[Nt]
    e = et[Nt]

    #=
    u = ut[Nt]
    v = vt[Nt]
    s = compute!(Field(sqrt(∂z(u)^2 + ∂z(v)^2)))
    zf = znodes(s)
    =#

    label = replace(case, "_" => " ")
    bline = lines!(axb, interior(b, 1, 1, :) .- b[1, 1, 1], z; label)
            #lines!(axu, interior(u, 1, 1, :), z; label)
            #lines!(axu, interior(v, 1, 1, :), z; label, linestyle=:dash)
            #lines!(axu, interior(s, 1, 1, :), zf; label)
            #case != "free_convection" && lines!(axu, interior(s, 1, 1, :), zf; label)
            lines!(axe, interior(e, 1, 1, :), z; label)
    push!(blines, bline)
end

obs = first(observations)
b0 = obs.field_time_serieses.b[1]
bline = lines!(axb, interior(b0, 1, 1, :) .- b0[1, 1, 1], z;
               color = (:black, 0.2),
               #linestyle = :dot,
               linewidth = 6,
               label = "Initial condition")

labels = [replace(case, "_" => " ") for case in cases]
push!(labels, "Initial condition")

push!(blines, bline)
Legend(fig[1, 2], blines, labels,
       framevisible = false,
       tellwidth = false,
       tellheight = false,
       margin = (10, 10, 10, 10),
       halign = :right,
       valign = :bottom)

#axislegend(axb, position=:rb)
#axislegend(axe, position=:rb)

display(fig)

save("two_day_suite.png", fig)

