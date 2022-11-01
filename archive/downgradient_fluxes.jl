using ParameterEstimocean
using Oceananigans
using GLMakie

cases = [#"free_convection",
         #"strong_wind_weak_cooling",
         #"med_wind_med_cooling",
         "weak_wind_strong_cooling",
         #"strong_wind",
         #"strong_wind_no_rotation",
         ]

suite = "one_day_suite"
case_path(case) = joinpath("data", suite, "384x384x384", case * "_time_averaged_statistics.jld2")

fontsize_theme = Theme(fontsize=28, linewidth=4)
set_theme!(fontsize_theme)

field_names = (:u, :v, :b, :e, :wb, :T, :wT)
regrid = RectilinearGrid(size=64, z=(-256, 0), topology=(Flat, Flat, Bounded))
observations = [SyntheticObservations(case_path(case); field_names, regrid) for case in cases]

obs1 = first(observations)
z = znodes(Center, obs1.grid)
zf = znodes(Face, obs1.grid)
t = obs1.times
Nt = length(t)

fig = Figure(resolution=(1600, 600))
axT = Axis(fig[1, 1], xlabel="Temperature (ᵒC)", ylabel="z (m)", xticks=LinearTicks(3))
#axe = Axis(fig[1, 2], xlabel="Turbulent kinetic energy (m² s⁻²)", ylabel="z (m)")
axN = Axis(fig[1, 2], xlabel="Temperature gradient (ᵒC m⁻¹)", ylabel="z (m)", xticks=LinearTicks(2))
axq = Axis(fig[1, 3], xlabel="Temperature flux (ᵒC m s⁻¹)", ylabel="z (m)", xticks=LinearTicks(2))
axK = Axis(fig[1, 4], xlabel="Diffusivity (m² s⁻¹)", ylabel="z (m)", yaxisposition=:right)#, xscale=log10)
# axu = Axis(fig[1, 3], xlabel="Velocities (m s⁻¹)", ylabel="z (m)", yaxisposition=:right, xscale=log10)

hidespines!(axT, :r, :t)
hidespines!(axe, :r, :t, :l)
hidespines!(axN, :r, :t, :l)
hidespines!(axq, :r, :t, :l)
hidespines!(axK, :t, :l)

#xlims!(axN, -1e-6, 5e-6)
#xlims!(axK, -1e-3, 1e1)

blines = []

for (case, obs) in zip(cases, observations)
    Tt = obs.field_time_serieses.T
    ut = obs.field_time_serieses.u
    vt = obs.field_time_serieses.v
    et = obs.field_time_serieses.e
    wTt = obs.field_time_serieses.wT

    T = Tt[Nt]
    e = et[Nt]
    wT = wTt[Nt]

    dTdz = compute!(Field(∂z(T)))
    K = compute!(Field(- wT / dTdz))

    #negatives = interior(K) .<= 0
    #interior(K)[negatives] .= NaN

    Nz = T.grid.Nz
    kf = 2:Nz

    label = replace(case, "_" => " ")
    bline = lines!(axT, interior(T, 1, 1, :), z; label)
            lines!(axN, interior(dTdz, 1, 1, kf), zf[kf]; label)
            lines!(axq, interior(wT, 1, 1, kf), zf[kf]; label)
            lines!(axK, interior(K, 1, 1, kf), zf[kf]; label)

    push!(blines, bline)
end

lines!(axK, [0, 0], z[[1, end]];
       color = (:black, 0.2),
       linewidth = 6)

#=
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
=#

#axislegend(axb, position=:rb)
#axislegend(axe, position=:rb)

display(fig)

#save(suite * ".png", fig)

