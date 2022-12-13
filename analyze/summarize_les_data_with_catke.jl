using Oceananigans
using JLD2
using GLMakie

using SingleColumnModelCalibration:
    dependent_parameter_sets,
    build_batched_inverse_problem,
    prior_library

set_theme!(Theme(fontsize=32))

suite = "24_hour_suite"
n = 70
ulims = (-0.15, 0.35)
elims = (-1e-4, 4e-3)
linewidth = 12
slicesdir = "/Users/gregorywagner/Projects/SingleColumnModelCalibration/data/slices/$suite/"
profilesdir = "/Users/gregorywagner/Projects/SingleColumnModelCalibration/data/profiles/$suite/1m"

#prefix = "weak_wind_strong_cooling"
prefix = "free_convection"
c1 = 1
xy_filepath = joinpath(slicesdir, prefix * "_xy_slice.jld2")
yz_filepath = joinpath(slicesdir, prefix * "_yz_slice.jld2")
xz_filepath = joinpath(slicesdir, prefix * "_xz_slice.jld2")
statistics_filepath = joinpath(profilesdir, prefix * "_instantaneous_statistics.jld2")

w_xy_t = FieldTimeSeries(xy_filepath, "w")
w_xz_t = FieldTimeSeries(xz_filepath, "w")
w_yz_t = FieldTimeSeries(yz_filepath, "w")

Tt = FieldTimeSeries(statistics_filepath, "T")
bt = FieldTimeSeries(statistics_filepath, "b")
Ut = FieldTimeSeries(statistics_filepath, "u")
Vt = FieldTimeSeries(statistics_filepath, "v")
Et = FieldTimeSeries(statistics_filepath, "e")

times = w_xy_t.times
Nt = length(times)
grid = w_xy_t.grid
Nx, Ny, Nz = size(grid)
x, y, z = nodes(w_xy_t)
Lx = grid.Lx
Ly = grid.Ly
Lz = grid.Lz

Nz += 1

x_xz = repeat(x, 1, Nz)
z_xz = repeat(reshape(z, 1, Nz), Nx, 1)
y_xz = 0.995 * Ly * ones(Nx, Nz)

y_yz = repeat(y, 1, Nz)
z_yz = repeat(reshape(z, 1, Nz), grid.Ny, 1)
x_yz = 0.995 * Lx * ones(Ny, Nz)

# Slight displacements to "stitch" the cube together
x_xy = x
y_xy = y
z_xy = - 0.001 * Lz * ones(Nx, Ny)

fig = Figure(resolution=(2200, 1600))

azimuth = 6.8
elevation = 0.50
perspectiveness = 0.5
Tticks = [19.65, 19.72, 19.8]
eticks = [0.0, 2e-3, 4e-3]

xlabel = "x (m)"
ylabel = "y (m)"
zlabel = "z (m)"

xlabeloffset = 100
ylabeloffset = 80
zlabeloffset = 100

ylabelpadding = 20.0

w = 5
h = 9
aspect = :data
ax_w1 = fig[h+1:2h, 1:w] = Axis3(fig; aspect, xlabel, ylabel, zlabel, azimuth, elevation=0.6, perspectiveness,
                                 xlabeloffset, ylabeloffset, zlabeloffset)
ax_w2 = fig[1:h, 1:w]    = Axis3(fig; aspect, xlabel, ylabel, zlabel, azimuth, elevation, perspectiveness,
                                 xlabeloffset, ylabeloffset, zlabeloffset)

k1 = k = h + 1 + Int((h-1)/2)
ax_b1 = Axis(fig[k-1:k+3, w+1]; xlabel="Buoyancy (m s⁻²)", ylabel="z (m)", xticks=[0.0386, 0.0390], ylabelpadding)
ax_u1 = Axis(fig[k-1:k+3, w+2]; xlabel="Velocities (m s⁻¹)", xticks=[-0.05, 0.0, 0.05])
ax_e1 = Axis(fig[k-1:k+3, w+3]; xlabel="Turbulent kinetic \n energy (m² s⁻²)", ylabel="z (m)", yaxisposition=:right,
                                xticks=[0.0, 2e-4, 4e-4], ylabelpadding)
Label(fig[k-2, w+1:w+3], text="Free convection", textsize=30)

k2 = k = 1 + Int((h-1)/2)
ax_b2 = Axis(fig[k-1:k+3, w+1]; xlabel="Buoyancy (m s⁻²)", ylabel="z (m)", ylabelpadding)
ax_u2 = Axis(fig[k-1:k+3, w+2]; xlabel="Velocities (m s⁻¹)", xticks=[-0.1, 0.0, 0.1, 0.2, 0.3])
ax_e2 = Axis(fig[k-1:k+3, w+3]; xlabel="Turbulent kinetic \n energy (m² s⁻²)", ylabel="z (m)", yaxisposition=:right,
             xticks=eticks, ylabelpadding)
Label(fig[k-2, w+1:w+3], text="Strong wind, weak cooling", textsize=30)

colsize!(fig.layout, 1, Relative(1/2))
colgap!(fig.layout, 1, -2200)
colgap!(fig.layout, w, -80)
colgap!(fig.layout, w+1, 40)
colgap!(fig.layout, w+2, 40)

rowgap!(fig.layout, h, -100)
rowgap!(fig.layout, k1+2, 200)
rowgap!(fig.layout, k2+2, 200)
rowgap!(fig.layout, k1+3, -50)
rowgap!(fig.layout, k2+3, -50)


w_xy = interior(w_xy_t[n], :, :, 1)
w_xz = interior(w_xz_t[n], :, 1, :)
w_yz = interior(w_yz_t[n], 1, :, :)

wlim = maximum(abs, w_xy) / 2
colorrange_w = (-wlim, wlim)
colormap_w = :balance

pl = surface!(ax_w1, x_xz, y_xz, z_xz; color=w_xz, colormap=colormap_w, colorrange=colorrange_w)
     surface!(ax_w1, x_yz, y_yz, z_yz; color=w_yz, colormap=colormap_w, colorrange=colorrange_w)
     surface!(ax_w1, x_xy, y_xy, z_xy; color=w_xy, colormap=colormap_w, colorrange=colorrange_w)

Colorbar(fig[0, 1:w], pl, vertical=false, flipaxis=true, width=700,
         label="Vertical velocity (m s⁻¹)", tellwidth=false, tellheight=false)

rowgap!(fig.layout, 1, -100)

#xlims!(ax_b1, 0.0385, 0.0392)
xlims!(ax_u1, -0.1, 0.1) #ulims...)
#xlims!(ax_e1, elims...)
#xlims!(ax_b2, 0.0385, 0.0392)
xlims!(ax_u2, ulims...)
#xlims!(ax_e2, elims...)

hidespines!(ax_b1, :t, :r)
hidespines!(ax_u1, :t, :l, :r)
hidespines!(ax_e1, :t, :l)

hidespines!(ax_b2, :t, :r)
hidespines!(ax_u2, :t, :l, :r)
hidespines!(ax_e2, :t, :l)

hideydecorations!(ax_u1, grid=false)
hideydecorations!(ax_u2, grid=false)

T = interior(Tt[Nt], 1, 1, :)
b = interior(bt[Nt], 1, 1, :)
U = interior(Ut[Nt], 1, 1, :)
V = interior(Vt[Nt], 1, 1, :)
E = interior(Et[Nt], 1, 1, :)
z = znodes(Ut)

lines!(ax_b1, b, z; color=(:royalblue1, 0.6), linewidth, label="LES")
lines!(ax_u1, U, z; color=(:royalblue1, 0.6), linewidth, label="u, LES")
lines!(ax_u1, V, z; color=(:orange1, 0.6),    linewidth, label="v, LES")
lines!(ax_e1, E, z; color=(:royalblue1, 0.6), linewidth, label="LES")

#####
##### Second plot
#####

prefix = "strong_wind"
c2 = 5
xy_filepath = joinpath(slicesdir, prefix * "_xy_slice.jld2")
yz_filepath = joinpath(slicesdir, prefix * "_yz_slice.jld2")
xz_filepath = joinpath(slicesdir, prefix * "_xz_slice.jld2")
statistics_filepath = joinpath(profilesdir, prefix * "_instantaneous_statistics.jld2")

w_xy_t = FieldTimeSeries(xy_filepath, "w")
w_xz_t = FieldTimeSeries(xz_filepath, "w")
w_yz_t = FieldTimeSeries(yz_filepath, "w")

Tt = FieldTimeSeries(statistics_filepath, "T")
bt = FieldTimeSeries(statistics_filepath, "b")
Ut = FieldTimeSeries(statistics_filepath, "u")
Vt = FieldTimeSeries(statistics_filepath, "v")
Et = FieldTimeSeries(statistics_filepath, "e")

w_xy = interior(w_xy_t[n], :, :, 1)
w_xz = interior(w_xz_t[n], :, 1, :)
w_yz = interior(w_yz_t[n], 1, :, :)

wlim = maximum(abs, w_xy) / 2
colorrange_w = (-wlim, wlim)
colormap_w = :balance

pl = surface!(ax_w2, x_xz, y_xz, z_xz; color=w_xz, colormap=colormap_w, colorrange=colorrange_w)
     surface!(ax_w2, x_yz, y_yz, z_yz; color=w_yz, colormap=colormap_w, colorrange=colorrange_w)
     surface!(ax_w2, x_xy, y_xy, z_xy; color=w_xy, colormap=colormap_w, colorrange=colorrange_w)

T = interior(Tt[Nt], 1, 1, :)
b = interior(bt[Nt], 1, 1, :)
U = interior(Ut[Nt], 1, 1, :)
V = interior(Vt[Nt], 1, 1, :)
E = interior(Et[Nt], 1, 1, :)
z = znodes(Ut)

lines!(ax_b2, b, z; color=(:royalblue1, 0.6), linewidth, label="LES")
lines!(ax_u2, U, z; color=(:royalblue1, 0.6), linewidth, label="u, LES")
lines!(ax_u2, V, z; color=(:orange, 0.6),     linewidth, label="v, LES")
lines!(ax_e2, E, z; color=(:royalblue1, 0.6), linewidth, label="LES")

#####
##### Run CATKE
#####

dir = "../parameters"
#name = "constant_Pr"
#name = "variable_Pr"
#name = "constant_Pr_conv_adj"
name = "variable_Pr_conv_adj"
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

grid_parameters = [(size=32, z=(-256, 0))]
suite_parameters = [(name = suite, resolution="1m", stop_time=24hours)]

batched_ip = build_batched_inverse_problem(closure, name;
                                           Nensemble = 1,
                                           Δt = 1minutes,
                                           grid_parameters,
                                           suite_parameters)

forward_run!(batched_ip, [optimal_parameters])

ip = batched_ip[1]

times = observation_times(ip.observations)
Nt = length(times)
b1 = interior(ip.time_series_collector.field_time_serieses.b[Nt], 1, c1, :)
u1 = interior(ip.time_series_collector.field_time_serieses.u[Nt], 1, c1, :)
v1 = interior(ip.time_series_collector.field_time_serieses.v[Nt], 1, c1, :)
e1 = interior(ip.time_series_collector.field_time_serieses.e[Nt], 1, c1, :)

z = znodes(ip.time_series_collector.field_time_serieses.b)

linestyle=:solid
color=:darkgreen
lines!(ax_b1, b1, z; color,   linestyle, linewidth=4, label="CATKE")
lines!(ax_u1, u1, z; color,   linestyle, linewidth=4, label="u, CATKE")
lines!(ax_u1, v1, z; color=:darkred, linestyle, linewidth=4, label="v, CATKE")
lines!(ax_e1, e1, z; color,   linestyle, linewidth=4, label="CATKE")

b2 = interior(ip.time_series_collector.field_time_serieses.b[Nt], 1, c2, :)
u2 = interior(ip.time_series_collector.field_time_serieses.u[Nt], 1, c2, :)
v2 = interior(ip.time_series_collector.field_time_serieses.v[Nt], 1, c2, :)
e2 = interior(ip.time_series_collector.field_time_serieses.e[Nt], 1, c2, :)

lines!(ax_b2, b2, z; color,   linestyle, linewidth=4, label="CATKE")
lines!(ax_u2, u2, z; color,   linestyle, linewidth=4, label="u, CATKE")
lines!(ax_u2, v2, z; color=:darkred, linestyle, linewidth=4, label="v, CATKE")
lines!(ax_e2, e2, z; color,   linestyle, linewidth=4, label="CATKE")

axislegend(ax_b2, position=:rb)
axislegend(ax_b1, position=:rb)
axislegend(ax_u2, position=:rb)

display(fig)

save("les_data_summary_with_catke.png", fig)
