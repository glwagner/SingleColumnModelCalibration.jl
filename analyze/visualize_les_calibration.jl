using Oceananigans
using JLD2
using GLMakie

set_theme!(Theme(fontsize=24))

dir = "/Users/gregorywagner/Projects/LocalOceanClosureCalibration/data/one_day_suite/384x384x384"
prefix = "weak_wind_strong_cooling"

xy_filepath = joinpath(dir, prefix * "_xy_slice.jld2")
yz_filepath = joinpath(dir, prefix * "_yz_slice.jld2")
xz_filepath = joinpath(dir, prefix * "_xz_slice.jld2")

w_xy_t = FieldTimeSeries(xy_filepath, "w")
w_xz_t = FieldTimeSeries(xz_filepath, "w")
w_yz_t = FieldTimeSeries(yz_filepath, "w")

statistics_filepath = joinpath(dir, "weak_wind_strong_cooling_time_averaged_statistics.jld2")

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

fig = Figure(resolution=(1600, 800))

azimuth = 6.8
elevation = 0.50
perspectiveness = 1
xlabel = "x (m)"
ylabel = "y (m)"
zlabel = "z (m)"
aspect = :data
fig = Figure(resolution=(2400, 1200))
ax_w = fig[1:8, 1:8] = Axis3(fig; aspect, xlabel, ylabel, zlabel, azimuth, elevation, perspectiveness)

n = 140
w_xy = interior(w_xy_t[n], :, :, 1)
w_xz = interior(w_xz_t[n], :, 1, :)
w_yz = interior(w_yz_t[n], 1, :, :)

wlim = maximum(abs, w_xy) / 2
colorrange_w = (-wlim, wlim)
colormap_w = :balance

pl = surface!(ax_w, x_xz, y_xz, z_xz; color=w_xz, colormap=colormap_w, colorrange=colorrange_w)
     surface!(ax_w, x_yz, y_yz, z_yz; color=w_yz, colormap=colormap_w, colorrange=colorrange_w)
     surface!(ax_w, x_xy, y_xy, z_xy; color=w_xy, colormap=colormap_w, colorrange=colorrange_w)

cp = fig[2, 4:5] = Colorbar(fig, pl, vertical=false, flipaxis=true, label="Vertical velocity (m s⁻¹)", tellwidth=false, tellheight=false)

Bt = FieldTimeSeries(statistics_filepath, "b")
Tt = FieldTimeSeries(statistics_filepath, "T")
Ut = FieldTimeSeries(statistics_filepath, "u")
Vt = FieldTimeSeries(statistics_filepath, "v")

z = znodes(Bt)

ax_b = Axis(fig[3:7, 7], xlabel="Temperature (ᵒC)", ylabel="z (m)")
ax_u = Axis(fig[3:7, 8], xlabel="Velocities (m s⁻¹)", ylabel="z (m)", yaxisposition=:right)

hidespines!(ax_b, :t, :r)
hidespines!(ax_u, :t, :l)

B = interior(Bt[7], 1, 1, :)
T = interior(Tt[7], 1, 1, :)
U = interior(Ut[7], 1, 1, :)
V = interior(Vt[7], 1, 1, :)

#lines!(ax_b, B, z, linewidth=4)
lines!(ax_b, T, z, linewidth=4)
lines!(ax_u, U, z, linewidth=4, label="u")
lines!(ax_u, V, z, linewidth=4, label="v")

axislegend(ax_u, position=:rb)

display(fig)

save("LES_summary.png", fig)
