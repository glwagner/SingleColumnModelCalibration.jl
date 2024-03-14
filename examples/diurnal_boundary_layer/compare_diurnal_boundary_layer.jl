using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity

#using CairoMakie
using GLMakie
using Printf
using Statistics

#####
##### Validation against LES
#####

set_theme!(Theme(fontsize=22))

les_filename = "diurnal_boundary_layer_les_averages.jld2"
Bt = FieldTimeSeries(les_filename, "B")
t_les = Bt.times
z_les = znodes(Bt)

bzt = hcat(bt...)'
zc = znodes(grid, Center())
t_scm = 0:simulation.Δt:simulation.stop_time

fig = Figure(size=(1200, 1000))
yticks = [-60, -30, 0]

axb_sfc = Axis(fig[1, 1]; xlabel="Time (hr)", ylabel="b(z=0) (m s⁻²)", xaxisposition=:top, yaxisposition=:right)
axb_les = Axis(fig[2, 1]; xlabel="Time (hr)", ylabel="z (m)", xaxisposition=:top, yaxisposition=:right, yticks)
axb_scm = Axis(fig[3, 1]; xlabel="Time (hr)", ylabel="z (m)", xaxisposition=:top, yaxisposition=:right, yticks)

Nz_les = size(Bt, 3)
Nz_scm = Nz
lines!(axb_sfc, t_les ./ hours, interior(Bt, 1, 1, Nz_les, :))
lines!(axb_sfc, t_scm ./ hours, view(bzt, :, Nz_scm))

cr = contourf!(axb_les, t_les ./ hours, z_les, interior(Bt, 1, 1, :, :)', levels=12, colormap=:viridis)
cr = contourf!(axb_scm, t_scm ./ hours, zc, bzt .- bzt[1, Nz], levels=12, colormap=:viridis)

Colorbar(fig[2:3, 2], cr, label="b (m s⁻²)") #, ticks=([0, 3e-5, 6e-5], ["0", "3×10⁻⁵", "6×10⁻⁵"]))

ylims!(axb_les, -80, 0)
ylims!(axb_scm, -80, 0)

