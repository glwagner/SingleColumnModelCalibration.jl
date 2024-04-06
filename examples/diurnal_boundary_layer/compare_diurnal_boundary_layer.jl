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

sc2_filename = "catke_diurnal_boundary_layer_dz2.jld2"
sc4_filename = "catke_diurnal_boundary_layer_dz4.jld2"
sc8_filename = "catke_diurnal_boundary_layer_dz8.jld2"

#les_filename = "wave_averaged_diurnal_boundary_layer_les_Nx128_Nz128_averages.jld2"
#les_filename = "diurnal_boundary_layer_les_Nx128_Nz128_averages.jld2"
les_filename = "wave_averaged_diurnal_boundary_layer_les_Nx64_Nz64_averages.jld2"
#les_filename = "diurnal_boundary_layer_les_averages.jld2"

Bt_les = FieldTimeSeries(les_filename, "B")
Bt_sc2 = FieldTimeSeries(sc2_filename, "b")
Bt_sc4 = FieldTimeSeries(sc4_filename, "b")
Bt_sc8 = FieldTimeSeries(sc8_filename, "b")

t_les = Bt_les.times
t_sc2 = Bt_sc2.times
t_sc4 = Bt_sc4.times
t_sc8 = Bt_sc8.times
z_les = znodes(Bt_les)
z_sc2 = znodes(Bt_sc2)
z_sc4 = znodes(Bt_sc4)
z_sc8 = znodes(Bt_sc8)

fig = Figure(size=(1200, 1000))
yticks = [-60, -30, 0]

axb_sfc = Axis(fig[1, 1]; xlabel="Time (hr)", ylabel="b(z=0) (m s⁻²)", xaxisposition=:top, yaxisposition=:right)
axb_les = Axis(fig[2, 1]; xlabel="Time (hr)", ylabel="z (m)", xaxisposition=:top, yaxisposition=:right, yticks)
axb_sc4 = Axis(fig[3, 1]; xlabel="Time (hr)", ylabel="z (m)", xaxisposition=:top, yaxisposition=:right, yticks)

Nz_les = size(Bt_les, 3)
Nz_sc2 = size(Bt_sc2, 3)
Nz_sc4 = size(Bt_sc4, 3)
Nz_sc8 = size(Bt_sc8, 3)

coarse_b = CenterField(Bt_sc8.grid)

function compute_surface_average(Bt)
    Nt = length(Bt)
    b₀ = zeros(Nt)
    for n = 1:Nt
        regrid!(coarse_b, Bt[n])
        Nz = size(coarse_b, 3)
        @inbounds b₀[n] = coarse_b[1, 1, Nz]
    end
    return b₀
end

# Compute surface-averaged buoyancy
b₀_les = compute_surface_average(Bt_les)
b₀_sc2 = compute_surface_average(Bt_sc2)
b₀_sc4 = compute_surface_average(Bt_sc4)

α = 2e-4
g = 9.81
T₀ = 30

b₀_sc8 = interior(Bt_sc8, 1, 1, Nz_sc8, :)
T₀_les = @. T₀ + b₀_les / (α * g)  
T₀_sc2 = @. T₀ + b₀_sc2 / (α * g)
T₀_sc4 = @. T₀ + b₀_sc4 / (α * g)
T₀_sc8 = @. T₀ + b₀_sc8 / (α * g)

#lines!(axb_sfc, t_les ./ hours, interior(Bt_les, 1, 1, Nz_les, :), label="LES")
#lines!(axb_sfc, t_sc4 ./ hours, interior(Bt_sc2, 1, 1, Nz_sc2, :), label="CATKE, Δz = 2 m")
#lines!(axb_sfc, t_sc4 ./ hours, interior(Bt_sc4, 1, 1, Nz_sc4, :), label="CATKE, Δz = 4 m")

lines!(axb_sfc, t_les ./ hours, T₀_les, label="LES")
lines!(axb_sfc, t_sc2 ./ hours, T₀_sc2, label="CATKE, Δz = 2 m")
lines!(axb_sfc, t_sc4 ./ hours, T₀_sc4, label="CATKE, Δz = 4 m")
lines!(axb_sfc, t_sc8 ./ hours, T₀_sc8, label="CATKE, Δz = 8 m")
axislegend(axb_sfc, position=:lt)


byz_les = interior(Bt_les, 1, 1, :, :)'
byz_sc4 = interior(Bt_sc4, 1, 1, :, :)'
Tt_les = @. T₀ + byz_les / (α * g)
Tt_sc4 = @. T₀ + byz_sc4 / (α * g)

Tmin = minimum(Tt_les)
Tmax = maximum(Tt_les)
levels = Tmin:0.1:Tmax

cr = contourf!(axb_les, t_les ./ hours, z_les, Tt_les, levels=levels, colormap=:thermal)
cr = contourf!(axb_sc4, t_sc4 ./ hours, z_sc4, Tt_sc4, levels=levels, colormap=:thermal)

#Colorbar(fig[2:3, 2], cr, label="b (m s⁻²)") #, ticks=([0, 3e-5, 6e-5], ["0", "3×10⁻⁵", "6×10⁻⁵"]))
Colorbar(fig[2:3, 2], cr, label="T (ᵒC)") #, ticks=([0, 3e-5, 6e-5], ["0", "3×10⁻⁵", "6×10⁻⁵"]))

ylims!(axb_les, -120, 0)
ylims!(axb_sc4, -120, 0)

