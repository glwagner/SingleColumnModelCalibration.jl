using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using GLMakie
using Printf

#=
Lz = 128
Δz = 8
zb = -(Lz + Δz) / 2
zt = +(Lz + Δz) / 2
z = collect(zb:Δz:zt)
=#

Lz = 128
Δz = 2
zL = (Lz + Δz) / 2
z = [-zL]
while z[end] < -Δz
    global Δz
    Δz = 1.01 * Δz
    push!(z, z[end] + Δz)
end

append!(z, -reverse(z))

@show z

Nz = length(z) - 1

grid = RectilinearGrid(size=Nz, z=z, topology=(Flat, Flat, Bounded))

closure = CATKEVerticalDiffusivity()
model = HydrostaticFreeSurfaceModel(; grid, buoyancy=BuoyancyTracer(), tracers=(:b, :e), closure)

Δ = 16
N² = 1e-6
Ri = 0.1 # N² / S²
dUdz = sqrt(N² / Ri)
@show U = dUdz * Δ
uᵢ(z) = U * tanh(z / Δ)
bᵢ(z) = N² * z
set!(model, u=uᵢ, b=bᵢ, e=1e-6)

u = model.velocities.u
for k = 1:Nz
    @info @sprintf("u[%d] = %.2f", k, u[1, 1, k])
end

simulation = Simulation(model, Δt=10minutes, stop_time=24hours)

bt = []
ut = []
et = []
κct = []
κut = []
κet = []
Rit = []
N²t = []
S²t = []

b = model.tracers.b
u = model.velocities.u
e = model.tracers.e
κc = model.diffusivity_fields.κᶜ
κu = model.diffusivity_fields.κᵘ
κe = model.diffusivity_fields.κᵉ
Ri = Field(∂z(b) / ∂z(u)^2)
N² = Field(∂z(b))
S² = Field(∂z(u)^2)

function collect_data(sim)
    compute!(Ri)
    compute!(N²)
    compute!(S²)
    push!(Rit, deepcopy(interior(Ri, 1, 1, :)))
    push!(N²t, deepcopy(interior(N², 1, 1, :)))
    push!(S²t, deepcopy(interior(S², 1, 1, :)))
    push!(bt,  deepcopy(interior(b, 1, 1, :)))
    push!(ut,  deepcopy(interior(u, 1, 1, :)))
    push!(et,  deepcopy(interior(e, 1, 1, :)))
    push!(κct, deepcopy(interior(κc, 1, 1, :)))
    push!(κut, deepcopy(interior(κu, 1, 1, :)))
    push!(κet, deepcopy(interior(κe, 1, 1, :)))

    return nothing
end

t = 0:10minutes:simulation.stop_time
Nt = length(t)

simulation.callbacks[:dc] = Callback(collect_data, SpecifiedTimes(t))

run!(simulation)

fig = Figure(resolution=(1400, 800))

axN0 = Axis(fig[1, 1], xlabel="N² (s⁻²)", ylabel="z (m)", xticks=[0.0, 2e-6])
axR0 = Axis(fig[1, 2], xlabel="U (m s⁻¹)", xticks=[-0.06, 0.0, 0.06])
axe0 = Axis(fig[2, 2], xlabel="E (m² s⁻²)", ylabel="z (m)", xticks=[0.0, 2e-5])
axκ0 = Axis(fig[3, 2], xlabel="κ (m² s⁻¹)", ylabel="z (m)", xticks=[0.0, 2e-3, 4e-3])

axR = Axis(fig[1, 3], xlabel="Time (hr)", ylabel="z (m)", xaxisposition=:top, yaxisposition=:right)
axe = Axis(fig[2, 3], ylabel="z (m)", yaxisposition=:right)
axκ = Axis(fig[3, 3], xlabel="Time (hr)", ylabel="z (m)", yaxisposition=:right)

xlims!(axN0, 0.0, 3e-6)
xlims!(axR0, -0.06, 0.06)
xlims!(axe0, -5e-6, 3e-5)

xlims!(axR, 0, 18)
xlims!(axe, 0, 18)
xlims!(axκ, 0, 18)

hidexdecorations!(axe)
hideydecorations!(axR0, grid=false)

hidespines!(axR0, :l, :r, :t)
hidespines!(axN0, :r, :t)
hidespines!(axe0, :r, :t)
hidespines!(axκ0, :r, :t)

Prt = map((κu, κc) -> κu ./ κc, κut, κct)

bzt = hcat(bt...)'
Rzt = hcat(Rit...)'
ezt = hcat(et...)'
κzt = hcat(κct...)'
Pzt = hcat(Prt...)'

zc = znodes(grid, Center())
zf = znodes(grid, Face())

for n = (7, 49, 24*4+1)
    @show tn = t[n] / hour
    label = @sprintf("t = %d hr", tn)
    #lines!(axR0, 1 ./ Rit[n], zf; label)
    #lines!(axR0, S²t[n][2:Nz], zf[2:Nz]; label)
    lines!(axR0, ut[n], zc; label)
    lines!(axN0, N²t[n][2:Nz], zf[2:Nz]; label)
    #lines!(axN0, S²t[n][2:Nz], zf[2:Nz]; label, linewidth=4, linestyle=:dash)
    ln = lines!(axe0, et[n], zc)
    lines!(axκ0, κct[n], zf, color = ln.color.val, label="κᶜ")
    lines!(axκ0, κut[n], zf, linestyle=:dash, color=ln.color.val, label="κᵘ")
    #lines!(axκ0, κet[n], zf, linestyle=:dot, color=ln.color.val, label="κᵉ")
end

Legend(fig[2, 1], axR0)
Legend(fig[3, 1], axκ0, merge=true)

vlines!(axR0, 4, color=(:gray, 0.5), linewidth=4)

cr = contourf!(axR, t ./ hour, zf, 1 ./ Rzt, colormap=:viridis)
Colorbar(fig[1, 4], cr, label="Ri⁻¹ = ∂z(U) / N²")

cr = contourf!(axe, t ./ hour, zc, ezt, colormap=:heat)
Colorbar(fig[2, 4], cr, label="E (m² s⁻²)")

cr = contourf!(axκ, t ./ hour, zf, κzt, colormap=:solar)
Colorbar(fig[3, 4], cr, label="κᶜ (m² s⁻¹)")

colsize!(fig.layout, 1, Relative(0.1))
colsize!(fig.layout, 2, Relative(0.1))

display(fig)

save("shear_layer.png", fig)
