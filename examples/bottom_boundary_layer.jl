using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    MixingLength,
    TurbulentKineticEnergyEquation

using GLMakie

Lz = 500
Nz = 10

grid = RectilinearGrid(size=Nz, z=(0, Lz), topology=(Flat, Flat, Bounded))

closure = CATKEVerticalDiffusivity()

@inline u_drag(x, y, t, u, v, cᴰ) = - cᴰ * u * sqrt(u^2 + v^2)
@inline v_drag(x, y, t, u, v, cᴰ) = - cᴰ * v * sqrt(u^2 + v^2)

cᴰ = 2e-3
u_bottom_bc = FluxBoundaryCondition(u_drag, field_dependencies=(:u, :v), parameters=cᴰ)
v_bottom_bc = FluxBoundaryCondition(v_drag, field_dependencies=(:u, :v), parameters=cᴰ)
u_bcs = FieldBoundaryConditions(bottom=u_bottom_bc)
v_bcs = FieldBoundaryConditions(bottom=v_bottom_bc)

model = HydrostaticFreeSurfaceModel(; grid, closure,
                                    buoyancy = BuoyancyTracer(),
                                    #coriolis = FPlane(f=1e-4),
                                    tracers = (:b, :e),
                                    boundary_conditions = (; u=u_bcs))

N² = 1e-6
U = 1
uᵢ(x, y, z) = U
bᵢ(x, y, z) = N² * z
set!(model, u=uᵢ, b=bᵢ, e=1e-4)

simulation = Simulation(model, Δt=1minutes, stop_time=2days)

bt = []
ut = []
vt = []
et = []
κct = []
κut = []
Rit = []

b = model.tracers.b
u, v, w = model.velocities
e = model.tracers.e
κc = model.diffusivity_fields.κᶜ
κu = model.diffusivity_fields.κᵘ
Ri = Field(∂z(b) / ∂z(u)^2)

function collect_data(sim)
    compute!(Ri)
    push!(Rit, deepcopy(interior(Ri, 1, 1, :)))
    push!(bt,  deepcopy(interior(b, 1, 1, :)))
    push!(ut,  deepcopy(interior(u, 1, 1, :)))
    push!(vt,  deepcopy(interior(v, 1, 1, :)))
    push!(et,  deepcopy(interior(e, 1, 1, :)))
    push!(κct, deepcopy(interior(κc, 1, 1, :)))
    push!(κut, deepcopy(interior(κu, 1, 1, :)))

    return nothing
end

t = 0:20minutes:simulation.stop_time
Nt = length(t)

simulation.callbacks[:dc] = Callback(collect_data, SpecifiedTimes(t))

run!(simulation)

fig = Figure(resolution=(1600, 600))

axb = Axis(fig[1, 1], xlabel="b", ylabel="z (m)")
axu = Axis(fig[1, 2], xlabel="u", ylabel="z (m)")
axR = Axis(fig[1, 3], xlabel="Ri⁻¹", ylabel="z (m)")
axe = Axis(fig[1, 4], xlabel="e", ylabel="z (m)")
axκ = Axis(fig[1, 5], xlabel="κ", ylabel="z (m)")

xlims!(axb, 0, Lz * N²)
xlims!(axu, 0, 1.1)
xlims!(axR, -5, 100)
xlims!(axe, -1e-5, 2e-4)
xlims!(axκ, -1e-1, 0.2)

slider = Slider(fig[2, 1:4], range=1:Nt, startvalue=Nt)
n = slider.value

bn  = @lift bt[$n]
un  = @lift ut[$n]
vn  = @lift vt[$n]
en  = @lift et[$n]
Rn  = @lift 1 ./ Rit[$n][2:Nz]
κcn = @lift κct[$n][2:Nz]
κun = @lift κut[$n][2:Nz]

zc = znodes(grid, Center())
zf = znodes(grid, Face())

lines!(axb, bn, zc)
lines!(axu, un, zc) #, label="u")
lines!(axe, en, zc)
lines!(axR, Rn, zf[2:Nz])
lines!(axκ, κcn, zf[2:Nz], label="κc")
lines!(axκ, κun, zf[2:Nz], label="κc")

display(fig)

