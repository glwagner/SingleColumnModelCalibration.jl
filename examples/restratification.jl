using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using GLMakie

Lz = 256
Nz = 32

grid = RectilinearGrid(size=Nz, z=(-Lz, 0), topology=(Flat, Flat, Bounded))

const c = Center()

@inline function b_restoring(i, j, k, grid, clock, fields, parameters)
    N² = parameters.N²
    τ = parameters.τ
    z = znode(i, j, k, grid, c, c, c)
    b★ = N² * z
    b = @inbounds fields.b[i, j, k]
    return (b★ - b) / τ
    #return -b / τ
end

Qᵇ = 1e-7
parameters = (N² = 1e-6, τ = 12hours)

b_forcing = Forcing(b_restoring; discrete_form=true, parameters)
top_b_bc = FluxBoundaryCondition(Qᵇ)
b_bcs = FieldBoundaryConditions(top=top_b_bc)

closure = CATKEVerticalDiffusivity()

model = HydrostaticFreeSurfaceModel(; grid,
                                    buoyancy = BuoyancyTracer(),
                                    tracers = (:b, :e),
                                    closure,
                                    forcing = (; b=b_forcing),
                                    boundary_conditions = (; b=b_bcs))

bᵢ(x, y, z) = parameters.N² * z
set!(model, b=bᵢ, e=1e-6)
simulation = Simulation(model, Δt=10minutes, stop_time=7days)

bt = []
et = []
κct = []

b = model.tracers.b
e = model.tracers.e
κc = model.diffusivity_fields.κᶜ

function collect_data(sim)
    push!(bt,  deepcopy(interior(b, 1, 1, :)))
    push!(et,  deepcopy(interior(e, 1, 1, :)))
    push!(κct, deepcopy(interior(κc, 1, 1, :)))

    return nothing
end

t = 0:10minutes:simulation.stop_time
Nt = length(t)

simulation.callbacks[:dc] = Callback(collect_data, SpecifiedTimes(t))

run!(simulation)

fig = Figure(resolution=(1600, 600))

axb = Axis(fig[1, 1])
axe = Axis(fig[1, 2])
axκ = Axis(fig[1, 3])

xlims!(axb, -Lz * parameters.N², 0)
xlims!(axe, -1e-4, 1e-3)
xlims!(axκ, -1e-1, 1e1)

slider = Slider(fig[2, 1:3], range=1:Nt, startvalue=Nt)
n = slider.value

bn  = @lift bt[$n]
en  = @lift et[$n]
κcn = @lift κct[$n][2:Nz]

zc = znodes(grid, Center())
zf = znodes(grid, Face())

lines!(axb, bn, zc)
lines!(axe, en, zc)
lines!(axκ, κcn, zf[2:Nz])

display(fig)

