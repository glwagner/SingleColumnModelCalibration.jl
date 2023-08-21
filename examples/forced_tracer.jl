using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode
using Oceananigans.Operators: Δzᶜᶜᶜ
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using GLMakie

Lz = 256
Nz = 32

grid = RectilinearGrid(size=Nz, z=(-Lz, 0), topology=(Flat, Flat, Bounded))

@inline function c_forcing_func(i, j, k, grid, clock, fields, parameters)
    k★ = parameters.k★
    Δ = parameters.Δ

    # Δ - δ * (Nz - k★) = 0
    δ = Δ / (Nz - k★)

    surface_layer = k > k★
    ϵ = ifelse(surface_layer, -δ, ifelse(k == k★, Δ, zero(grid)))

    return ϵ
end

Qᵇ = 1e-7
parameters = (k★ = Nz-8, Δ=1/Lz)
c_forcing = Forcing(c_forcing_func; discrete_form=true, parameters)
top_b_bc = FluxBoundaryCondition(Qᵇ)
b_bcs = FieldBoundaryConditions(top=top_b_bc)

closure = CATKEVerticalDiffusivity()

model = HydrostaticFreeSurfaceModel(; grid,
                                    buoyancy = BuoyancyTracer(),
                                    tracers = (:b, :c, :e),
                                    closure,
                                    forcing = (; c=c_forcing),
                                    boundary_conditions = (; b=b_bcs))

N² = 1e-5
bᵢ(x, y, z) = N² * z
set!(model, b=bᵢ, e=1e-6)
simulation = Simulation(model, Δt=10minutes, stop_time=7days)

bt = []
ct = []
et = []
κct = []

b = model.tracers.b
c = model.tracers.c
e = model.tracers.e
κc = model.diffusivity_fields.κᶜ

function collect_data(sim)
    push!(bt,  deepcopy(interior(b, 1, 1, :)))
    push!(ct,  deepcopy(interior(c, 1, 1, :)))
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
axc = Axis(fig[1, 2])
axe = Axis(fig[1, 3])
axκ = Axis(fig[1, 4])

xlims!(axb, -Lz * N², 0)
xlims!(axe, -1e-4, 1e-3)
xlims!(axκ, -1e-1, 1e1)

slider = Slider(fig[2, 1:4], range=1:Nt, startvalue=Nt)
n = slider.value

bn  = @lift bt[$n]
cn  = @lift ct[$n]
en  = @lift et[$n]
κcn = @lift κct[$n][2:Nz]

zc = znodes(grid, Center())
zf = znodes(grid, Face())

lines!(axb, bn, zc)
lines!(axc, cn, zc)
lines!(axe, en, zc)
lines!(axκ, κcn, zf[2:Nz])

xlims!(axc, -1, 1)

display(fig)

