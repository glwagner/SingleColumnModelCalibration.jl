using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode
using Printf
using Statistics

Lx = 256
Lz = 512
Nx = 256
Nz = 512

grid = RectilinearGrid(size = (Nx, Nz),
                       x = (0, Lx),
                       z = (-Lz/2, Lz/2),
                       topology=(Periodic, Flat, Bounded))

# Ri = N² / (U/d)^2
# → U² = d² * N² / Ri
#
@inline U₀(p) = sqrt(p.N² * p.d^2 / p.Ri)
@inline u★(z, p) = U₀(p) * tanh(z / p.d)
@inline b★(z, p) = p.N² * z # p.δ * p.N² * tanh(z / p.δ)
@inline u_restoring(x, y, z, t, u, p) = (u★(z, p) - u) / p.τ

const c = Center()

@inline function b_restoring(i, j, k, grid, clock, fields, p)
    z = znode(i, j, k, grid, c, c, c)
    b = @inbounds fields.b[i, j, k]
    return (b★(z, p) - b) / p.τ
end

parameters = (; Ri=0.1, τ = 6hour, N² = 1e-6, d=50.0)
@show U₀(parameters)

u_forcing = Forcing(u_restoring; field_dependencies=:u, parameters)
b_forcing = Forcing(b_restoring; discrete_form=true, parameters)

model = HydrostaticFreeSurfaceModel(; grid,
                                    advection = WENO(),
                                    buoyancy = BuoyancyTracer(),
                                    tracers = :b,
                                    forcing = (; u=u_forcing, b=b_forcing))


model.clock.time = 0
model.clock.iteration = 0

bᵢ(x, y, z) = b★(z, parameters) * (1 + 1e-2 * randn())
uᵢ(x, y, z) = u★(z, parameters)
set!(model, b=bᵢ, u=uᵢ)
simulation = Simulation(model, Δt=10.0, stop_time=1day)

wizard = TimeStepWizard(cfl=0.5, max_change=1.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

b = model.tracers.b
u, v, w = model.velocities
B = Average(b, dims=(1, 2))
U = Average(u, dims=(1, 2))

start_time = Ref(time_ns())
function progress(sim)
    elapsed = time_ns() - start_time[]
    msg = @sprintf("Iter: %d, time: %s, wall time: %s, max|w|: %.2e",
                   iteration(sim), prettytime(sim), prettytime(1e-9 * elapsed), maximum(abs, w))
    start_time[] = time_ns()
    @info msg
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

simulation.output_writers[:avg] = JLD2OutputWriter(model, (; U, B),
                                                   schedule = TimeInterval(5minutes),
                                                   filename = "forced_shear_les_averages.jld2",
                                                   overwrite_existing = true)

simulation.output_writers[:fields] = JLD2OutputWriter(model, (; u, v, w, b),
                                                      schedule = TimeInterval(30minutes),
                                                      filename = "forced_shear_les_fields.jld2",
                                                      overwrite_existing = true)

run!(simulation)

