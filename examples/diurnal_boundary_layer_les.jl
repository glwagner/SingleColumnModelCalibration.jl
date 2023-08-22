using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode

using Printf
using Statistics

Lx = Ly = 256
Lz = 128
Nx = Ny = 256
Nz = 128

grid = RectilinearGrid(GPU(),
                       size = (Nx, Ny, Nz),
                       x = (0, Lx),
                       y = (0, Ly),
                       z = (-Lz, 0),
                       topology = (Periodic, Periodic, Bounded))

Qᵘ = -4e-5
const ω = 2π / 1day
Q₀ = 4e-7
ϕ₀ = π
@inline Qᵇ(x, y, t, p) = p.Q₀ * min(sin(p.ω * t + p.ϕ₀), 1/4)

top_b_bc = FluxBoundaryCondition(Qᵇ, parameters=(; ω, ϕ₀, Q₀))
b_bcs = FieldBoundaryConditions(top=top_b_bc)

top_u_bc = FluxBoundaryCondition(Qᵘ)
u_bcs = FieldBoundaryConditions(top=top_u_bc)

closure = AnisotropicMinimumDissipation()

model = NonhydrostaticModel(; grid, closure,
                            timestepper = :RungeKutta3,
                            buoyancy = BuoyancyTracer(),
                            advection = WENO(),
                            tracers = :b,
                            boundary_conditions = (; u=u_bcs, b=b_bcs))

N² = 1e-5
bᵢ(x, y, z) = N² * z + 1e-3 * rand() * N² * Lz
set!(model, b=bᵢ)
simulation = Simulation(model, Δt=1minutes, stop_time=4.5days)

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
                                                   schedule = TimeInterval(10minutes),
                                                   filename = "diurnal_boundary_layer_les_averages.jld2",
                                                   overwrite_existing = true)

simulation.output_writers[:fields] = JLD2OutputWriter(model, (; u, v, w, b),
                                                      schedule = TimeInterval(2hours),
                                                      filename = "diurnal_boundary_layer_les_fields.jld2",
                                                      overwrite_existing = true)

run!(simulation)

#using GLMakie

