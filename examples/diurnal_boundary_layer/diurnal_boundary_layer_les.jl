using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode

using Printf
using Statistics
using CUDA

using LESbrary.IdealizedExperiments: ConstantFluxStokesDrift

Lx = Ly = 256
Lz = 256
Nx = Ny = 256
Nz = 256

prefix = "wave_averaged_diurnal_boundary_layer_les_Nx$(Nx)_Nz$(Nz)"

grid = RectilinearGrid(GPU(),
                       size = (Nx, Ny, Nz),
                       halo = (5, 5, 5),
                       x = (0, Lx),
                       y = (0, Ly),
                       z = (-Lz, 0),
                       topology = (Periodic, Periodic, Bounded))

Jᵘ = -1e-4
const ω = 2π / 1day
J₀ = 4e-7
ϕ₀ = π
@inline Jᵇ(x, y, t, p) = p.J₀ * min(sin(p.ω * t + p.ϕ₀), 1/4)

top_b_bc = FluxBoundaryCondition(Jᵇ, parameters=(; ω, ϕ₀, J₀))
b_bcs = FieldBoundaryConditions(top=top_b_bc)

top_u_bc = FluxBoundaryCondition(Jᵘ)
u_bcs = FieldBoundaryConditions(top=top_u_bc)

g = 9.81
kᵖ = 1e-6 * g / abs(Jᵘ)
stokes_drift = ConstantFluxStokesDrift(grid, Jᵘ, kᵖ)
#closure = AnisotropicMinimumDissipation()

@info "Whipping up the Stokes drift..."

uˢ₀ = CUDA.@allowscalar stokes_drift.uˢ[1, 1, grid.Nz]
a★ = stokes_drift.air_friction_velocity
ρʷ = stokes_drift.water_density
ρᵃ = stokes_drift.air_density
u★ = a★ * sqrt(ρᵃ / ρʷ)
La = sqrt(u★ / uˢ₀)
@info @sprintf("Air u★: %.4f, water u★: %.4f, λᵖ: %.4f, La: %.3f, Surface Stokes drift: %.4f m s⁻¹",
                a★, u★, 2π/kᵖ, La, uˢ₀)

model = NonhydrostaticModel(; grid, stokes_drift, # closure,
                            timestepper = :RungeKutta3,
                            buoyancy = BuoyancyTracer(),
                            advection = WENO(order=9),
                            tracers = :b,
                            boundary_conditions = (; u=u_bcs, b=b_bcs))

N² = 1e-5
bᵢ(x, y, z) = N² * z + 1e-3 * rand() * N² * Lz
set!(model, b=bᵢ)
simulation = Simulation(model, Δt=1.0, stop_time=4.5days)

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
                   iteration(sim), prettytime(sim), prettytime(1e-9 * elapsed), maximum(abs, interior(w)))
    start_time[] = time_ns()
    @info msg
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

simulation.output_writers[:avg] = JLD2OutputWriter(model, (; U, B),
                                                   schedule = TimeInterval(10minutes),
                                                   filename = prefix * "_averages.jld2",
                                                   overwrite_existing = true)

simulation.output_writers[:fields] = JLD2OutputWriter(model, (; u, v, w, b),
                                                      schedule = TimeInterval(2hours),
                                                      filename = prefix * "_fields.jld2",
                                                      overwrite_existing = true)

run!(simulation)

