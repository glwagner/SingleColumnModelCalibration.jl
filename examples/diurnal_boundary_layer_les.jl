using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode

using Printf
using Statistics

Lx = Ly = 256
Lz = 256
Nx = Ny = 128
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
                            advection = WENO(order=9),
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

#=
using GLMakie
using Oceananigans.Units

filename = "diurnal_boundary_layer_les_averages.jld2"
Bt = FieldTimeSeries(filename, "B")

t = Bt.times
Nt = length(t)
x, y, z = nodes(Bt)
Nx, Ny, Nz = size(Bt.grid)

# Compute mixed layer depth
h = zeros(Nt)
Δb = 1e-4
for n = 1:Nt
    # Surface buoyancy
    b₀ = Bt[1, 1, Nz, n]
    for k = Nz-1:-1:1
        bk = Bt[1, 1, k, n]
        # b is _decreasing_ downward
        if (b₀ - bk) >= Δb || k == 1
            h[n] = - z[k]
            break
        end
    end
end

fig = Figure()
axh = Axis(fig[1, 1])
axb0 = Axis(fig[2, 1])
axb = Axis(fig[3, 1])
lines!(axh, t ./ hours, h)
lines!(axb0, t ./ hours, interior(Bt, 1, 1, Nz, :))
contourf!(axb, t ./ hours, z, interior(Bt, 1, 1, :, :)')
ylims!(ax, -65, 0)
display(fig)
=#
