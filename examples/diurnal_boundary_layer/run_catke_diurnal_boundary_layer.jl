using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity

#using CairoMakie
using GLMakie
using Printf
using Statistics

Δz = 4
Lz = 256
Nz = round(Int, Lz/Δz)
Jᵘ = -4e-5
ω = 2π / 1day
J₀ = 4e-7
ϕ₀ = π
N² = 1e-5

filename = "catke_diurnal_boundary_layer_dz$(Δz).jld2"

grid = RectilinearGrid(size=Nz, z=(-Lz, 0), topology=(Flat, Flat, Bounded))

@inline Jᵇ(t, p) = p.J₀ * min(sin(p.ω * t + p.ϕ₀), 1/4)
parameters = (; J₀, ω, ϕ₀)
top_b_bc = FluxBoundaryCondition(Jᵇ; parameters)
b_bcs = FieldBoundaryConditions(top=top_b_bc)

top_u_bc = FluxBoundaryCondition(Jᵘ)
u_bcs = FieldBoundaryConditions(top=top_u_bc)

closure = CATKEVerticalDiffusivity()

model = HydrostaticFreeSurfaceModel(; grid, closure,
                                    buoyancy = BuoyancyTracer(),
                                    tracers = (:b, :e),
                                    boundary_conditions = (; u=u_bcs, b=b_bcs))

bᵢ(z) = N² * z
set!(model, b=bᵢ, e=1e-6)
simulation = Simulation(model, Δt=1minutes, stop_time=4.5days)

include("tracer_length_scale_operations.jl")

bt = []
ut = []
et = []
ℓᶜsheart = []
ℓᶜconvt = []
ℓᶜtotalt = []
κct = []
N²t = []

b = model.tracers.b
e = model.tracers.e
u = model.velocities.u
κc = model.diffusivity_fields.κᶜ

N² = Field(∂z(b))
ℓᶜshear = tracer_stable_length_scale_operation(model)
ℓᶜconv = tracer_convective_length_scale_operation(model)
ℓᶜtotal = tracer_mixing_length_operation(model)

outputs = (; b, u, e, κc, N², ℓᶜshear, ℓᶜconv, ℓᶜtotal)
simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs; filename,
                                                    schedule = IterationInterval(1),
                                                    overwrite_existing = true)

#=
ℓᶜshear = Field(ℓᶜshear)
ℓᶜconv = Field(ℓᶜconv)
ℓᶜtotal = Field(ℓᶜtotal)

function collect_data(sim)
    compute!(ℓᶜshear)
    compute!(ℓᶜconv)
    compute!(ℓᶜtotal)
    compute!(N²)

    push!(bt,  deepcopy(interior(b, 1, 1, :)))
    push!(ut,  deepcopy(interior(u, 1, 1, :)))
    push!(et,  deepcopy(interior(e, 1, 1, :)))
    push!(κct, deepcopy(interior(κc, 1, 1, 2:Nz)))
    push!(N²t, deepcopy(interior(N², 1, 1, :)))

    push!(ℓᶜsheart, deepcopy(interior(ℓᶜshear, 1, 1, 2:Nz)))
    push!(ℓᶜconvt,  deepcopy(interior(ℓᶜconv,  1, 1, 2:Nz)))
    push!(ℓᶜtotalt, deepcopy(interior(ℓᶜtotal, 1, 1, 2:Nz)))

    return nothing
end

add_callback!(simulation, collect_data, IterationInterval(1))
=#

t = 0:simulation.Δt:simulation.stop_time
Nt = length(t)

run!(simulation)

