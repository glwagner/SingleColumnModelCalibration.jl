using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using GLMakie
using JLD2

function run_simple_mixing!(axs, grid, closure;
                            τˣ = 1e-4,
                            Jᵇ = 1e-7,
                            N² = 1e-5,
                            Δt = 1minute,
                            stop_time = 24hours)

    tracers = (:b, :e)
    buoyancy = BuoyancyTracer()
    b_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵇ))
    u_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(τˣ))
    boundary_conditions = (b=b_bcs, u=u_bcs)

    model = HydrostaticFreeSurfaceModel(; grid, tracers, buoyancy, closure,
                                        boundary_conditions)

    bᵢ(z) = N² * z
    set!(model, b=bᵢ)

    simulation = Simulation(model; Δt, stop_time)
    run!(simulation)

    un = model.velocities.u
    bn = model.tracers.b
    en = model.tracers.e
    z = znodes(bn)

    axb, axu, axe = axs
    lines!(axb, interior(bn, 1, 1, :), z)
    lines!(axu, interior(un, 1, 1, :), z)
    lines!(axe, interior(en, 1, 1, :), z)

    return model
end

fig = Figure(size=(800, 400), title="0.91.1")
axb = Axis(fig[1, 1], xlabel="Buoyancy", ylabel="z (m)")
axu = Axis(fig[1, 2], xlabel="x-velocity", ylabel="z (m)")
axe = Axis(fig[1, 3], xlabel="TKE", ylabel="z (m)")
Label(fig[0, 1:3], "0.91.1")
axs = (axb, axu, axe)

grid = RectilinearGrid(size=128, z=(-256, 0), topology=(Flat, Flat, Bounded))

@load "optimal_catke.jld2" optimal_catke
closure = optimal_catke
#closure = CATKEVerticalDiffusivity()

run_simple_mixing!(axs, grid, closure, τˣ=0.0, Jᵇ=1e-6)
run_simple_mixing!(axs, grid, closure, τˣ=5e-4, Jᵇ=2e-7)
run_simple_mixing!(axs, grid, closure, τˣ=1e-3, Jᵇ=0.0)

display(fig)

