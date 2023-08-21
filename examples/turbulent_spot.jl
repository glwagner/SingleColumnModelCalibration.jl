using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: MixingLength
using GLMakie
using Printf

Δz = 5
Lz = 500 + Δz/2
z = -Lz:Δz:Lz
Nz = length(z) - 1

grid = RectilinearGrid(size=Nz, z=z, topology=(Flat, Flat, Bounded))

closure = CATKEVerticalDiffusivity()

model = HydrostaticFreeSurfaceModel(; grid,
                                    buoyancy = BuoyancyTracer(),
                                    tracers = (:b, :e),
                                    closure)

function evolve_turbulent_spot!(model, N²ᵢ = 1e-5)
    model.clock.time = 0
    model.clock.iteration = 0

    δ = 20
    bᵢ(x, y, z) = N²ᵢ * z
    #eᵢ(x, y, z) = 1e-6 + 1e-2 * ifelse(z > 0, 0, 1)
    eᵢ(x, y, z) = 1e-6 + 1e-2 * exp(-z^2 / 2δ^2)
    set!(model, b=bᵢ, e=eᵢ)
    simulation = Simulation(model, Δt=1minute, stop_time=1hour)

    bt = []
    et = []
    κct = []
    κet = []
    N²t = []

    b = model.tracers.b
    u = model.velocities.u
    e = model.tracers.e
    κc = model.diffusivity_fields.κᶜ
    κe = model.diffusivity_fields.κᵉ
    N² = Field(∂z(b))

    function collect_data(sim)
        compute!(N²)
        push!(N²t, deepcopy(interior(N², 1, 1, :)))
        push!(bt,  deepcopy(interior(b, 1, 1, :)))
        push!(et,  deepcopy(interior(e, 1, 1, :)))
        push!(κct, deepcopy(interior(κc, 1, 1, :)))
        push!(κet, deepcopy(interior(κe, 1, 1, :)))

        return nothing
    end

    t = 0:simulation.Δt:simulation.stop_time
    Nt = length(t)

    simulation.callbacks[:dc] = Callback(collect_data, SpecifiedTimes(t))

    run!(simulation)

    timeseries = (; bt, N²t, et, κct, κet)

    return timeseries
end

fig = Figure(resolution=(1200, 600))

axN0 = Axis(fig[1, 2], xlabel="N² (s⁻²)", ylabel="z (m)")
axe0 = Axis(fig[2, 2], xlabel="E (m² s⁻²)", ylabel="z (m)")

axN = Axis(fig[1, 3], xlabel="Time (minute)", ylabel="z (m)", xaxisposition=:top, yaxisposition=:right)
axe = Axis(fig[2, 3], ylabel="z (m)", yaxisposition=:right)

hidexdecorations!(axe)
hidespines!(axN0, :r, :t)
hidespines!(axe0, :r, :t)

zc = znodes(grid, Center())
zf = znodes(grid, Face())

bt, N²t, et, κct, κet = evolve_turbulent_spot!(model, 1e-6)

for n = (1, 11, 61)
    @show tn = t[n] / minute
    label = @sprintf("t = %d min", tn)
    lines!(axN0, N²t[n][2:Nz], zf[2:Nz]; label)
    ln = lines!(axe0, et[n], zc)
end

Legend(fig[2, 1], axN0)

bzt = hcat(bt...)'
Nzt = hcat(N²t...)'
ezt = hcat(et...)'

cr = heatmap!(axN, t ./ minute, zf, Nzt, colormap=:viridis, colorrange=(8e-7, 1.5e-6))
Colorbar(fig[1, 4], cr, label="N²")

cr = contourf!(axe, t ./ minute, zc, ezt, colormap=:heat, colorrange=(0, 1e-3))
Colorbar(fig[2, 4], cr, label="E (m² s⁻²)")

for ax in (axN0, axe0, axN, axe)
    ylims!(ax, -250, 250)
end

for ax in (axN, axe)
    xlims!(ax, 0, 60)
end


display(fig)

# save("turbulent_spot.png", fig)
