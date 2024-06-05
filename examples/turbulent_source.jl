using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity, ExplicitTimeDiscretization
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: MixingLength, TurbulentKineticEnergyEquation
using CairoMakie
using Printf

Δz = 5
Lz = 500 + Δz/2
z = -Lz:Δz:Lz
Nz = length(z) - 1

grid = RectilinearGrid(size=Nz, z=z, topology=(Flat, Flat, Bounded))
closure = CATKEVerticalDiffusivity()

@inline tke_source(z, t, p) = p.q * exp(-z^2 / (2 * p.d^2)) * exp(-(t - p.τ)^2 / (2 * p.τ^2))
parameters = (q = 1e-2 / hour, d = 16.0, τ = 6hours)
e_forcing = Forcing(tke_source; parameters)

model = HydrostaticFreeSurfaceModel(; grid,
                                    buoyancy = BuoyancyTracer(),
                                    tracers = (:b, :e),
                                    forcing = (; e=e_forcing),
                                    closure)

include("tracer_length_scale_operations.jl")         

function evolve_turbulent_source!(model, N²ᵢ = 1e-5)
    model.clock.time = 0
    model.clock.iteration = 0

    δ = 20
    bᵢ(z) = N²ᵢ * z
    set!(model, b=bᵢ, e=1e-6)
    simulation = Simulation(model, Δt=10minute, stop_time=1day)

    bt = []
    et = []
    κct = []
    κet = []
    N²t = []
    wbt = []
    ϵt = []

    b = model.tracers.b
    u = model.velocities.u
    e = model.tracers.e
    κc = model.diffusivity_fields.κc
    κe = model.diffusivity_fields.κe
    N² = Field(∂z(b))

    etd = ExplicitTimeDiscretization()
    #explicit_catke = CATKEVerticalDiffusivity(etd; mixing_length, turbulent_kinetic_energy_equation)
    explicit_catke = CATKEVerticalDiffusivity(etd)
    wb_op = buoyancy_flux_operation(model, explicit_catke)
    ϵ_op = dissipation_operation(model, explicit_catke)
    wb = Field(wb_op)
    ϵ = Field(ϵ_op)

    function collect_data(sim)
        compute!(N²)
        compute!(wb)
        compute!(ϵ)

        push!(N²t, deepcopy(interior(N², 1, 1, :)))
        push!(bt,  deepcopy(interior(b, 1, 1, :)))
        push!(et,  deepcopy(interior(e, 1, 1, :)))
        push!(κct, deepcopy(interior(κc, 1, 1, :)))
        push!(κet, deepcopy(interior(κe, 1, 1, :)))
        push!(wbt, deepcopy(interior(wb, 1, 1, :)))
        push!(ϵt, deepcopy(interior(ϵ, 1, 1, :)))

        return nothing
    end

    t = 0:simulation.Δt:simulation.stop_time

    simulation.callbacks[:dc] = Callback(collect_data, SpecifiedTimes(t))

    run!(simulation)

    timeseries = (; t, bt, wbt, ϵt, N²t, et, κct, κet)

    return timeseries
end

set_theme!(Theme(fontsize=24))

fig = Figure(size=(1200, 800))

axN0 = Axis(fig[1:2, 2], xlabel="N² (s⁻²)", ylabel="z (m)", xticks=([1e-6, 1e-5], ["10⁻⁶", "10⁻⁵"]))
axe0 = Axis(fig[3:4, 2], xlabel="e (m² s⁻²)", ylabel="z (m)", xticks=[0, 0.001])

axN = Axis(fig[1:2, 3], xlabel="Time (hours)", ylabel="z (m)", xaxisposition=:top, yaxisposition=:right)
axe = Axis(fig[3:4, 3], xlabel="Time (hours)", ylabel="z (m)", yaxisposition=:right)

axb = Axis(fig[3:4, 1], xlabel="TKE budget terms (m² s⁻³)", ylabel="z (m)", xticks=([-1e-6, 0, 1e-6], ["-10⁻⁶", "0", "10⁻⁶"]))

#hidexdecorations!(axe)
hideydecorations!(axe0)
hidespines!(axN0, :r, :t)
hidespines!(axe0, :l, :r, :t)
hidespines!(axb, :r, :t)

zc = znodes(grid, Center())
zf = znodes(grid, Face())

t, bt, wbt, ϵt, N²t, et, κct, κet = evolve_turbulent_source!(model, 1e-6)
Nt = length(t)

# For diagnostics
tke_source(z) = tke_source(z, t[73], parameters)
P = CenterField(grid)
set!(P, tke_source)
Pi = interior(P, 1, 1, :)

colors = Makie.wong_colors()
lines!(axb, wbt[37], zc, color=colors[4], linewidth=3, label="Buoyancy flux")
lines!(axb, -ϵt[37], zc, color=colors[5], linewidth=3, label="Dissipation")
lines!(axb, Pi,       zc, color=colors[6], linewidth=3, label="Production")
Legend(fig[1, 1], axb, framevisible=false)

for n = (4, 37, 145)
    @show tn = t[n] / hour
    label = if tn < 1
        @sprintf("t = %.1f hr", tn)
    else
        @sprintf("t = %d hr", tn)
    end

    lines!(axN0, N²t[n][2:Nz], zf[2:Nz]; linewidth=3, label)
    ln = lines!(axe0, et[n], zc, linewidth=3)
end

Legend(fig[2, 1], axN0, tellheight=false, tellwidth=false)

bzt = hcat(bt...)'
Nzt = hcat(N²t...)'
ezt = hcat(et...)'

cr = heatmap!(axN, t ./ hour, zf, Nzt, colormap=:viridis, colorrange=(8e-7, 1.5e-6))
Colorbar(fig[1:2, 4], cr, label="N² (s⁻²)") #, vertical=false)

cr = contourf!(axe, t ./ hour, zc, ezt, colormap=:heat, colorrange=(0, 1e-3))
Colorbar(fig[3:4, 4], cr, label="e (m² s⁻²)", tellwidth=false)

for ax in (axN0, axe0, axN, axe, axb)
    ylims!(ax, -250, 250)
end

for ax in (axN, axe)
    xlims!(ax, 0, 1 * 24)
end

xlims!(axN0, -2e-7, 1e-5)
xlims!(axb, -1e-6, 2e-6)

#colgap!(fig.layout, 1, -20)
colsize!(fig.layout, 1, Relative(0.2))
colsize!(fig.layout, 2, Relative(0.2))
#rowsize!(fig.layout, 3, Relative(0.2))

text!(axN0, 0.98, 0.02, align=(:right, :bottom), space=:relative, text="(a)", color=:gray)
text!(axN,  0.02, 0.02, align=(:left, :bottom), space=:relative, text="(b)", color=:lightgray)
text!(axb,  0.98, 0.02, align=(:right, :bottom), space=:relative, text="(c)", color=:gray)
text!(axe0, 0.98, 0.02, align=(:right, :bottom), space=:relative, text="(d)", color=:gray)
text!(axe,  0.02, 0.02, align=(:left, :bottom), space=:relative, text="(e)", color=:gray)

display(fig)

save("turbulent_source.pdf", fig)

