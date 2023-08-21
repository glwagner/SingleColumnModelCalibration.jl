#=
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity, ExplicitTimeDiscretization
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: MixingLength, TurbulentKineticEnergyEquation
using CairoMakie
using Printf
using Statistics

Δz = 5
Lz = 500 + Δz/2
z = -Lz:Δz:Lz
Nz = length(z) - 1

grid = RectilinearGrid(size=Nz, z=z, topology=(Flat, Flat, Bounded))
closure = CATKEVerticalDiffusivity()

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

include("tracer_length_scale_operations.jl")         

function shear_instability!(grid, Ri)
    parameters = (; Ri, τ = 6hour, N² = 1e-6, d=50.0)
    @show U₀(parameters)

    u_forcing = Forcing(u_restoring; field_dependencies=:u, parameters)
    b_forcing = Forcing(b_restoring; discrete_form=true, parameters)

    model = HydrostaticFreeSurfaceModel(; grid,
                                        buoyancy = BuoyancyTracer(),
                                        tracers = (:b, :e),
                                        forcing = (; u=u_forcing, b=b_forcing),
                                        closure)


    model.clock.time = 0
    model.clock.iteration = 0

    bᵢ(x, y, z) = b★(z, parameters)
    uᵢ(x, y, z) = u★(z, parameters)
    set!(model, b=bᵢ, u=uᵢ, e=1e-6)
    simulation = Simulation(model, Δt=1minute, stop_time=2day)

    bt = []
    ut = []
    et = []
    κct = []
    κet = []
    N²t = []
    Rit = []
    St = []
    wbt = []
    ϵt = []
    Pt = []

    b = model.tracers.b
    u = model.velocities.u
    e = model.tracers.e
    κc = model.diffusivity_fields.κᶜ
    κe = model.diffusivity_fields.κᵉ
    N² = Field(∂z(b))
    Ri = Field(∂z(b) / ∂z(u)^2)
    S = Field(∂z(u))

    etd = ExplicitTimeDiscretization()
    explicit_catke = CATKEVerticalDiffusivity(etd)
    wb_op = buoyancy_flux_operation(model, explicit_catke)
    ϵ_op = dissipation_operation(model, explicit_catke)
    P_op = shear_production_operation(model, explicit_catke)
    wb = Field(wb_op)
    ϵ = Field(ϵ_op)
    P = Field(P_op)

    function collect_data(sim)
        compute!(N²)
        compute!(Ri)
        compute!(S)
        compute!(wb)
        compute!(ϵ)
        compute!(P)

        push!(N²t, deepcopy(interior(N², 1, 1, :)))
        push!(Rit, deepcopy(interior(Ri, 1, 1, :)))
        push!(St, deepcopy(interior(S, 1, 1, :)))
        push!(bt,  deepcopy(interior(b, 1, 1, :)))
        push!(ut,  deepcopy(interior(u, 1, 1, :)))
        push!(et,  deepcopy(interior(e, 1, 1, :)))
        push!(κct, deepcopy(interior(κc, 1, 1, :)))
        push!(κet, deepcopy(interior(κe, 1, 1, :)))
        push!(wbt, deepcopy(interior(wb, 1, 1, :)))
        push!(ϵt, deepcopy(interior(ϵ, 1, 1, :)))
        push!(Pt, deepcopy(interior(P, 1, 1, :)))

        return nothing
    end

    t = 0:simulation.Δt:simulation.stop_time

    simulation.callbacks[:dc] = Callback(collect_data, SpecifiedTimes(t))

    run!(simulation)

    timeseries = (; t, bt, ut, et, wbt, ϵt, Pt, N²t, St, Rit, κct, κet)

    return timeseries
end

#####
##### Generate data
#####

timeserieses = []

Ris = 0.05:0.01:0.25

for Ri in Ris
    timeseries = shear_instability!(grid, Ri)
    push!(timeserieses, timeseries)
end

#####
##### Visualize results
#####

timeserieses = []

Ris = 0.05:0.01:0.25

for Ri in Ris
    timeseries = shear_instability!(grid, Ri)
    push!(timeserieses, timeseries)
end
=#

set_theme!(Theme(fontsize=24))
fig = Figure(resolution=(1200, 600))

labelpadding = 10

axe = Axis(fig[1, 1],
           xlabel = "Time (hours)",
           ylabel = "E(t) = ∫ e dz (m³ s⁻²)",
           ylabelpadding = labelpadding,
           xlabelpadding = labelpadding,
           xaxisposition = :top,
           yscale = log10)

axR = Axis(fig[2, 1], xlabel="Time (hours)", ylabel="min(Ri)",
           ylabelpadding = 20,
           xlabelpadding = 20,
           yticks = 0.0:0.05:0.25)

axm = Axis(fig[1, 2], xlabel="Ri⋆", ylabel="max(E) (m³ s⁻²)",
           ylabelpadding = labelpadding,
           xlabelpadding = labelpadding,
           yscale = log10,
           xaxisposition = :top,
           yaxisposition = :right,
           yticks = ([1e-3, 1e-2, 1e-1, 1e0], ["10⁻³", "10⁻²", "10⁻¹", "10⁰"]),
           ytrimspine = true)

ylabel = "Equilibrium \n Ri"

axf = Axis(fig[2, 2], xlabel="Ri⋆", ylabel=ylabel,
           yaxisposition = :right,
           ylabelpadding = labelpadding,
           xlabelpadding = labelpadding,
           yticks = 0.0:0.05:0.25)

colors = cgrad(:blues, rev=true, length(Ris), categorical=true)

for (n, (timeseries, Ri)) in enumerate(zip(timeserieses, Ris))
    t, bt, ut, et, wbt, ϵt, Pt, N²t, St, Rit, κct, κet = timeseries

    min_Ri = [minimum(Ri[2:Nz]) for Ri in Rit]
    int_e  = [Δz * sum(e) for e in et]

    lines!(axR, t ./ hour, min_Ri, label="Ri = $Ri", color=colors[n])
    lines!(axe, t ./ hour, int_e, label="Ri = $Ri", color=colors[n])

    scatter!(axf, Ri, min_Ri[end], color=colors[n])
    scatter!(axm, Ri, maximum(int_e), color=colors[n])
end

guide = 0.18:0.01:0.25
lines!(axf, guide, guide, color=:black)

xlims!(axR, 0, 24)
xlims!(axe, 0, 24)

ylims!(axR, 0.05, 0.25)
ylims!(axf, 0.05, 0.25)

ylims!(axR, 0.05, 0.25)
ylims!(axf, 0.05, 0.25)

ylims!(axe, 1e-4, 1e0)
ylims!(axm, 1e-4, 1e0)

hidespines!(axe, :b, :r)
hidespines!(axR, :t, :r)
hidespines!(axm, :l, :b)
hidespines!(axf, :t, :l)

#hidexdecorations!(axf)
#axm.xgridvisible[] = false

rowgap!(fig.layout, 1, 30)

text!(axe, 0.98, 0.98, text="(a)", color=:gray, align=(:right, :top), space=:relative) 
text!(axR, 0.98, 0.03, text="(b)", color=:gray, align=(:right, :bottom), space=:relative) 
text!(axm, 0.98, 0.98, text="(c)", color=:gray, align=(:right, :top), space=:relative) 
text!(axf, 0.98, 0.03, text="(d)", color=:gray, align=(:right, :bottom), space=:relative) 

display(fig)

save("forced_stratified_turbulence_Ri_study.pdf", fig)

