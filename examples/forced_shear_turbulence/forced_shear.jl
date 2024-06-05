using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity, ExplicitTimeDiscretization
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: MixingLength, TurbulentKineticEnergyEquation
using CairoMakie
using Printf

set_theme!(Theme(fontsize=22))

Δz = 5
Lz = 500 + Δz/2
z = -Lz:Δz:Lz
Nz = length(z) - 1

grid = RectilinearGrid(size=Nz, z=z, topology=(Flat, Flat, Bounded))

closure = CATKEVerticalDiffusivity()

# Ri = N² / (U/d)^2
# → U² = d² * N² / Ri
@inline U₀(p) = sqrt(p.N² * p.d^2 / p.Ri)
@inline u★(z, p) = U₀(p) * tanh(z / p.d)
@inline b★(z, p) = p.N² * z # p.δ * p.N² * tanh(z / p.δ)
@inline u_restoring(z, t, u, p) = (u★(z, p) - u) / p.τ

const c = Center()

@inline function b_restoring(i, j, k, grid, clock, fields, p)
    z = znode(i, j, k, grid, c, c, c)
    b = @inbounds fields.b[i, j, k]
    return (b★(z, p) - b) / p.τ
end

parameters = (Ri=0.1, τ=6hour, N²=1e-6, d=50.0, δ=80.0)
u_forcing = Forcing(u_restoring; field_dependencies=:u, parameters)
b_forcing = Forcing(b_restoring; discrete_form=true, parameters)

model = HydrostaticFreeSurfaceModel(; grid, closure,
                                    buoyancy = BuoyancyTracer(),
                                    tracers = (:b, :e),
                                    forcing = (; u=u_forcing, b=b_forcing))

include("../tracer_length_scale_operations.jl")         

function evolve_forced_shear!(model)
    model.clock.time = 0
    model.clock.iteration = 0

    bᵢ(z) = b★(z, parameters)
    uᵢ(z) = u★(z, parameters)
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
    κc = model.diffusivity_fields.κc
    κe = model.diffusivity_fields.κe
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
        push!(St,  deepcopy(interior(S, 1, 1, :)))
        push!(bt,  deepcopy(interior(b, 1, 1, :)))
        push!(ut,  deepcopy(interior(u, 1, 1, :)))
        push!(et,  deepcopy(interior(e, 1, 1, :)))
        push!(κct, deepcopy(interior(κc, 1, 1, :)))
        push!(κet, deepcopy(interior(κe, 1, 1, :)))
        push!(wbt, deepcopy(interior(wb, 1, 1, :)))
        push!(ϵt,  deepcopy(interior(ϵ, 1, 1, :)))
        push!(Pt,  deepcopy(interior(P, 1, 1, :)))

        return nothing
    end

    t = 0:simulation.Δt:simulation.stop_time
    simulation.callbacks[:dc] = Callback(collect_data, SpecifiedTimes(t))
    run!(simulation)

    timeseries = (; t, bt, ut, et, wbt, ϵt, Pt, N²t, St, Rit, κct, κet)
    return timeseries
end

fig = Figure(size=(1100, 900))

axR0 = Axis(fig[1, 1], xlabel="Ri", ylabel="z (m)", xticks=([0, 0.25, 0.5, 1.0], ["0", "0.25", "0.5", "1"]))
axu0 = Axis(fig[2, 1], xlabel="u (m s⁻¹)", ylabel="z (m)")
axN0 = Axis(fig[1, 2], xlabel="N² (s⁻²)", ylabel="z (m)", xticks=([0, 1e-6, 2e-6], ["0", "10⁻⁶", "2×10⁻⁶"]))
axS0 = Axis(fig[2, 2], xlabel="∂u/∂z (s⁻¹)", ylabel="z (m)")
axe0 = Axis(fig[3, 2], xlabel="e (m² s⁻²)", ylabel="z (m)", xticks=([0, 0.0001, 0.0002], ["0", "10⁻⁴", "2×10⁻⁴"]))

axN = Axis(fig[1, 3], xlabel="Time (hours)", ylabel="z (m)", xaxisposition=:top, yaxisposition=:right)
axS = Axis(fig[2, 3], xlabel="Time (hours)", ylabel="z (m)", xaxisposition=:top, yaxisposition=:right)
axe = Axis(fig[3, 3], xlabel="Time (hours)", ylabel="z (m)", yaxisposition=:right)

axb = Axis(fig[3, 1], xlabel="TKE budget terms \n (m² s⁻³)", ylabel="z (m)", xticks=([-1e-7, 0, 1e-7], ["-10⁻⁷", "0", "10⁻⁷"]))
axt = Axis(fig[4, 3], xlabel="Time (hours)", ylabel="Γ = - wb / ϵ", yaxisposition=:right, yticks=[0.14, 0.18, 0.22])

#hidexdecorations!(axe)
hidexdecorations!(axS)
hideydecorations!(axN0, grid=false)
hideydecorations!(axS0, grid=false)
hideydecorations!(axe0, grid=false)
hidespines!(axR0, :r, :t)
hidespines!(axN0, :l, :r, :t)
hidespines!(axS0, :l, :r, :t)
hidespines!(axu0, :r, :t)
hidespines!(axe0, :l, :r, :t)
hidespines!(axb, :r, :t)
hidespines!(axt, :l, :t)

zc = znodes(grid, Center())
zf = znodes(grid, Face())

t, bt, ut, et, wbt, ϵt, Pt, N²t, St, Rit, κct, κet = evolve_forced_shear!(model)
Nt = length(t)

colors = Makie.wong_colors()
lines!(axb, wbt[end], zc, color=colors[4], linewidth=2, label="Buoyancy flux")
lines!(axb, -ϵt[end], zc, color=colors[5], linewidth=2, label="Dissipation")
lines!(axb,  Pt[end], zc, color=colors[6], linewidth=2, label="Shear production")
Legend(fig[4, 1], axb)

# Mixing ratio:
Γ = [-sum(wbt[n]) / sum(ϵt[n]) for n = 1:Nt]
lines!(axt, t ./ hour, Γ)

for n = (1, 19, 145)
    @show tn = t[n] / hour
    label = @sprintf("t = %d hr", tn)
    lines!(axN0, N²t[n][2:Nz], zf[2:Nz]; linewidth=2, label)
    lines!(axS0, St[n][2:Nz], zf[2:Nz], linewidth=2)
    lines!(axu0, ut[n], zc, linewidth=2)
    lines!(axR0, Rit[n][2:Nz], zf[2:Nz], linewidth=2)
    lines!(axe0, et[n], zc, linewidth=2)
end

Legend(fig[4, 2], axN0, tellheight=false, tellwidth=false)

bzt = hcat(bt...)'
Szt = hcat(St...)'
Nzt = hcat(N²t...)'
ezt = hcat(et...)'

cr = heatmap!(axN, t ./ hour, zf, Nzt, colormap=:viridis, colorrange=(8e-7, 1.5e-6))
Colorbar(fig[1, 4], cr, label="N² (s⁻²)") #, vertical=false)

cr = contourf!(axS, t ./ hour, zf, Szt, colormap=:solar)#, colorrange=(0, 1e-3))
Colorbar(fig[2, 4], cr, label="∂u/∂z (s⁻¹)", tellwidth=false)

cr = contourf!(axe, t ./ hour, zc, ezt, colormap=:heat, colorrange=(0, 1e-3))
Colorbar(fig[3, 4], cr, label="e (m² s⁻²)", tellwidth=false)

for ax in (axN0, axu0, axS0, axR0, axe0, axN, axe, axS, axb)
    ylims!(ax, -150, 150)
end

xlims!(axN0, 0, 2e-6)

xlims!(axR0, 0, 1)
vlines!(axR0, 0.25, color=(:black, 0.5), linewidth=4)

for ax in (axN, axe, axt)
    xlims!(ax, 0, 24)
end

colsize!(fig.layout, 1, Relative(0.25))
colsize!(fig.layout, 2, Relative(0.25))
rowsize!(fig.layout, 4, Relative(0.15))

lb = (:left, :bottom)
rb = (:right, :bottom)

text!(axR0, 0.02, 0.03, space=:relative, align=lb, text="(a)", color=:gray)
text!(axN0, 0.02, 0.03, space=:relative, align=lb, text="(b)", color=:gray)
text!(axN,  0.02, 0.03, space=:relative, align=lb, text="(c)", color=:lightgray)
text!(axu0, 0.98, 0.03, space=:relative, align=rb, text="(d)", color=:gray)
text!(axS0, 0.98, 0.03, space=:relative, align=rb, text="(e)", color=:gray)
text!(axS,  0.02, 0.03, space=:relative, align=lb, text="(f)", color=:lightgray)

text!(axb,  0.02, 0.03, space=:relative, align=lb, text="(g)", color=:gray)
text!(axe0, 0.98, 0.03, space=:relative, align=rb, text="(h)", color=:gray)
text!(axe,  0.02, 0.03, space=:relative, align=lb, text="(i)", color=:gray)
text!(axt,  0.98, 0.03, space=:relative, align=rb, text="(j)", color=:gray)

display(fig)

save("forced_shear.pdf", fig)

