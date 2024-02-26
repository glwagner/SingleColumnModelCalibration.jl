using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    CATKEVerticalDiffusivity,
    MixingLength

#using CairoMakie
using GLMakie
using Printf
using Statistics

Lz = 80
Δz = 4
Nz = round(Int, Lz/Δz)

grid = RectilinearGrid(size=Nz, z=(-Lz, 0), topology=(Flat, Flat, Bounded))

Qᵘ = -4e-5
const ω = 2π / 1day
Q₀ = 4e-7
ϕ₀ = π
Qᵇ(x, y, t) = Q₀ * min(sin(ω * t + ϕ₀), 1/4)

top_b_bc = FluxBoundaryCondition(Qᵇ)
b_bcs = FieldBoundaryConditions(top=top_b_bc)

top_u_bc = FluxBoundaryCondition(Qᵘ)
u_bcs = FieldBoundaryConditions(top=top_u_bc)

closure = CATKEVerticalDiffusivity()

model = HydrostaticFreeSurfaceModel(; grid, closure,
                                    buoyancy = BuoyancyTracer(),
                                    tracers = (:b, :e),
                                    boundary_conditions = (; u=u_bcs, b=b_bcs))

N² = 1e-5
h₀ = 10 # initial mixed layer depth
#bᵢ(x, y, z) = min(N² * h₀, N² * z)
bᵢ(x, y, z) = N² * z
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
ℓᶜshear_op = tracer_stable_length_scale_operation(model)
ℓᶜconv_op = tracer_convective_length_scale_operation(model)
ℓᶜtotal_op = tracer_mixing_length_operation(model)

ℓᶜshear = Field(ℓᶜshear_op)
ℓᶜconv = Field(ℓᶜconv_op)
ℓᶜtotal = Field(ℓᶜtotal_op)

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

simulation.callbacks[:dc] = Callback(collect_data, IterationInterval(1))

t = 0:simulation.Δt:simulation.stop_time
Nt = length(t)

run!(simulation)

#####
##### Diurnal boundary layer
#####

set_theme!(Theme(fontsize=22))

#####
##### Validation against LES
#####

fig = Figure(resolution=(1400, 1200))

#####
##### Diagnostics plot
#####

#=
fig = Figure(resolution=(1400, 1200))

yticks = [-60, -30, 0]
axb0 = Axis(fig[1, 2]; xlabel="b (m s⁻²)", ylabel="z (m)", yticks, xticks=([-5e-4, 2e-3], ["-5×10⁻⁴", "2×10⁻³"]))
axe0 = Axis(fig[2, 2]; xlabel="e (m² s⁻²)", ylabel="z (m)", yticks, xticks=([0.0, 2e-4], ["0", "2×10⁻⁴"]))
axκ0 = Axis(fig[3, 2]; xlabel="log₁₀(κ)", ylabel="z (m)", yticks, xticks=[-4, -2, 0])
axℓ0 = Axis(fig[4, 2]; xlabel="ℓᶜ (m)", ylabel="z (m)", xscale=log10, yticks, xticks=([0.1, 1, 20], ["0.1", "1", "20"]))

axN = Axis(fig[1, 3]; xlabel="Time (hr)", ylabel="z (m)", xaxisposition=:top, yaxisposition=:right, yticks)
axe = Axis(fig[2, 3]; ylabel="z (m)", yaxisposition=:right, yticks)
axκ = Axis(fig[3, 3]; xlabel="Time (hours)", ylabel="z (m)", yaxisposition=:right, yticks)
axℓ = Axis(fig[4, 3]; xlabel="Time (hours)", ylabel="z (m)", yaxisposition=:right, yticks)

axℓt = Axis(fig[5, 3], xlabel="Time (hours)", yscale=log10, ylabel="max ℓᶜ \n (m)", yaxisposition=:right,
            yticks=([0.1, 1, 20], ["0.1", "1", "20"]))

axQ1 = Axis(fig[0, 3], xlabel="Time (hours)", ylabel="Qᵇ \n (m² s⁻³)",
            yticks=([-4e-7, 0, 1e-7], ["-4×10⁻⁷", "\n 0", "10⁻⁷"]),
            xticks=0:24:(6*24),
            xaxisposition=:top, yaxisposition=:right)

lines!(axQ1, t ./ hour, Qᵇ.(0, 0, t))

ylims!(axℓt, 7e-2, 100)

for ax in (axN, axe, axκ, axℓ, axℓt, axQ1)
    xlims!(ax, 0, simulation.stop_time/hour)
end

hidexdecorations!(axN)
hidexdecorations!(axe)
hidexdecorations!(axκ)
hidexdecorations!(axℓ)

hidespines!(axQ1, :l, :b)
hidespines!(axb0, :r, :t)
hidespines!(axe0, :r, :t)
hidespines!(axκ0, :r, :t)
hidespines!(axℓ0, :r, :t)
hidespines!(axℓt, :l, :t)

bzt = hcat(bt...)'
Nzt = hcat(N²t...)'
ezt = hcat(et...)'
κzt = hcat(κct...)'
ℓzt = hcat(ℓᶜtotalt...)'

zc = znodes(grid, Center())
zf = znodes(grid, Face())

#cr = contourf!(axN, t ./ hour, zf, Nzt, levels=8, colormap=:viridis)
#Colorbar(fig[1, 4], cr, label="N² (s⁻²)", ticks=([0, 3e-5, 6e-5], ["0", "3×10⁻⁵", "6×10⁻⁵"]))
cr = contourf!(axN, t ./ hour, zc, bzt, levels=12, colormap=:viridis)
Colorbar(fig[1, 4], cr, label="b (m s⁻²)") #, ticks=([0, 3e-5, 6e-5], ["0", "3×10⁻⁵", "6×10⁻⁵"]))

elevels = 2e-6:2e-5:2e-4
cr = contourf!(axe, t ./ hour, zc, ezt, levels=elevels, extendhigh=:auto, colormap=:heat)
Colorbar(fig[2, 4], cr, label="e (m² s⁻²)")

cr = contourf!(axκ, t ./ hour, zf[2:Nz], log10.(κzt .+ 1e-6), levels=6, colorrange=(-6, 0), colormap=:solar)
Colorbar(fig[3, 4], cr, label="log10(κᶜ)", ticks=[-6, -4, -2, 0])

cr = contourf!(axℓ, t ./ hour, zf[2:Nz], clamp.(ℓzt, 0, 15), levels=[0.25, 0.5, 1, 2, 4, 8, 16], colormap=:blues)
Colorbar(fig[4, 4], cr, label="ℓ (m)", ticks=[1, 4, 16])

cr = lines!(axℓt, t ./ hour, [maximum(ℓ .+ 1e-2) for ℓ in ℓᶜsheart], linewidth=3, color=(:purple, 0.6), label="ℓᶜ shear")
cr = lines!(axℓt, t ./ hour, [maximum(ℓ .+ 1e-2) for ℓ in ℓᶜconvt],  linewidth=3, color=(:black,  0.6), label="ℓᶜ convective")
Legend(fig[5, 2], axℓt, tellwidth=false, framevisible=false)

Nday  = round(Int, 24hour / simulation.Δt)
Nhour = round(Int, hour / simulation.Δt)

for n = (18Nhour+1, 56Nhour+1, 94Nhour+1)
    @show tn = t[n] / hour
    label = @sprintf("t = %d hr", tn)
    ln = lines!(axb0, bt[n], zc; label, linewidth=3)
    lines!(axe0, et[n], zc, linewidth=3)
    lines!(axκ0, log10.(κct[n] .+ 1e-6), zf[2:Nz], color = ln.color.val, label="κᶜ", linewidth=3)
    lines!(axℓ0, ℓᶜtotalt[n], zf[2:Nz], color = ln.color.val, linewidth=3)

    vlines!(axN, t[n] / hour,  color=(ln.color.val, 0.6), linewidth=4)
    vlines!(axe, t[n] / hour,  color=(ln.color.val, 0.6), linewidth=4)
    vlines!(axκ, t[n] / hour,  color=(ln.color.val, 0.6), linewidth=4)
    vlines!(axℓ, t[n] / hour,  color=(ln.color.val, 0.6), linewidth=4)
    vlines!(axQ1, t[n] / hour,  color=(ln.color.val, 0.6), linewidth=4)
    vlines!(axℓt, t[n] / hour, color=(ln.color.val, 0.6), linewidth=4)
end

for ax in (axb0, axe0, axκ0, axℓ0, axN, axe, axκ, axℓ)
    ylims!(ax, -65, 0)
end

xlims!(axb0, -7e-4, 2e-3)
xlims!(axℓ0, 1e-2, 40)

Legend(fig[0, 2], axb0, framevisible=false)

rowsize!(fig.layout, 0, Relative(0.1))
colsize!(fig.layout, 1, Relative(0.1))
colsize!(fig.layout, 2, Relative(0.1))
colgap!(fig.layout, 3, Relative(0.0))
rowsize!(fig.layout, 5, Relative(0.1))

ϵ = 0.05
δ = 0.01
lb = (:left, :bottom)
lt = (:left, :top)
rb = (:right, :bottom)
space = :relative
text!(axQ1, ϵ/2, 1-ϵ; space, align=lt, text="(a)", color=:gray)
text!(axb0, 1-ϵ,   ϵ; space, align=rb, text="(b)", color=:gray)
text!(axN,    δ,   ϵ; space, align=lb, text="(c)", color=:lightgray)
text!(axe0, 1-ϵ,   ϵ; space, align=rb, text="(d)", color=:gray)
text!(axe,    δ,   ϵ; space, align=lb, text="(e)", color=:gray)
text!(axκ0, 1-ϵ,   ϵ; space, align=rb, text="(f)", color=:gray)
text!(axκ,    δ,   ϵ; space, align=lb, text="(g)", color=:lightgray)
text!(axℓ0, 1-ϵ,   ϵ; space, align=rb, text="(h)", color=:gray)
text!(axℓ,    δ,   ϵ; space, align=lb, text="(i)", color=:gray)
text!(axℓt,   δ, 1-ϵ; space, align=lt, text="(j)", color=:gray)

display(fig)

save("diurnal_boundary_layer.pdf", fig)
=#

#####
##### Diurnal structure
#####

#=
fig = Figure(resolution=(1200, 400))

color1 = Makie.wong_colors()[1]
color2 = Makie.wong_colors()[2]
color3 = Makie.wong_colors()[3]

kws = []
for c in (color1, color2, color3)
    kw = Dict()
    for p in (:leftspinecolor, :rightspinecolor, :ytickcolor, :yticklabelcolor, :ylabelcolor) 
        kw[p] = c
    end
    kw[:ytrimspine] = true
    kw[:xlabel] = "Time (hours)"
    push!(kws, kw)
end

axq = Axis(fig[1, 1]; ylabel="Qᵇ (m² s⁻³)", yticks=([0, 1e-7], ["0", "10⁻⁷"]), ylabelpadding=20, xaxisposition=:top, kws[1]...)
axe = Axis(fig[2, 1]; ylabel="⟨e⟩ (m² s⁻²)", yticks=([0, 1e-4], ["0", "10⁻⁴"]), yaxisposition=:right, kws[2]...)
axκ = Axis(fig[3, 1]; ylabel="⟨κᶜ⟩ (m² s⁻¹)", yticks=[0.0, 0.1], xaxisposition=:bottom, kws[3]...)

e_mean = [mean(e) for e in et]
κ_mean = [mean(κᶜ) for κᶜ in κct]

κ_max = maximum(κ_mean)

lines!(axq, t ./ hour, Qᵇ.(0, 0, t), linewidth=4, color=color1)
lines!(axe, t ./ hour, e_mean,       linewidth=4, color=color2)
lines!(axκ, t ./ hour, κ_mean,       linewidth=4, color=color3)

hidespines!(axq, :b, :r)
hidespines!(axe, :l, :t, :b)
hidespines!(axκ, :t, :r)

hideydecorations!(axq, label=false, ticklabels=false, ticks=false)
hideydecorations!(axe, label=false, ticklabels=false, ticks=false)
hideydecorations!(axκ, label=false, ticklabels=false, ticks=false)

hidexdecorations!(axq, label=false, ticklabels=false, ticks=false)
hidexdecorations!(axe)
hidexdecorations!(axκ, label=false, ticklabels=false, ticks=false)

#lines!(axt, t ./ hour, e_mean ./ maximum(e_mean))
#lines!(axt, t ./ hour, κ_mean ./ maximum(κ_mean))
#lines!(axt, t ./ hour, max.(0.0, Qᵇ.(0, 0, t) ./ 1e-7))

axb = Axis(fig[3, 1])
hidespines!(axb)
hidedecorations!(axb)

Ndays = 4

for n in 1:Ndays
    @show t_hr = 24 * (n-1) .+ [12, 24]
    upper = 0.6 .* [κ_max, κ_max]
    lower = [0.0, 0.0] 
    band!(axκ, t_hr, lower, upper, color=(:black, 0.1))

    upper = +1e-7 .* [1, 1]
    lower = -1e-7 .* [1, 1]
    band!(axq, t_hr, lower, upper, color=(:black, 0.1))
end

ylims!(axe, 0, 1e-4)

for ax in (axq, axe, axκ, axb)
    xlims!(ax, 50, 100) #-6, 162)
end

rowgap!(fig.layout, 1, -100)
rowgap!(fig.layout, 2, -120)

display(fig)

save("diurnal_mixing_length_memory.pdf", fig)

=#
