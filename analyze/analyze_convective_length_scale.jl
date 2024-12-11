using Oceananigans
using JLD2
using CairoMakie
using ParameterEstimocean
using Oceananigans.Units
using MathTeXEngine

using ParameterEstimocean.Parameters: build_parameters_named_tuple

using SingleColumnModelCalibration:
    dependent_parameter_sets,
    build_batched_inverse_problem,
    prior_library

using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity, tracer_mixing_lengthᶜᶜᶠ

fonts = (; regular=texfont())
set_theme!(Theme(fontsize=32, linewidth=4; fonts))

n = 70
ulims = (-0.15, 0.35)
elims = (-1e-4, 4e-4)
linewidth = 12
prefix = "free_convection"

fig = Figure(size=(1600, 600))

xticks = ([0, 1e-4, 2e-4, 3e-4], ["0", "1", "2", "3"])
ax_b = Axis(fig[1, 1]; xlabel="Buoyancy \n (10⁻⁴ × m s⁻²)", ylabel="z (m)", xticks)
ax_e = Axis(fig[1, 2]; xlabel="Turbulent kinetic \n energy (m² s⁻²)", xticks=[0, 1e-3, 2e-3]) #, yaxisposition=:right)
ax_l = Axis(fig[1, 3]; xlabel="Tracer \n mixing length (m)")
ax_κ = Axis(fig[1, 4]; xlabel="Eddy \n diffusivity (m² s⁻¹)", ylabel="z (m)", yaxisposition=:right)

xlims!(ax_b, 1e-4, 3.4e-4)
xlims!(ax_κ, -0.1, 8.8)

ylims!(ax_b, -192, 0)
ylims!(ax_e, -192, 0)
ylims!(ax_l, -192, 0)
ylims!(ax_κ, -192, 0)

text!(ax_b, 0.03, 0.98, text="(a)", align=(:left, :top), space=:relative)
text!(ax_e, 0.03, 0.98, text="(b)", align=(:left, :top), space=:relative)
text!(ax_l, 0.03, 0.98, text="(c)", align=(:left, :top), space=:relative)
text!(ax_κ, 0.03, 0.98, text="(d)", align=(:left, :top), space=:relative)

hidespines!(ax_b, :t, :r)
hidespines!(ax_e, :t, :l, :r)
hidespines!(ax_l, :t, :l, :r)
hidespines!(ax_κ, :t, :l)

hideydecorations!(ax_e, grid=false)
hideydecorations!(ax_l, grid=false)

suite_parameters = [
    (name = "6_hour_suite",  resolution="1m", stop_time=6hours),
    (name = "24_hour_suite", resolution="1m", stop_time=24hours),
    (name = "72_hour_suite", resolution="1m", stop_time=66hours),
]   

suite_names = [L"J_b = 9.6 × 10^{-7} \, \mathrm{m^2 \, s^{-3}}",
               L"J_b = 2.4 × 10^{-7} \, \mathrm{m^2 \, s^{-3}}",
               L"J_b = 8.8 × 10^{-8} \, \mathrm{m^2 \, s^{-3}}"]

name = "extended_stability_conv_adj"
dir = "../parameters/"
filepath = joinpath(dir, string(name) * "_best_parameters.jld2")
file = jldopen(filepath)
optimal_parameters = file["optimal_parameters"]
close(file)

#=
dependent_parameters = dependent_parameter_sets[string(name)]
parameter_names = keys(optimal_parameters)
free_parameters = FreeParameters(prior_library; names=parameter_names, dependent_parameters)
optimal_parameters = build_parameters_named_tuple(free_parameters, optimal_parameters)

optimal_parameters = Dict(name => value for (name, value) in pairs(optimal_parameters))
optimal_parameters[:Cᶜc] = 1.0
@show optimal_parameters[:Cᵉc]
optimal_parameters[:Cᶜe] = 1.0
optimal_parameters = NamedTuple(name => value for (name, value) in pairs(optimal_parameters))

#@show optimal_parameters

optimal_parameters = (
    Cˢ   = 2.4, # 0.41,
    Cᵇ   = Inf,
    Cᶜc  = 1.5,
    Cᶜe  = 1.2,
    Cᵉc  = 0.2,
    Cᵉe  = 0.0,
    Cˢᵖ  = 0.14,
    Cˡᵒu = 0.19, # 0.46,
    Cʰⁱu = 0.086, # 0.21,
    Cˡᵒc = 0.20, # 0.49,
    Cʰⁱc = 0.045, # 0.11,
    Cˡᵒe = 1.9, # 4.5,
    Cʰⁱe = 0.57, # 1.4,
    CRiᵟ = 0.45,
    CRi⁰ = 0.47,
    # CˡᵒD = 2.3,
    # CʰⁱD = 6.7,
    CˡᵒD = 1.1, # 0.43,
    CʰⁱD = 0.37, # 0.15,
    CᶜD  = 0.88,
    Cᵂu★ = 1.1,
    CᵂwΔ = 4.0,
)
=#

parameters = (;
#    Cᶜc = 0.1,
#    Cᵘⁿc = 0.1,
#    Cʰⁱe = 0.01,
)

colors = [
    :royalblue1,
    :darkred,
    :black,
]

for (n, p) in enumerate(suite_parameters)

    grid_parameters = [
        (size=32, z=(-256, 0)),
        (size=128, z=(-256, 0)),
    ]

    closure = CATKEVerticalDiffusivity()
    
    batched_ip = build_batched_inverse_problem(closure, name;
                                               Nensemble = 1,
                                               Δt = 5minutes,
                                               grid_parameters,
                                               suite_parameters = [p])
    
    forward_run!(batched_ip, [optimal_parameters])
    #forward_run!(batched_ip, [parameters])

    @show s = p.name

    # Low resolution
    ip = batched_ip[1]
    c = 1 #free convection
    obs = ip.observations[c]
    bt = obs.field_time_serieses.b
    #et = obs.field_time_serieses.e
    t = bt.times
    grid = bt.grid
    Nt = length(t)
    Nx, Ny, Nz = size(grid)

    model = ip.simulation.model
    b, e = model.tracers
    b.boundary_conditions.top.condition[1, c]

    z = znodes(b)
    b_CATKE = interior(b, 1, c, :)
    e_CATKE = interior(e, 1, c, :)
    b₀ = b_CATKE[1]

    ss = suite_names[n]
    lines!(ax_b, b_CATKE .- b₀, z; linewidth=4, color=(colors[n], 1), label="$ss")
    lines!(ax_e, e_CATKE, z; linewidth=4, color=(colors[n], 1), label="$ss")

    κ = model.diffusivity_fields.κc
    κ_CATKE = interior(κ, 1, c, :)
    z = znodes(κ)
    κ_CATKE[Nz+1] = NaN
    lines!(ax_κ, κ_CATKE, z; linewidth=4, color=(colors[n], 0.8), label=s)

    # Mixing length
    grid = model.grid
    closure = first(model.closure)
    velocities = model.velocities
    tracers = model.tracers
    buoyancy = model.buoyancy
    diffusivities = model.diffusivity_fields
    clock = model.clock
    
    Cᵇ = closure.mixing_length.Cᵇ
    Cᶜ = closure.mixing_length.Cᶜc
    Cᵉ = closure.mixing_length.Cᵉc
    Cˢᵖ = closure.mixing_length.Cˢᵖ
    
    top_tracer_bcs = NamedTuple(c => tracers[c].boundary_conditions.top for c in propertynames(tracers))
   
    ℓᶜ = KernelFunctionOperation{Center, Center, Face}(tracer_mixing_lengthᶜᶜᶠ, grid,
                                                       closure, velocities, tracers, buoyancy, diffusivities.Jᵇ)

    l_CATKE = interior(compute!(Field(ℓᶜ)), 1, c, :)
    l_CATKE[Nz+1] = NaN
    @show maximum(l_CATKE)

    lines!(ax_l, l_CATKE, z; linewidth=4, color=(colors[n], 0.8), label=s)

    # LES data
    z = znodes(bt)
    b_LES = interior(bt[Nt], 1, 1, :)
    # e_LES = interior(et[Nt], 1, 1, :)

    #=
    if n == 3
        lines!(ax_b, b_LES, z; linewidth=8, color=(colors[n], 0.4), label="LES, $ss")
        lines!(ax_e, e_LES, z; linewidth=8, color=(colors[n], 0.4), label="LES, $ss")
    end
    =#
end

Legend(fig[0, 1:4], ax_b, nbanks=3, tellheight=true, framevisible=false)
# axislegend(ax_e)
# axislegend(ax_κ)

display(fig)

save("convective_length_scale.pdf", fig)

