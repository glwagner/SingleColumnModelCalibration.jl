using Oceananigans
using JLD2
using GLMakie
using ParameterEstimocean
using Oceananigans.Units

using ParameterEstimocean.Parameters: build_parameters_named_tuple

using SingleColumnModelCalibration:
    dependent_parameter_sets,
    build_batched_inverse_problem,
    prior_library

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity, tracer_mixing_lengthᶜᶜᶠ

set_theme!(Theme(fontsize=32))

n = 70
ulims = (-0.15, 0.35)
elims = (-1e-4, 4e-4)
linewidth = 12
prefix = "free_convection"

fig = Figure(resolution=(1600, 600))

ax_b = Axis(fig[1, 1]; xlabel="Buoyancy (m s⁻²)", ylabel="z (m)", xticks=[0.0386, 0.0388])
ax_e = Axis(fig[1, 2]; xlabel="Turbulent kinetic \n energy (m² s⁻²)", xticks=[0, 1e-3, 2e-3]) #, yaxisposition=:right)
ax_l = Axis(fig[1, 3]; xlabel="Tracer \n mixing length (m)")
ax_κ = Axis(fig[1, 4]; xlabel="Eddy \n diffusivity (m² s⁻¹)", ylabel="z (m)", yaxisposition=:right)

xlims!(ax_b, 0.0386, 0.03885)
xlims!(ax_κ, -0.1, 6)

ylims!(ax_b, -192, 0)
ylims!(ax_e, -192, 0)
ylims!(ax_l, -192, 0)
ylims!(ax_κ, -192, 0)

hidespines!(ax_b, :t, :r)
hidespines!(ax_e, :t, :l, :r)
hidespines!(ax_l, :t, :l, :r)
hidespines!(ax_κ, :t, :l)

hideydecorations!(ax_e, grid=false)
hideydecorations!(ax_l, grid=false)

suite_parameters = [
    (name = "6_hour_suite", resolution="0.75m", stop_time=6hours),
    (name = "24_hour_suite", resolution="1m", stop_time=24hours),
    (name = "72_hour_suite", resolution="1m", stop_time=66hours),
]   

suite_names = [
    "Qᵇ = 9.6 × 10⁻⁷",
    "Qᵇ = 2.4 × 10⁻⁷",
    "Qᵇ = 8.8 × 10⁻⁸",
]

name = "variable_Pr_conv_adj"
dir = "../parameters/"
filepath = joinpath(dir, string(name) * "_best_parameters.jld2")
file = jldopen(filepath)
optimal_parameters = file["optimal_parameters"]
close(file)

dependent_parameters = dependent_parameter_sets[string(name)]
parameter_names = keys(optimal_parameters)
free_parameters = FreeParameters(prior_library; names=parameter_names, dependent_parameters)
optimal_parameters = build_parameters_named_tuple(free_parameters, optimal_parameters)

optimal_parameters = Dict(name => value for (name, value) in pairs(optimal_parameters))
optimal_parameters[:Cᶜc] = 1.0
@show optimal_parameters[:Cᵉc]
optimal_parameters[:Cᶜe] = 1.0
optimal_parameters = NamedTuple(name => value for (name, value) in pairs(optimal_parameters))

@show optimal_parameters

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
                                               Δt = 1minutes,
                                               grid_parameters,
                                               suite_parameters = [p])
    
    forward_run!(batched_ip, [optimal_parameters])

    @show s = p.name

    # Low resolution
    ip = batched_ip[1]
    c = 1 #free convection
    obs = ip.observations[c]
    bt = obs.field_time_serieses.b
    et = obs.field_time_serieses.e
    t = bt.times
    grid = bt.grid
    Nt = length(t)
    Nx, Ny, Nz = size(grid)

    model = ip.simulation.model
    b, e = model.tracers
    @show b.boundary_conditions.top.condition[1, c]

    z = znodes(b)
    b_CATKE = interior(b, 1, c, :)
    e_CATKE = interior(e, 1, c, :)

    ss = suite_names[n]
    lines!(ax_b, b_CATKE, z; linewidth=4, color=(colors[n], 1), label="$ss")
    lines!(ax_e, e_CATKE, z; linewidth=4, color=(colors[n], 1), label="$ss")

    κ = model.diffusivity_fields.κᶜ
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
    clock = model.clock
    
    Cᵇ = closure.mixing_length.Cᵇ
    Cᶜ = closure.mixing_length.Cᶜc
    Cᵉ = closure.mixing_length.Cᵉc
    Cˢᶜ = closure.mixing_length.Cˢᶜ
    
    top_tracer_bcs = NamedTuple(c => tracers[c].boundary_conditions.top for c in propertynames(tracers))
   
    ℓᶜ = KernelFunctionOperation{Center, Center, Face}(tracer_mixing_lengthᶜᶜᶠ, grid,
                                                       closure, velocities, tracers, buoyancy, clock, top_tracer_bcs)

    l_CATKE = interior(compute!(Field(ℓᶜ)), 1, c, :)
    l_CATKE[Nz+1] = NaN
    @show maximum(l_CATKE)

    lines!(ax_l, l_CATKE, z; linewidth=4, color=(colors[n], 0.8), label=s)

    # LES data
    z = znodes(bt)
    b_LES = interior(bt[Nt], 1, 1, :)
    e_LES = interior(et[Nt], 1, 1, :)

    #=
    if n == 3
        lines!(ax_b, b_LES, z; linewidth=8, color=(colors[n], 0.4), label="LES, $ss")
        lines!(ax_e, e_LES, z; linewidth=8, color=(colors[n], 0.4), label="LES, $ss")
    end
    =#
end

axislegend(ax_b, position=:lt, labelsize=30)
# axislegend(ax_e)
# axislegend(ax_κ)

display(fig)

save("convective_length_scale.png", fig)

#=
#####
##### Run the parameterization
#####

#=
closure = first(model.closure)
velocities = model.velocities
tracers = model.tracers
buoyancy = model.buoyancy
clock = model.clock

Cᵇ = closure.mixing_length.Cᵇ
Cᶜ = closure.mixing_length.Cᶜc
Cᵉ = closure.mixing_length.Cᵉc
Cˢᶜ = closure.mixing_length.Cˢᶜ

top_tracer_bcs = NamedTuple(c => tracers[c].boundary_conditions.top for c in propertynames(tracers))
convective_length_scale_args = (first(model.closure), Cᶜ, Cᵉ, Cˢᶜ,
                                velocities,
                                tracers,
                                buoyancy,
                                clock,
                                top_tracer_bcs)

ℓᶜconv = KernelFunctionOperation{Center, Center, Face}(convective_length_scaleᶜᶜᶠ, grid, convective_length_scale_args...)

C⁻ = closure.mixing_length.C⁻c
C⁺ = closure.mixing_length.C⁺c
σ = KernelFunctionOperation{Center, Center, Face}(stability_functionᶜᶜᶠ, grid, closure, C⁻, C⁺, velocities, tracers, buoyancy)
ℓ_stable = KernelFunctionOperation{Center, Center, Face}(stable_length_scaleᶜᶜᶠ, grid, Cᵇ, tracers.e, velocities, tracers, buoyancy)
ℓᶜshear = σ * ℓ_stable

ℓᶜ = KernelFunctionOperation{Center, Center, Face}(tracer_mixing_lengthᶜᶜᶠ, grid,
                                                   closure, velocities, tracers, buoyancy, clock, top_tracer_bcs)

=#


dir = "../parameters"
name = "variable_Pr_conv_adj"
closure_label = "CATKE"
closure = CATKEVerticalDiffusivity()

filepath = joinpath(dir, string(name) * "_best_parameters.jld2")
file = jldopen(filepath)
optimal_parameters = file["optimal_parameters"]
close(file)

# Overwrite...
optimal_parameters = NamedTuple()
dependent_parameters = dependent_parameter_sets[string(name)]
parameter_names = keys(optimal_parameters)
free_parameters = FreeParameters(prior_library; names=parameter_names, dependent_parameters)
optimal_parameters = build_parameters_named_tuple(free_parameters, optimal_parameters)
@show optimal_parameters


ip = batched_ip[1]
times = observation_times(ip.observations)
Nt = length(times)

c1 = 1
b1 = @lift interior(ip.time_series_collector.field_time_serieses.b[$n], 1, c1, :)
u1 = @lift interior(ip.time_series_collector.field_time_serieses.u[$n], 1, c1, :)
v1 = @lift interior(ip.time_series_collector.field_time_serieses.v[$n], 1, c1, :)
e1 = @lift interior(ip.time_series_collector.field_time_serieses.e[$n], 1, c1, :)

z = znodes(ip.time_series_collector.field_time_serieses.b)

linestyle = :solid
color = :darkgreen
lines!(ax_b1, b1, z; color,   linestyle, linewidth=4, label="$closure_label")
lines!(ax_u1, u1, z; color,   linestyle, linewidth=4, label="u, $closure_label")
lines!(ax_u1, v1, z; color=:darkred, linestyle, linewidth=4, label="v, $closure_label")
lines!(ax_e1, e1, z; color,   linestyle, linewidth=4, label="$closure_label")

axislegend(ax_b1, position=:rb)

display(fig)

# record(fig, "animate_les_data_with_catke.mp4", 1:Nt, framerate=12) do nn
#     n[] = nn
# end
=#
