using Oceananigans
using Oceananigans.Units

using Printf
using JLD2
using LinearAlgebra

using CairoMakie
using ElectronDisplay

set_theme!(Theme(fontsize=24))

case_path(case, suite, resolution="1m") = joinpath("..", "data",
                                                   suite,
                                                   resolution,
                                                   case * "_instantaneous_statistics.jld2")

momentum_fluxes = Dict()
buoyancy_fluxes = Dict()
heat_fluxes = Dict()
u₁₀s = Dict()

suites = [
    "12_hour_suite",
    "24_hour_suite",
    "48_hour_suite",
]

cases = [
    "free_convection",
    "weak_wind_strong_cooling",
    "med_wind_med_cooling",
    "strong_wind_weak_cooling",
    "strong_wind",
    "strong_wind_no_rotation",
]

ρᵃ = 1.225
ρʷ = 1024
cᴰ = 1e-3
α = 2e-4
cᴾ = 3991
g = 9.81

for suite in suites
    suite_u₁₀ = Float64[]
    suite_momentum_fluxes = Float64[]
    suite_buoyancy_fluxes = Float64[]
    suite_heat_fluxes = Float64[]
    for case in cases
        file = jldopen(case_path(case, suite))
        momentum_flux = abs(file["parameters/momentum_flux"])
        buoyancy_flux = file["parameters/buoyancy_flux"]
        close(file)

        push!(suite_momentum_fluxes, momentum_flux)
        push!(suite_buoyancy_fluxes, buoyancy_flux)

        # u10
        τ = ρʷ * momentum_flux
        u₁₀ = sqrt(τ / ρᵃ / cᴰ)
        push!(suite_u₁₀, u₁₀)

        # Heat flux
        # Q = ρ * cᴾ * Qᵀ = ρ * cᴾ / (α * g)
        Q = ρʷ * cᴾ / (α * g) * buoyancy_flux
        push!(suite_heat_fluxes, Q)
    end

    u₁₀s[suite] = suite_u₁₀
    momentum_fluxes[suite] = suite_momentum_fluxes
    buoyancy_fluxes[suite] = suite_buoyancy_fluxes
    heat_fluxes[suite] = suite_heat_fluxes
end

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Momentum flux (N m⁻²)", ylabel="Heat flux (W m⁻²)")
markersize=20

for suite in suites
    #scatter!(ax, momentum_fluxes[suite], buoyancy_fluxes[suite])
    #scatter!(ax, u₁₀s[suite], buoyancy_fluxes[suite])
    sc = scatter!(ax, ρʷ .* momentum_fluxes[suite][1:end-1], heat_fluxes[suite][1:end-1]; markersize)
    scatter!(ax, ρʷ .* momentum_fluxes[suite][end], heat_fluxes[suite][end]; color=sc.color, markersize, marker=:utriangle)
end

display(fig)

