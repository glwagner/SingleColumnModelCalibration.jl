using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Printf
using JLD2
using CairoMakie

set_theme!(Theme(fontsize=24, linewidth=4, alpha=0.6))

d = 4
dash = Linestyle([0.0, d, 1.6d, 2.6d])

@load "optimal_catke.jld2"
@show optimal_catke

Cˡᵒc = optimal_catke.mixing_length.Cˡᵒc
Cˡᵒu = optimal_catke.mixing_length.Cˡᵒu
Cˡᵒe = optimal_catke.mixing_length.Cˡᵒe
Cʰⁱc = optimal_catke.mixing_length.Cʰⁱc
Cʰⁱu = optimal_catke.mixing_length.Cʰⁱu
Cʰⁱe = optimal_catke.mixing_length.Cʰⁱe
Cᶜc  = optimal_catke.mixing_length.Cᶜc
Cᶜu  = optimal_catke.mixing_length.Cᶜu
Cᶜe  = optimal_catke.mixing_length.Cᶜe
Cᵘⁿc = optimal_catke.mixing_length.Cᵘⁿc
Cᵘⁿu = optimal_catke.mixing_length.Cᵘⁿu
Cᵘⁿe = optimal_catke.mixing_length.Cᵘⁿe
Cˢ   = optimal_catke.mixing_length.Cˢ
CRiᵟ = optimal_catke.mixing_length.CRiᵟ
CRi⁰ = optimal_catke.mixing_length.CRi⁰
CˡᵒD = optimal_catke.turbulent_kinetic_energy_equation.CˡᵒD
CʰⁱD = optimal_catke.turbulent_kinetic_energy_equation.CʰⁱD
CᵘⁿD = optimal_catke.turbulent_kinetic_energy_equation.CᵘⁿD

optimal_parameters = (;
    Cˡᵒc,
    Cˡᵒu,
    Cˡᵒe, 
    Cʰⁱc, 
    Cʰⁱu, 
    Cʰⁱe, 
    Cᶜc,  
    Cᶜu, 
    Cᶜe,  
    Cᵘⁿc, 
    Cᵘⁿu, 
    Cᵘⁿe, 
    Cˢ,   
    CRiᵟ,
    CRi⁰, 
    CˡᵒD,
    CʰⁱD,
    CᵘⁿD)

derived_parameters = (
    Pr₀ = Cˡᵒu / Cˡᵒc,
    Pr∞ = Cʰⁱu / Cʰⁱc,
    Pr⁻ = Cᵘⁿu / Cᵘⁿc,
    Prᶜ = Cᶜu / Cᶜc,
    Sc₀ = Cˡᵒu / Cˡᵒe,
    Sc∞ = Cʰⁱu / Cʰⁱe,
    Sc⁻ = Cᵘⁿu / Cᵘⁿe,
    Scᶜ = Cᶜu / Cᶜe,
    Ri★ = Cˡᵒu / (Cˡᵒc + CˡᵒD),
    ϰvk = Cˢ * (Cˡᵒu^3 / CˡᵒD)^(1/4),
    Γ₀ = Cˡᵒc / CˡᵒD,
    Γ∞ = Cʰⁱc / CʰⁱD,
    Ri⁺ = CRi⁰ + CRiᵟ
)

for name in keys(optimal_parameters)
    println(name, ": ", optimal_parameters[name])
end

println()

for name in keys(derived_parameters)
    println(name, ": ", derived_parameters[name])
end

fig = Figure(size=(1200, 500))

#xlabel = L"Richardson number, $N^2 / | \partial_z u |^2$"
xlabel = "Richardson number"

axSt = Axis(fig[1, 1]; xlabel, ylabel="Stability \n functions")
axSte = Axis(fig[2, 1]; xlabel, ylabel="TKE stability \n function")
axPr = Axis(fig[1:2, 2]; xlabel, yaxisposition=:right, ylabel="Prandtl and TKE Schmidt numbers")


Ri = -1:1e-3:2

@inline step(x, c, w) = max(zero(x), min(one(x), (x - c) / w))

@inline function 𝕊(Ri, σ⁻, σ⁰, σ∞, c=CRi⁰, w=CRiᵟ)
    σ⁺ = σ⁰ + (σ∞ - σ⁰) * step(Ri, c, w)
    σ = σ⁻ * (Ri < 0) + σ⁺ * (Ri ≥ 0)
    return σ
end

alpha = 0.6
tkecolor = :tomato

lines!(axSt, Ri, 𝕊.(Ri, Cᵘⁿc, Cˡᵒc, Cʰⁱc); alpha, label="Tracers") 
lines!(axSt, Ri, 𝕊.(Ri, Cᵘⁿu, Cˡᵒu, Cʰⁱu); alpha, label="Momentum") 
lines!(axSt, Ri, 1 ./ 𝕊.(Ri, CᵘⁿD, CˡᵒD, CʰⁱD); alpha, label="Dissipation (1/𝕊ᴰ)") 

lne = lines!(axSt, Ri, NaN .* 𝕊.(Ri, Cᵘⁿe, Cˡᵒe, Cʰⁱe), color=tkecolor, label="TKE") 
lines!(axSte, Ri, 𝕊.(Ri, Cᵘⁿe, Cˡᵒe, Cʰⁱe), color=lne.color, label="TKE") 

lines!(axPr, Ri, 𝕊.(Ri, Cᵘⁿu, Cˡᵒu, Cʰⁱu) ./ 𝕊.(Ri, Cᵘⁿc, Cˡᵒc, Cʰⁱc),
       color=:black, label="Prandtl number")
lines!(axPr, Ri, 𝕊.(Ri, Cᵘⁿu, Cˡᵒu, Cʰⁱu) ./ 𝕊.(Ri, Cᵘⁿe, Cˡᵒe, Cʰⁱe),
       color=:black, linestyle=:dot, label="Schmidt number for TKE")

#axislegend(axSt, position=:lt)
Legend(fig[0, 1:2], axSt, framevisible=false, nbanks=4, tellheight=true)
axislegend(axPr, position=:lt)

hidexdecorations!(axSt, grid=false)
hidespines!(axSt, :t, :r, :b)
hidespines!(axPr, :t, :l)
hidespines!(axSte, :t, :r)

#xlims!(axSt, -1, 2)
ylims!(axSt, 0, 2.0)

display(fig)

save("stability_Prandtl_Schmidt_numbers.pdf", fig)
