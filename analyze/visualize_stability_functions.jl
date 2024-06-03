using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Printf
using JLD2
using CairoMakie

set_theme!(Theme(fontsize=20, linewidth=3))

@load "optimal_catke.jld2"

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
    Ri⁺ = CRi⁰ + CRiᵟ
)

for name in keys(optimal_parameters)
    println(name, ": ", optimal_parameters[name])
end

println()

for name in keys(derived_parameters)
    println(name, ": ", derived_parameters[name])
end

fig = Figure(size=(1200, 400))
axSt = Axis(fig[1, 1], xlabel="Richardson number, N² / |∂z u|²", ylabel="Stability functions")
axPr = Axis(fig[1, 2], xlabel="Richardson number, N² / |∂z u|²", ylabel="Prandtl and TKE Schmidt numbers")

Ri = -1:0.0001:2

function 𝕊(Ri, un, lo, hi)
    ϵ = (Ri - CRi⁰) / CRiᵟ
    if Ri < 0
        return un
    elseif Ri < CRi⁰
        return lo
    elseif Ri < CRi⁰ + CRiᵟ
        return lo + (hi - lo) * ϵ * Ri
    elseif Ri >= CRi⁰ + CRiᵟ
        return hi
    end
end

lines!(axSt, Ri, 𝕊.(Ri, Cᵘⁿc, Cˡᵒc, Cʰⁱc), label="Tracers") 
lines!(axSt, Ri, 𝕊.(Ri, Cᵘⁿu, Cˡᵒu, Cʰⁱu), label="Momentum") 
lines!(axSt, Ri, 𝕊.(Ri, Cᵘⁿe, Cˡᵒe, Cʰⁱe), label="TKE") 
lines!(axSt, Ri, 1 ./ 𝕊.(Ri, CᵘⁿD, CˡᵒD, CʰⁱD), label="Dissipation (1/𝕊ᴰ)") 

lines!(axPr, Ri, 𝕊.(Ri, Cᵘⁿu, Cˡᵒu, Cʰⁱu) ./ 𝕊.(Ri, Cᵘⁿc, Cˡᵒc, Cʰⁱc), label="Prandtl number")
lines!(axPr, Ri, 𝕊.(Ri, Cᵘⁿu, Cˡᵒu, Cʰⁱu) ./ 𝕊.(Ri, Cᵘⁿe, Cˡᵒe, Cʰⁱe), label="Schmidt number for TKE")

axislegend(axSt, position=:rt)
axislegend(axPr, position=:lt)

save("stability_Prandtl_Schmidt_numbers.pdf", fig)
