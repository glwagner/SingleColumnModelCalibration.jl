
parameter_sets = Dict()
dependent_parameter_sets = Dict()
bounds_library = Dict()

#####
##### CATKE
#####

parameter_sets["constant_Pr"] = (:CᵂwΔ, :Cᵂu★, :Cʰⁱc, :Cʰⁱu, :Cʰⁱe, :CʰⁱD, :Cˢ)
parameter_sets["variable_Pr"] = (:CᵂwΔ, :Cᵂu★, :Cʰⁱc, :Cʰⁱu, :Cʰⁱe, :CʰⁱD, :Cˢ,
                                 :Cˡᵒc, :Cˡᵒu, :Cˡᵒe, :CˡᵒD, :CRi⁰, :CRiᵟ)

parameter_sets["extended_stability"] = (:CᵂwΔ, :Cᵂu★, :Cʰⁱc, :Cʰⁱu, :Cʰⁱe, :CʰⁱD, :Cˢ,
                                        :Cˡᵒc, :Cˡᵒu, :Cˡᵒe, :CˡᵒD, :CRi⁰, :CRiᵟ,
                                        :Cᵘⁿc, :Cᵘⁿu, :Cᵘⁿe, :CᵘⁿD)
                                        
conv_adj_names = (:Cᶜc, :Cᶜu, :Cᶜe, :CᶜD, :Cᵉc, :Cˢᵖ)

for set in ["constant_Pr", "variable_Pr", "extended_stability"]
    names = parameter_sets[set]
    conv_adj_set = set * "_conv_adj"
    parameter_sets[conv_adj_set] =  tuple(names..., conv_adj_names...)
    dependent_parameter_sets[set] = NamedTuple()
    dependent_parameter_sets[conv_adj_set] = NamedTuple()
end

# Neutralize convective mixing length
CᶜD(θ) = 0.0 
Cᶜu(θ) = 0.0
Cᶜc(θ) = 0.0
Cᶜe(θ) = 0.0
Cᵉc(θ) = 0.0
Cˢᵖ(θ) = 0.0

dependent_parameter_sets["extended_stability"] = (; CᶜD, Cᶜu, Cᶜc, Cᶜe, Cᵉc, Cˢᵖ)

CˡᵒD(θ) = θ.CʰⁱD
Cˡᵒu(θ) = θ.Cʰⁱu
Cˡᵒc(θ) = θ.Cʰⁱc
Cˡᵒe(θ) = θ.Cʰⁱe

# Functions for the "unstable" branch of the stability functions
CᵘⁿD(θ) = θ.CʰⁱD
Cᵘⁿu(θ) = θ.Cʰⁱu
Cᵘⁿc(θ) = θ.Cʰⁱc
Cᵘⁿe(θ) = θ.Cʰⁱe

for set in ["constant_Pr", "constant_Pr_conv_adj"]
    dependent_parameter_sets[set] = (; Cˡᵒu, Cˡᵒc, Cˡᵒe, CˡᵒD, Cᵘⁿu, Cᵘⁿc, Cᵘⁿe, CᵘⁿD) 
end

for set in ["variable_Pr", "variable_Pr_conv_adj"]
    dependent_parameter_sets[set] = (; Cᵘⁿu, Cᵘⁿc, Cᵘⁿe, CᵘⁿD)
end

parameter_sets["fixed_Ric"] = (:CᵂwΔ, :Cᵂu★, :Cˢ,
                               :Cʰⁱc, :Cʰⁱu, :Cʰⁱe, :CʰⁱD,
                               :Cˡᵒc, :Cˡᵒu, :Cˡᵒe,
                               :CRi⁰, :CRiᵟ, :Cᶜc,
                               :Cᶜe, :CᶜD, :Cᵉc, :Cˢᵖ)

function CˡᵒD_fixed_Riᶜ(θ)
    Riᶜ = 0.25
    return max(0.0, θ.Cˡᵒu / Riᶜ - θ.Cˡᵒc)
end

dependent_parameter_sets["fixed_Ric"] = (; CˡᵒD = CˡᵒD_fixed_Riᶜ)

# Turbulent kinetic energy parameters
bounds_library[:CᵂwΔ] = (0.0, 8.0)
bounds_library[:Cᵂu★] = (0.0, 8.0)

# Mixing length parameters
bounds_library[:Cˢ]   = (0.0, 2.0)
bounds_library[:Cᵇ]   = (0.0, 2.0)

bounds_library[:Cˡᵒu] = (0.0, 2.0)
bounds_library[:Cˡᵒc] = (0.0, 2.0)
bounds_library[:Cˡᵒe] = (0.0, 10.0)
bounds_library[:CˡᵒD] = (0.0, 10.0)
bounds_library[:Cʰⁱc] = (0.0, 2.0)
bounds_library[:Cʰⁱu] = (0.0, 2.0)
bounds_library[:Cʰⁱe] = (0.0, 10.0)
bounds_library[:CʰⁱD] = (0.0, 10.0)
bounds_library[:Cᵘⁿc] = (0.0, 2.0)
bounds_library[:Cᵘⁿu] = (0.0, 2.0)
bounds_library[:Cᵘⁿe] = (0.0, 10.0)
bounds_library[:CᵘⁿD] = (0.0, 10.0)

bounds_library[:CRi⁰] = (0.0, 2.0)
bounds_library[:CRiᵟ] = (0.0, 2.0)

# Convective adjustment parameters
bounds_library[:Cᶜu]  = (0.0, 8.0)
bounds_library[:Cᶜc]  = (0.0, 8.0)
bounds_library[:Cᶜe]  = (0.0, 8.0)
bounds_library[:CᶜD]  = (0.0, 8.0)
bounds_library[:Cᵉc]  = (0.0, 1.0)
bounds_library[:Cᵉe]  = (0.0, 1.0)
bounds_library[:CᵉD]  = (0.0, 1.0)
bounds_library[:Cˢᵖ]  = (0.0, 1.0)

#####
##### Ri-based
#####

parameter_sets["ri_based"] = (:ν₀, :κ₀, :κᶜᵃ, :Cᵉⁿ, :Cᵃᵛ, :Ri₀, :Riᵟ)

bounds_library[:ν₀]  = (0.0, 1.0)
bounds_library[:κ₀]  = (0.0, 1.0)
bounds_library[:κᶜᵃ] = (0.0, 2.0)
bounds_library[:Cᵉⁿ] = (0.0, 2.0)
bounds_library[:Cᵃᵛ] = (0.0, 2.0)
bounds_library[:Ri₀] = (0.0, 2.0)
bounds_library[:Riᵟ] = (0.0, 2.0)

#####
##### TKEDissipation
#####

parameter_sets["constant_stabilities"] = (:Cσe, :Cσϵ, :Cu₀, :Cc₀, :Cᵋϵ, :Cᴾϵ, :Cᵇϵ, :Cᵂu★, :CᵂwΔ, :𝕊u₀)
dependent_parameter_sets["constant_stabilities"] = NamedTuple()

parameter_sets["variable_stabilities"] = (:Cσe, :Cσϵ,
                                          :Cu₀, :Cu₁, :Cu₂,  
                                          :Cc₀, :Cc₁, :Cc₂, 
                                          :Cd₁, :Cd₂, :Cd₃, :Cd₄, :Cd₅, 
                                          :Cᵋϵ, :Cᴾϵ, :Cᵇϵ, :Cᵂu★, :CᵂwΔ)

parameter_sets["dissipation_and_transport"] = (:Cσe, :Cσϵ, :Cᵋϵ, :Cᴾϵ, :Cᵇϵ⁺, :Cᵇϵ⁻, :Cᵂu★, :CᵂwΔ)
dependent_parameter_sets["dissipation_and_transport"] = NamedTuple()


# Neutralize convective mixing length
Cd₃(θ) = 0.0 
Cd₄(θ) = 0.0 
Cd₅(θ) = 0.0 
Cu₁(θ) = 0.0
Cu₂(θ) = 0.0
Cc₁(θ) = 0.0
Cc₂(θ) = 0.0

parameter_sets["very_simple_stabilities"] = (:Cσe, :Cσϵ, :Cᵋϵ, :Cᴾϵ, :Cᵇϵ⁺, :Cᵇϵ⁻, :Cᵂu★, :CᵂwΔ,
                                             :Cu₀, :Cc₀, :Cd₁, :Cd₂)
dependent_parameter_sets["very_simple_stabilities"] = (; Cd₃, Cd₄, Cd₅, Cu₁, Cu₂, Cc₁, Cc₂)

parameter_sets["pretty_simple_stabilities"] = (:Cσe, :Cσϵ, :Cᵋϵ, :Cᴾϵ, :Cᵇϵ⁺, :Cᵇϵ⁻, :Cᵂu★, :CᵂwΔ,
                                               :Cu₀, :Cu₁, :Cu₂, :Cc₀, :Cc₁, :Cc₂, :Cd₁, :Cd₂)
dependent_parameter_sets["pretty_simple_stabilities"] = (; Cd₃, Cd₄, Cd₅, Cu₁, Cu₂, Cc₁, Cc₂)

bounds_library[:Cσe] = (0.0, 10.0)      
bounds_library[:Cσϵ] = (0.0, 10.0)
bounds_library[:Cᵋϵ] = (0.0, 10.0)
bounds_library[:Cᴾϵ] = (0.0, 10.0)
bounds_library[:Cᵇϵ] = (-10.0, 0.0)
bounds_library[:Cᵇϵ⁻] = (-10.0, 0.0)
bounds_library[:Cᵇϵ⁺] = (0.0, 10.0)

bounds_library[:Cu₀] = (0.0, 0.2)
bounds_library[:Cu₁] = (0.0, 0.04)
bounds_library[:Cu₂] = (-0.001, 0.001)
bounds_library[:Cc₀] = (0.0, 0.2)
bounds_library[:Cc₁] = (0.0, 0.008)
bounds_library[:Cc₂] = (0.0, 0.01)
bounds_library[:Cd₁] = (0.0, 0.4)
bounds_library[:Cd₂] = (0.0, 0.05)
bounds_library[:Cd₃] = (-0.1, 0.1)
bounds_library[:Cd₄] = (-0.1, 0.1)
bounds_library[:Cd₅] = (-0.1, 0.1)
bounds_library[:𝕊u₀] = (0.0, 2.0)

prior_library = Dict()
for p in keys(bounds_library)
    prior_library[p] = ScaledLogitNormal(; bounds=bounds_library[p])
end

get_free_parameters(name) = FreeParameters(prior_library, names = parameter_sets[name],
                                           dependent_parameters = dependent_parameter_sets[name])

