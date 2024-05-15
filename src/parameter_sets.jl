parameter_sets = Dict(
    "ri_based"    => (:ν₀, :κ₀, :κᶜᵃ, :Cᵉⁿ, :Cᵃᵛ, :Ri₀, :Riᵟ),
    "constant_Pr" => (:CᵂwΔ, :Cᵂu★, :Cʰⁱc, :Cʰⁱu, :Cʰⁱe, :CʰⁱD, :Cˢ),
    "variable_Pr" => (:CᵂwΔ, :Cᵂu★, :Cʰⁱc, :Cʰⁱu, :Cʰⁱe, :CʰⁱD, :Cˢ, :Cˡᵒc, :Cˡᵒu, :Cˡᵒe, :CˡᵒD, :CRi⁰, :CRiᵟ),
    "extended_stability" => (:CᵂwΔ, :Cᵂu★, :Cʰⁱc, :Cʰⁱu, :Cʰⁱe, :CʰⁱD, :Cˢ, :Cˡᵒc, :Cˡᵒu, :Cˡᵒe, :CˡᵒD, :CRi⁰, :CRiᵟ, :Cᵘⁿc, :Cᵘⁿu, :Cᵘⁿe, :CᵘⁿD),
                             
)

conv_adj_names = (:Cᶜc, :Cᶜu, :Cᶜe, :CᶜD, :Cᵉc, :Cˢᵖ)
#conv_adj_names = (:Cᶜc, :Cᶜe, :CᶜD, :Cˢᵖ)

for set in ["constant_Pr", "variable_Pr", "extended_stability"]
    names = parameter_sets[set]
    conv_adj_set = set * "_conv_adj"
    parameter_sets[conv_adj_set] =  tuple(names..., conv_adj_names...)
end

# Some dependent parameters
dependent_parameter_sets = Dict()
for (set, names) in parameter_sets
    dependent_parameter_sets[set] = NamedTuple()
end

CˡᵒD(θ) = θ.CʰⁱD
Cˡᵒu(θ) = θ.Cʰⁱu
Cˡᵒc(θ) = θ.Cʰⁱc
Cˡᵒe(θ) = θ.Cʰⁱe

CᵘⁿD(θ) = θ.CˡᵒD
Cᵘⁿu(θ) = θ.Cˡᵒu
Cᵘⁿc(θ) = θ.Cˡᵒc
Cᵘⁿe(θ) = θ.Cˡᵒe

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

#####
##### Bounds and priors
#####

bounds_library = Dict()

# Turbulent kinetic energy parameters
bounds_library[:CᵂwΔ] = (0.0, 8.0)
bounds_library[:Cᵂu★] = (0.0, 8.0)

# Mixing length parameters
bounds_library[:Cˢ]   = (0.0, 2.0)
bounds_library[:Cᵇ]   = (0.0, 2.0)

bounds_library[:Cˡᵒu] = (0.0, 1.0)
bounds_library[:Cˡᵒc] = (0.0, 1.0)
bounds_library[:Cˡᵒe] = (0.0, 8.0)
bounds_library[:CˡᵒD] = (0.0, 8.0)
bounds_library[:Cʰⁱc] = (0.0, 1.0)
bounds_library[:Cʰⁱu] = (0.0, 1.0)
bounds_library[:Cʰⁱe] = (0.0, 8.0)
bounds_library[:CʰⁱD] = (0.0, 8.0)
bounds_library[:Cᵘⁿc] = (0.0, 2.0)
bounds_library[:Cᵘⁿu] = (0.0, 2.0)
bounds_library[:Cᵘⁿe] = (0.0, 8.0)
bounds_library[:CᵘⁿD] = (0.0, 8.0)

bounds_library[:CRi⁰] = (0.0, 1.0)
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

# Ri-based
bounds_library[:ν₀]  = (0.0, 1.0)
bounds_library[:κ₀]  = (0.0, 1.0)
bounds_library[:κᶜᵃ] = (0.0, 2.0)
bounds_library[:Cᵉⁿ] = (0.0, 2.0)
bounds_library[:Cᵃᵛ] = (0.0, 2.0)
bounds_library[:Ri₀] = (0.0, 2.0)
bounds_library[:Riᵟ] = (0.0, 2.0)

prior_library = Dict()

for p in keys(bounds_library)
    prior_library[p] = ScaledLogitNormal(; bounds=bounds_library[p])
end

get_free_parameters(name) = FreeParameters(prior_library, names = parameter_sets[name],
                                           dependent_parameters = dependent_parameter_sets[name])

