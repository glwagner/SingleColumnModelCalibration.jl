parameter_sets = Dict(
    "ri_based"    => (:ν₀, :κ₀, :κᶜᵃ, :Cᵉⁿ, :Cᵃᵛ, :Ri₀, :Riᵟ),
    "constant_Pr" => (:CᵂwΔ, :Cᵂu★, :Cʰⁱc, :Cʰⁱu, :Cʰⁱe, :CʰⁱD, :Cˢ),
    "variable_Pr" => (:CᵂwΔ, :Cᵂu★, :Cʰⁱc, :Cʰⁱu, :Cʰⁱe, :CʰⁱD, :Cˢ, :Cˡᵒc, :Cˡᵒu, :Cˡᵒe, :CˡᵒD, :CRi⁰, :CRiᵟ),
)

conv_adj_names = (:Cᶜc, :Cᶜe, :CᶜD, :Cᵉc, :Cˢᵖ)

for set in ["constant_Pr", "variable_Pr"]
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

for set in ["constant_Pr", "constant_Pr_conv_adj"]
    dependent_parameter_sets[set] = (; Cˡᵒu, Cˡᵒc, Cˡᵒe, CˡᵒD) 
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

parameter_guide = Dict(:C⁻D    => (name = "Dissipation parameter (TKE equation)", latex = L"C^D", bounds = (0.1, 10.0)), 
                        :C⁺D    => (name = "Dissipation parameter (TKE equation)", latex = L"C^D", bounds = (0.1, 10.0)), 
                        # :Cᵟu   => (name = "Ratio of mixing length to grid spacing", latex = L"C^{\delta}u" bounds = (0.0, 4.0)), 
                        # :Cᵟc   => (name = "Ratio of mixing length to grid spacing", latex = L"C^{\delta}c" bounds = (0.0, 4.0)), 
                        # :Cᵟe   => (name = "Ratio of mixing length to grid spacing", latex = L"C^{\delta}e" bounds = (0.0, 4.0)), 
                        :Cᵂu★  => (name = "TKE subgrid flux parameter", latex = L"C^W_{u\star}", bounds = (0.1,  10.0)), 
                        :CᵂwΔ  => (name = "TKE subgrid flux parameter", latex = L"C^Ww\Delta", bounds = (0.1, 10.0)), 
                        :CRiʷ => (name = "Stability function parameter", latex = L"CRi^w", bounds = (0.1, 2.0)), 
                        :CRiᶜ => (name = "Stability function parameter", latex = L"CRi^c", bounds = (0.1, 1.0)), 
                        :C⁻u  => (name = "Velocity diffusivity LB", latex = L"C^Ku^-", bounds = (0.1, 1.0)), 
                        :C⁺u  => (name = "Velocity diffusivity (UB-LB)/LB", latex = L"C^Ku^r", bounds = (0.1, 1.0)), 
                        :C⁻c  => (name = "Tracer diffusivity LB", latex = L"C^Kc^-", bounds = (0.1, 1.0)), 
                        :C⁺c  => (name = "Tracer diffusivity (UB-LB)/LB", latex = L"C^Kc^r", bounds = (0.1, 1.0)), 
                        :C⁻e  => (name = "TKE diffusivity LB", latex = L"C^Ke^-", bounds = (0.1, 10.0)), 
                        :C⁺e  => (name = "TKE diffusivity (UB-LB)/LB", latex = L"C^Ke^r", bounds = (0.1, 10.0)), 
                        # :Cᴬu   => (name = "Convective mixing length parameter", latex = L"C^A_U", bounds = (0.0, 10.0)), 
                        # :Cᴬc   => (name = "Convective mixing length parameter", latex = L"C^A_C", bounds = (0.0, 5.0)), 
                        # :Cᴬe   => (name = "Convective mixing length parameter", latex = L"C^A_E", bounds = (0.0, 30.0)),
                        :Cᵇ   => (name = "Stratified mixing length parameter", latex = L"C^b_U", bounds = (0.2, 2.0)),
                        # :Cᵇu   => (name = "Stratified mixing length parameter", latex = L"C^b_U", bounds = (0.0, 2.0)),
                        # :Cᵇc   => (name = "Stratified mixing length parameter", latex = L"C^b_C", bounds = (0.0, 2.0)),
                        # :Cᵇe   => (name = "Stratified mixing length parameter", latex = L"C^b_E", bounds = (0.0, 2.0)),
                        :Cˢ   => (name = "Shear mixing length coefficient", latex = L"C^s_U", bounds = (0.2, 2.0)),
                        # :Cˢu   => (name = "Shear mixing length coefficient", latex = L"C^s_U", bounds = (0.0, 3.0)),
                        # :Cˢc   => (name = "Shear mixing length coefficient", latex = L"C^s_C", bounds = (0.0, 3.0)),
                        # :Cˢe   => (name = "Shear mixing length coefficient", latex = L"C^s_E", bounds = (0.0, 3.0)),
                        # :Cʷ★   => (name = "Softmin strength parameter", latex = L"C^w_\star", bounds = (1.0, 4.0)),
                        # :Cʷℓ   => (name = "Softmax strength parameter", latex = L"C^w_\ell", bounds = (0.0, 4.0)),
                        :Cᶜc => (latex = L"C^c_C", bounds = (0.01, 10.0)),
                        :Cᶜe => (latex = L"C^c_E", bounds = (0.01, 10.0)),
                        :CᶜD => (latex = L"C^c_D", bounds = (0.01, 10.0)),
                        :Cᵉc => (latex = L"C^e_C", bounds = (0.01, 2.0)),
                        :Cᵉe => (latex = L"C^e_E", bounds = (0.01, 2.0)),
                        :CᵉD => (latex = L"C^e_D", bounds = (0.01, 2.0)),
                        :Cˢᶜ => (latex = L"C^{sc}", bounds = (0.01, 10.0)),
)

#####
##### Bounds and priors
#####

bounds_library = Dict()

# Turbulent kinetic energy parameters
bounds_library[:CᵂwΔ] = (0.0, 4.0)
bounds_library[:Cᵂu★] = (0.0, 4.0)

# Mixing length parameters
bounds_library[:Cˢ]   = (0.0, 1.0)
bounds_library[:Cᵇ]   = (0.0, 1.0)

bounds_library[:CˡᵒD] = (0.0, 1.0)
bounds_library[:Cˡᵒu] = (0.0, 1.0)
bounds_library[:Cˡᵒc] = (0.0, 1.0)
bounds_library[:Cˡᵒe] = (0.0, 4.0)

bounds_library[:Cʰⁱc] = (0.0, 1.0)
bounds_library[:Cʰⁱu] = (0.0, 1.0)
bounds_library[:CʰⁱD] = (0.0, 1.0)
bounds_library[:Cʰⁱe] = (0.0, 4.0)

bounds_library[:CRi⁰] = (0.0, 1.0)
bounds_library[:CRiᵟ] = (0.0, 1.0)

# Convective adjustment parameters
bounds_library[:Cᶜc]  = (0.0, 1.0)
bounds_library[:Cᶜe]  = (0.0, 1.0)
bounds_library[:CᶜD]  = (0.0, 2.0)
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

bounds(name) = parameter_guide[name].bounds

function named_tuple_map(names, f)
    names = Tuple(names)
    return NamedTuple{names}(f.(names))
end

f(p) = ScaledLogitNormal(; bounds=parameter_guide[p].bounds)

function get_free_parameters(name; f = f)

    names = parameter_sets[name]
    dependent_parameters = dependent_parameter_sets[name]

    prior_library = named_tuple_map(names, f)

    return FreeParameters(prior_library; names, dependent_parameters)
end