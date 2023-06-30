using LaTeXStrings

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


parameter_guide = Dict(:CˡᵒD    => (name = "Dissipation parameter (TKE equation)", latex = L"C^{lo}_D", bounds = (0.1, 10.0)), 
                        :CʰⁱD    => (name = "Dissipation parameter (TKE equation)", latex = L"C^{hi}_D", bounds = (0.1, 10.0)), 

                        :Cᵂu★  => (name = "TKE subgrid flux parameter", latex = L"C^W_{u\star}", bounds = (0.0,  4.0)), 
                        :CᵂwΔ  => (name = "TKE subgrid flux parameter", latex = L"C^Ww\Delta", bounds = (0.0, 4.0)), 

                        :CRiᵟ => (name = "Stability function parameter", latex = L"CRi^{\delta}", bounds = (0.1, 1.0)), 
                        :CRi⁰ => (name = "Stability function parameter", latex = L"CRi^{\theta}", bounds = (0.1, 1.0)), 

                        :Cˡᵒu  => (name = "Velocity diffusivity LB", latex = L"C^{lo}u", bounds = (0.1, 1.0)), 
                        :Cʰⁱu  => (name = "Velocity diffusivity (UB-LB)/LB", latex = L"C^{hi}_u", bounds = (0.1, 1.0)), 
                        :Cˡᵒc  => (name = "Tracer diffusivity LB", latex = L"C^{lo}_c", bounds = (0.1, 1.0)), 
                        :Cʰⁱc  => (name = "Tracer diffusivity (UB-LB)/LB", latex = L"C^{hi}_c", bounds = (0.1, 1.0)), 
                        :Cˡᵒe  => (name = "TKE diffusivity LB", latex = L"C^{lo}_e", bounds = (0.1, 4.0)), 
                        :Cʰⁱe  => (name = "TKE diffusivity (UB-LB)/LB", latex = L"C^{hi}_e", bounds = (0.1, 4.0)), 

                        :Cᵇ   => (name = "Stratified mixing length parameter", latex = L"C^b", bounds = (0.0, 1.0)),
                        :Cˢ   => (name = "Shear mixing length coefficient", latex = L"C^s", bounds = (0.0, 1.0)),

                        :Cᶜc => (name = "Convective adjustment parameter", latex = L"C^c_c", bounds = (0.01, 1.0)),
                        :Cᶜe => (name = "Convective adjustment parameter", latex = L"C^c_e", bounds = (0.01, 1.0)),
                        :CᶜD => (name = "Convective adjustment parameter", latex = L"C^c_D", bounds = (0.01, 2.0)),
                        :Cᵉc => (name = "Convective adjustment parameter", latex = L"C^e_c", bounds = (0.01, 1.0)),
                        :Cᵉe => (name = "Convective adjustment parameter", latex = L"C^e_e", bounds = (0.01, 1.0)),
                        :CᵉD => (name = "Convective adjustment parameter", latex = L"C^e_D", bounds = (0.01, 1.0)),
                        :Cˢᵖ => (name = "Convective adjustment parameter", latex = L"C^{sp}", bounds = (0.01, 1.0)),
)

# Ri-based
# bounds_library[:ν₀]  = (0.0, 1.0)
# bounds_library[:κ₀]  = (0.0, 1.0)
# bounds_library[:κᶜᵃ] = (0.0, 2.0)
# bounds_library[:Cᵉⁿ] = (0.0, 2.0)
# bounds_library[:Cᵃᵛ] = (0.0, 2.0)
# bounds_library[:Ri₀] = (0.0, 2.0)
# bounds_library[:Riᵟ] = (0.0, 2.0)

# prior_library = Dict()

# for p in keys(bounds_library)
#     prior_library[p] = ScaledLogitNormal(; bounds=bounds_library[p])
# end

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