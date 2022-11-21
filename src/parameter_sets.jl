parameter_sets = Dict(
    "ri_based"    => (:ν₀, :κ₀, :κᶜ, :Cᵉ, :Ri₀, :Riᵟ),
    "constant_Pr" => (:CᵂwΔ, :Cᵂu★, :Cᴷc⁺, :Cᴷu⁺, :Cᴷe⁺, :Cᵇ, :Cˢ, :Cᴰ⁺),
    "variable_Pr" => (:CᵂwΔ, :Cᵂu★, :Cᴷc⁺, :Cᴷu⁺, :Cᴷe⁺, :Cᵇ, :Cˢ, :Cᴰ⁺, :Cᴰ⁻,
                      :Cᴷc⁻, :Cᴷu⁻, :Cᴷe⁻, :CRiᶜ, :CRiʷ),
)

conv_adj_names = (:Cᶜc, :Cᶜe, :Cᵉc, :Cᵉe, :Cᶜˢ)

for (set, names) in parameter_sets
    conv_adj_set = set * "_conv_adj"
    parameter_sets[conv_adj_set] =  tuple(names..., conv_adj_names...)
end

# Some dependent parameters
dependent_parameter_sets = Dict()
for (set, names) in parameter_sets
    dependent_parameter_sets[set] = NamedTuple()
end

Cᴰ⁻(θ)  = θ.Cᴰ⁺
Cᴷu⁻(θ) = θ.Cᴷu⁺
Cᴷc⁻(θ) = θ.Cᴷc⁺
Cᴷe⁻(θ) = θ.Cᴷe⁺

dependent_parameter_sets["constant_Pr"]          = (; Cᴷu⁻, Cᴷc⁻, Cᴷe⁻, Cᴰ⁻) 
dependent_parameter_sets["constant_Pr_conv_adj"] = (; Cᴷu⁻, Cᴷc⁻, Cᴷe⁻, Cᴰ⁻) 

#####
##### Bounds and priors
#####

bounds_library = Dict()

# Ri-based
bounds_library[:ν₀]  = (0.0,  1.0)
bounds_library[:κ₀]  = (0.0,  1.0)
bounds_library[:κᶜ]  = (0.0, 10.0)
bounds_library[:Cᵉ]  = (0.0, 10.0)
bounds_library[:Ri₀] = (0.0,  4.0)
bounds_library[:Riᵟ] = (0.0,  4.0)

# Turbulent kinetic energy parameters
bounds_library[:CᵂwΔ] = (0.1, 10.0)
bounds_library[:Cᵂu★] = (0.1, 10.0)
bounds_library[:Cᴰ⁻]  = (0.1,  1.0)
bounds_library[:Cᴰ⁺]  = (0.1,  1.0)

# Mixing length parameters
bounds_library[:Cᵇ]   = (0.1, 1.0)
bounds_library[:Cˢ]   = (0.1, 4.0)

bounds_library[:Cᴷu⁻] = (0.1, 1.0)
bounds_library[:Cᴷc⁻] = (0.1, 1.0)
bounds_library[:Cᴷe⁻] = (1.0, 4.0)

bounds_library[:Cᴷu⁺] = (0.1, 1.0)
bounds_library[:Cᴷc⁺] = (0.1, 1.0)
bounds_library[:Cᴷe⁺] = (1.0, 4.0)

bounds_library[:CRiᶜ] = (0.1, 2.0)
bounds_library[:CRiʷ] = (0.1, 2.0)

# Convective adjustment parameters
bounds_library[:Cᶜc]  = (0.1, 2.0)
bounds_library[:Cᵉc]  = (0.1, 2.0)
bounds_library[:Cᶜe]  = (0.1, 2.0)
bounds_library[:Cᵉe]  = (0.1, 2.0)
bounds_library[:Cᶜˢ]  = (0.1, 5.0)

#bounds_library[:Cᶜe]  = ( 1.0,  30.0)
#bounds_library[:Cᵉe]  = ( 1.0,  30.0)

prior_library = Dict()

for p in keys(bounds_library)
    prior_library[p] = ScaledLogitNormal(; bounds=bounds_library[p])
end

get_free_parameters(name) = FreeParameters(prior_library, names = parameter_sets[name],
                                           dependent_parameters = dependent_parameter_sets[name])

