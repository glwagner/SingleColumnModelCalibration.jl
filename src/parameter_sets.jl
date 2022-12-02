parameter_sets = Dict(
    "ri_based"    => (:ν₀, :κ₀, :κᶜ, :Cᵉ, :Ri₀, :Riᵟ),
    "constant_Pr" => (:CᵂwΔ, :Cᵂu★, :C⁺c, :C⁺u, :C⁺e, :C⁺D, :Cᵇ, :Cˢ),
    "constant_Pr_no_shear" => (:CᵂwΔ, :Cᵂu★, :C⁺c, :C⁺u, :C⁺e, :C⁺D, :Cᵇ),
    "variable_Pr" => (:CᵂwΔ, :Cᵂu★, :C⁺c, :C⁺u, :C⁺e, :C⁺D, :Cᵇ, #:Cˢ, 
                                    :C⁻c, :C⁻u, :C⁻e, :C⁻D, :CRiᶜ, :CRiʷ),
)

conv_adj_names = (:Cᶜc, :Cᶜe, :CᶜD, :Cᵉc, :Cᵉe, :CᵉD, :Cˢᶜ)

for (set, names) in parameter_sets
    conv_adj_set = set * "_conv_adj"
    parameter_sets[conv_adj_set] =  tuple(names..., conv_adj_names...)
end

# Some dependent parameters
dependent_parameter_sets = Dict()
for (set, names) in parameter_sets
    dependent_parameter_sets[set] = NamedTuple()
end

C⁻D(θ)  = θ.C⁺D
C⁻u(θ) = θ.C⁺u
C⁻c(θ) = θ.C⁺c
C⁻e(θ) = θ.C⁺e

dependent_parameter_sets["constant_Pr"]          = (; C⁻u, C⁻c, C⁻e, C⁻D) 
dependent_parameter_sets["constant_Pr_conv_adj"] = (; C⁻u, C⁻c, C⁻e, C⁻D) 

#####
##### Bounds and priors
#####

bounds_library = Dict()

# Turbulent kinetic energy parameters
bounds_library[:CᵂwΔ] = (0.1, 10.0)
bounds_library[:Cᵂu★] = (0.1, 10.0)
bounds_library[:C⁻D]  = (0.1, 10.0)
bounds_library[:C⁺D]  = (0.1, 10.0)

# Mixing length parameters
bounds_library[:Cᵇ]   = (0.2, 2.0)
bounds_library[:Cˢ]   = (0.2, 2.0)

bounds_library[:C⁻u] = (0.1, 1.0)
bounds_library[:C⁺u] = (0.1, 1.0)
bounds_library[:C⁻c] = (0.1, 1.0)
bounds_library[:C⁺c] = (0.1, 1.0)
bounds_library[:C⁻e] = (0.1, 10.0)
bounds_library[:C⁺e] = (0.1, 10.0)

bounds_library[:CRiᶜ] = (0.1, 1.0)
bounds_library[:CRiʷ] = (0.1, 2.0)

# Convective adjustment parameters
bounds_library[:Cᶜc]  = (0.01, 10.0)
bounds_library[:Cᶜe]  = (0.01, 10.0)
bounds_library[:CᶜD]  = (0.01, 10.0)
bounds_library[:Cᵉc]  = (0.01,  2.0)
bounds_library[:Cᵉe]  = (0.01,  2.0)
bounds_library[:CᵉD]  = (0.01,  2.0)
bounds_library[:Cˢᶜ]  = (0.01, 10.0)

# Ri-based
bounds_library[:ν₀]  = (0.0,  1.0)
bounds_library[:κ₀]  = (0.0,  1.0)
bounds_library[:κᶜ]  = (0.0, 10.0)
bounds_library[:Cᵉ]  = (0.0, 10.0)
bounds_library[:Ri₀] = (0.0,  4.0)
bounds_library[:Riᵟ] = (0.0,  4.0)


prior_library = Dict()

for p in keys(bounds_library)
    prior_library[p] = ScaledLogitNormal(; bounds=bounds_library[p])
end

get_free_parameters(name) = FreeParameters(prior_library, names = parameter_sets[name],
                                           dependent_parameters = dependent_parameter_sets[name])

