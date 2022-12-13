parameter_sets = Dict(
    "ri_based"             => (:ν₀, :κ₀, :κᶜ, :Cᵉ, :Ri₀, :Riᵟ),
    "constant_Pr"          => (:CᵂwΔ, :Cᵂu★, :C⁺c, :C⁺u, :C⁺e, :C⁺D, :Cᵇ, :Cˢ),
    "constant_Pr_no_shear" => (:CᵂwΔ, :Cᵂu★, :C⁺c, :C⁺u, :C⁺e, :C⁺D, :Cᵇ),
    "variable_Pr"          => (:CᵂwΔ, :Cᵂu★, :C⁺c, :C⁺u, :C⁺e, :C⁺D, :Cᵇ, :C⁻c, :C⁻u, :C⁻e, :C⁻D, :CRiᶜ, :CRiʷ),
)

conv_adj_names = (:Cᶜc, :Cᶜe, :CᶜD, :Cᵉc, :Cᵉe, :CᵉD, :Cˢᶜ)
tracer_conv_adj_names = (:Cᶜc, :Cᵉc, :Cˢᶜ)
simple_conv_adj_names = (:Cᶜc, :Cᵉc, :CᶜD, :Cˢᶜ)

sets = ["constant_Pr",
        "constant_Pr_no_shear",
        "variable_Pr"]

for set in sets
    names = parameter_sets[set]
        
    conv_adj_set = set * "_conv_adj"
    parameter_sets[conv_adj_set] =  tuple(names..., conv_adj_names...)

    tracer_conv_adj_set = set * "_tracer_conv_adj"
    parameter_sets[tracer_conv_adj_set] =  tuple(names..., tracer_conv_adj_names...)

    simple_conv_adj_set = set * "_simple_conv_adj"
    parameter_sets[simple_conv_adj_set] =  tuple(names..., simple_conv_adj_names...)
end

# Some dependent parameters
dependent_parameter_sets = Dict()
for (set, names) in parameter_sets
    dependent_parameter_sets[set] = NamedTuple()
end

C⁻D(θ) = θ.C⁺D
C⁻u(θ) = θ.C⁺u
C⁻c(θ) = θ.C⁺c
C⁻e(θ) = θ.C⁺e

for name in ["constant_Pr",
             "constant_Pr_no_shear",
             "constant_Pr_no_shear_simple_conv_adj",
             "constant_Pr_simple_conv_adj",
             "constant_Pr_no_shear_tracer_conv_adj",
             "constant_Pr_tracer_conv_adj",
             "constant_Pr_no_shear_conv_adj",
             "constant_Pr_conv_adj"]

    dependent_parameter_sets[name] = (; C⁻u, C⁻c, C⁻e, C⁻D) 
end

#####
##### Bounds and priors
#####

bounds_library = Dict()

# Turbulent kinetic energy parameters
bounds_library[:CᵂwΔ] = (0.0, 10.0)
bounds_library[:Cᵂu★] = (0.0, 10.0)
bounds_library[:C⁻D]  = (0.0, 10.0)
bounds_library[:C⁺D]  = (0.0, 10.0)

# Mixing length parameters
bounds_library[:Cᵇ]   = (0.0, 2.0)
bounds_library[:Cˢ]   = (0.0, 2.0)

bounds_library[:C⁻u] = (0.0, 1.0)
bounds_library[:C⁺u] = (0.0, 1.0)
bounds_library[:C⁻c] = (0.0, 1.0)
bounds_library[:C⁺c] = (0.0, 1.0)
bounds_library[:C⁻e] = (0.0, 10.0)
bounds_library[:C⁺e] = (0.0, 10.0)

bounds_library[:CRiᶜ] = (0.0, 1.0)
bounds_library[:CRiʷ] = (0.0, 2.0)

# Convective adjustment parameters
bounds_library[:Cᶜc]  = (0.0, 10.0)
bounds_library[:Cᶜe]  = (0.0, 10.0)
bounds_library[:CᶜD]  = (0.0, 10.0)
bounds_library[:Cᵉc]  = (0.0, 10.0)
bounds_library[:Cᵉe]  = (0.0, 10.0)
bounds_library[:CᵉD]  = (0.0, 10.0)
bounds_library[:Cˢᶜ]  = (0.0, 10.0)

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

