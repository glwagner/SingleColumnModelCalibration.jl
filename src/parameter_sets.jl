parameter_sets = Dict(
    "ri_based"    => (:ν₀, :κ₀, :κᶜ, :Cᵉ, :Ri₀, :Riᵟ),
    "constant_Pr" => (:CᵂwΔ, :Cᵂu★, :Cᴷc⁺, :Cᴷu⁺, :Cᴷe⁺, :Cᵇc, :Cˢc, :Cᴰ⁺),
    "variable_Pr" => (:CᵂwΔ, :Cᵂu★, :Cᴷc⁺, :Cᴷu⁺, :Cᴷe⁺, :Cᵇc, :Cˢc, :Cᴰ⁺, :Cᴰ⁻, :Cᴷc⁻, :Cᴷu⁻, :Cᴷe⁻, :CRiᶜ, :CRiʷ),
    "complex"     => (:CᵂwΔ, :Cᵂu★, :Cᴷc⁺, :Cᴷu⁺, :Cᴷe⁺, :Cᵇc, :Cˢc, :Cᴰ⁺, :Cᴰ⁻, :Cᴷc⁻, :Cᴷu⁻, :Cᴷe⁻, :CRiᶜ, :CRiʷ, :Cᵇu, :Cᵇe, :Cˢu, :Cˢe),
)

conv_adj_names = (:Cᴬc, :Cᴬe, :Cʰˢ)

for (set, names) in parameter_sets
    conv_adj_set = set * "_conv_adj"
    parameter_sets[conv_adj_set] =  tuple(names..., conv_adj_names...)
end

# Some dependent parameters
dependent_parameter_sets = Dict()
for (set, names) in parameter_sets
    dependent_parameter_sets[set] = NamedTuple()
end

Cᵟu(θ)  = θ.Cᵟc
Cᵟe(θ)  = θ.Cᵟc
Cᵇ(θ)   = θ.Cᵇc
Cˢ(θ)   = θ.Cˢc
Cᴷu⁻(θ) = θ.Cᴷu⁺
Cᴷc⁻(θ) = θ.Cᴷc⁺
Cᴷe⁻(θ) = θ.Cᴷe⁺
Cᴰ⁻(θ)  = θ.Cᴰ⁺

dependent_parameter_sets["constant_Pr"]              = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ, Cᴷu⁻, Cᴷc⁻, Cᴷe⁻, Cᴰ⁻) 
dependent_parameter_sets["constant_Pr_conv_adj"]     = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ, Cᴷu⁻, Cᴷc⁻, Cᴷe⁻, Cᴰ⁻) 
dependent_parameter_sets["variable_Pr"]              = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ)
dependent_parameter_sets["variable_Pr_conv_adj"]     = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ)

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
bounds_library[:CᵂwΔ] = ( 1.0,  10.0)
bounds_library[:Cᵂu★] = ( 1.0,  10.0)
bounds_library[:Cᴰ⁻]  = ( 1.0,  4.0)
bounds_library[:Cᴰ⁺]  = ( 1.0,  4.0)

# Mixing length parameters
bounds_library[:Cᴷu⁻] = ( 0.1, 2.0)
bounds_library[:Cᴷc⁻] = ( 0.1, 2.0)
bounds_library[:Cᴷe⁻] = ( 0.1, 4.0)

bounds_library[:Cᴷu⁺] = ( 0.1, 2.0)
bounds_library[:Cᴷc⁺] = ( 0.1, 4.0)
bounds_library[:Cᴷe⁺] = ( 0.1, 4.0)

bounds_library[:CRiᶜ] = ( 0.1,  2.0)
bounds_library[:CRiʷ] = ( 0.1,  2.0)

bounds_library[:Cᵇu]  = ( 0.1, 2.0)
bounds_library[:Cᵇc]  = ( 0.1, 2.0)
bounds_library[:Cᵇe]  = ( 0.1, 2.0)

bounds_library[:Cˢu]  = ( 0.1, 4.0)
bounds_library[:Cˢc]  = ( 0.1, 4.0)
bounds_library[:Cˢe]  = ( 0.1, 4.0)

bounds_library[:Cᴬu]  = ( 0.0,  10.0)
bounds_library[:Cᴬc]  = ( 0.1,   2.0)
bounds_library[:Cᴬe]  = ( 1.0,  30.0)

bounds_library[:Cʰ]   = ( 0.0,   0.1)
bounds_library[:Cʰˢ]  = ( 0.1,   5.0)

prior_library = Dict()

for p in keys(bounds_library)
    prior_library[p] = ScaledLogitNormal(; bounds=bounds_library[p])
end

get_free_parameters(name) = FreeParameters(prior_library, names = parameter_sets[name],
                                           dependent_parameters = dependent_parameter_sets[name])

