parameter_sets = Dict(
    "ri_based"                    => (:ν₀, :κ₀, :κᶜ, :Cᵉ, :Ri₀, :Riᵟ),
    "nemo_like"                   => (:CᵂwΔ, :Cᵂu★, :Cᴷc⁺, :Cᴷu⁺, :Cᴷe⁺, :Cᵇc, :Cᴰ⁺),
    "ri_dependent"                => (:CᵂwΔ, :Cᵂu★, :Cᴷc⁺, :Cᴷu⁺, :Cᴷe⁺, :Cᵇc, :Cˢc, :Cᴰ⁺, :Cᴷc⁻, :Cᴷu⁻, :Cᴷe⁻, :CᴷRiᶜ, :CᴷRiʷ),
    "ri_dependent_dissipation"    => (:CᵂwΔ, :Cᵂu★, :Cᴷc⁺, :Cᴷu⁺, :Cᴷe⁺, :Cᵇc, :Cˢc, :Cᴰ⁺, :Cᴷc⁻, :Cᴷu⁻, :Cᴷe⁻, :CᴷRiᶜ, :CᴷRiʷ, :Cᴰ⁻, :CᴰRiᶜ, :CᴰRiʷ),
    "complex"                     => (:CᵂwΔ, :Cᵂu★, :Cᴷc⁺, :Cᴷu⁺, :Cᴷe⁺, :Cᵇc, :Cᵇu, :Cᵇe, :Cˢc, :Cˢu, :Cˢe, :Cᴰ⁺, :Cᴷc⁻, :Cᴷu⁻, :Cᴷe⁻, :CᴷRiᶜ, :CᴷRiʷ),
    "complex_dissipation"         => (:CᵂwΔ, :Cᵂu★, :Cᴷc⁺, :Cᴷu⁺, :Cᴷe⁺, :Cᵇc, :Cᵇu, :Cᵇe, :Cˢc, :Cˢu, :Cˢe, :Cᴰ⁺, :Cᴷc⁻, :Cᴷu⁻, :Cᴷe⁻, :CᴷRiᶜ, :CᴷRiʷ, :Cᴰ⁻, :CᴰRiᶜ, :CᴰRiʷ),
)

parameter_sets["shear_nemo_like"] = tuple(parameter_sets["nemo_like"]..., :Cˢc)

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

Cᵟu(θ) = θ.Cᵟc
Cᵟe(θ) = θ.Cᵟc
Cᵇ(θ) = θ.Cᵇc
Cˢ(θ) = θ.Cˢc
Cᴷu⁻(θ) = θ.Cᴷu⁺
Cᴷc⁻(θ) = θ.Cᴷc⁺
Cᴷe⁻(θ) = θ.Cᴷe⁺
Cᴰ⁻(θ) = θ.Cᴰ⁺

dependent_parameter_sets["basic"]             = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ, Cᴰ⁻)
dependent_parameter_sets["basic_conv_adj"]    = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ, Cᴰ⁻)

dependent_parameter_sets["basic_dissipation_length"]             = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ)
dependent_parameter_sets["basic_dissipation_length_conv_adj"]    = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ)

dependent_parameter_sets["nemo_like"]             = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cᴷu⁻, Cᴷc⁻, Cᴷe⁻, Cᴰ⁻) 
dependent_parameter_sets["nemo_like_conv_adj"]    = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cᴷu⁻, Cᴷc⁻, Cᴷe⁻, Cᴰ⁻) 

dependent_parameter_sets["ri_dependent"]             = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ, Cᴰ⁻) 
dependent_parameter_sets["ri_dependent_conv_adj"]    = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ, Cᴰ⁻) 
dependent_parameter_sets["ri_dependent_dissipation"]           = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ)
dependent_parameter_sets["ri_dependent_dissipation_conv_adj"]  = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ)

dependent_parameter_sets["shear_nemo_like"]             = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ, Cᴷu⁻, Cᴷc⁻, Cᴷe⁻, Cᴰ⁻) 
dependent_parameter_sets["shear_nemo_like_conv_adj"]    = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ, Cᴷu⁻, Cᴷc⁻, Cᴷe⁻, Cᴰ⁻) 

dependent_parameter_sets["complex_dissipation_length"]             = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ, Cᴷu⁻, Cᴷc⁻, Cᴷe⁻)
dependent_parameter_sets["complex_dissipation_length_conv_adj"]    = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ, Cᴷu⁻, Cᴷc⁻, Cᴷe⁻)

dependent_parameter_sets["goldilocks"]             = (; Cᴷu⁻, Cᴷc⁻, Cᴷe⁻) 
dependent_parameter_sets["goldilocks_conv_adj"]    = (; Cᴷu⁻, Cᴷc⁻, Cᴷe⁻) 

dependent_parameter_sets["goldilocks"]             = (; Cᴰ⁻) 
dependent_parameter_sets["goldilocks_conv_adj"]    = (; Cᴰ⁻) 

#####
##### Bounds and priors
#####

bounds_library = Dict()

# Ri-based
bounds_library[:ν₀]  = (0.0,  1.0)
bounds_library[:κ₀]  = (0.0,  1.0)
bounds_library[:κᶜ]  = (0.0, 10.0)
bounds_library[:Cᵉ]  = (0.0, 10.0)
bounds_library[:Ri₀] = (0.0,  1.0)
bounds_library[:Riᵟ] = (0.0,  1.0)

# Turbulent kinetic energy parameters
bounds_library[:CᵂwΔ]  = ( 0.0,  10.0)
bounds_library[:Cᵂu★]  = ( 0.0,  10.0)
bounds_library[:Cᴰ⁻]   = ( 1.0,  6.0)
bounds_library[:Cᴰ⁺]   = ( 1.0,  6.0)
bounds_library[:CᴰRiᶜ] = ( 0.1,  1.0)
bounds_library[:CᴰRiʷ] = ( 0.1,  1.0)

# Mixing length parameters
bounds_library[:Cᴷu⁻]  = ( 0.1, 1.0)
bounds_library[:Cᴷc⁻]  = ( 0.1, 4.0)
bounds_library[:Cᴷe⁻]  = ( 0.1, 4.0)

bounds_library[:Cᴷu⁺]  = ( 0.1, 2.0)
bounds_library[:Cᴷc⁺]  = ( 0.1, 4.0)
bounds_library[:Cᴷe⁺]  = ( 0.1, 4.0)

bounds_library[:CᴷRiᶜ] = ( 0.1, 1.0)
bounds_library[:CᴷRiʷ] = ( 0.1, 1.0)

bounds_library[:Cᵇu]   = ( 0.1, 2.0)
bounds_library[:Cᵇc]   = ( 0.1, 2.0)
bounds_library[:Cᵇe]   = ( 0.1, 4.0)

bounds_library[:Cˢu]   = ( 0.1, 4.0)
bounds_library[:Cˢc]   = ( 0.1, 4.0)
bounds_library[:Cˢe]   = ( 0.1, 4.0)

bounds_library[:Cᵟu]   = ( 0.0, 4.0)
bounds_library[:Cᵟc]   = ( 0.0, 4.0)
bounds_library[:Cᵟe]   = ( 0.0, 4.0)

bounds_library[:Cᴬu]   = ( 0.0,  10.0)
bounds_library[:Cᴬc]   = ( 0.0,   5.0)
bounds_library[:Cᴬe]   = ( 0.0,  30.0)

bounds_library[:Cʰ]    = ( 0.0,   0.1)
bounds_library[:Cʰˢ]   = ( 0.0,   4.0)

bounds_library[:Cʷ★]   = ( 1.0,  4.0)
bounds_library[:Cʷℓ]   = ( 0.0,  4.0)

prior_library = Dict()

for p in keys(bounds_library)
    prior_library[p] = ScaledLogitNormal(; bounds=bounds_library[p])
end

get_free_parameters(name) = FreeParameters(prior_library, names = parameter_sets[name],
                                           dependent_parameters = dependent_parameter_sets[name])

