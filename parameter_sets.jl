parameter_sets = Dict(
    #"goldilocks" => (:Cʷ★, :Cʷℓ, :CᵂwΔ, :Cᵂu★, :Cᴷc⁻, :Cᴷu⁻, :Cᴷe⁻, :CᴷRiᶜ, :CᴷRiʷ, :Cᴰ⁻, :Cᵇc, :Cᵇu, :Cᵇe, :Cˢc, :Cˢu, :Cˢe),
    #"nemo_like"  => (:Cʷ★, :Cʷℓ, :CᵂwΔ, :Cᵂu★, :Cᴷc⁻, :Cᴷu⁻, :Cᴷe⁻, :Cᴰ⁻, :Cᵇc),
    "goldilocks" => (:CᵂwΔ, :Cᵂu★, :Cᴷc⁻, :Cᴷu⁻, :Cᴷe⁻, :Cᴷc⁺, :Cᴷu⁺, :Cᴷe⁺, :CᴷRiᶜ, :CᴷRiʷ, :Cᴰ⁻, :Cᵇc, :Cᵇu, :Cᵇe, :Cˢc, :Cˢu, :Cˢe),
    "nemo_like"  => (:CᵂwΔ, :Cᵂu★, :Cᴷc⁻, :Cᴷu⁻, :Cᴷe⁻, :Cᴰ⁻, :Cᵇc),
    "basic"  => (:CᵂwΔ, :Cᵂu★, :Cᴷc⁻, :Cᴷu⁻, :Cᴷc⁺, :Cᴷu⁺, :Cᴷe⁺, :CᴷRiᶜ, :CᴷRiʷ, :Cᴷe⁻, :Cᴰ⁻, :Cᵇc, :Cˢc),
)

parameter_sets["complex"] = tuple(parameter_sets["goldilocks"]..., :Cᴰ⁺, :CᴰRiᶜ, :CᴰRiʷ)
parameter_sets["shear_nemo_like"] = tuple(parameter_sets["nemo_like"]..., :Cˢc)

#conv_adj_names = (:Cᴬc, :Cᴬu, :Cᴬe, :Cʰˢ)
conv_adj_names = (:Cᴬc, :Cᴬe, :Cʰˢ)
grid_length_names = (:Cᵟc, :Cᵟu, :Cᵟe)

for (set, names) in parameter_sets
    conv_adj_set = set * "_conv_adj"
    parameter_sets[conv_adj_set] =  tuple(names..., conv_adj_names...)

    grid_length_set = set * "_grid_length"
    parameter_sets[grid_length_set] =  tuple(names..., grid_length_names...)
end

# Some dependent parameters
dependent_parameter_sets = Dict()
for (set, names) in parameter_sets
    dependent_parameter_sets[set] = NamedTuple()
end

Cᵇ(θ) = θ.Cᵇc
Cˢ(θ) = θ.Cˢc
Cᴷu⁺(θ) = θ.Cᴷu⁻
Cᴷc⁺(θ) = θ.Cᴷc⁻
Cᴷe⁺(θ) = θ.Cᴷe⁻
Cᴰ⁺(θ) = θ.Cᴰ⁻

dependent_parameter_sets["basic"]             = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ)
dependent_parameter_sets["basic_conv_adj"]    = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ)
dependent_parameter_sets["basic_grid_length"] = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ)

dependent_parameter_sets["nemo_like"]             = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ, Cᴷu⁺, Cᴷc⁺, Cᴷe⁺, Cᴰ⁺) 
dependent_parameter_sets["nemo_like_conv_adj"]    = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ, Cᴷu⁺, Cᴷc⁺, Cᴷe⁺, Cᴰ⁺) 
dependent_parameter_sets["nemo_like_grid_length"] = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cˢu=Cˢ, Cˢe=Cˢ, Cᴷu⁺, Cᴷc⁺, Cᴷe⁺, Cᴰ⁺) 

dependent_parameter_sets["shear_nemo_like"]             = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cᴷu⁺, Cᴷc⁺, Cᴷe⁺, Cᴰ⁺) 
dependent_parameter_sets["shear_nemo_like_conv_adj"]    = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cᴷu⁺, Cᴷc⁺, Cᴷe⁺, Cᴰ⁺) 
dependent_parameter_sets["shear_nemo_like_grid_length"] = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ, Cᴷu⁺, Cᴷc⁺, Cᴷe⁺, Cᴰ⁺) 

dependent_parameter_sets["goldilocks"]             = (; Cᴷu⁺, Cᴷc⁺, Cᴷe⁺) 
dependent_parameter_sets["goldilocks_conv_adj"]    = (; Cᴷu⁺, Cᴷc⁺, Cᴷe⁺) 
dependent_parameter_sets["goldilocks_grid_length"] = (; Cᴷu⁺, Cᴷc⁺, Cᴷe⁺) 

dependent_parameter_sets["goldilocks2"]             = (; Cᴰ⁺) 
dependent_parameter_sets["goldilocks2_conv_adj"]    = (; Cᴰ⁺) 
dependent_parameter_sets["goldilocks2_grid_length"] = (; Cᴰ⁺) 

# Some defaults
neutral_default_mixing_length_parameters = Dict(
    :Cᵟc => 0.5,
    :Cᵟu => 0.5,
    :Cᵟe => 0.5,
    :Cᵇu => Inf,
    :Cᵇc => Inf,
    :Cᵇe => Inf,
    :Cˢu => Inf,
    :Cˢc => Inf,
    :Cˢe => Inf,
    :CᴷRiᶜ => Inf,
    :CᴷRiʷ => 0.0,
)

neutral_default_tke_parameters = Dict(
    :CᵂwΔ => 0.0,
    :Cᵂu★ => 0.0,
    :Cᴰ⁻ => 0.0,
    :Cᴰ⁺ => 0.0,
    :CᴰRiᶜ => Inf,
    :CᴰRiʷ => 0.0,
)

