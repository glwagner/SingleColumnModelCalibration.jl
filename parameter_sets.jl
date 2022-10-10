parameter_sets = Dict(
    "goldilocks" => (:CᵂwΔ, :Cᵂu★, :Cᴷc⁻, :Cᴷu⁻, :Cᴷe⁻, :Cᴰ⁻, :Cᴰʳ, :CᴰRiᶜ, :CᴰRiʷ, :Cᵇc, :Cᵇu, :Cᵇe, :Cˢc, :Cˢu, :Cˢe),
    "nemo_like"  => (:CᵂwΔ, :Cᵂu★, :Cᴷc⁻, :Cᴷu⁻, :Cᴷe⁻, :Cᴰ⁻, :Cᵇc),
)

parameter_sets["complex"] = tuple(parameter_sets["goldilocks"]..., :CᴷRiᶜ, :CᴷRiʷ)
parameter_sets["shear_nemo_like"] = tuple(parameter_sets["nemo_like"]..., :Cˢc)

conv_adj_names = (:Cᴬc, :Cᴬu, :Cᴬe)
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

dependent_parameter_sets["nemo_like"]             = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ) 
dependent_parameter_sets["nemo_like_conv_adj"]    = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ) 
dependent_parameter_sets["nemo_like_grid_length"] = (; Cᵇu=Cᵇ, Cᵇe=Cᵇ) 

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
    :Cᴷcʳ => 0.0,
    :Cᴷuʳ => 0.0,
    :Cᴷeʳ => 0.0,
)

neutral_default_tke_parameters = Dict(
    :CᵂwΔ => 0.0,
    :Cᵂu★ => 0.0,
    :Cᴰ⁻ => 0.0,
    :Cᴰʳ => 0.0,
    :CᴰRiᶜ => 0.0,
    :CᴰRiʷ => 0.0,
)

mixing_length = MixingLength(; neutral_default_mixing_length_parameters...)
turbulent_kinetic_energy_equation = TurbulentKineticEnergyEquation(; neutral_default_tke_parameters...)
closure = CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation)

