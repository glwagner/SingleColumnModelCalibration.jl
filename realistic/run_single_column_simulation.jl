using Oceananigans

include("single_column_omip_simulation.jl")

# locations:
#    :eastern_mediterranean
#    :ocean_station_papa
#    :north_atlantic
#    :drake_passage
#    :weddell_sea
#    :tasman_southern_ocean

simulation = single_column_simulation(:north_atlantic)

jld2_output_writer = simulation.output_writers[:jld2]
jld2_filepath = jld2_output_writer.filepath

@info "Writing data to $jld2_filepath"

run!(simulation)

