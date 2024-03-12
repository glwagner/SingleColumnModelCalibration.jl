using Oceananigans

include("single_column_omip_simulation.jl")

# locations:
#    :eastern_mediterranean
#    :ocean_station_papa
#    :north_atlantic
#    :drake_passage
#    :weddell_sea
#    :tasman_southern_ocean

#for location in (:north_atlantic, :ocean_station_papa)
for location in (:ocean_station_papa,)
    simulation = single_column_simulation(location)
    jld2_output_writer = simulation.output_writers[:jld2]
    jld2_filepath = jld2_output_writer.filepath

    @info "Writing data to $jld2_filepath"

    run!(simulation)
end

