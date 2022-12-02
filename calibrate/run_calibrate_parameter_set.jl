using Oceananigans
using Oceananigans.Units
using JLD2

using Oceananigans.TurbulenceClosures:
    RiBasedVerticalDiffusivity,
    CATKEVerticalDiffusivity

using ParameterEstimocean: iterate!

using SingleColumnModelCalibration:
    build_ensemble_kalman_inversion,
    generate_filepath,
    parameter_sets

#name = "ri_based"
#closure = RiBasedVerticalDiffusivity()

names = [
    "constant_Pr",
    "variable_Pr",
    "constant_Pr_conv_adj",
    "variable_Pr_conv_adj",
]

grid_parameters = [
    (size=32, z=(-256, 0)),
    (size=64, z=(-256, 0)),
    #(size=128, z=(-256, 0)),
]

suite_parameters = [
    (name = "12_hour_suite", stop_time=12hours),
    (name = "24_hour_suite", stop_time=24hours),
    (name = "48_hour_suite", stop_time=48hours),
]

resultsdir = "../results"

closure = CATKEVerticalDiffusivity()
architecture = CPU()
Nensemble = 100
resample_failure_fraction = 0.05
stop_pseudotime = 1000.0
Δt = 20minutes
Nrepeats = 10

for name in names
    start_time = time_ns()
    
    repeat_ekis = []
    for i = 1:Nrepeats
        eki = build_ensemble_kalman_inversion(closure, name;
                                              architecture,
                                              Nensemble,
                                              Δt,
                                              grid_parameters,
                                              suite_parameters,
                                              resample_failure_fraction)
        push!(repeat_ekis, eki)
    end

    asyncmap(1:Nrepeats) do i
        eki = repeat_ekis[i]

        logname = string(name, "_", i, ".txt")

        while eki.pseudotime < stop_pseudotime
            iterate!(eki)

            if eki.iteration % 10 == 0
                open(logname, "a") do io
                    show(io, "text/plain", eki.iteration_summaries[end])
                    write(io, '\n')
                    write(io, '\n')
                end

                @show eki.iteration_summaries[end]
            end
        end

        filename = string(name, "_", i)
        filepath = generate_filepath(; Δt,
                                     dir = resultsdir,
                                     suite_parameters,
                                     grid_parameters,
                                     stop_pseudotime,
                                     Nensemble,
                                     filename)


        rm(filepath; force=true)

        @info "Saving data to $filepath..."
        file = jldopen(filepath, "a+")
        file["resample_failure_fraction"] = resample_failure_fraction
        file["stop_pseudotime"] = stop_pseudotime
        file["iteration_summaries"] = eki.iteration_summaries
        close(file)
    end

    elapsed = 1e-9 * (time_ns() - start_time)
    
    @info "Calibrating $name parameters took " * prettytime(elapsed)
end

