using Oceananigans
using Oceananigans.Units
using Oceananigans.Simulations: reset!
using ParameterEstimocean
using ParameterEstimocean: Transformation
using ParameterEstimocean.InverseProblems: initialize_forward_run!
using LinearAlgebra
using DataDeps
using Distributions
using GLMakie

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity

#####
##### Compile LESbrary
#####

field_names = (:b, :e, :u, :v)
Δt = 20minutes
closure = CATKEVerticalDiffusivity()

z = [-256.0,
     -238.3,
     -207.2,
     -182.3,
     -162.5,
     -146.5,
     -133.0,
     -121.3,
     -110.5,
     -100.2,
     - 90.0,
     - 80.0,
     - 70.0,
     - 60.0,
     - 50.0,
     - 40.0,
     - 30.0,
     - 20.0,
     - 10.0,
        0.0]

regrid = RectilinearGrid(size=length(z)-1; z, topology=(Flat, Flat, Bounded))
#regrid = RectilinearGrid(size=128; z=(-256, 0), topology=(Flat, Flat, Bounded))

cases = ["free_convection",
         "strong_wind_weak_cooling",
         "med_wind_med_cooling",
         "weak_wind_strong_cooling",
         "strong_wind",
         "strong_wind_no_rotation"]

suite = "two_day_suite"
case_path(case) = joinpath("data", suite, case * "_instantaneous_statistics.jld2")
test_obs = SyntheticObservations(case_path("free_convection"); field_names, regrid)
times = test_obs.times[100:end]

observation_library = Dict()

# The rest are standard
for case in cases #["strong_wind", "med_wind_med_cooling", "strong_wind_weak_cooling", "weak_wind_strong_cooling"]
    observation_library[case] = SyntheticObservations(case_path(case); field_names, regrid, times)
end

observations = [observation_library[case] for case in cases]
     
#####
##### Calibration
#####

bounds_library = Dict()

# Turbulent kinetic energy parameters
bounds_library[:CᵂwΔ]  = ( 2.0,  20.0)
bounds_library[:Cᵂu★]  = ( 2.0,  20.0)
bounds_library[:Cᴰ⁻]   = ( 0.0,   2.0)
bounds_library[:Cᴰʳ]   = (-1.0,  10.0)
bounds_library[:CᴰRiᶜ] = (-4.0,   4.0)
bounds_library[:CᴰRiʷ] = ( 0.0,   4.0)

# Mixing length parameters
#
#   Recall σ = σ⁻ (1 + σʳ * step(x, c, w))
#
bounds_library[:Cᴷu⁻]  = ( 0.0,  10.0)
bounds_library[:Cᴷc⁻]  = ( 0.0,  10.0)
bounds_library[:Cᴷe⁻]  = ( 0.0,  10.0)
bounds_library[:Cᴷuʳ]  = (-1.0,   2.0)
bounds_library[:Cᴷcʳ]  = (-1.0,   2.0)
bounds_library[:Cᴷeʳ]  = (-1.0,   2.0)
bounds_library[:CᴷRiᶜ] = (-4.0,   4.0)
bounds_library[:CᴷRiʷ] = ( 0.0,   4.0)
bounds_library[:Cᵇu]   = ( 0.0,   4.0)
bounds_library[:Cᵇc]   = ( 0.0,   4.0)
bounds_library[:Cᵇe]   = ( 0.0,   4.0)
bounds_library[:Cˢu]   = ( 0.0,   4.0)
bounds_library[:Cˢc]   = ( 0.0,   4.0)
bounds_library[:Cˢe]   = ( 0.0,   4.0)

bounds_library[:Cᴬˢu]  = ( 0.0,  2.0)
bounds_library[:Cᴬˢc]  = ( 0.0,  2.0)
bounds_library[:Cᴬˢe]  = ( 0.0,  2.0)
bounds_library[:Cᴬu]   = ( 0.0,  0.01)
bounds_library[:Cᴬc]   = (1e-3,  0.1)
bounds_library[:Cᴬe]   = ( 0.0,  0.01)

# Extras
bounds_library[:Cᵇ]    = ( 0.0,   4.0)
bounds_library[:Cˢ]    = ( 0.0,   4.0)
bounds_library[:Cᴰ]    = ( 0.0,   2.0)
bounds_library[:Cᵟu]   = ( 0.0,  10.0)
bounds_library[:Cᵟc]   = ( 0.0,  10.0)
bounds_library[:Cᵟe]   = ( 0.0,  10.0)

prior_library = Dict()

for p in keys(bounds_library)
    prior_library[p] = ScaledLogitNormal(; bounds=bounds_library[p])
end

# :Cᴸˢ
parameter_names = (:CᵂwΔ,  :Cᵂu★, :Cᴰ⁻, :Cᴰʳ,
                   :Cˢc,   :Cˢu,  :Cˢe,
                   :Cᵇc,   :Cᵇu,  :Cᵇe,
                   :Cᴷc⁻,  :Cᴷu⁻, :Cᴷe⁻,
                   :Cᴷcʳ,  :Cᴷuʳ, :Cᴷeʳ,
                   :CᴰRiᶜ, :CᴰRiʷ,
                   :CᴷRiᶜ, :CᴷRiʷ)

free_parameters = FreeParameters(prior_library, names=parameter_names)

function build_simulation()
    simulation = ensemble_column_model_simulation(observations;
                                                  closure = CATKEVerticalDiffusivity(),
                                                  Nensemble = 1,
                                                  architecture = CPU(),
                                                  tracers = (:b, :e))

    simulation.Δt = Δt    

    Qᵘ = simulation.model.velocities.u.boundary_conditions.top.condition
    Qᵇ = simulation.model.tracers.b.boundary_conditions.top.condition
    N² = simulation.model.tracers.b.boundary_conditions.bottom.condition
    
    for (case, obs) in enumerate(observations)
        f = obs.metadata.parameters.coriolis_parameter
        view(Qᵘ, :, case) .= obs.metadata.parameters.momentum_flux
        view(Qᵇ, :, case) .= obs.metadata.parameters.buoyancy_flux
        view(N², :, case) .= obs.metadata.parameters.N²_deep
        view(simulation.model.coriolis, :, case) .= Ref(FPlane(f=f))
    end

    return simulation
end

simulation = build_simulation()
evaluation = InverseProblem(observations, simulation, free_parameters)

velocities = simulation.model.velocities
tracers = simulation.model.tracers
outputs = merge(velocities, tracers)
filename = "catke_evalution_$suite.jld2"
obstimes = first(observations).times

simulation.output_writers[:jld] = JLD2OutputWriter(simulation.model, outputs; filename,
                                                   #schedule = TimeInterval(10minutes),
                                                   schedule = SpecifiedTimes(obstimes),
                                                   overwrite_existing = true)

progress(sim) = @info string("Iteration: ", iteration(sim), ", time: ", prettytime(sim))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

initialize_forward_run!(simulation, observations, evaluation.time_series_collector, evaluation.initialize_simulation)
run!(simulation)

filename = "catke_evalution_$suite.jld2"
btc = FieldTimeSeries(filename, "b")
utc = FieldTimeSeries(filename, "u")
vtc = FieldTimeSeries(filename, "v")
etc = FieldTimeSeries(filename, "e")
Nt = length(btc.times)
z = znodes(btc)

bto, uto, vto, eto = [], [], [], []
for obs in observations
    push!(bto, obs.field_time_serieses.b)
    push!(uto, obs.field_time_serieses.u)
    push!(vto, obs.field_time_serieses.v)
    push!(eto, obs.field_time_serieses.e)
end

fig = Figure(resolution=(1800, 300))
axs = [Axis(fig[1, i]) for i=1:length(observations)]
slider = Slider(fig[2, :], range=1:Nt, startvalue=1)
n = slider.value

for (i, ax) in enumerate(axs)
    bc = @lift interior(btc[$n], 1, i, :)
    bo = @lift interior(bto[i][$n], 1, 1, :)
    lines!(ax, bo, z, linewidth=4, color=(:black, 0.3))
    lines!(ax, bc, z, linewidth=2, color=:blue)
end

display(fig)

