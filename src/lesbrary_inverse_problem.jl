function batched_lesbrary_observations(regrid; times, suite,
                                       resolution = "1m",
                                       field_names = (:b, :e, :u, :v),
                                       tke_weight = 0.0,
                                       cases = default_cases)

    normalizations = (b = ZScore(),
                      u = ZScore(),
                      v = ZScore(),
                      e = RescaledZScore(tke_weight))

    Nz = regrid.Nz
    zf = znodes(Face, regrid)
    kb = round(Int, Nz/3) # exclude the bottom third, ie ≈ - 170m
    kt = findfirst(z -> z > -16, zf) - 1 # exclude the top 16 meters
    space = SpaceIndices(z=kb:kt)

    transformation = NamedTuple(n => Transformation(; space, normalization=normalizations[n])
                                for n in keys(normalizations))

    case_path(case) = joinpath("/Users/andresouza/Desktop/Repositories/SingleColumnModelCalibration.jl/data", suite, resolution, case * "_instantaneous_statistics.jld2")

    observation_library = Dict()

    # Don't optimize u, v for free_convection
    free_convection_names = filter(n -> n ∈ (:b, :e), field_names)
    observation_library["free_convection"] =
        SyntheticObservations(case_path("free_convection"); transformation, times, regrid,
                              field_names = free_convection_names)
                                                                    
    # Don't optimize v for non-rotating cases
    strong_wind_no_rotation_names = filter(n -> n ∈ (:b, :e, :u), field_names)
    observation_library["strong_wind_no_rotation"] =
        SyntheticObservations(case_path("strong_wind_no_rotation"); transformation, times, regrid,
                              field_names = strong_wind_no_rotation_names)

    # The rest are standard
    for case in ["strong_wind", "med_wind_med_cooling", "strong_wind_weak_cooling", "weak_wind_strong_cooling"]
        observation_library[case] = SyntheticObservations(case_path(case); transformation, times, regrid, field_names)
    end

    observations = [observation_library[case] for case in cases]
    weights = [observation_weights[name] for name in cases]
    batched_observations = BatchedSyntheticObservations(observations; weights)

    return batched_observations
end

function estimate_noise_covariance(grid; kwargs...)
    # Estimate noise covariance
    obs_1m = batched_lesbrary_observations(grid; resolution="1m", kwargs...)
    obs_2m = batched_lesbrary_observations(grid; resolution="2m", kwargs...)
    obs_4m = batched_lesbrary_observations(grid; resolution="4m", kwargs...)
    Γ = cov([obs_1m, obs_2m, obs_4m])
    ϵ = 1e-2 #* mean(abs, [Γ[n, n] for n=1:size(Γ, 1)])
    Γ .+= ϵ * Diagonal(I, size(Γ, 1))

    # Remove off-diagonal elements
    Γ = diagm(diag(Γ)) .* 2
    # Γ = Γ + I(size(Γ, 1))
    
    return Γ
end
    
function lesbrary_inverse_problem(regrid;
                                  free_parameters,
                                  times = [48hours - 10minutes, 48hours],
                                  Nensemble = 500,
                                  observations_resolution = "1m",
                                  Δt = 10minutes,
                                  field_names = (:b, :e, :u, :v),
                                  closure = CATKEVerticalDiffusivity(),
                                  non_ensemble_closure = nothing,
                                  suite = "one_day_suite",
                                  tke_weight = 0.0,
                                  cases = default_cases,
                                  architecture = CPU())

    batched_observations = batched_lesbrary_observations(regrid; resolution=observations_resolution,
                                                         times, field_names, suite, tke_weight, cases)

    observations = batched_observations.observations

    simulation = ensemble_column_model_simulation(observations;
                                                  closure,
                                                  Nensemble,
                                                  architecture,
                                                  non_ensemble_closure,
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

    ip = InverseProblem(batched_observations, simulation, free_parameters)

    return ip
end

