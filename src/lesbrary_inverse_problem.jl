function batched_lesbrary_observations(regrid; times, suite,
                                       resolution = "1m",
                                       field_names = (:b, :c, :u, :v),
                                       tke_weight = 0.0)

    normalizations = (b = ZScore(),
                      c = ZScore(),
                      u = ZScore(),
                      v = ZScore(),
                      e = RescaledZScore(tke_weight))

    Nz = regrid.Nz
    zf = znodes(regrid, Face())
    kb = round(Int, Nz/3) # exclude the bottom third, ie ≈ - 170m
    k4 = findfirst(z -> z > -4, zf) - 1 # exclude the top 8 meters
    kt = min(Nz-1, k4) # exclude the top 4 meters, or top grid point --- whichever is larger
    space = SpaceIndices(z=kb:kt)

    transformation = NamedTuple(n => Transformation(; space, normalization=normalizations[n])
                                for n in keys(normalizations))

    case_path(case) = joinpath("..", "data", "profiles", suite, resolution,
                               case * "_with_tracer_instantaneous_statistics.jld2")

    observation_library = Dict()

    field_names = (:b, :c, :u, :v)

    # Don't optimize u, v for free_convection
    free_convection_names = filter(n -> n ∈ (:b, :c, :e), field_names)
    observation_library["free_convection"] =
        SyntheticObservations(case_path("free_convection"); transformation, times, regrid,
                              field_names = free_convection_names)
                                                                    
    # Don't optimize v for non-rotating cases
    no_rotation_names = filter(n -> n ∈ (:b, :c, :e, :u), field_names)
    for case in ["strong_wind_no_rotation", "strong_wind_and_sunny"]
        observation_library[case] =
            SyntheticObservations(case_path(case); transformation, times, regrid,
                                  field_names = no_rotation_names)
    end

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

    # Add constant component to diagonal to "regularize" Γ
    diag = [Γ[n, n] for n=1:size(Γ, 1)]
    ϵ = 1e0 * mean(abs, diag)
    Γ .+= ϵ * Diagonal(I, size(Γ, 1))

    return Γ
end
    
function lesbrary_inverse_problem(regrid;
                                  free_parameters,
                                  times = [48hours - 10minutes, 48hours],
                                  Nensemble = 100,
                                  observations_resolution = "1m",
                                  Δt = 5minutes,
                                  field_names = (:b, :e, :u, :v),
                                  closure = CATKEVerticalDiffusivity(),
                                  non_ensemble_closure = nothing,
                                  suite = "one_day_suite",
                                  tke_weight = 0.0,
                                  architecture = CPU())

    batched_observations = batched_lesbrary_observations(regrid; resolution=observations_resolution,
                                                         times, field_names, suite, tke_weight)

    observations = batched_observations.observations

    simulation = ensemble_column_model_simulation(observations;
                                                  closure,
                                                  Nensemble,
                                                  architecture,
                                                  non_ensemble_closure,
                                                  verbose = false,
                                                  forced_fields = (:b, :c),
                                                  tracers = (:b, :e, :ϵ, :c))

    simulation.Δt = Δt    

    τˣ = simulation.model.velocities.u.boundary_conditions.top.condition
    Jᵇ = simulation.model.tracers.b.boundary_conditions.top.condition
    N² = simulation.model.tracers.b.boundary_conditions.bottom.condition

    grid = simulation.model.grid
    I_field = ZFaceField(grid)
    Fᵇ = simulation.model.forcing.b.parameters
    Fᶜ = simulation.model.forcing.c.parameters

    for (case, obs) in enumerate(observations)
        f = obs.metadata.parameters.coriolis_parameter
        view(τˣ, :, case) .= obs.metadata.parameters.momentum_flux
        view(Jᵇ, :, case) .= obs.metadata.parameters.buoyancy_flux
        view(N², :, case) .= obs.metadata.parameters.N²_deep
        view(simulation.model.coriolis, :, case) .= Ref(FPlane(f=f))

        τ = obs.metadata.parameters.tracer_forcing_timescale
        z₀ = - obs.metadata.parameters.tracer_forcing_depth
        λ = obs.metadata.parameters.tracer_forcing_width
        μ⁺ = 1 / τ
        μ⁻ = √(2π) * λ / grid.Lz * μ⁺
        c_forcing_func(z) = μ⁺ * exp(-(z - z₀)^2 / (2 * λ^2)) - μ⁻
        c_forcing = view(Fᶜ, :, case, :)
        set!(c_forcing, c_forcing_func)

        I₀ = obs.metadata.parameters.penetrating_buoyancy_flux
        # Not saved in file now (but they will / might be)
        ϵ₁ = 0.6
        λ₁ = 1.0
        λ₂ = 20.0

        if !isnothing(I₀)
            I_case = view(I_field, :, case, :)
            I(z) = I₀ * (ϵ₁ * exp(z / λ₁) + (1 - ϵ₁) * exp(z / λ₂))
            set!(I_case, I)
            dIdz = Field(-1 * ∂z(I_case))
            compute!(dIdz)
            interior(Fᵇ, :, case, :) .= interior(dIdz, :, 1, :)
        end
    end

    ip = InverseProblem(batched_observations, simulation, free_parameters)

    return ip
end

