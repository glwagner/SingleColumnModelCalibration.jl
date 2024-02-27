using Oceananigans
using Oceananigans.Units
using Oceananigans.BuoyancyModels: buoyancy_frequency
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: MixingLength
using Oceananigans.Fields: interpolate!

using ClimaOcean
using ClimaOcean.DataWrangling.ECCO2: ecco2_column

using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

using GLMakie
using Printf
using Dates

# Determine the simulation location
locations = (
    eastern_mediterranean = (λ =  30, φ = 32), 
    ocean_station_papa    = (λ = 215, φ = 50), 
    north_atlantic        = (λ = 325, φ = 50), 
    drake_passage         = (λ = 300, φ = -60), 
    weddell_sea           = (λ = 325, φ = -70), 
    tasman_southern_ocean = (λ = 145, φ = -55), 
)

location = :ocean_station_papa

epoch = Date(1992, 1, 1)

function single_column_simulation(location;
                                  date = Date(1992, 10, 1),
                                  Nz = 100,
                                  Lz = 400)

    # Build the initial condition
    start_seconds = Second(date - epoch).value
    Tᵢ = ecco2_field(:temperature, date)
    Sᵢ = ecco2_field(:salinity, date)

    # Build the grid
    arch = CPU()
    λ★, φ★ = locations[location]
    i★, j★, longitude, latitude = ecco2_column(λ★, φ★)

    grid = LatitudeLongitudeGrid(arch; longitude, latitude,
                                 size = (1, 1, Nz),
                                 z = (-Lz, 0),
                                 topology = (Periodic, Periodic, Bounded))

    # Build the ocean model
    top_ocean_heat_flux          = Qᵀ = Field{Center, Center, Nothing}(grid)
    top_salt_flux                = Fˢ = Field{Center, Center, Nothing}(grid)
    top_zonal_momentum_flux      = τˣ = Field{Face, Center, Nothing}(grid)
    top_meridional_momentum_flux = τʸ = Field{Center, Face, Nothing}(grid)

    ocean_boundary_conditions = (u = FieldBoundaryConditions(top=FluxBoundaryCondition(τˣ)),
                                 v = FieldBoundaryConditions(top=FluxBoundaryCondition(τʸ)),
                                 T = FieldBoundaryConditions(top=FluxBoundaryCondition(Qᵀ)),
                                 S = FieldBoundaryConditions(top=FluxBoundaryCondition(Fˢ)))

    teos10 = TEOS10EquationOfState()
    buoyancy = SeawaterBuoyancy(equation_of_state=teos10)

    mixing_length = MixingLength(Cᵇ=0.01)
    closure = CATKEVerticalDiffusivity(; mixing_length,
                                       minimum_turbulent_kinetic_energy = 1e-6,
                                       minimum_convective_buoyancy_flux = 1e-11)

    ocean_model = HydrostaticFreeSurfaceModel(; grid, buoyancy, closure,
                                              tracers = (:T, :S, :e),
                                              tracer_advection = nothing,
                                              momentum_advection = nothing,
                                              free_surface = SplitExplicitFreeSurface(cfl=0.7; grid),
                                              boundary_conditions = ocean_boundary_conditions,
                                              coriolis = HydrostaticSphericalCoriolis())

    ocean = Simulation(ocean_model; Δt=2minutes, verbose=false)

    forcing_days = 60
    backend = JRA55NetCDFBackend(8 * forcing_days)
    atmosphere = JRA55_prescribed_atmosphere(Colon(); longitude, latitude, backend)

    ocean.model.clock.time = start_seconds
    ocean.model.clock.iteration = 0
    interpolate!(ocean.model.tracers.T, Tᵢ)
    interpolate!(ocean.model.tracers.S, Sᵢ)
    set!(ocean.model, e=1e-6)

    radiation = Radiation()
    coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)
    coupled_simulation = Simulation(coupled_model, Δt=2minutes, stop_time=start_seconds + forcing_days * days)

    wall_clock = Ref(time_ns())

    function progress(sim)
        msg = string("(", location, ")")
        msg *= string(", iter: ", iteration(sim), ", time: ", prettytime(sim))

        elapsed = 1e-9 * (time_ns() - wall_clock[])
        msg *= string(", wall time: ", prettytime(elapsed))
        wall_clock[] = time_ns()

        u, v, w = sim.model.ocean.model.velocities
        msg *= @sprintf(", max|u|: (%.2e, %.2e)", maximum(abs, u), maximum(abs, v))

        T = sim.model.ocean.model.tracers.T
        S = sim.model.ocean.model.tracers.S
        e = sim.model.ocean.model.tracers.e

        τˣ = first(sim.model.fluxes.total.ocean.momentum.τˣ)
        τʸ = first(sim.model.fluxes.total.ocean.momentum.τʸ)
        u★ = sqrt(sqrt(τˣ^2 + τʸ^2))
        Q = first(sim.model.fluxes.total.ocean.heat)

        Nz = size(T, 3)
        msg *= @sprintf(", u★: %.2f m s⁻¹", u★)
        msg *= @sprintf(", Q: %.2f W m⁻²", Q)
        msg *= @sprintf(", T₀: %.2f ᵒC",     first(interior(T, 1, 1, Nz)))
        msg *= @sprintf(", extrema(T): (%.2f, %.2f) ᵒC", minimum(T), maximum(T))
        msg *= @sprintf(", S₀: %.2f g/kg",   first(interior(S, 1, 1, Nz)))
        msg *= @sprintf(", e₀: %.2e m² s⁻²", first(interior(e, 1, 1, Nz)))

        @info msg
    end

    coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

    # Build flux outputs
    Ju = coupled_model.fluxes.total.ocean.momentum.u
    Jv = coupled_model.fluxes.total.ocean.momentum.v
    JT = coupled_model.fluxes.total.ocean.tracers.T
    JS = coupled_model.fluxes.total.ocean.tracers.S
    E  = coupled_model.fluxes.turbulent.fields.water_vapor
    Qc = coupled_model.fluxes.turbulent.fields.sensible_heat
    Qv = coupled_model.fluxes.turbulent.fields.latent_heat
    ρₒ = coupled_model.fluxes.ocean_reference_density
    cₚ = coupled_model.fluxes.ocean_heat_capacity

    Q = ρₒ * cₚ * JT
    τx = ρₒ * Ju
    τy = ρₒ * Jv
    N² = buoyancy_frequency(ocean.model)
    κc = ocean.model.diffusivity_fields.κᶜ

    fluxes = (; τx, τy, E, JS, Q, Qc, Qv)

    auxiliary_fields = (; N², κc)
    fields = merge(ocean.model.velocities, ocean.model.tracers, auxiliary_fields)

    # Slice fields at the surface
    outputs = merge(fields, fluxes)

    output_attributes = Dict{String, Any}(
        "κc"  => Dict("long_name" => "Tracer diffusivity",          "units" => "m^2 / s"),
        "Q"   => Dict("long_name" => "Net heat flux",               "units" => "W / m^2", "convention" => "positive upwards"),
        "Qv"  => Dict("long_name" => "Latent heat flux",            "units" => "W / m^2", "convention" => "positive upwards"),
        "Qc"  => Dict("long_name" => "Sensible heat flux",          "units" => "W / m^2", "convention" => "positive upwards"),
        "Js"  => Dict("long_name" => "Salt flux",                   "units" => "g kg⁻¹ m s⁻¹", "convention" => "positive upwards"),
        "E"   => Dict("long_name" => "Freshwater evaporation flux", "units" => "m s⁻¹", "convention" => "positive upwards"),
        "e"   => Dict("long_name" => "Turbulent kinetic energy",    "units" => "m^2 / s^2"),
        "τx"  => Dict("long_name" => "Zonal momentum flux",         "units" => "m^2 / s^2"),
        "τx"  => Dict("long_name" => "Meridional momentum flux",    "units" => "m^2 / s^2"),
    )

    filename = "single_column_omip_$location"

    coupled_simulation.output_writers[:jld2] = JLD2OutputWriter(ocean.model, outputs; filename,
                                                                schedule = TimeInterval(1hours),
                                                                overwrite_existing = true)

    #=
    coupled_simulation.output_writers[:nc] = NetCDFOutputWriter(ocean.model, outputs; filename,
                                                                schedule = AveragedTimeInterval(1days),
                                                                output_attributes,
                                                                overwrite_existing = true)
    =#

    return coupled_simulation
end

