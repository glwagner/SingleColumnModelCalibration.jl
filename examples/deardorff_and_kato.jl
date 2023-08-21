using Oceananigans
using Oceananigans.Units
using Oceananigans.Simulations: reset!
using Oceananigans.Grids: znode

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    CATKEVerticalDiffusivity,
    TurbulentKineticEnergyEquation,
    MixingLength

using GLMakie
using Printf
using Statistics

const c = Center()
const f = Face()

function mixed_layer_depth(dz_b, N²★=5e-6)
    grid = dz_b.grid
    dz_b_i = interior(dz_b, 1, 1, :)

    # Integrating downwards from top
    k⁻ = findlast(N² -> N² >= N²★, dz_b_i)
    isnothing(k⁻) && return grid.Lz

    N²⁻ = @inbounds dz_b_i[k⁻]
    N²⁺ = @inbounds dz_b_i[k⁻ + 1]
    z⁻ = znode(1, 1, k⁻, grid, c, c, f)
    z⁺ = znode(1, 1, k⁻ + 1, grid, c, c, f)

    dN²dz = (N²⁺ - N²⁻) / (z⁺ - z⁻)
    z★ = z⁻ + (N²★ - N²⁻) / dN²dz

    return - z★
end

Qᵘ = zeros(1, 1)
Qᵇ = zeros(1, 1)

top_b_bc = FluxBoundaryCondition(Qᵇ)
b_bcs = FieldBoundaryConditions(top=top_b_bc)

top_u_bc = FluxBoundaryCondition(Qᵘ)
u_bcs = FieldBoundaryConditions(top=top_u_bc)

Cᵉc = 1.0
CᵉD = 0.0
Cᵉe = 0.0
Cʰⁱc = 0.1
mixing_length = MixingLength(; Cᵉc, Cᵉe, Cʰⁱc)
turbulent_kinetic_energy_equation = TurbulentKineticEnergyEquation(; CᵉD)
minimum_convective_buoyancy_flux = 0.0
closure = CATKEVerticalDiffusivity(; minimum_convective_buoyancy_flux,
                                   turbulent_kinetic_energy_equation,
                                   mixing_length)

#closure = CATKEVerticalDiffusivity(; minimum_convective_buoyancy_flux)

function simulate_boundary_layer_deepening(simulation; N²=1e-5, Qᵇ=1e-7, Qᵘ=0.0, stop_time)
    reset!(simulation)
    Nt = ceil(Int, stop_time / simulation.Δt)
    simulation.stop_time = stop_time

    model = simulation.model
    u = model.velocities.u
    b = model.tracers.b
    u.boundary_conditions.top.condition[1, 1] = Qᵘ
    b.boundary_conditions.top.condition[1, 1] = Qᵇ

    bᵢ(x, y, z) = N² * z
    set!(model, b=bᵢ, e=1e-4)

    ∂z_b = Field(∂z(b))

    t = Float64[]
    h = Float64[]

    function collect_mixed_layer_depth(sim)
        compute!(∂z_b)
        hn = mixed_layer_depth(∂z_b, N²/2)
        push!(h, hn)
        push!(t, time(sim))
        return nothing
    end

    simulation.callbacks[:mld] = Callback(collect_mixed_layer_depth)

    run!(simulation)

    return h, t
end

fig = Figure()
axw = Axis(fig[1, 1])
axu = Axis(fig[1, 2])

Lz = 200
h_stop = Lz / 2 # meters
Qᵇ_min = 1e-8
N²_max = 1e-5

τ = h_stop^2 * N²_max / Qᵇ_min

markers = [:rect, :circle, :utriangle]
colors = Makie.wong_colors(0.5)

for (n, Δz) = enumerate((3, 5, 10))
    @show Nz = round(Int, Lz/Δz)
    grid = RectilinearGrid(size=Nz, z=(-Lz, 0), topology=(Flat, Flat, Bounded))
    model = HydrostaticFreeSurfaceModel(; grid, closure,
                                        buoyancy = BuoyancyTracer(),
                                        tracers = (:b, :e),
                                        boundary_conditions = (; u=u_bcs, b=b_bcs))

    simulation = Simulation(model, Δt=10minutes, verbose=false)
    
    for N² = (1e-5, 1e-6)
        for Qᵇ = (2e-8, 5e-8, 1e-7)
            @info "N²: $N², Qᵇ: $Qᵇ"
            # h² = 3 * Qᵇ t / N² → t ~ h² * N² / Qᵇ
            stop_time = h_stop^2 * N² / Qᵇ
            h, t = simulate_boundary_layer_deepening(simulation; Qᵇ, Qᵘ=0.0, N², stop_time)

            # h² = 3 * Qᵇ t / N² → t = h² * N² / Qᵇ
            # √t ~ (h N) / √Qᵇ
            h′ = @. h * sqrt(N² / (t * Qᵇ))
            lines!(axw, t ./ stop_time, h′, color=colors[n])
        end

        for Qᵘ = -[5e-5, 1e-4, 2e-4, 5e-4]
            @info "N²: $N², Qᵘ: $Qᵘ"
            # h² = Qᵘ t / N → t ~ h² * N / Qᵘ
            stop_time = h_stop^2 * sqrt(N²) / abs(Qᵘ)
            h, t = simulate_boundary_layer_deepening(simulation; Qᵇ=0.0, Qᵘ, N², stop_time)

            # h² = 3 * Qᵇ t / N² → t = h² * N² / Qᵇ
            # √t ~ (h N) / √Qᵇ
            h′ = @. h * sqrt(sqrt(N²) / (t * abs(Qᵘ)))
            lines!(axu, t ./ stop_time, h′, color=colors[n])
        end

    end
end

xlims!(axu, 0, 1)
xlims!(axw, 0, 1)
#ylims!(axu, 1, 3)
hlines!(axu, 1, color=:black)
hlines!(axw, sqrt(3), color=:black)

display(fig)
