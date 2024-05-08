using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode
using ParameterEstimocean
using GLMakie

Qᵇ = 1e-7
N² = 1e-6
filename = "free_convection_Qb2.0e-08_Nsq1.0e-05_Nh256_Nz256_statistics.jld2" 

bt = FieldTimeSeries(filename, "b")
wbt = FieldTimeSeries(filename, "wb")

#t = bt.times
#
t = collect(0:hour:14days)
Nt = length(t)
Nz = 128 #bt.grid.Nz
Lz = 1024 #bt.grid.Lz

bts = []
qts = []
hts = []
Nts = []

grid = RectilinearGrid(size=Nz, halo=4, z=(-Lz, 0), topology=(Flat, Flat, Bounded))
top_b_bc = FluxBoundaryCondition(Qᵇ)
b_bcs = FieldBoundaryConditions(top=top_b_bc)

const c = Center()
const f = Face()

function mixed_layer_depth(dz_b, N²★=N²/2)
    dz_b_i = interior(dz_b, 1, 1, :)

    k⁻ = findlast(N² -> N² >= N²★, dz_b_i)
    isnothing(k⁻) && return dz_b.grid.Lz

    N²⁻ = @inbounds dz_b_i[k⁻]
    N²⁺ = @inbounds dz_b_i[k⁻ + 1]
    z⁻ = znode(1, 1, k⁻, grid, c, c, f)
    z⁺ = znode(1, 1, k⁻ + 1, grid, c, c, f)

    dN²dz = (N²⁺ - N²⁻) / (z⁺ - z⁻)
    z★ = z⁻ + (N²★ - N²⁻) / dN²dz

    return - z★
end

for (κᶜᵃ, Cᵉⁿ) in zip((2.0, 2.0), (0.0, 20.0))
    # Single column model
    closure = RiBasedVerticalDiffusivity(κ₀ = 0.0, κᶜᵃ = κᶜᵃ, Cᵉⁿ = Cᵉⁿ, Cᵃᵛ=0.0)
    model = HydrostaticFreeSurfaceModel(; grid, closure,
                                        tracers = :b,
                                        buoyancy = BuoyancyTracer(),
                                        boundary_conditions = (; b=b_bcs))
    
    bᵢ(x, y, z) = N² * z
    set!(model, b=bᵢ)
    simulation = Simulation(model, Δt=10minutes, stop_time=t[end])
    
    bt_model = []
    qt_model = []
    Nt_model = []
    ht_model = Float64[]
    
    κ = model.diffusivity_fields.κᶜ
    b = model.tracers.b
    dz_b = Field(∂z(b))
    q = Field(- κ * ∂z(b))
    
    function collect_data(sim)
        compute!(q)
        compute!(dz_b)
        push!(bt_model, deepcopy(interior(b, 1, 1, :)))
        push!(qt_model, deepcopy(interior(q, 1, 1, :)))
        push!(Nt_model, deepcopy(interior(dz_b, 1, 1, :)))

        h = mixed_layer_depth(dz_b)
        push!(ht_model, h)

        return nothing
    end
    
    simulation.callbacks[:collector] = Callback(collect_data, SpecifiedTimes(t))
    
    run!(simulation)

    push!(bts, bt_model)
    push!(qts, qt_model)
    push!(hts, ht_model)
    push!(Nts, Nt_model)
end

set_theme!(Theme(fontsize=16))

fig = Figure()

axb = Axis(fig[1, 1])
axq = Axis(fig[1, 2])
axN = Axis(fig[1, 3])
axh = Axis(fig[2, 1:3])

slider = Slider(fig[3, 1:3], range=1:Nt, startvalue=1)
n = slider.value

zc = znodes(grid, Center())
zf = znodes(grid, Face())

#=
b = @lift interior(bt[$n], 1, 1, :)
wb = @lift interior(wbt[$n], 1, 1, :)

lines!(axb, b, zc)
lines!(axq, wb, zf)
=#

bm1 = @lift bts[1][$n]
qm1 = @lift qts[1][$n][2:Nz]
Nm1 = @lift Nts[1][$n][2:Nz]

bm2 = @lift bts[2][$n]
qm2 = @lift qts[2][$n][2:Nz]
Nm2 = @lift Nts[2][$n][2:Nz]

zm1 = @lift -hts[1][$n]
zm2 = @lift -hts[2][$n]

lines!(axb, bm1, zc)
lines!(axq, qm1, zf[2:Nz])
lines!(axN, Nm1, zf[2:Nz])
hlines!(axb, zm1)
hlines!(axq, zm1)
hlines!(axN, zm1)

lines!(axb, bm2, zc)
lines!(axq, qm2, zf[2:Nz])
lines!(axN, Nm2, zf[2:Nz])
hlines!(axb, zm2)
hlines!(axq, zm2)
hlines!(axN, zm2)

h_LES = @. sqrt(3 * Qᵇ * t / N²)
h_CA  = @. sqrt(2 * Qᵇ * t / N²)

h1 = hts[1] #.+ h_CA[end] .- hts[1][end]
h2 = hts[2] #.+ h_LES[end] .- hts[2][end]

lines!(axh, t, h1, label="Convective adjustment")
lines!(axh, t, h2, label="Adjustment + entrainment")
lines!(axh, t, h_LES, label="√ 3 Q t / N²")
lines!(axh, t, h_CA, label="√ 2 Q t / N²")

axislegend(axh, position=:rb)

xlims!(axq, -1e-7, 1e-7)
xlims!(axN, -1e-5, 1e-5)

display(fig)

