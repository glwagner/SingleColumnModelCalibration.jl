using Oceananigans
using Oceananigans.Units
using ParameterEstimocean
using GLMakie

Qᵇ = 1e-7
N² = 1e-5
filename = "free_convection_Qb2.0e-08_Nsq1.0e-05_Nh256_Nz256_statistics.jld2" 

bt = FieldTimeSeries(filename, "b")
wbt = FieldTimeSeries(filename, "wb")

#t = bt.times
#
t = collect(0:hour:7days)
Nt = length(t)
Nz = 32  #bt.grid.Nz
Lz = 256 #bt.grid.Lz

bts = []
qts = []

grid = RectilinearGrid(size=Nz, z=(-Lz, 0), topology=(Flat, Flat, Bounded))
top_b_bc = FluxBoundaryCondition(Qᵇ)
b_bcs = FieldBoundaryConditions(top=top_b_bc)

for Cᵉⁿ = (0.0, 0.5)
    # Single column model
    closure = RiBasedVerticalDiffusivity(κ₀ = 0.0, κᶜᵃ = 1.0, Cᵉⁿ=Cᵉⁿ, Cᵃᵛ=0.0)
    model = HydrostaticFreeSurfaceModel(; grid, closure,
                                        tracers = :b,
                                        buoyancy = BuoyancyTracer(),
                                        boundary_conditions = (; b=b_bcs))
    
    bᵢ(x, y, z) = N² * z
    set!(model, b=bᵢ)
    simulation = Simulation(model, Δt=5minutes, stop_time=t[end])
    
    bt_model = []
    qt_model = []
    
    κ = model.diffusivity_fields.κᶜ
    b = model.tracers.b
    q = Field(- κ * ∂z(b))
    
    function collect_data(sim)
        compute!(q)
        push!(bt_model, deepcopy(interior(b, 1, 1, :)))
        push!(qt_model, deepcopy(interior(q, 1, 1, :)))
        return nothing
    end
    
    simulation.callbacks[:collector] = Callback(collect_data, SpecifiedTimes(t))
    
    run!(simulation)

    push!(bts, bt_model)
    push!(qts, qt_model)
end

set_theme!(Theme(fontsize=16))

fig = Figure()

axb = Axis(fig[1, 1])
axq = Axis(fig[1, 2])

slider = Slider(fig[2, 1:2], range=1:Nt, startvalue=1)
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

bm2 = @lift bts[2][$n]
qm2 = @lift qts[2][$n][2:Nz]

lines!(axb, bm1, zc)
lines!(axq, qm1, zf[2:Nz])

lines!(axb, bm2, zc)
lines!(axq, qm2, zf[2:Nz])

xlims!(axq, -1e-7, 1e-7)

display(fig)

