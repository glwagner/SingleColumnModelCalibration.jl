using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using JLD2

filename = "single_column_omip_north_atlantic.jld2"

n₀ = 1 # starting time index
arch = CPU()

Lx = Ly = 256
Lz = 256

Nx = Ny = 64
Nz = 64

grid = RectilinearGrid(arch,
                       size = (Nx, Ny, Nz),
                       halo = (5, 5, 5),
                       x = (0, Lx),
                       y = (0, Ly),
                       z = (-Lz, 0),
                       topology = (Periodic, Periodic, Bounded))

# Boundary conditions
Ju = FieldTimeSeries(filename, "Ju")
Jv = FieldTimeSeries(filename, "Jv")
JS = FieldTimeSeries(filename, "JS")
JT = FieldTimeSeries(filename, "JT")

u_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(Ju))
v_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(Jv))
T_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(JT))
S_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(JS))

file = jldopen(filename)
f₀ = file["f₀"]
close(file)

teos10 = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(equation_of_state=teos10)
coriolis = FPlane(f=f₀)

model = NonhydrostaticModel(; grid, buoyancy, coriolis,
                            tracers = (:T, :S),
                            boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs))

# Initial condition
u = FieldTimeSeries(filename, "u")
v = FieldTimeSeries(filename, "v")
T = FieldTimeSeries(filename, "T")
S = FieldTimeSeries(filename, "S")

data_grid = u.grid
ic_grid = RectilinearGrid(arch,
                          size = data_grid.Nz,
                          halo = data_grid.Hz,
                          z = (-data_grid.Lz, 0),
                          topology = (Flat, Flat, Bounded))

t₀ = u.times[n₀]

uᵢ = Field{Nothing, Nothing, Center}(ic_grid)
vᵢ = Field{Nothing, Nothing, Center}(ic_grid)
Tᵢ = Field{Nothing, Nothing, Center}(ic_grid)
Sᵢ = Field{Nothing, Nothing, Center}(ic_grid)

copyto!(parent(uᵢ), parent(u[n₀]))
copyto!(parent(vᵢ), parent(v[n₀]))
copyto!(parent(Tᵢ), parent(T[n₀]))
copyto!(parent(Sᵢ), parent(S[n₀]))

interpolate!(model.velocities.u, uᵢ) 
interpolate!(model.velocities.v, vᵢ) 
interpolate!(model.tracers.T, Tᵢ) 
interpolate!(model.tracers.S, Sᵢ) 


