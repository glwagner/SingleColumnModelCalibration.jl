using Oceananigans
using ParameterEstimocean
using GLMakie

suite = "12_hour_suite"
resolution = "1m"
case_path(case) = joinpath("..", "data", "profiles", suite, resolution, case * "_instantaneous_statistics.jld2")
obs = SyntheticObservations(case_path("free_convection"); field_names = (:b, :wb))

t = obs.times
Nt = length(t)

set_theme!(Theme(fontsize=16))

fig = Figure()

axb = Axis(fig[1, 1])
axq = Axis(fig[1, 2])

slider = Slider(fig[2, 1:2], range=1:Nt, startvalue=1)
n = slider.value

bt = obs.field_time_serieses.b
wbt = obs.field_time_serieses.wb

b = @lift interior(bt[$n], 1, 1, :)
wb = @lift interior(wbt[$n], 1, 1, :)

grid = obs.grid
zc = znodes(grid, Center())
zf = znodes(grid, Face())

lines!(axb, b, zc)
lines!(axq, wb, zf)

display(fig)

