using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: Time
using ClimaOcean
using ClimaOcean.DataWrangling.ECCO2: ecco2_column
using GLMakie
using Printf
using Dates
using JLD2

# Determine the simulation location
locations = (
    eastern_mediterranean = (λ =  30, φ = 32),
    ocean_station_papa    = (λ = 215, φ = 50),
    north_atlantic        = (λ = 325, φ = 50),
    drake_passage         = (λ = 300, φ = -60),
    weddell_sea           = (λ = 325, φ = -70),
    tasman_southern_ocean = (λ = 145, φ = -55),
)

#location = :ocean_station_papa
location = :north_atlantic
filename = "single_column_omip_$location.jld2"

file = jldopen(filename)
ρₒ = file["ρₒ"]
close(file)

ut  = FieldTimeSeries(filename, "u")
vt  = FieldTimeSeries(filename, "v")
Tt  = FieldTimeSeries(filename, "T")
St  = FieldTimeSeries(filename, "S")
et  = FieldTimeSeries(filename, "e")
N²t = FieldTimeSeries(filename, "N²")
κt  = FieldTimeSeries(filename, "κc")

Qt  = FieldTimeSeries(filename, "Q")
Qct = FieldTimeSeries(filename, "Qc")
Qvt = FieldTimeSeries(filename, "Qv")
Jˢt = FieldTimeSeries(filename, "JS")
Et  = FieldTimeSeries(filename, "E")
τxt = FieldTimeSeries(filename, "τx")
τyt = FieldTimeSeries(filename, "τy")

Nz = size(Tt, 3)
times = Qt.times

λ★, φ★ = locations[location]
i★, j★, longitude, latitude = ecco2_column(λ★, φ★)
forcing_days = 60
backend = JRA55NetCDFBackend(8 * forcing_days)
atmosphere = JRA55_prescribed_atmosphere(Colon(); longitude, latitude, backend)

ua = atmosphere.velocities.u
va = atmosphere.velocities.v
Ta = atmosphere.tracers.T
qa = atmosphere.tracers.q
Ql = atmosphere.downwelling_radiation.longwave
Qs = atmosphere.downwelling_radiation.shortwave
Fr = atmosphere.freshwater_flux.rain
Fs = atmosphere.freshwater_flux.snow

Nt = length(times)
uat = zeros(Nt)
vat = zeros(Nt)
Tat = zeros(Nt)
qat = zeros(Nt)
Qst = zeros(Nt)
Qlt = zeros(Nt)
Ft = zeros(Nt)

for n = 1:Nt
    t = times[n]
    uat[n] = ua[1, 1, 1, Time(t)]
    vat[n] = va[1, 1, 1, Time(t)]
    Tat[n] = Ta[1, 1, 1, Time(t)]
    qat[n] = qa[1, 1, 1, Time(t)]
    Qst[n] = Qs[1, 1, 1, Time(t)]
    Qlt[n] = Ql[1, 1, 1, Time(t)]
    Ft[n]  = Fr[1, 1, 1, Time(t)] + Fs[1, 1, 1, Time(t)]
end

set_theme!(Theme(linewidth=3))

fig = Figure(size=(2400, 1800))

axτ = Axis(fig[1, 1:2], xlabel="Days since Oct 1 1992", ylabel="Wind stress (N m⁻²)")
axu = Axis(fig[2, 1:2], xlabel="Days since Oct 1 1992", ylabel="Velocities (m s⁻¹)")
axQ = Axis(fig[1, 3:4], xlabel="Days since Oct 1 1992", ylabel="Heat flux (W m⁻²)")
axT = Axis(fig[2, 3:4], xlabel="Days since Oct 1 1992", ylabel="Surface temperature (ᵒC)")
axF = Axis(fig[1, 5:6], xlabel="Days since Oct 1 1992", ylabel="Freshwater volume flux (m s⁻¹)")
axS = Axis(fig[2, 5:6], xlabel="Days since Oct 1 1992", ylabel="Surface salinity (g kg⁻¹)")

axuz = Axis(fig[3, 1], xlabel="Velocities (m s⁻¹)",                ylabel="z (m)")
axTz = Axis(fig[3, 2], xlabel="Temperature (ᵒC)",                  ylabel="z (m)")
axSz = Axis(fig[3, 3], xlabel="Salinity (g kg⁻¹)",                 ylabel="z (m)")
axNz = Axis(fig[3, 4], xlabel="Buoyancy frequency (s⁻²)",          ylabel="z (m)")
axκz = Axis(fig[3, 5], xlabel="Eddy diffusivity (m² s⁻¹)",         ylabel="z (m)", xscale=log10)
axez = Axis(fig[3, 6], xlabel="Turbulent kinetic energy (m² s⁻²)", ylabel="z (m)", xscale=log10)

title = @sprintf("Single column simulation at %.2f, %.2f", φ★, λ★)
Label(fig[0, 1:6], title)

slider = Slider(fig[4, 1:6], range=1:Nt, startvalue=1)
n = slider.value

times = (times .- times[1]) ./days
tn = @lift times[$n]

colors = Makie.wong_colors()

#lines!(axu, times, uat, color=colors[1])
#lines!(axu, times, vat, color=colors[2])

Jut = interior(τxt, 1, 1, 1, :) ./ ρₒ
Jvt = interior(τyt, 1, 1, 1, :) ./ ρₒ
u★ = @. (Jut^2 + Jvt^2)^(1/4)

lines!(axu, times, interior(ut, 1, 1, Nz, :), color=colors[1], label="Zonal")
lines!(axu, times, interior(vt, 1, 1, Nz, :), color=colors[2], label="Meridional")
lines!(axu, times, u★, color=colors[3], label="Ocean-side u★") 
vlines!(axu, tn, linewidth=4, color=(:black, 0.5))
axislegend(axu)

lines!(axτ, times, interior(τxt, 1, 1, 1, :), label="Zonal")
lines!(axτ, times, interior(τyt, 1, 1, 1, :), label="Meridional")
vlines!(axτ, tn, linewidth=4, color=(:black, 0.5))
axislegend(axτ)

lines!(axT, times, Tat .- 273.15,             color=colors[1], linewidth=2, linestyle=:dash, label="Atmosphere temperature")
lines!(axT, times, interior(Tt, 1, 1, Nz, :), color=colors[2], linewidth=4, label="Ocean surface temperature")
vlines!(axT, tn, linewidth=4, color=(:black, 0.5))
axislegend(axT)

lines!(axQ, times, interior(Qt, 1, 1, 1, :),   color=colors[1], label="Total",     linewidth=6)
lines!(axQ, times, interior(Qct, 1, 1, 1, :), color=colors[2], label="Sensible",  linewidth=2)
lines!(axQ, times, interior(Qvt, 1, 1, 1, :), color=colors[3], label="Latent",    linewidth=2)
lines!(axQ, times, - Qst,                     color=colors[4], label="Shortwave", linewidth=2)
lines!(axQ, times, - Qlt,                     color=colors[5], label="Longwave",  linewidth=2)
vlines!(axQ, tn, linewidth=4, color=(:black, 0.5))
axislegend(axQ)

#lines!(axF, times, interior(Jˢt, 1, 1, 1, :), label="Net freshwater flux")
lines!(axF, times, Ft, label="Prescribed freshwater flux")
lines!(axF, times, - interior(Et, 1, 1, 1, :), label="Evaporation")
vlines!(axF, tn, linewidth=4, color=(:black, 0.5))
axislegend(axF)

lines!(axS, times, interior(St, 1, 1, Nz, :))
vlines!(axS, tn, linewidth=4, color=(:black, 0.5))

zc = znodes(Tt)
zf = znodes(κt)
un = @lift interior(ut[$n], 1, 1, :)
vn = @lift interior(vt[$n], 1, 1, :)
Tn = @lift interior(Tt[$n], 1, 1, :)
Sn = @lift interior(St[$n], 1, 1, :)
κn = @lift interior(κt[$n], 1, 1, :)
en = @lift max.(1e-6, interior(et[$n], 1, 1, :))
N²n = @lift interior(N²t[$n], 1, 1, :)

scatterlines!(axuz, un,  zc, label="u") 
scatterlines!(axuz, vn,  zc, label="v") 
scatterlines!(axTz, Tn,  zc) 
scatterlines!(axSz, Sn,  zc) 
scatterlines!(axez, en,  zc) 
scatterlines!(axNz, N²n, zf) 
scatterlines!(axκz, κn,  zf) 

axislegend(axuz)

Tmax = maximum(interior(Tt))
Tmin = minimum(interior(Tt))
xlims!(axTz, Tmin - 0.1, Tmax + 0.1)

Nmax = maximum(interior(N²t))
Nmin = minimum(interior(N²t))
xlims!(axNz, Nmin / 2, Nmin * 1.1)

emax = maximum(interior(et))
xlims!(axez, 8e-7, emax * 1.1)
xlims!(axκz, 1e-7, 10)

Smax = maximum(interior(St))
Smin = minimum(interior(St))
xlims!(axSz, Smin - 0.2, Smax + 0.2)

display(fig)

record(fig, "$(location)_single_column_simulation.mp4", 1:Nt, framerate=24) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end

