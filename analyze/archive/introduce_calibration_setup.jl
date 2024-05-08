using Oceananigans
using Oceananigans.Units
using GLMakie
using JLD2

using ParameterEstimocean
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using SingleColumnModelCalibration: build_ensemble_kalman_inversion

set_theme!(Theme(fontsize=32))

grid_parameters = [
    (size=32, z=(-256, 0)),
    (size=64, z=(-256, 0)),
]

suite_parameters = [
    (name = "12_hour_suite", stop_time=12hours),
    (name = "24_hour_suite", stop_time=24hours),
    (name = "48_hour_suite", stop_time=48hours),
]

profilesdir = "../data/profiles"
slicesdir = "../data/slices"
parametersdir = "../parameters"
name =  "variable_Pr_conv_adj"

filepath = joinpath(parametersdir, string(name) * "_best_parameters.jld2")
file = jldopen(filepath)
optimal_parameters = file["optimal_parameters"]
first_best_parameters = file["first_best_parameters"]
first_mean_parameters = file["example_first_mean_parameters"]
close(file)

closure = CATKEVerticalDiffusivity()
architecture = CPU()
Nensemble = 2
resample_failure_fraction = 0.05
stop_pseudotime = 1000.0
Δt = 20minutes
linewidth = 6

start_time = time_ns()

suites = [p.name for p in suite_parameters]

suitelabels = [
    "12 hour suite",
    "24 hour suite",
    "48 hour suite",
]

cases = [
    "free_convection",
    "weak_wind_strong_cooling",
    "med_wind_med_cooling",
    "strong_wind_weak_cooling",
    "strong_wind",
    "strong_wind_no_rotation",
]

fig = Figure(resolution=(2000, 1400))
wlim = 0.1

for (s, suite) in enumerate(suites)
    l = 2s-1
    Label(fig[1, l:l+1], suitelabels[s])
    for (c, case) in enumerate(cases)
        # For doing x-y plots instead:
        #dataname = case * "_xy_slice.jld2"
        #ax = Axis(fig[i, j], aspect=1)
        #heatmap!(ax, interior(wn, :, :, 1), colormap=:balance, colorrange=(-wlim, wlim))
        
        datadir = joinpath(slicesdir, suite)
        dataname = case * "_xz_slice.jld2"
        datapath = joinpath(datadir, dataname)
        w = FieldTimeSeries(datapath, "w")
        wn = w[end]
        # @show wlim = maximum(abs, w) / 2
        
        # (row, column)
        i = c <= 3 ? c : c-3
        j = c <= 3 ? l : l+1

        ax = Axis(fig[i+1, j], aspect=2)
        global hm = heatmap!(ax, interior(wn, :, 1, :), colormap=:balance, colorrange=(-wlim, wlim))
        hidedecorations!(ax)
    end
end

Nsuites = length(suites)
Ncases = length(cases)

for (j, suite) in enumerate(suites)
    Np = Nsuites+2

    axT = Axis(fig[Np:Np+1, 2j - 1];
               xlabel = "Buoyancy (m s⁻²)",
               ylabel = "z (m)",
               xticks = [0.0386, 0.039])

    if j == 1
        hidespines!(axT, :t, :r)
    else
        hidespines!(axT, :t, :r, :l)
        hideydecorations!(axT; grid=false)
    end

    axU = Axis(fig[Np:Np+1, 2j],
               yaxisposition = :right,
               xlabel = "Velocities (m s⁻¹)",
               ylabel = "z (m)",
               xticks = [-0.1, 0.1, 0.3])

    xlims!(axU, -0.2, 0.5)

    if j == 3
        hidespines!(axU, :t, :l)
    else
        hideydecorations!(axU; grid=false)
        hidespines!(axU, :t, :l, :r)
    end

    if suite == "24_hour_suite"
        text!(axU, +0.3, -80, text="u", textsize=32)
        text!(axU, -0.17, -110, text="v", textsize=32)
    end

    for (j, case) in enumerate(cases)
        datadir = joinpath(profilesdir, suite, "1m")
        dataname = case * "_instantaneous_statistics.jld2"
        datapath = joinpath(datadir, dataname)
        T = FieldTimeSeries(datapath, "T")
        b = FieldTimeSeries(datapath, "b")
        u = FieldTimeSeries(datapath, "u")
        v = FieldTimeSeries(datapath, "v")
        z = znodes(T)

        Tn = T[end]
        bn = b[end]
        un = u[end]
        vn = v[end]

        #ln = lines!(axT, interior(Tn, 1, 1, :), z; linewidth)
        ln = lines!(axT, interior(bn, 1, 1, :), z; linewidth)
        ln.color = (ln.color.val, 0.7)

        if case != "free_convection"
            lines!(axU, interior(un, 1, 1, :), z; color=ln.color, linewidth)
            if case != "strong_wind_no_rotation"
                lines!(axU, interior(vn, 1, 1, :), z; color=ln.color, linestyle=:dash, linewidth=4)
            end
        end
    end
end

Ny = Nsuites+4
ax = Axis(fig[Ny:Ny+2, :], ylabel=L"\mathcal{Y}", xlabel="Data index")
hidespines!(ax, :t, :r)

#=
eki = build_ensemble_kalman_inversion(closure, name;
                                      architecture,
                                      Nensemble = 3,
                                      Δt,
                                      grid_parameters, # = grid_parameters[2:2],
                                      suite_parameters, # = suite_parameters[3:3],
                                      resample_failure_fraction)


G = forward_map(eki.inverse_problem, [first_mean_parameters, optimal_parameters])
G₁ = G[:, 1]
G★ = G[:, 2]
=#

Y = observation_map(eki.inverse_problem)[:]

lines!(ax, Y, linewidth=5, color=(:black, 0.5), label=L"\mathcal{Y}")

display(fig)

save("introduce_calibration_setup.png", fig)

