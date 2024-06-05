using Oceananigans
using Oceananigans.Operators: Δzᶜᶜᶜ
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity

using ParameterEstimocean
using ParameterEstimocean.InverseProblems: BatchedInverseProblem, inverting_forward_map
using ParameterEstimocean.PseudoSteppingSchemes: Kovachki2018InitialConvergenceRatio
using ParameterEstimocean.Parameters: build_parameters_named_tuple

using Printf
using JLD2
using NCDatasets
using LinearAlgebra
using CairoMakie

using SingleColumnModelCalibration:
    dependent_parameter_sets,
    build_batched_inverse_problem,
    prior_library

set_theme!(Theme(fontsize=19))

@load "optimal_catke.jld2" optimal_catke
closure = optimal_catke
#closure = CATKEVerticalDiffusivity()
closure_label = "CATKE"

cases = ["free_convection",
         "weak_wind_strong_cooling",
         "med_wind_med_cooling",
         "strong_wind_weak_cooling",
         "strong_wind",
         "strong_wind_no_rotation"]

titles = ["6 hour simulation \n (extreme forcing)",
          "12 hour simulation \n (strong forcing)",
          "24 hour simulation \n (medium forcing)",
          "48 hour simulation \n (weak forcing)",
          "72 hour simulation \n (very weak forcing)"]

#####
##### Figure
#####

resolutions = ["1m" for c = 1:5]
#for cc = 1:6

    cc = 2
    fig = Figure(size=(1200, 500))

    Nz = 128
    k₀ = 28
    ax_b = []

    suites = [6, 12, 24, 48, 72]

    zlims = [-200 for c = 1:6]
    zlims[1] = -170
    zlims[2] = -170

    xticks = [0:2e-4:4e-4 for c = 1:6]
    xticks[1] = 0:1e-4:4e-4
    xticks[2] = 0:1e-4:4e-4

    zlim = zlims[cc]

    grid_parameters = [(size=Nz, z=(-256, 0))]
    time_steps = [10.0, 20.0, 30.0, 1minute]

    for c = 1:5
        resolution = resolutions[c]
        suitehrs = suites[c]

        # Batch the inverse problems
        suite = "$(suitehrs)_hour_suite"
        resolution = "1m"
        #catkecolor = (:royalblue1, 0.8)
        catkecolor = (:black, 0.9)
        LEScolor = (:seagreen, 0.6)
        suite_parameters = (; name=suite, resolution, stop_time=suitehrs*hours)

        inverse_problems = []
        
        for Δt in time_steps
            batched_ip = build_batched_inverse_problem(closure, "variable_Pr_conv_adj";
                                                       Nensemble = 1,
                                                       Δt = Δt,
                                                       start_time = 10minutes,
                                                       grid_parameters,
                                                       suite_parameters)
        
            parameters = (; minimum_turbulent_kinetic_energy=1e-9)
            forward_run!(batched_ip, [parameters])
            push!(inverse_problems, batched_ip[1])
        end

        ip = first(inverse_problems)
        times = observation_times(ip.observations)
        Nt = length(times)
        LES_str = string("Large eddy \n simulation")

        grid = ip.simulation.model.grid
        z = znodes(grid, Center())
        Δz = Δzᶜᶜᶜ(1, 1, grid.Nz, grid)
        Δz_str = @sprintf("%d", Δz)

        case = cases[cc]
        yaxisposition = c < 5 ? :left : :right
        Label(fig[2, c], titles[c], tellwidth=false)

        if c == 1
            ax_bc = Axis(fig[3, c]; ylabel="z (m)", xlabel="Buoyancy \n (m s⁻²)", yaxisposition, xticks=[1e-4, 2e-4])
        else
            ax_bc = Axis(fig[3, c]; ylabel="z (m)", xlabel="Buoyancy \n (m s⁻²)", yaxisposition, xticks=xticks[cc])
        end
        push!(ax_b, ax_bc)
        ylims!(ax_b[c], zlim, 5)

        b_init = interior(ip.observations[cc].field_time_serieses.b[1], 1, 1, :)
        b_obs  = interior(ip.observations[cc].field_time_serieses.b[Nt], 1, 1, :)
        b₀ = b_init[k₀]

        lines!(ax_b[c], b_obs .- b₀, z, linewidth=8, label=LES_str, color=LEScolor)

        for (m, ip) = enumerate(inverse_problems)
            Δt = time_steps[m]
            b = interior(ip.time_series_collector.field_time_serieses.b[Nt], 1, cc, :)
            label = string("Δt = ", prettytime(Δt))
            lines!(ax_b[c], b .- b₀, z, linewidth=2; label) #, color=catkecolor)
        end

        hidespines!(ax_b[c], :t)
        c != 1 && hidespines!(ax_b[c], :l)
        c != 1 && c != 5 && hideydecorations!(ax_b[c], grid=false)
        c != 5 && hidespines!(ax_b[c], :r)

        Legend(fig[1, 1:5], ax_b[1], nbanks=5, framevisible=false)
    end

    if cc == 1 # free_convection
        xlims!(ax_b[1], 5e-5, 2.99e-4)
        xlims!(ax_b[2], 5e-5, 2.8e-4)
        xlims!(ax_b[3], 5e-5, 2.6e-4)
        xlims!(ax_b[4], 5e-5, 2.5e-4)
        xlims!(ax_b[5], 5e-5, 2.5e-4)
        xtxt, ztxt = 6e-5, -10
    elseif cc == 2 # weak wind, strong cooling
        xlims!(ax_b[1], 5e-5, 2.9e-4)
        xlims!(ax_b[2], 5e-5, 2.9e-4)
        xlims!(ax_b[3], 5e-5, 2.8e-4)
        xlims!(ax_b[4], 5e-5, 2.8e-4)
        xlims!(ax_b[5], 5e-5, 2.6e-4)
        xtxt, ztxt = 6e-5, -5
    elseif cc == 3 # medium wind, medium cooling
        xlims!(ax_b[1], 1e-5, 3.4e-4)
        xlims!(ax_b[2], 1e-5, 3.3e-4)
        xlims!(ax_b[3], 1e-5, 3.2e-4)
        xlims!(ax_b[4], 1e-5, 3.2e-4)
        xlims!(ax_b[5], 1e-5, 3.3e-4)
        xtxt, ztxt = 2e-5, -5
    elseif cc == 4 # strong wind, weak cooling
        xlims!(ax_b[1], 1e-5, 3.8e-4)
        xlims!(ax_b[2], 1e-5, 3.8e-4)
        xlims!(ax_b[3], 1e-5, 3.8e-4)
        xlims!(ax_b[4], 1e-5, 3.8e-4)
        xlims!(ax_b[5], 1e-5, 3.8e-4)
        xtxt, ztxt = 2e-5, -5
    elseif cc == 5 # strong wind
        xlims!(ax_b[1], 1e-5, 5.4e-4)
        xlims!(ax_b[2], 1e-5, 5.3e-4)
        xlims!(ax_b[3], 1e-5, 5.2e-4)
        xlims!(ax_b[4], 1e-5, 5.2e-4)
        xlims!(ax_b[5], 1e-5, 5.3e-4)
        xtxt, ztxt = 2e-5, -3
    elseif cc == 6 # strong wind, no rotation
        xlims!(ax_b[1], 1e-5, 5.4e-4)
        xlims!(ax_b[2], 1e-5, 5.3e-4)
        xlims!(ax_b[3], 1e-5, 5.2e-4)
        xlims!(ax_b[4], 1e-5, 5.2e-4)
        xlims!(ax_b[5], 1e-5, 5.3e-4)
        xtxt, ztxt = 2e-5, -3
    end

    txts = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    for a = 1:5
        ax = ax_b[a]
        txt = txts[a]
        text!(ax, xtxt, ztxt, text=txt, align=(:left, :top))
    end

    rowsize!(fig.layout, 1, Relative(0.1))

    display(fig)
    case = cases[cc]
#    save("catke_time_step_dependence_$(case)_$Nz.pdf", fig)
#end
