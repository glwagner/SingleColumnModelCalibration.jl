using JLD2
using ParameterEstimocean
using LinearAlgebra
using Statistics
using Printf
using GLMakie

fileprefix = "progress_"

diags = []
times = []

for i = 1:50
    filename = @sprintf("%s%02d.jld2", fileprefix, i)

    file = jldopen(filename)
    summary = file["latest_summary"]
    close(file)

    Nensemble = length(summary.parameters)

    θ = hcat(Tuple(collect(summary.parameters[k]) for k in 1:Nensemble)...)'

    push!(diags, diag(cov(θ)))
    push!(times, summary.pseudotime)
end
