using GLMakie
using Statistics

n, m = 100, 101
t = range(0, 1, length=m)
X = cumsum(randn(n, m), dims=2)
X = X .- X[:, 1]
μ = vec(mean(X, dims=1))
σ = vec(std(X, dims=1))

left  = [Point2f(μ[i] - σ[i], t[i]) for i=1:n]
right = [Point2f(μ[i] + σ[i], t[i]) for i=1:n]

fig = Figure()
ax = Axis(fig[1, 1])
band!(ax, left, right)

display(fig)
