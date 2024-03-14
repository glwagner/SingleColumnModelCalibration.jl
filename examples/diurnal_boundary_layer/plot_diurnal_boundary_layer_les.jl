using GLMakie
using Oceananigans
using Oceananigans.Units

filename = "diurnal_boundary_layer_les_averages.jld2"
Bt = FieldTimeSeries(filename, "B")

t = Bt.times
Nt = length(t)
x, y, z = nodes(Bt)
Nx, Ny, Nz = size(Bt.grid)

# Compute mixed layer depth
h = zeros(Nt)
Δb = 1e-4
for n = 1:Nt
    # Surface buoyancy
    b₀ = Bt[1, 1, Nz, n]
    for k = Nz-1:-1:1
        bk = Bt[1, 1, k, n]
        # b is _decreasing_ downward
        if (b₀ - bk) >= Δb || k == 1
            h[n] = - z[k]
            break
        end
    end
end

fig = Figure()

axh = Axis(fig[1, 1])
axb0 = Axis(fig[2, 1])
axb = Axis(fig[3, 1])

lines!(axh, t ./ hours, h)
lines!(axb0, t ./ hours, interior(Bt, 1, 1, Nz, :))
contourf!(axb, t ./ hours, z, interior(Bt, 1, 1, :, :)')

ylims!(axb, -65, 0)

display(fig)
