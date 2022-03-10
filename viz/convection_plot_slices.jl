using GLMakie, JLD2, Statistics

function get_edges(state)
    # since periodic in first two we take different points
    s_bottom = state[:, :, 1]
    s_top = state[:, :, end]
    s_west = state[1, :, :]
    s_east = state[64, :, :]
    s_south = state[:, 1, :]
    s_north = state[:, 64, :]
    return (; s_top, s_bottom, s_west, s_east, s_south, s_north)
end


# jl_file = jldopen("/home/sandre/Repositories/ClimaAtmos.jl/interpolated_convection_3.jld2", "a+")
# jl_file = jldopen("/home/sandre/Repositories/ClimaAtmos.jl/new_interpolated_convection_1.jld2", "a+")
# cjl_file = jldopen("/home/sandre/Repositories/ClimaAtmos.jl/new_interpolated_convection_more_gridpoints.jld2", "a+")
# jl_file = jldopen("/home/sandre/Repositories/ClimaAtmos.jl/new_interpolated_convection_just_tryin.jld2", "a+")
tkeys = keys(jl_file["state"])

state = jl_file["state"][tkeys[end]]
ρ = state[:, :, :, 1]
w = state[:, :, :, 4]
T = state[:, :, :, 7]
p = state[:, :, :, 6]
θ = T .* (1e5 ./ p) .^ (0.286)
θ̅ = θ .- mean(θ, dims = (1, 2))

x = jl_file["grid"]["x"]
y = jl_file["grid"]["y"]
z = jl_file["grid"]["z"]

# x = x .- (x[1] + x[end]) * 0.5
# y = y .- (y[1] + y[end]) * 0.5
# z = z .- (z[1] + z[end]) * 0.5

xsurf = range(x[1], x[end], length = length(x))
ysurf = range(y[1], y[end], length = length(y))
zsurf = range(z[1], z[end], length = length(z))

fig = Figure(resolution = (1000, 500));
ax = LScene(fig[1:3, 1:3], title = "Cooling", show_axis = false);

plot_state = θ
clims = quantile.(Ref(plot_state[:]), (0.05, 0.6))


s_top, s_bottom, s_west, s_east, s_south, s_north = get_edges(plot_state)

colormap = :balance


xsurf2 = [xsurf[i] for i in eachindex(xsurf), j in eachindex(zsurf)]
ysurf2 = [ysurf[j] for i in eachindex(xsurf), j in eachindex(zsurf)]
zsurf21 = [zsurf[1] for i in eachindex(xsurf), j in eachindex(zsurf)]
zsurf2end = [zsurf[end] for i in eachindex(xsurf), j in eachindex(zsurf)]

GLMakie.surface!(ax, xsurf2, ysurf2, zsurf21, color = s_bottom, colorrange = clims, colormap = colormap)
GLMakie.surface!(ax, xsurf2, ysurf2, zsurf2end, color = s_top, colorrange = clims, colormap = colormap)


ysurf2 = [ysurf[i] for i in eachindex(ysurf), j in eachindex(zsurf)]
zsurf2 = [zsurf[j] for i in eachindex(ysurf), j in eachindex(zsurf)]
xsurf21 = [x[1] for i in eachindex(ysurf), j in eachindex(zsurf)]
xsurf2end = [x[end] for i in eachindex(ysurf), j in eachindex(zsurf)]

GLMakie.surface!(ax, xsurf21, ysurf2, zsurf2, color = s_west, colorrange = clims, colormap = colormap)
GLMakie.surface!(ax, xsurf2end, ysurf2, zsurf2, color = s_east, colorrange = clims, colormap = colormap)

zsurf2 = [zsurf[j] for i in eachindex(zsurf), j in eachindex(xsurf)]
xsurf2 = [xsurf[i] for i in eachindex(zsurf), j in eachindex(xsurf)]
ysurf21 = [y[1] for i in eachindex(zsurf), j in eachindex(xsurf)]
ysurf2end = [y[end] for i in eachindex(zsurf), j in eachindex(xsurf)]

GLMakie.surface!(ax, xsurf2, ysurf21, zsurf2, color = s_south, colorrange = clims, colormap = colormap)
GLMakie.surface!(ax, xsurf2, ysurf2end, zsurf2, color = s_north, colorrange = clims, colormap = colormap)

#=
#edge1
GLMakie.surface!(ax, xsurf, zsurf, s_south, transformation = (:xz, y[1]), colorrange = clims, colormap = colormap)

# edge 2
GLMakie.surface!(ax, xsurf, zsurf, s_north, transformation = (:xz, y[end]), colorrange = clims, colormap = colormap)

# edge 3
GLMakie.surface!(ax, ysurf, zsurf, s_west, transformation = (:yz, x[1]), colorrange = clims, colormap = colormap)

# edge 4
GLMakie.surface!(ax, ysurf, zsurf, s_east, transformation = (:yz, x[end]), colorrange = clims, colormap = colormap)

# edge 5
GLMakie.surface!(ax, xsurf, ysurf, s_bottom, transformation = (:xy, z[1]), colorrange = clims, colormap = colormap)
# edge 6
GLMakie.surface!(ax, xsurf, ysurf, s_top, transformation = (:xy, z[end]), colorrange = clims, colormap = colormap)
=#
