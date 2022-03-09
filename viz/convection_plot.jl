using GLMakie, JLD2, Statistics

# jl_file = jldopen("/home/sandre/Repositories/ClimaAtmos.jl/interpolated_convection_3.jld2", "a+")
# jl_file = jldopen("/home/sandre/Repositories/ClimaAtmos.jl/new_interpolated_convection_1.jld2", "a+")
# jl_file = jldopen("/home/sandre/Repositories/ClimaAtmos.jl/new_interpolated_convection_more_gridpoints.jld2", "a+")
jl_file = jldopen("/home/sandre/Repositories/ClimaAtmos.jl/new_interpolated_convection_just_tryin.jld2", "a+")
tkeys = keys(jl_file["state"])
# state = jl_file["state"][tkeys[end]]

fig = Figure(resolution = (600, 600));
ax = LScene(fig[1, 1], scenekw = (camera = cam3d!, show_axis = true));


tindex = Observable(1)
state = @lift(jl_file["state"][tkeys[$tindex]])

ρ = @lift($state[:, :, :, 1])
w = @lift($state[:, :, :, 4])
T = @lift($state[:, :, :, 7])
p = @lift($state[:, :, :, 6])
θ = @lift($T .* (1e5 ./ $p) .^ (0.286))
θ̅ = @lift($θ .- mean($θ, dims = (1, 2)))

cmap = :balance # :Blues_9
cmapa = RGBAf.(to_colormap(cmap), 1);
cmap = vcat(cmapa[1:15], fill(RGBAf(0, 0, 0, 0), 10), cmapa[25:end])
# cmap = vcat(fill(RGBAf(0, 0, 0, 0), 15), cmapa[25:end])
#=
cmap = vcat(cmapa[1:5], fill(RGBAf(0, 0, 0, 0), 5))
for i in 1:7
    cmap = vcat(cmap..., cmapa[5*i+1:5*i+5], fill(RGBAf(0, 0, 0, 0), 5))
end
=#
# wmax = maximum(abs.(w))
clims = (-7, 7) # w
plot_state = w
# clims = (351.71910286663194, 424.3662951232513)# @lift(extrema($plot_state))

v1 = volume!(ax, plot_state,
    colorrange = clims, algorithm = :absorption, absorption = 10.0f0,
    colormap = cmap)
tindex[] = length(tkeys)

# scatter(mean(θ[], dims = (1,2))[1,1,:], collect(1:128))

#=
framerate = 30
record(fig, "theta_convection.mp4", 1:length(tkeys), framerate = framerate) do t
    if t == 1
        # zoom!(ax.scene, 0.8)
        # zoom!(ax.scene, 0.8)
    end
    tindex[] = t
    println("at time t= ", t)
end
=#