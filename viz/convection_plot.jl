using GLMakie, JLD2, Statistics

jl_file = jldopen("/home/sandre/Repositories/ClimaAtmos.jl/interpolated_convection_3.jld2", "a+")
tkeys = keys(jl_file["state"])
# state = jl_file["state"][tkeys[end]]

fig = Figure(resolution = (600, 600));
ax = LScene(fig[1, 1], scenekw = (camera = cam3d!, show_axis = true));


tindex = Observable(1)
state = @lift(jl_file["state"][tkeys[$tindex]])

ρ = @lift($state[:, :, :, 1])
ρw = @lift($state[:, :, :, 4])
w = @lift($ρw ./ $ρ)

cmap = :balance # :Blues_9
cmapa = RGBAf.(to_colormap(cmap), 1);
cmap = vcat(cmapa[1:15], fill(RGBAf(0, 0, 0, 0), 10), cmapa[25:end])
# wmax = maximum(abs.(w))
clims = (-6, 6)

v1 = volume!(ax, w,
    colorrange = clims, algorithm = :absorption, absorption = 10.0f0,
    colormap = cmap)
tindex[] = length(tkeys)


framerate = 10
record(fig, "convection.mp4", 1:101, framerate = framerate) do t
    if t == 1
        # zoom!(ax.scene, 0.8)
        # zoom!(ax.scene, 0.8)
    end
    tindex[] = t
    println("at time t= ", t)
end