using JLD2, GLMakie
ib_file = jldopen("interpolated_rising_bubble_5.jld2", "a+")

time_keys = keys(ib_file["state"])


time_node = Node(1)

_ρ = 4
state_base = sum(ib_file["state"][time_keys[1]][:, :, :, _ρ], dims = 2)[1:1, 1, :]  # should be ρ
x = ib_file["grid"]["x"]
z = ib_file["grid"]["z"]

# density anamoly
# state = @lift(state_base .- sum(ib_file["state"][time_keys[$time_node]][:, :, :, _ρ], dims = 2)[:, 1, :])
state = @lift(sum(ib_file["state"][time_keys[$time_node]][:, :, :, _ρ], dims = 2)[:, 1, :])

fig, hm, ax = heatmap(x, z, state, colormap = :balance)

display(fig)

for i in 1:4:length(time_keys)
    sleep(0.1)
    time_node[] = i
end

#=
framerate = 30

record(fig, "rising_bubble_animation.mp4", 1:200,
    framerate = framerate) do t
    time_node[] = t
end

=#