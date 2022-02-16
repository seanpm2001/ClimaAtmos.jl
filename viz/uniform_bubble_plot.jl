using JLD2, GLMakie
ib_file = jldopen("interpolated_rising_bubble_4.jld2", "a+")

time_keys = keys(ib_file["state"])


time_node = Node(1)

_ρ = 1
state_base = ib_file["state"][time_keys[1]][:, 1, :, _ρ]  # should be ρ
x = ib_file["grid"]["x"]
z = ib_file["grid"]["y"]

state = @lift(ib_file["state"][time_keys[$time_node]][:, 1, :, _ρ] - state_base)


fig, hm, ax = heatmap(x, z, state, colormap = :balance )

display(fig)

for i in 1:length(time_keys)
    sleep(0.1)
    time_node[] = i
end