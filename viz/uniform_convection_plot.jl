using JLD2, GLMakie
c_file = jldopen("interpolated_convection_3.jld2", "a+")
c_file = jldopen("/home/sandre/Repositories/ClimaAtmos.jl/new_interpolated_convection_1.jld2", "a+")
c_file = jldopen("/home/sandre/Repositories/ClimaAtmos.jl/new_interpolated_convection_1.jld2", "a+")

time_keys = keys(c_file["state"])


time_node = Node(1)

_ρ = 1
state_base = c_file["state"][time_keys[1]][:, :, :, _ρ]  # should be ρ
x = ib_file["grid"]["x"]
y = ib_file["grid"]["y"]
z = ib_file["grid"]["z"]



state = @lift(c_file["state"][time_keys[$time_node]][:, 1, :, 4])
# state = @lift(sum(ib_file["state"][time_keys[$time_node]][:, :, :, _ρ], dims = 2)[:, 1, :])


fig, hm, ax = heatmap(x, z, state, colormap = :balance)

display(fig)


for i in 1:1:length(time_keys)
    sleep(0.1)
    time_node[] = i
end
