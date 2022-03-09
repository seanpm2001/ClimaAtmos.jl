bubble_file = jldopen(filename, "a+")

ib = InterpolationBrick(
    simulation,
    xlength = 128,
    ylength = 128,
    zlength = 128
)
geopotential = simulation.rhs.state_auxiliary.Φ

ib_file = jldopen("new_" * "interpolated_" * filename, "a+")
JLD2.Group(ib_file, "state")
JLD2.Group(ib_file, "time")
JLD2.Group(ib_file, "grid")


ib_file["grid"]["x"] = ib.x1g
ib_file["grid"]["y"] = ib.x2g
ib_file["grid"]["z"] = ib.x3g

time_keys = keys(bubble_file["state"])
state = ClimateMachine.CUDA.CuArray(bubble_file["state"][time_keys[1]])
istate = ClimateMachine.CUDA.CuArray(similar(state, ib.Npl, 7)) # 7 because, ρ, ρu, ρv, ρw, ρe, p, T
m1, vnames = get_state(state, geopotential, parameters)
ib_file["names"] = vnames


for time_key in time_keys
    println("time = ", time_key)
    state = ClimateMachine.CUDA.CuArray(bubble_file["state"][time_key])
    moment_1, names = get_state(state, geopotential, parameters)
    # istate = similar(state, ib.Npl, size(state)[2])
    interpolate_local!(ib, ClimateMachine.CUDA.CuArray(moment_1), istate) #interpolate_local!(ib, state, istate)
    new_istate = Array(accumulate_interpolated_data(MPI.COMM_WORLD, ib, istate))
    ib_file["state"][time_key] = new_istate
    ib_file["time"][time_key] = bubble_file["time"][time_key]
end

close(ib_file)