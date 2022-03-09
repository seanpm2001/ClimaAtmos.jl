bubble_file = jldopen(filename, "a+")

ib = InterpolationBrick(
    simulation,
    xlength = 64,
    ylength = 64,
    zlength = 64
)

ib_file = jldopen("interpolated_" * filename, "a+")
JLD2.Group(ib_file, "state")
JLD2.Group(ib_file, "time")
JLD2.Group(ib_file, "grid")

ib_file["grid"]["x"] = ib.x1g
ib_file["grid"]["y"] = ib.x2g
ib_file["grid"]["z"] = ib.x3g

time_keys = keys(bubble_file["state"])
for time_key in time_keys
    println("time = ", time_key)
    state = ClimateMachine.CUDA.CuArray(bubble_file["state"][time_key])
    istate = similar(state, ib.Npl, size(state)[2])
    interpolate_local!(ib, state, istate)
    new_istate = Array(accumulate_interpolated_data(MPI.COMM_WORLD, ib, istate))
    ib_file["state"][time_key] = new_istate
    ib_file["time"][time_key] = bubble_file["time"][time_key]
end

close(ib_file)