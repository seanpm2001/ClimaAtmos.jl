using JLD2

filenames = []
newfilenames = []

push!(filenames, "earth_hs_he_12_hp_4_ve_12_vp_4_roefanov_lat_lon.jld2")
push!(newfilenames, "earth_snapshots.jld2")

# push!(filenames, "small_earth_lat_lon.jld2")
# push!(newfilenames, "small_earth_snapshots.jld2")

for i in 1:1
    filename = filenames[i]
    newfilename = newfilenames[i]

    jl_file = jldopen(filename, "r+")

    jl_groups = keys(jl_file)
    state_keys = jl_groups[1:end-1]
    grid_key = jl_groups[end]

    file = jldopen(newfilename, "a+")
    for jl_group in jl_groups
        JLD2.Group(file, jl_group)
    end

    grid_names = keys(jl_file[grid_key])
    for grid_name in grid_names
        file[grid_key][grid_name] = jl_file[grid_key][grid_name]
    end

    time_keys = keys(jl_file[state_keys[1]])
    for state in state_keys 
        file[state]["0"] = jl_file[state][time_keys[end-1]]
        file[state]["1"] = jl_file[state][time_keys[end]]
    end

    close(file)
end
