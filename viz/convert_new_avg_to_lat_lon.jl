using JLD2
filename = jld_filepath

prepend_name = "avg_" 
jl_file = jldopen(prepend_name * filename, "r+")

moment_1 = ClimateMachine.CUDA.CuArray(jl_file["moment_1"])
moment_2 = ClimateMachine.CUDA.CuArray(jl_file["moment_2"])
test_state = ClimateMachine.CUDA.CuArray(moment_1 ./ jl_file["times"] )  # compute average
# Create Interpolation Structure
using ClimateMachine.Mesh.Interpolation

# grab old grid info
grid = simulation.rhs[1].grid

vert_range = grid1d(
        domain.radius, 
        domain.radius + domain.height, 
        nothing,
        nelem = discretized_domain.discretization.vertical.elements,
)
# create new grid info
lat_grd = collect(-89:1:89) .* 1.0
long_grd = collect(-180:1:180) .* 1.0
rad_grd = collect(domain.radius:1e3:(domain.radius + domain.height)) .* 1.0
tic = time()
println("creating interpolation object")
# first three arguments are from the old grid 
# next three arguments are fdor the new grid
interpol = InterpolationCubedSphere(
        grid,
        vert_range,
        discretized_domain.discretization.horizontal.elements,
        lat_grd,
        long_grd,
        rad_grd
) 
toc = time()
println("done creating interpolation object: took ", toc - tic, " seconds")

bl = simulation.rhs[1]
n_states_1 = size(moment_1)[2]
istate_1 = similar(test_state, interpol.Npl, n_states_1)

n_states_2 = size(moment_2)[2]
istate_2 = similar(test_state, interpol.Npl, n_states_2)

split_filename = split(filename, ".")
lat_lon_filename = prepend_name * split_filename[1] * "_lat_lon.jld2"

lat_lon_file = jldopen(lat_lon_filename, "a+")

# convention for lat-lon filename
moment_1_names = jl_file["moment_1_names"]
moment_2_names = jl_file["moment_2_names"]

for jld_state in moment_1_names
    JLD2.Group(lat_lon_file, jld_state)
end

for jld_state in moment_2_names
    JLD2.Group(lat_lon_file, jld_state)
end

JLD2.Group(lat_lon_file, "grid")
lat_lon_file["grid"]["latitude"] = lat_grd
lat_lon_file["grid"]["longitude"] = long_grd
lat_lon_file["grid"]["radius"] = rad_grd


# interpolate data
it = "0"

# first moment
test_state = ClimateMachine.CUDA.CuArray(moment_1 ./ jl_file["times"] )  # compute average
interpolate_local!(interpol, test_state, istate_1) 
# istate has linear indexing
# project to usual velocity components 
_ρu, _ρv, _ρw = 2, 3, 4
project_cubed_sphere!(interpol, istate_1, (_ρu, _ρv, _ρw))
# longitude latitude radial state
all_state_data = Array(accumulate_interpolated_data(MPI.COMM_WORLD, interpol, istate_1))
for state_index in 1:length(moment_1_names)
    lat_lon_file[moment_1_names[state_index]][it] = all_state_data[:,:,:,state_index]
    println("extrema of ", moment_1_names[state_index], " is ", extrema(all_state_data[:,:,:,state_index]))
end

# second moment 
test_state = ClimateMachine.CUDA.CuArray(moment_2 ./ jl_file["times"] )  # compute average
interpolate_local!(interpol, test_state, istate_2) 
# istate has linear indexing
# DONT project to usual velocity components 
# project_cubed_sphere!(interpol, istate, (_ρu, _ρv, _ρw))
# longitude latitude radial state
all_state_data_2 = Array(accumulate_interpolated_data(MPI.COMM_WORLD, interpol, istate_2))
for state_index in 1:length(moment_2_names)
    lat_lon_file[moment_2_names[state_index]][it] = all_state_data_2[:,:,:,state_index]
    println("extrema of ", moment_2_names[state_index], " is ", extrema(all_state_data_2[:,:,:,state_index]))
end

close(lat_lon_file)