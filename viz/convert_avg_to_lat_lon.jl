using JLD2
filename = jld_filepath

prepend_name = "avg_" 
jl_file = jldopen(prepend_name * filename, "r+")
test_state = ClimateMachine.CUDA.CuArray(jl_file["state"])
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
n_states = size(test_state)[2]
istate = similar(test_state, interpol.Npl, n_states)

istate_1 = similar(test_state, interpol.Npl, 1)

split_filename = split(filename, ".")
lat_lon_filename = prepend_name * split_filename[1] * "_lat_lon.jld2"

lat_lon_file = jldopen(lat_lon_filename, "a+")

# pressure 
function coordinates(grid::DiscontinuousSpectralElementGrid)
    x = view(grid.vgeo, :, grid.x1id, :)   # x-direction	
    y = view(grid.vgeo, :, grid.x2id, :)   # y-direction	
    z = view(grid.vgeo, :, grid.x3id, :)   # z-direction
    return x, y, z
end

x,y,z = coordinates(grid)

r = sqrt.(x .^2 + y .^2 + z .^2)

function pressure(Q, radius; g = parameters.g, gamma = 1.4)
    ρ  = Q[:,1,:]
    ρu = Q[:,2,:]
    ρv = Q[:,3,:]
    ρw = Q[:,4,:]
    ρe = Q[:,5,:]

    γ = 1.4
    ϕ = radius .* g
    p = (γ-1) .* (ρe - 0.5 * (ρu .^2 + ρv .^2 + ρw .^2) ./ ρ .- ρ .* ϕ)
    return p
end 

# convention for lat-lon filename
states = ["ρ", "ρu", "ρv", "ρw", "ρe"]
_ρ , _ρu, _ρv, _ρw, _ρe = 1, 2, 3, 4, 5
for jld_state in states
    JLD2.Group(lat_lon_file, jld_state)
end
JLD2.Group(lat_lon_file, "p")
JLD2.Group(lat_lon_file, "T")
JLD2.Group(lat_lon_file, "grid")
lat_lon_file["grid"]["latitude"] = lat_grd
lat_lon_file["grid"]["longitude"] = long_grd
lat_lon_file["grid"]["radius"] = rad_grd


# interpolate data
it = "0"
test_state = ClimateMachine.CUDA.CuArray(jl_file["state"] ./ jl_file["times"] )  # compute average
interpolate_local!(interpol, test_state, istate) 
# istate has linear indexing
# project to usual velocity components 
project_cubed_sphere!(interpol, istate, (_ρu, _ρv, _ρw))
# longitude latitude radial state
all_state_data = Array(accumulate_interpolated_data(MPI.COMM_WORLD, interpol, istate))
for state_index in 1:5
    lat_lon_file[states[state_index]][it] = all_state_data[:,:,:,state_index]
end

# now calculate pressure. Do so on DG grid, then interpolate
t_state = Array(test_state)
t_r = Array(r)
press = pressure(t_state, t_r)
Temp = press ./ t_state[:, _ρ ,:] ./ parameters.R_d

r_p = reshape(press, (size(press)[1], 1, size(press)[2]))
r_T = reshape(Temp, (size(press)[1], 1, size(press)[2]))

interpolate_local!(interpol, ClimateMachine.CUDA.CuArray(r_p), istate_1)
i_p = Array(accumulate_interpolated_data(MPI.COMM_WORLD, interpol, istate_1))
lat_lon_file["p"][it] = i_p
interpolate_local!(interpol, ClimateMachine.CUDA.CuArray(r_T), istate_1)
i_T = Array(accumulate_interpolated_data(MPI.COMM_WORLD, interpol, istate_1))
lat_lon_file["T"][it] = i_T

close(lat_lon_file)