filename = "viz_he_18_hp_4_ve_6_vp_4_roefanov.jld2"
# write function to reconstruct grid
function reconstruct_info(filename)
    split_string = split(split(filename, "he_")[2], "_")
    he = parse(Int, split_string[1])
    hp = parse(Int, split_string[3])
    ve = parse(Int, split_string[5])
    vp = parse(Int, split_string[7])
    return (;he=he, hp = hp, ve = ve, vp = vp)
end
gridinfo = reconstruct_info(filename)

filename = "small_earth.jld2"
fhdof = 90
vdof = 30
hp = 4
vp = 4
he = ceil(Int, fhdof / (hp+1))
ve = ceil(Int, vdof / (vp+1))
gridinfo = (;he=he, hp = hp, ve = ve, vp = vp)
include("../smallearth/smallearth.jl")
#=
include("../src/interface/domains.jl")
include("../src/interface/models.jl")
include("../src/interface/physics.jl")
include("../src/interface/boundary_conditions.jl")
include("../src/interface/grids.jl")
include("../src/backends/backends.jl")
include("../src/interface/timestepper_abstractions.jl")
include("../src/backends/dg_model_backends/backend_hook.jl")
include("../src/interface/simulations.jl")
include("../src/interface/callbacks.jl")
include("../src/backends/dg_model_backends/boilerplate.jl")
include("../src/utils/sphere_utils.jl")
=#
# to be removed

#=
using CLIMAParameters #: AbstractEarthParameterSet
struct PlanetParameterSet <: AbstractEarthParameterSet end
get_planet_parameter(p::Symbol) = getproperty(CLIMAParameters.Planet, p)(PlanetParameterSet())

parameters = (
    a    = get_planet_parameter(:planet_radius),
    Ω    = get_planet_parameter(:Omega),
    g    = get_planet_parameter(:grav),
    κ    = get_planet_parameter(:kappa_d),
    R_d  = get_planet_parameter(:R_d),
    cv_d = get_planet_parameter(:cv_d),
    cp_d = get_planet_parameter(:cp_d),
    γ    = get_planet_parameter(:cp_d)/get_planet_parameter(:cv_d),
    H    = 30e3,
    pₒ   = 1.0e5,
    k    = 3.0,
    Γ    = 0.005,
    T_0  = 0.0,
    T_E  = 310.0,
    T_P  = 240.0,
    b    = 2.0,
    z_t  = 15e3,
    λ_c  = π / 9,
    ϕ_c  = 2 * π / 9,
    V_p  = 1.0,
    day = 86400,
    p0 = 1e5,
    T_ref = 255,
)
=#

domain = SphericalShell(radius = parameters.a, height = parameters.H)
discretized_domain = DiscretizedDomain(
    domain = domain,
    discretization = (
	    horizontal = SpectralElementGrid(elements = gridinfo.he, polynomial_order = gridinfo.hp),
	    vertical = SpectralElementGrid(elements = gridinfo.ve, polynomial_order = gridinfo.vp)
	),
)

backend = DiscontinuousGalerkinBackend(numerics = (flux = :roefanov,),)
grid = create_grid(backend, discretized_domain)


jl_file = jldopen(filename, "r+")
it_keys = keys(jl_file["state"])

test_state = ClimateMachine.CUDA.CuArray(jl_file["state"][it_keys[1]])
# Create Interpolation Structure

vert_range = grid1d(
        domain.radius, 
        domain.radius + domain.height, 
        nothing,
        nelem = discretized_domain.discretization.vertical.elements,
)
# create new grid info
lat_grd = collect(-90:1:90) .* 1.0
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

n_states = size(test_state)[2]
istate = similar(test_state, interpol.Npl, n_states)

istate_1 = similar(test_state, interpol.Npl, 1)

split_filename = split(filename, ".")
lat_lon_filename = split_filename[1] * "_lat_lon.jld2"

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

    γ = gamma
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
JLD2.Group(lat_lon_file, "time")
JLD2.Group(lat_lon_file, "p")
JLD2.Group(lat_lon_file, "T")
JLD2.Group(lat_lon_file, "grid")
lat_lon_file["grid"]["latitude"] = lat_grd
lat_lon_file["grid"]["longitude"] = long_grd
lat_lon_file["grid"]["radius"] = rad_grd

for it in it_keys
    println("starting ", it)
    println("done with  ", parse(Int,it) / parse(Int,it_keys[end]) * 100, "%")
    # interpolate data
    local test_state = ClimateMachine.CUDA.CuArray(jl_file["state"][it])
    interpolate_local!(interpol, test_state, istate) 
    # istate has linear indexing
    # project to usual velocity components 
    project_cubed_sphere!(interpol, istate, (_ρu, _ρv, _ρw))
    # longitude latitude radial state
    local all_state_data = Array(accumulate_interpolated_data(MPI.COMM_WORLD, interpol, istate))
    for state_index in 1:5
        lat_lon_file[states[state_index]][it] = all_state_data[:,:,:,state_index]
    end
    lat_lon_file["time"][it] = jl_file["time"][it]

    # now calculate pressure. Do so on DG grid, then interpolate
    t_state = Array(test_state)
    t_r = Array(r)
    local p = pressure(t_state, t_r)
    local T = p ./ t_state[:, _ρ ,:] ./ parameters.R_d
    
    r_p = reshape(p, (size(p)[1], 1, size(p)[2]))
    r_T = reshape(T, (size(p)[1], 1, size(p)[2]))

    interpolate_local!(interpol, ClimateMachine.CUDA.CuArray(r_p), istate_1)
    i_p = Array(accumulate_interpolated_data(MPI.COMM_WORLD, interpol, istate_1))
    lat_lon_file["p"][it] = i_p
    interpolate_local!(interpol, ClimateMachine.CUDA.CuArray(r_T), istate_1)
    i_T = Array(accumulate_interpolated_data(MPI.COMM_WORLD, interpol, istate_1))
    lat_lon_file["T"][it] = i_T
end
close(lat_lon_file)
