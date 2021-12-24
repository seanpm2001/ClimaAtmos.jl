using JLD2

filename = jld_filepath

jl_file = jldopen(filename, "r+")

avg_jl_file = jldopen("avg_" * filename, "r+")

t_keys = keys(jl_file["state"])
Q = jl_file["state"][t_keys[1]] .* 0.0

for t_key in t_keys
    Q .+= jl_file["state"][t_key]
end

Q_avg = avg_jl_file["state"]
Q - Q_avg # should be zeroish

times = avg_jl_file["times"]



#=
it_keys = keys(jl_file["state"])

test_state = jl_file["state"][it_keys[end]]
# Create Interpolation Structure
using ClimateMachine.Mesh.Interpolation

# grab old grid info
grid = simulation.rhs[1].grid

function coordinates(grid::DiscontinuousSpectralElementGrid)
    x = view(grid.vgeo, :, grid.x1id, :)   # x-direction	
    y = view(grid.vgeo, :, grid.x2id, :)   # y-direction	
    z = view(grid.vgeo, :, grid.x3id, :)   # z-direction
    return x, y, z
end
x,y,z = Array.(coordinates(grid))

r = sqrt.(x .^2 + y .^2 + z .^2)
g = parameters.g

function pressure(ρ, ρu, ρv, ρw, ρe, radius; g = parameters.g, gamma = 1.4)
    γ = 1.4
    ϕ = radius * g
    p = (γ-1) * (ρe - 0.5 * (ρu^2 + ρv^2 + ρw^2) / ρ - ρ*ϕ)
    return p
end 
Q = test_state
Q_ρ  = Q[:,1,:]
Q_ρu = Q[:,2,:]
Q_ρv = Q[:,3,:]
Q_ρw = Q[:,4,:]
Q_ρe = Q[:,5,:]

Q_p = pressure.(Q_ρ, Q_ρu, Q_ρv, Q_ρw, Q_ρe, r)
println("pressure extrema " , extrema(Q_p))

Q_T = Q_p ./ Q_ρ ./ parameters.R_d .- 273 # convert to celcius
println("temperature extrema in celcius ", extrema(Q_T))

Q_ρ0 = ρ₀ᶜᵃʳᵗ.(Ref(parameters), x,y,z)

Q_T0 = T.(Ref(parameters), lat.(x,y,z) , rad.(x,y,z)) 
Q_p0 = p.(Ref(parameters), lat.(x,y,z) , rad.(x,y,z))
Q_ρe0 = ρe.(Ref(parameters), lon.(x,y,z), lat.(x,y,z) , rad.(x,y,z))
Q_ϕ =  e_pot.(Ref(parameters), lon.(x,y,z), lat.(x,y,z) , rad.(x,y,z))
Q_u20 = e_kin.(Ref(parameters), lon.(x,y,z), lat.(x,y,z) , rad.(x,y,z))

Q_ρ0 .* parameters.R_d .* Q_T0 .- Q_p0

Q_ρe0 - Q_ρe

Q_u20 .* Q_ρ - 0.5 .* (Q_ρu .^2 + Q_ρv .^2 + Q_ρw .^2) ./ Q_ρ

Q_T1 = (Q_ρe0  ./ Q_ρ0 .- Q_u20 .- Q_ϕ ) ./ (parameters.R_d / parameters.κ - parameters.R_d)

##
filename = "hs_he_30_hp_2_ve_10_vp_2_lat_lon.jld2"
jl_file = jldopen(filename, "r+")
ρ_file = jl_file["ρ"]
ρu_file = jl_file["ρu"]
ρv_file = jl_file["ρv"]
ρw_file = jl_file["ρw"]
ρe_file = jl_file["ρe"]
t_keys = keys(ρ_file)

lat_grid = jl_file["grid"]["latitude"]
lon_grid = jl_file["grid"]["latitude"]
rad_grid = jl_file["grid"]["radius"]

t_key = t_keys[end]
ll_ρ = ρ_file[t_key]
ll_ρu = ρu_file[t_key]
ll_ρv = ρv_file[t_key]
ll_ρw = ρw_file[t_key]
ll_ρe = ρe_file[t_key]

ll_r = ll_ρ .* 0
for i in eachindex(rad_grid)
    ll_r[:,:,i] .= rad_grid[i]
end

Q_p = pressure.(ll_ρ, ll_ρu, ll_ρv, ll_ρw, ll_ρe, ll_r)
=#