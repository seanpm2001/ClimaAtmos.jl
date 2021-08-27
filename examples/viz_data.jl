using JLD2
#=
filename = "hs_data.jld2"
jl_file = jldopen(filename, "r+")
all_state_data = jl_file["all_state_data"]
heatmap(all_state_data[:,:,1,end], colormap = :balance, interpolate = true)
=#

filename = "hs_lat_lon.jld2"
jl_file = jldopen(filename, "r+")
ρ_file = jl_file["ρ"]
ρu_file = jl_file["ρu"]
ρv_file = jl_file["ρv"]
ρw_file = jl_file["ρw"]
ρe_file = jl_file["ρe"]
t_keys = keys(ρ_file)

lat_grd = collect(-89:1:89) .* 1.0
long_grd = collect(-180:1:180) .* 1.0

using GLMakie

#=
t_index = Node(1)
t_key = @lift(t_keys[$t_index])
state = @lift(ρe_file[$t_key][:,:,1])
fig = heatmap(long_grd,lat_grd, state, colormap = :balance, interpolate = true)

# movietime
iterations = 1:length(t_keys)
record(fig.figure, "makiehs.mp4", iterations, framerate=30) do i
    t_index[] = i
    println("finishing ", i)
end
=#
λ = long_grd 
ϕ = lat_grd
x = [cos(λ[j]) * cos(ϕ[k]) for j in 1:360, k in 1:178] # eachindex(ϕ)]
y = [sin(λ[j]) * cos(ϕ[k]) for j in 1:360, k in 1:178] # eachindex(ϕ)]
z = [sin(ϕ[k])             for j in 1:360, k in 1:178] # eachindex(ϕ)]
t_index = Node(48)
t_key = @lift(t_keys[$t_index])
# state = @lift(ρe_file[$t_key][:,:,1])
state = x .^2 + y .^2
fig = surface(x, y, z, color = state, colormap = :balance, interpolate = true, shading = false, show_axis=false)

# movietime
#=
iterations = 1:length(t_keys)
record(fig.figure, "makiehs_sphere.mp4", iterations, framerate=30) do i
    t_index[] = i
    println("finishing ", i)
end
=#