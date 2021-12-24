using GLMakie
using JLD2


# filename = "earth_lat_lon.jld2"
filename = "small_earth_lat_lon.jld2"


jl_file = jldopen(filename, "r+")

t_keys = keys(jl_file["ρ"])

lat_grid = jl_file["grid"]["latitude"]
lon_grid = jl_file["grid"]["longitude"]
rad_grid = jl_file["grid"]["radius"]

using GLMakie
using Statistics

λ = lon_grid 
ϕ = lat_grid
r = rad_grid 

# surface of sphere
x = [cosd(λ[j]) * cosd(ϕ[k]) for j in eachindex(λ), k in eachindex(ϕ)]
y = [sind(λ[j]) * cosd(ϕ[k]) for j in eachindex(λ), k in eachindex(ϕ)]
z = [sind(ϕ[k])              for j in eachindex(λ), k in eachindex(ϕ)]

# annulus 
m_r = [0.5 + 0.5 * (i-1) / (length(r) -1) for i in 1:length(r)] # modified radius for plotting
a_x = [cosd(λ[j]) * m_r[k] for j in eachindex(λ), k in eachindex(r)]
a_y = [sind(λ[j]) * m_r[k] for j in eachindex(λ), k in eachindex(r)]
a_z = [m_r[k]*eps(1.0)                 for j in eachindex(λ), k in eachindex(r)]

# half annulus
ha_x = [cosd(ϕ[j]) * m_r[k] for j in eachindex(ϕ), k in eachindex(r)]
ha_y = [sind(ϕ[j]) * m_r[k] for j in eachindex(ϕ), k in eachindex(r)]
ha_z = [m_r[k]*eps(1.0)*0   for j in eachindex(ϕ), k in eachindex(r)]

# state_file[t_keys[1]]

t_index = Node(length(t_keys)) #length(t_keys) is last value
t_key = @lift(t_keys[$t_index])

s_string = "ρw" # state string
state_file =  jl_file[s_string]
# ρ  = @lift(jl_file["ee"][$t_key] - jl_file["e"][$t_key] .^2)
ρ  = @lift(state_file[$t_key])
global_clims = quantile.(Ref(state_file[t_keys[end]][:]), [0.1, 0.9])
# symmetrize about 0
# maxc = maximum(abs.(global_clims))
# global_clims = (-maxc, maxc)
# global_clims = (-13.5, 13.5)
use_global_clims = false
# global_clims = (-20, 20)
# global_clims = (84464, 90103)

fig = Figure(resolution = (1100, 800))
colormap = :balance # :thermal # :balance # :Blues_9 #:thermal
# slices
axρ1 = fig[2,1:3] = LScene(fig)
ϕ_eq = argmin(abs.(ϕ .- 45)) 

slice1 = @lift($ρ[:,ϕ_eq,:])
clims1 = @lift(quantile.(Ref($slice1[:]), [0.05,0.95]))
clims1 = use_global_clims ? global_clims : clims1 
surface!(axρ1, a_x, a_y, a_z, color = slice1, colorrange = clims1, colormap = colormap, shading = false, show_axis=false)
fig[1,2] = Label(fig, s_string * ": lat slice at 45ᵒ", textsize = 30) 
rotate_cam!(fig.scene.children[1], (2*π/3, 0, 0))
update!(fig.scene)

axρ2 = fig[2,1+3:3+3] = LScene(fig)
λ_eq = argmin(abs.(λ .- 25)) 

slice2 = @lift($ρ[λ_eq,:,:])
clims2 = @lift(quantile.(Ref($slice2[:]), [0.05,0.95]))
clims2 = use_global_clims ? global_clims : clims2 
surface!(axρ2, ha_x, ha_y, ha_z, color = slice2, colorrange = clims2, colormap = colormap, shading = false, show_axis=false)
fig[1,2+3] = Label(fig, s_string * ": lon slice at 25ᵒ", textsize = 30) 
# rotate_cam!(fig.scene.children[2], (2*π/3, 0, 0))
# update!(fig.scene)

axρ3 = fig[2,1+3*2:3+3*2] = LScene(fig)
λ_eq = argmin(abs.(λ .- 0)) 
slice3 = @lift(mean($ρ[1:end-1,:,:], dims=1)[1,:,:])
clims3 = @lift(quantile.(Ref($slice3[:]), [0.05,0.95]))
clims3 = use_global_clims ? global_clims : clims3 
surface!(axρ3, ha_x, ha_y, ha_z, color = slice3, colorrange = clims3, colormap = colormap, shading = false, show_axis=false)
fig[1,2+3*2] = Label(fig, s_string * ": zonal avg", textsize = 30) 
# rotate_cam!(fig.scene.children[3], (2*π/3, 0, 0))
# update!(fig.scene)

# sphere 
height_index = 2
axρ4 = fig[4,1:3] = LScene(fig) 
fig[3,2] = Label(fig, s_string * ": sphere at 1km ", textsize = 30) 

# axρ4 = fig[3:4,1:3] = Axis3(fig, title = "ρv: sphere", titlesize = 30, show_axis=false) 
slice4 = @lift($ρ[:,:, height_index])
clims4 = @lift(quantile.(Ref($slice4[:]), [0.05,0.95]))
clims4 = use_global_clims ? global_clims : clims4 
surface!(axρ4, x, y, z, color = slice4, colormap = colormap, colorrange = clims4, shading = false, show_axis = false)


# left right bottom top for ax.padding = (0, 6, 16, 0)
# axρ5 = fig[4,1+3*2:3+3*2] = LScene(fig) # s_string * ":lat lon at 1km" * 
x = "an equation"
titlestring =  L"does this work now $\langle$T$\rangle$"
axρ5 = fig[3:4,1+3:3+3*2] = Axis(fig, title = titlestring, titlesize = 30) 
# fig[3,2+3*2] = Label(fig, "ρv: latlon ", textsize = 30) 
scene5 = heatmap!(axρ5, λ, ϕ, slice4, colormap = colormap, colorrange = clims4, interpolate = true, shading = false, show_axis=false)
# scene5.padding = (0, 0, 0, 00)
# update!(fig.scene)
update!(fig.scene)