using JLD2
#=
filename = "hs_data.jld2"
jl_file = jldopen(filename, "r+")
all_state_data = jl_file["all_state_data"]
heatmap(all_state_data[:,:,1,end], colormap = :balance, interpolate = true)
=#

filename = "hs_lat_lon.jld2"
filename = "hs_he_15_hp_5_ve_7_vp_2_lat_lon.jld2"
filename = "hs_he_30_hp_2_ve_7_vp_2_lat_lon.jld2"
filename = "hs_he_30_hp_3_ve_10_vp_2_lat_lon.jld2"
filename = "hs_he_11_hp_3_ve_7_vp_2_lat_lon.jld2"
filename = "hs_he_22_hp_1_ve_7_vp_2_lat_lon.jld2"
# filename = "hs_he_9_hp_4_ve_7_vp_2_lat_lon.jld2"
jl_file = jldopen(filename, "r+")
ρ_file = jl_file["ρ"]
ρu_file = jl_file["ρu"]
ρv_file = jl_file["ρv"]
ρw_file = jl_file["ρw"]
ρe_file = jl_file["ρe"]
p_file = jl_file["p"]
T_file = jl_file["T"]
t_keys = keys(ρ_file)

lat_grid = jl_file["grid"]["latitude"]
lon_grid = jl_file["grid"]["longitude"]
rad_grid = jl_file["grid"]["radius"]

# lat_grd = collect(-89:1:89) .* 1.0
# long_grd = collect(-180:1:180) .* 1.0

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
ha_z = [m_r[k]*eps(1.0)*0     for j in eachindex(ϕ), k in eachindex(r)]



t_index = Node(30)
t_key = @lift(t_keys[$t_index])
state_file = T_file
ρ  = @lift(state_file[$t_key])
global_clims = quantile.(Ref(state_file[t_keys[end]][:]), [0.05, 0.95])
# symmetrize about 0
maxc = maximum(abs.(global_clims))
global_clims = (-maxc, maxc)
# global clims

fig = Figure(resolution = (1100, 800)) 
# slices
axρ1 = fig[2,1:3] = LScene(fig)
ϕ_eq = argmin(abs.(ϕ .- 45)) 

slice1 = @lift($ρ[:,ϕ_eq,:])
clims1 = @lift(quantile.(Ref($slice1[:]), [0.05,0.95]))
# clims1 = global_clims
surface!(axρ1, a_x, a_y, a_z, color = slice1, colorrange = clims1, colormap = :balance, shading = false, show_axis=false)
fig[1,2] = Label(fig, "ρv: lat slice at 45ᵒ", textsize = 30) 
rotate_cam!(fig.scene.children[1], (2*π/3, 0, 0))
update!(fig.scene)

axρ2 = fig[2,1+3:3+3] = LScene(fig)
λ_eq = argmin(abs.(λ .- 25)) 

slice2 = @lift($ρ[λ_eq,:,:])
clims2 = @lift(quantile.(Ref($slice2[:]), [0.05,0.95]))
# clims2 = global_clims
surface!(axρ2, ha_x, ha_y, ha_z, color = slice2, colorrange = clims2, colormap = :balance, shading = false, show_axis=false)
fig[1,2+3] = Label(fig, "ρv: lon slice at 25ᵒ", textsize = 30) 
# rotate_cam!(fig.scene.children[2], (2*π/3, 0, 0))
# update!(fig.scene)

axρ3 = fig[2,1+3*2:3+3*2] = LScene(fig)
λ_eq = argmin(abs.(λ .- 0)) 
slice3 = @lift(mean($ρ[1:end-1,:,:], dims=1)[1,:,:])
clims3 = @lift(quantile.(Ref($slice3[:]), [0.05,0.95]))
# clims3 = global_clims
surface!(axρ3, ha_x, ha_y, ha_z, color = slice3, colorrange = clims3, colormap = :balance, shading = false, show_axis=false)
fig[1,2+3*2] = Label(fig, "ρv: zonal avg", textsize = 30) 
# rotate_cam!(fig.scene.children[3], (2*π/3, 0, 0))
# update!(fig.scene)

# sphere 
height_index = 1
axρ4 = fig[4,1:3] = LScene(fig) 
fig[3,2] = Label(fig, "ρv: sphere ", textsize = 30) 
# axρ4 = fig[3:4,1:3] = Axis3(fig, title = "ρv: sphere", titlesize = 30, show_axis=false) 
slice4 = @lift($ρ[:,:, height_index])
clims4 = @lift(quantile.(Ref($slice4[:]), [0.05,0.95]))
# clims4 = global_clims
surface!(axρ4, x, y, z, color = slice4, colormap = :balance, colorrange = clims4, shading = false, show_axis = false)


# left right bottom top for ax.padding = (0, 6, 16, 0)
# axρ5 = fig[4,1+3*2:3+3*2] = LScene(fig) 
axρ5 = fig[3:4,1+3:3+3*2] = Axis(fig, title = "ρv:lat lon", titlesize = 30) 
# fig[3,2+3*2] = Label(fig, "ρv: latlon ", textsize = 30) 
scene5 = heatmap!(axρ5, λ, ϕ, slice4, colormap = :balance, colorrange = clims4, interpolate = true, shading = false, show_axis=false)
# scene5.padding = (0, 0, 0, 00)
# update!(fig.scene)

#=
tic = time()
iterations = 1:length(t_keys)
record(fig, "makiehs_sphere.mp4", iterations, framerate=30) do i
    t_index[] = i
    println("finishing ", i)
end
toc = time()
println("the time for generating the figure is ", (toc - tic)/60, "minutes")
=#

#=
t_index = Node(1)
t_key = @lift(t_keys[$t_index])
state = @lift(ρu_file[$t_key][:,:,1])
clims = quantile.(Ref(ρu_file[t_keys[end]][:]), [0.1,0.90])
fig = surface(x, y, z, color = state, colormap = :balance, interpolate = true, shading = false, show_axis=false)
rotate_cam!(fig.figure.scene.children[1], (π/2, π/6, π/3))
# movietime

iterations = 1:length(t_keys)
record(fig.figure, "makiehs_sphere.mp4", iterations, framerate=30) do i
    t_index[] = i
    println("finishing ", i)
end
=#

#=
fig = Figure() # resolution = [750, 450]

t_index = Node(1800)
t_key = @lift(t_keys[$t_index])
height_index = 1 # 1:31
ρ  = @lift(ρ_file[$t_key][:,:,height_index])
ρu = @lift(ρu_file[$t_key][:,:,height_index])
ρv = @lift(ρv_file[$t_key][:,:,height_index])
ρw = @lift(ρw_file[$t_key][:,:,height_index])
ρe = @lift(ρe_file[$t_key][:,:,height_index])
e = @lift($ρ ./ $ρ)
# clims = quantile.(Ref(ρu_file[t_keys[end]][:]), [0.1,0.90])

ρclims = quantile.(Ref(ρ_file[t_keys[end]][:,:,height_index][:]), [0.05,0.95])
eclims = @lift(quantile.(Ref($ρe[:] ./ $ρ[:]), [0.05,0.95]))


uclims = quantile.(Ref(ρu_file[t_keys[end]][:,:,height_index][:]), [0.05,0.95])
vclims = quantile.(Ref(ρv_file[t_keys[end]][:,:,height_index][:]), [0.05,0.95])
wclims = quantile.(Ref(ρw_file[t_keys[end]][:,:,height_index][:]), [0.05,0.95])
# fig[0,:] = Label(fig, "Held-Suarez", textsize = 50) 


fig[1,2] = Label(fig, "ρ", textsize = 30) 
fig[1,5] = Label(fig, "e", textsize = 30) 

fig[3,2] = Label(fig, "ρu", textsize = 30) 
fig[3,5] = Label(fig, "ρv", textsize = 30) 
fig[3,8] = Label(fig, "ρw", textsize = 30) 

axρ = fig[2,1:3] = LScene(fig) 
surface!(axρ, x, y, z, color = ρ, colorrange = ρclims, colormap = :balance, interpolate = true, shading = false, show_axis=false)

axρe = fig[2,4:6] = LScene(fig)
surface!(axρe, x, y, z, color = e,  colormap = :balance, interpolate = true, shading = false, show_axis=false)

axρu = fig[4,1:3] = LScene(fig)
surface!(axρu, x, y, z, color = ρu, colorrange = uclims, colormap = :balance, interpolate = true, shading = false, show_axis=false)
axρv = fig[4,4:6] = LScene(fig)
surface!(axρv, x, y, z, color = ρv, colorrange = vclims, colormap = :balance, interpolate = true, shading = false, show_axis=false)
axρw = fig[4,7:9] = LScene(fig)
surface!(axρw, x, y, z, color = ρw, colorrange = wclims, colormap = :balance, interpolate = true, shading = false, show_axis=false)

for i in 1:5
    rotate_cam!(fig.scene.children[i], (-π/6, 0, 0))
    rotate_cam!(fig.scene.children[i], (0, -π/8, 0))
end
=#
#=
iterations = 1:length(t_keys)
record(fig, "makiehs_sphere.mp4", iterations, framerate=30) do i
    t_index[] = i
    println("finishing ", i)
    println("done with ", i/length(t_keys)*100, " percent")
end
=#
#=
t_key = t_keys[end]
cρ = ρ_file[t_key]
avg_cρ = mean(cρ[1:360, :, :], dims = 1)[1,:,:]

cu = ρu_file[t_key] ./ cρ
avg_cu = mean(cu[1:360, :, :], dims = 1)[1,:,:]
heatmap(avg_cu, interpolate = true, colormap = :balance)

ce = ρe_file[t_key] ./ cρ
avg_ce = mean(ce[1:360, :, :], dims = 1)[1,:,:]
heatmap(avg_ce, interpolate = true, colormap = :balance)

cw = ρw_file[t_key] ./ cρ
avg_cw = mean(cw[1:360, :, :], dims = 1)[1,:,:]
heatmap(avg_cw, interpolate = true, colormap = :balance)

# contour(avg_ce)
contour(avg_ce, levels = 10, color = :black)
=#
#=

function pressure(ρ, ρu, ρv, ρw, ρe, radius)
    γ = 1.4
    g = 9.8
    ϕ = radius * g
    p = (γ-1) * (ρe - 0.5 * (ρu^2 + ρv^2 + ρw^2) / ρ - ρ*ϕ)
    return p
end 

t_key = t_keys[end]
ρ = ρ_file[t_key]
ρu = ρu_file[t_key]
ρv = ρv_file[t_key]
ρw = ρw_file[t_key]
ρe = ρe_file[t_key]

avg_ρ = mean(ρ[1:360, :, :], dims = 1)[1,:,:]

p = pressure.(ρ, ρu, ρv, ρw, ρe, radius)

avg_p = mean(p[1:360, :, :], dims = 1)[1,:,:]
avgp = avg_p .- minimum(avg_p)
heatmap(avg_p , interpolate = true, colormap = :balance)


avg_T = mean((p ./ ρ)[1:360, :, :], dims = 1)[1,:,:]
# heatmap(avg_T , interpolate = true, colormap = :balance)
contour(avg_T , levels = 10)
=#
#=
filename = "hs.jld2"
ojl_file = jldopen(filename, "r+")
=#