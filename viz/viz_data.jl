using JLD2, GLMakie
using JLD2

filename = "earth_hs_he_5_hp_5_ve_5_vp_5_roefanov_lat_lon.jld2"
# filename = "earth_lat_lon.jld2"
filename = "small_earth_lat_lon.jld2"
##
#=
# save reduced data 
file = jldopen("reduced_data", "a+")
JLD2.Group(file, "grid")
JLD2.Group(file, "T")
for t_key in keys(jl_file["T"])
    file["T"][t_key] = jl_file["T"][t_key][:,:,1:4]
end
for grid_key in keys(jl_file["grid"])
    file["grid"][grid_key] = jl_file["grid"][grid_key]
end
=#
##

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
titlestring =  L" \langle T \rangle"
axρ5 = fig[3:4,1+3:3+3*2] = Axis(fig, title = titlestring, titlesize = 30) 
# fig[3,2+3*2] = Label(fig, "ρv: latlon ", textsize = 30) 
scene5 = heatmap!(axρ5, λ, ϕ, slice4, colormap = colormap, colorrange = clims4, interpolate = true, shading = false, show_axis=false)
# scene5.padding = (0, 0, 0, 00)
# update!(fig.scene)


#=
tic = time()
iterations = 1:length(t_keys)
record(fig, "temperature_smallearth_sphere.mp4", iterations, framerate=30) do i
    t_index[] = i
    println("finishing ", i)
end
toc = time()
println("the time for generating the figure is ", (toc - tic)/60, "minutes")
=#

#=
# usual height plot
s_string = "u"
slice_zonal = mean(jl_file[s_string]["0"][1:end-1, :, :], dims = 1)[1,:,:]
p_coord = mean(jl_file["p"]["0"][1:end-1, :, :], dims = (1,2))[1,1,:]
# , colorrange = [-28,28]
heatmap(ϕ, -p_coord, slice_zonal, colormap = :balance, interpolate = true, shading = false, show_axis=false)

fig, ax, hp = contour(ϕ, -p_coord, slice_zonal, levels= collect(4:4:28), colormap = :reds)
contour!(ax, ϕ, -p_coord, slice_zonal, levels = [-8, -4], colormap = :blues)
cplot = contour!(ax, ϕ, -p_coord, slice_zonal, levels = [0], color = :purple, linewidth = 3.0, visible = false)
ax.title = "Zonal Velocity [m/s]"
ax.titlesize = 40
ax.xlabel = "Latitude [ᵒ]"
ax.ylabel = "Average Pressure [hPa]"
ax.xlabelsize = 25
ax.ylabelsize = 25 
ax.xticks = ([-60, -30,0, 30, 60], ["60S", "30S", "0", "30N", "60N"])

pressure_levels = [1000, 850, 700, 550, 400, 250, 100, 10]
ax.yticks = (pressure_levels .* -1e2, string.(pressure_levels))


contour_levels = collect(-8:4:28)
list_o_stuff = []
for level in contour_levels
    fig_t, ax_t, cp_t = contour(ϕ, -p_coord, slice_zonal, levels= [level], linewidth = 0)
    segments = cp_t.plots[1][1][]
    index_vals = []
    for (i, p) in enumerate(segments)
        # the segments are separated by NaN, which signals that a new contour starts
        if isnan(p)
            push!(beginnings, segments[i-1])
            push!(index_vals, i)
        end
    end
    push!(list_o_stuff, (; segments, beginnings, index_vals))
end

for contour_index in 1:length(contour_levels)

    contour_val = contour_levels[contour_index]
    segments = list_o_stuff[contour_index].segments

    for index in list_o_stuff[contour_index].index_vals[1:end]
        location = Point3(segments[index-1]..., 2f0)
        # sc = scatter!(ax, location, markersize=30, color=(:white, 0.001), strokecolor=:white)
        anno = text!(ax, [("$contour_val", location)], textsize = 20)
        # translate!(sc, 0, 0, 1)
        # translate!(anno, 0, 0, 2)
    end
end

for i in 1:3:length(segments)
    if i in index_vals
    else
        sc = scatter!(ax, Point3(segments[i]..., 2f0), markersize=30, color=(:white, 0.001), strokecolor=:white)
        anno = text!(ax, [("$i", Point3(segments[i]..., 2f0))])
        translate!(sc, 0, 0, 1)
        translate!(anno, 0, 0, 2)
    end
end
beginnings = Point2f0[]; colors = RGBAf0[]
# First plot in contour is the line plot, first arguments are the points of the contour
segments = hp.plots[1][1][]
index_vals = []
for (i, p) in enumerate(segments)
    # the segments are separated by NaN, which signals that a new contour starts
    if isnan(p)
        push!(beginnings, segments[i-1])
        push!(index_vals, i)
    end
end
text!(ax, [("hello", Point3(segments[64]..., 2f0))])
sc = scatter!(ax, beginnings, markersize=30, color=(:white, 0.001), strokecolor=:white)
anno = text!(ax, [(string(float(i)), Point3(p..., 2f0)) for (i, p) in enumerate(beginnings)], align=(:center, :center), color=:black)
translate!(sc, 0, 0, 1)
translate!(anno, 0, 0, 2)
# Reshuffle the plot order, so that the scatter plot gets drawn before the line plot
delete!(ax, sc)
delete!(ax, cplot)
delete!(ax, anno)
push!(ax.scene, anno)
push!(ax.scene, sc)
push!(ax.scene, cplot)

display(fig)

=#

#=
x = range(-3, 3, length=200)
y = range(-2, 2, length=100)
z = @. x^2 + y'^2
fig, ax, hp = heatmap(x, y, z)
levels = 0:1:100
cplot = contour!(ax, x, y, z, color=:black, levels=levels)

beginnings = Point2f0[]; colors = RGBAf0[]
# First plot in contour is the line plot, first arguments are the points of the contour
segments = cplot.plots[1][1][]
for (i, p) in enumerate(segments)
    # the segments are separated by NaN, which signals that a new contour starts
    if isnan(p)
        push!(beginnings, segments[i-1])
    end
end
# sc = scatter!(ax, beginnings, markersize=30, color=(:white, 0.001), strokecolor=:white)
# anno = text!(ax, [(string(float(i)), Point3(p..., 2f0)) for (i, p) in enumerate(beginnings)], align=(:center, :center), color=:black)

sc = scatter!(ax, beginnings, markersize=30, color=(:white, 0.001), strokecolor=:white)
anno = text!(ax, [(string(float(i)), Point3(p..., 2f0)) for (i, p) in enumerate(beginnings)], align=(:center, :center), color=:black)


translate!(sc, 0, 0, 1)
translate!(anno, 0, 0, 2)
# Reshuffle the plot order, so that the scatter plot gets drawn before the line plot
delete!(ax, sc)
delete!(ax, cplot)
delete!(ax, anno)
push!(ax.scene, anno)
push!(ax.scene, sc)
push!(ax.scene, cplot)
# move, so that text is in front

fig
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