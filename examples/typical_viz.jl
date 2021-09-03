filename = "avg_long_hs_he_15_hp_2_ve_7_vp_2_lat_lon.jld2"
# filename = "avg_long_hs_he_12_hp_3_ve_7_vp_2_lat_lon.jld2"
# filename = "avg_long_hs_he_7_hp_6_ve_7_vp_2_lat_lon.jld2"
jl_file = jldopen(filename, "r+")
t_keys = keys(ρ_file)

lat_grid = jl_file["grid"]["latitude"]
lon_grid = jl_file["grid"]["longitude"]
rad_grid = jl_file["grid"]["radius"]

using GLMakie
using Statistics

λ = lon_grid 
ϕ = lat_grid
r = rad_grid 

s_string = "u"
slice_zonal = mean(jl_file[s_string]["0"][1:end-1, :, :], dims = 1)[1,:,:]
p_coord = mean(jl_file["p"]["0"][1:end-1, :, :], dims = (1,2))[1,1,:]
# , colorrange = [-28,28]
heatmap(ϕ, -p_coord, slice_zonal, colormap = :balance, interpolate = true, shading = false, show_axis=false)

fig, ax, cplot_p = contour(ϕ, -p_coord, slice_zonal, levels= collect(4:4:28), colormap = :reds)
cplot_n = contour!(ax, ϕ, -p_coord, slice_zonal, levels = [-8, -4], colormap = :blues)
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
    local fig_t, ax_t, cp_t = contour(ϕ, -p_coord, slice_zonal, levels= [level], linewidth = 0)
    local segments = cp_t.plots[1][1][]
    local index_vals = []
    local beginnings = []
    for (i, p) in enumerate(segments)
        # the segments are separated by NaN, which signals that a new contour starts
        if isnan(p)
            push!(beginnings, segments[i-1])
            push!(index_vals, i)
        end
    end
    push!(list_o_stuff, (; segments, beginnings, index_vals))
end

#=
contour_index = 3
contour_val = contour_levels[contour_index]
segments = list_o_stuff[contour_index].segments

index = list_o_stuff[contour_index].index_vals[1]
translate = Float32.([0,0]) # translate 2 degree to the left
location = segments[index-40]
t_location = Point2(location .- translate)
d_location = Point3(location..., 3f0)
sc = scatter!(ax, [d_location], markersize=30, color=(:white, 0.1), strokecolor=:white, align = (:center, :center))
anno = text!(ax, [("$contour_val", d_location)], textsize = 20, align = (:center, :center), color = :black)
# translate!(sc, 0, 0, 1)
# translate!(anno, 0, 0, 2)

delete!(ax, sc)
delete!(ax, cplot)
delete!(ax, anno)
push!(ax.scene, anno)
push!(ax.scene, sc)
push!(ax.scene, cplot)
=#
using Random
Random.seed!(300)
for contour_index in 1:length(contour_levels)

    local contour_val = contour_levels[contour_index]
    local segments = list_o_stuff[contour_index].segments

    local indices = [0, list_o_stuff[contour_index].index_vals[1:end]...]
    for i in 1:length(indices)-1
        local index = rand(indices[i]+1:indices[i+1]-1) # choose random point on segment
        local location = Point3(segments[index-1]..., 2f0)
        local sc = scatter!(ax, location, markersize=20, align = (:center, :center), color=(:white, 0.1), strokecolor=:white)
        local anno = text!(ax, [("$contour_val", location)], align = (:center, :center), textsize = 20)
        # translate!(sc, 0, 0, 1)
        # translate!(anno, 0, 0, 2)
        delete!(ax, sc)
        delete!(ax, cplot)
        delete!(ax, cplot_n)
        delete!(ax, cplot_p)
        delete!(ax, anno)

        push!(ax.scene, anno)
        push!(ax.scene, sc)
        push!(ax.scene, cplot)
        push!(ax.scene, cplot_n)
        push!(ax.scene, cplot_p)

    end
end
