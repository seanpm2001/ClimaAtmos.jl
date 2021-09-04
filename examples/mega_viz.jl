filename = []
push!(filename, "avg_long_hs_he_23_hp_1_ve_7_vp_2_lat_lon.jld2")
push!(filename, "avg_long_hs_he_15_hp_2_ve_7_vp_2_lat_lon.jld2")
push!(filename, "avg_long_hs_he_12_hp_3_ve_7_vp_2_lat_lon.jld2")
push!(filename, "avg_long_hs_he_9_hp_4_ve_7_vp_2_lat_lon.jld2")
push!(filename, "avg_long_hs_he_8_hp_5_ve_7_vp_2_lat_lon.jld2")
push!(filename, "avg_long_hs_he_7_hp_6_ve_7_vp_2_lat_lon.jld2")

jl_file = jldopen(filename[1], "r+")
t_keys = keys(ρ_file)

lat_grid = jl_file["grid"]["latitude"]
lon_grid = jl_file["grid"]["longitude"]
rad_grid = jl_file["grid"]["radius"]

using GLMakie
using Statistics

λ = lon_grid 
ϕ = lat_grid
r = rad_grid 

p0 = 1e5 # pressure bottom
pH = 800 # pressure top
H = rad_grid[end] - rad_grid[1] # height of domain
sH = -H / log(pH/p0) # scale height
# roughly p = p0 * exp(-z / sH), so log(p / p0) * sH = z
p_coord = p0 * exp.( -(rad_grid .- rad_grid[1]) ./ sH)

# makes the usual plot

function typical_plot!(ax, ϕ, p_coord, slice_zonal; title = "Zonal Velocity [m/s]", random_seed = 300, markersize = 10)
    cplot_p = contour!(ax, ϕ, -p_coord, slice_zonal, levels= collect(4:4:28), colormap = :reds)
    cplot_n = contour!(ax, ϕ, -p_coord, slice_zonal, levels = [-8, -4], colormap = :blues)
    cplot = contour!(ax, ϕ, -p_coord, slice_zonal, levels = [0], color = :purple, linewidth = 3.0, visible = false)
    ax.title = title
    ax.titlesize = 40
    ax.xlabel = "Latitude"
    ax.ylabel = "Stretched Height"
    ax.xlabelsize = 25
    ax.ylabelsize = 25 
    ax.xticks = ([-80, -60, -30, 0, 30, 60, 80], ["80S", "60S", "30S", "0", "30N", "60N", "80N"])

    pressure_levels = [1000, 850, 700, 550, 400, 250, 100, 10]
    ax.yticks = (pressure_levels .* -1e2, string.(pressure_levels))

    # don't judge me
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

    Random.seed!(random_seed)
    for contour_index in 1:length(contour_levels)
        local contour_val = contour_levels[contour_index]
        local segments = list_o_stuff[contour_index].segments
        local indices = [0, list_o_stuff[contour_index].index_vals[1:end]...]
        for i in 1:length(indices)-1
            local index1 = rand(indices[i]+1:indices[i+1]-1) # choose random point on segment
            local index2 = round(Int, 0.5 * indices[i] + 0.5 * indices[i+1]) # choose point in middle
            local index  = rand([index1, index2]) # choose between random point and point in middle
            local location = Point3(segments[index]..., 2f0)
            local sc = scatter!(ax, location, markersize= markersize, align = (:center, :center), color=(:white, 0.1), strokecolor=:white)
            local anno = text!(ax, [("$contour_val", location)], align = (:center, :center), textsize =  markersize)
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
    return nothing
end


fig = Figure(resolution = (1700, 1000))

s_string = "u"

i = 1
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[i], "r+")
slice_zonal = mean(jl_file[s_string]["0"][1:end-1, :, :], dims = 1)[1,:,:]
ax1 = fig[jj,ii] = Axis(fig)
typical_plot!(ax1, ϕ, p_coord, slice_zonal, title = "p = $i")
hidexdecorations!(ax2, grid = false)

i = 2
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[i], "r+")
slice_zonal = mean(jl_file[s_string]["0"][1:end-1, :, :], dims = 1)[1,:,:]
ax2 = fig[jj,ii] = Axis(fig)
typical_plot!(ax2, ϕ, p_coord, slice_zonal, title = "p = $i")
# hidexdecorations!(ax2, grid = false)
hideydecorations!(ax2, grid = false)

i = 3
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[i], "r+")
slice_zonal = mean(jl_file[s_string]["0"][1:end-1, :, :], dims = 1)[1,:,:]
ax3 = fig[jj,ii] = Axis(fig)
typical_plot!(ax3, ϕ, p_coord, slice_zonal, title = "p = $i")
# hidexdecorations!(ax3, grid = false)
hideydecorations!(ax3, grid = false)

i = 4
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[i], "r+")
slice_zonal = mean(jl_file[s_string]["0"][1:end-1, :, :], dims = 1)[1,:,:]
ax4 = fig[jj,ii] = Axis(fig)
typical_plot!(ax4, ϕ, p_coord, slice_zonal, title = "p = $i")

i = 5
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[i], "r+")
slice_zonal = mean(jl_file[s_string]["0"][1:end-1, :, :], dims = 1)[1,:,:]
ax5 = fig[jj,ii] = Axis(fig)
typical_plot!(ax5, ϕ, p_coord, slice_zonal, title = "p = $i")
hideydecorations!(ax5, grid = false)

i = 6
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[i], "r+")
slice_zonal = mean(jl_file[s_string]["0"][1:end-1, :, :], dims = 1)[1,:,:]
ax6 = fig[jj,ii] = Axis(fig)
typical_plot!(ax6, ϕ, p_coord, slice_zonal, title = "p = $i")
hideydecorations!(ax6, grid = false)

