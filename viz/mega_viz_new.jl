using Random
using GLMakie
using Statistics
using JLD2


function mega_viz!(ax, ϕ, p_coord, slice_zonal, contour_levels, colorrange; add_labels = false, title_string = "", colormap = :balance, heuristic = 1, random_seed = 12345)
    cplot = contour!(ax, ϕ, p_coord, slice_zonal, levels = contour_levels, color = :black, interpolate = true,  show_axis = false)
    hm = heatmap!(ax, ϕ, p_coord, slice_zonal, colorrange = colorrange, colormap = colormap, interpolate = true)

    ax.limits = (extrema(ϕ)..., extrema(p_coord)...)

    ax.title = title_string
    ax.titlesize = 40
    ax.xlabel = "Latitude [ᵒ]"
    ax.ylabel = "Stretched Height"
    ax.xlabelsize = 25
    ax.ylabelsize = 25 
    ax.xticks = ([-80, -60, -30, 0, 30, 60, 80], ["80S", "60S", "30S", "0", "30N", "60N", "80N"])
    pressure_levels = [1000, 850, 700, 550, 400, 250, 100, 10]
    ax.yticks = (pressure_levels .* 1e2, string.(pressure_levels))
    ax.yreversed = true

    # hack 
    Random.seed!(random_seed)
    if add_labels
        list_o_stuff = []
        labeled_contours = contour_levels[1:1:end]
        for level in labeled_contours
            local fig_t, ax_t, cp_t = contour(ϕ, p_coord, slice_zonal, levels= [level], linewidth = 0)
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

        for contour_index in 1:length(labeled_contours)

            local contour_val = labeled_contours[contour_index]
            local segments = list_o_stuff[contour_index].segments

            local indices = [0, list_o_stuff[contour_index].index_vals[1:end]...]
            for i in 1:length(indices)-1
                # heuristics for choosing where on line
                local index1 = rand(indices[i]+1:indices[i+1]-1) # choose random point on segment
                local index2 = round(Int, 0.5 * indices[i] + 0.5 * indices[i+1]) # choose point in middle
                β = (rand() - 0.5) * 0.9 + 0.5 # α ∈ [0,1]
                # α = (contour_index-1) / (length(labeled_contours)-1)
                α = contour_index % 2 == 0 ? 0.15 : 0.85
                α = rand([α,β])
                local index3 = round(Int, α * (indices[i]+1)+ (1-α) * (indices[i+1]-1)) # choose point in middle
                if heuristic == 3
                    local index  = index3 # rand([index1, index2]) # choose between random point and point in middle
                elseif heuristic == 1
                    local index  = index1 
                elseif heuristic == 2
                    local index  = index2
                end
                # end of heuristics
                local location = Point3(segments[index]..., 2f0)
                local sc = scatter!(ax, location, markersize=20, align = (:center, :center), color=(:white, 0.1), strokecolor=:white)
                local anno = text!(ax, [("$contour_val", location)], align = (:center, :center), textsize = 20, color = :black)

                delete!(ax, sc)
                delete!(ax, cplot)
                delete!(ax, anno)

                push!(ax.scene, anno)
                push!(ax.scene, sc)
                push!(ax.scene, cplot)
            end
        end
    end # end of adding labels
end


filename = "avg_long_hs_he_10_hp_2_ve_7_vp_2.jld2"
filename = "avg_long_hs_he_8_hp_2_ve_7_vp_2.jld2"

# new versions 
# filename = "avg_long_hs_he_45_hp_1_ve_7_vp_2.jld2"
# filename = "avg_long_hs_he_23_hp_3_ve_7_vp_2.jld2"
# filename = "avg_long_hs_he_15_hp_5_ve_7_vp_2.jld2"

# refanov check, 3
filename = []
push!(filename, "avg_long_hs_he_12_hp_3_ve_10_vp_2_refanov.jld2")
push!(filename, "avg_long_hs_he_12_hp_3_ve_8_vp_3_refanov.jld2")
push!(filename, "avg_long_hs_he_12_hp_3_ve_6_vp_4_refanov.jld2")
# roefanov, 3
push!(filename, "avg_long_hs_he_12_hp_3_ve_10_vp_2_roefanov.jld2")
push!(filename, "avg_long_hs_he_12_hp_3_ve_8_vp_3_roefanov.jld2")
push!(filename, "avg_long_hs_he_12_hp_3_ve_6_vp_4_roefanov.jld2")

# check higherrez, 3
push!(filename, "avg_long_hs_he_23_hp_3_ve_8_vp_3.jld2")
push!(filename, "avg_long_hs_he_7_hp_2_ve_8_vp_3.jld2")
push!(filename, "avg_long_hs_he_3_hp_6_ve_8_vp_3.jld2")

# checking super low rez, 24 dof horizontal and 30 dof vertical
# perhaps do 1, 4, and 5, (2 * 3 * 5), 
# 10
push!(filename, "avg_long_hs_he_6_hp_3_ve_6_vp_4_refanov.jld2")
push!(filename, "avg_long_hs_he_8_hp_2_ve_6_vp_4_refanov.jld2")
push!(filename, "avg_long_hs_he_8_hp_2_ve_6_vp_4_roefanov.jld2")
push!(filename, "avg_long_hs_he_4_hp_6_ve_6_vp_4_roefanov.jld2")

push!(filename, "avg_long_hs_he_3_hp_7_ve_6_vp_4_refanov.jld2") # 14
push!(filename, "avg_long_hs_he_3_hp_7_ve_6_vp_4_roefanov.jld2") # 15

push!(filename, "avg_long_hs_he_12_hp_3_ve_8_vp_4_refanov.jld2") # 16
push!(filename, "avg_long_hs_he_12_hp_3_ve_8_vp_4_roefanov.jld2") # 17

# checking highest rez 
push!(filename, "avg_hr_long_hs_he_13_hp_6_ve_8_vp_4_roefanov.jld2") # 18
push!(filename, "avg_hr_long_hs_he_9_hp_4_ve_16_vp_4_roefanov.jld2") # 19

# check new ic 
push!(filename, "avg_new_ic_hs_he_12_hp_3_ve_6_vp_4_roefanov.jld2") # 20

# small planet
push!(filename, "averages_small_earth_long_hs_he_12_hp_4_ve_12_vp_4_roefanov.jld2") # 21

title_names = []
# 3
push!(title_names, "Rusanov, p = 2")
push!(title_names, "Rusanov, p = 3")
push!(title_names, "Rusanov, p = 4")
# 3
push!(title_names, "Roe, p = 2")
push!(title_names, "Roe, p = 3")
push!(title_names, "Roe, p = 4")
# 4
push!(title_names, "higher rez")
push!(title_names, "lower rez")
push!(title_names, "lower rez")

# 6
push!(title_names, "he 6 hp 3, refanov, low rez")
push!(title_names, "he 8 hp 2, refanov, low rez")

push!(title_names, "he 8 hp 2, roefanov, low rez")
push!(title_names, "he 4 hp 6, roefanov, low rez")

push!(title_names, "he 3 hp 7, refanov, low rez")
push!(title_names, "he 3 hp 7, roefanov, low rez")

push!(title_names, "he 12 hp 3, refanov,  medium rez")
push!(title_names, "he 12 hp 3, roefanov,  medium rez")

# 
push!(title_names, "highest resolution")
push!(title_names, "jordan level vert")

# 
push!(title_names, "he 12 hp 3 ve 6 vp 4 new ic")
push!(title_names, "12 hp 4 ve 12 vp 4 roefanov small planet")

# filename = ["checkme.jld2" for i in 1:6]

# 22
push!(filename, "avg_earth_hs_he_12_hp_4_ve_12_vp_4_roefanov.jld2") 
push!(title_names, "12 hp 4 ve 12 vp 4 roefanov planet")

# comparison starts here
# 23
push!(filename, "averages_small_earth_long_hs_he_23_hp_1_ve_20_vp_1_roefanov.jld2")
push!(title_names, "hs_he_23_hp_1_ve_20_vp_1_roefanov")

# 24
push!(filename, "averages_small_earth_long_hs_he_23_hp_1_ve_8_vp_4_roefanov.jld2")
push!(title_names, "hs_he_23_hp_1_ve_8_vp_4_roefanov")

# 25
push!(filename, "averages_small_earth_long_hs_he_9_hp_4_ve_20_vp_1_roefanov.jld2")
push!(title_names, "hs_he_9_hp_4_ve_20_vp_1_roefanov")

# 26 
push!(filename, "averages_small_earth_long_hs_he_9_hp_4_ve_8_vp_4_roefanov.jld2")
push!(title_names, "hs_he_9_hp_4_ve_8_vp_4_roefanov")

# 27 
push!(filename, "averages_small_earth_long_hs_he_5_hp_8_ve_20_vp_1_roefanov.jld2")
push!(title_names, "hs_he_5_hp_8_ve_20_vp_1_roefanov")

# 28
push!(filename, "averages_small_earth_long_hs_he_5_hp_8_ve_10_vp_3_roefanov.jld2")
push!(title_names, "he_5_hp_8_ve_10_vp_3")

# 29 # Back to Earth
push!(filename, "avg_earth_hs_he_15_hp_1_ve_15_vp_1_roefanov.jld2")
push!(title_names, "he_15_hp_1_ve_15_vp_1")

# 30 
push!(filename, "avg_earth_hs_he_10_hp_2_ve_10_vp_2_roefanov.jld2")
push!(title_names, "he_10_hp_2_ve_10_vp_2")

# 31 
push!(filename, "avg_earth_hs_he_6_hp_4_ve_15_vp_1_roefanov.jld2")
push!(title_names, "avg_earth_hs_he_6_hp_4_ve_15_vp_1_roefanov")

# 32 
push!(filename, "avg_earth_hs_he_6_hp_4_ve_6_vp_4_roefanov.jld2")
push!(title_names, "he_6_hp_4_ve_6_vp_4")

# 33
push!(filename, "avg_earth_hs_he_5_hp_5_ve_5_vp_5_roefanov.jld2")
push!(title_names, "he_5_hp_5_ve_5_vp_5")

# convenience functions 
function grab_state(s_string, jl_file)
    m_string = "moment_" * "$(length(s_string))"
    m_names = jl_file[m_string * "_names"]
    m_ind = findall(x->x==s_string, m_names)[1]
    state = jl_file[m_string][1:end-1, :, :, m_ind] 
    return mean(state, dims = 1)[1,:,:] ./ jl_file["times"]
end

# convenience functions 
function grab_grid(jl_file)
    lat_grid = jl_file["grid"]["latitude"]
    lon_grid = jl_file["grid"]["longitude"]
    rad_grid = jl_file["grid"]["radius"]

    # 6 is pressure
    eq_ind = argmin(abs.(lat_grid)) # pressure at equator
    p2 = mean(jl_file["moment_1"][1:end-1, eq_ind, :, 6], dims = 1)[1,:] ./ jl_file["times"]

    λ = lon_grid 
    ϕ = lat_grid
    r = rad_grid 

    p0 = 1e5 # pressure bottom
    pH = 800 # pressure top
    H = rad_grid[end] - rad_grid[1] # height of domain
    sH = -H / log(pH/p0) # scale height
    # roughly p = p0 * exp(-z / sH), so log(p / p0) * sH = z
    p = p0 * exp.( -(rad_grid .- rad_grid[1]) ./ sH)
    return λ, ϕ, r, p2
end

function eddy_variance(m2_string, jl_file)
    @assert length(m2_string)==2
    string_1 = m2_string[1]
    string_2 = m2_string[end]

    state_1 = grab_state(string(string_1), jl_file)
    state_2 = grab_state(string(string_2), jl_file)
    state_12 = grab_state(m2_string, jl_file)
    state = state_12 .- state_1 .* state_2
    name = string_1 * "'" * string_2 * "'"
    return state, name
end

function plot_helper(s_string, slice_zonal)
    zonal_quantiles = extrema(slice_zonal[:])
    colorrange = zonal_quantiles
    # contour_levels = collect(0:20:120)
    if s_string == "u"
        contour_levels = collect(-8:4:28) # u-velocity 
    elseif s_string == "T"
        contour_levels = collect(190:10:300) # T
    elseif s_string == "T'T'"
        contour_levels = collect(0:8:40) # TT
    elseif s_string == "u'u'"
        # contour_levels = collect(0:40:160) # uu
        contour_levels = [0, 40, 80, 160, 240]
    elseif s_string =="v'T'"
        contour_levels = collect(-21:3:21) # vT
    elseif s_string == "u'v'"
        contour_levels = collect(-30:10:30) # uv
    else
        contour_levels = range(colorrange[1], colorrange[2], length = 10) # default
    end
    return colorrange, contour_levels
end


# Fixed File looking over states: 
# start here
fig = Figure(resolution = (1700+600, 1000+400))
add_label = true
# fi = 22
fi = 33
# fi = 14
println("looking at ", title_names[fi])
state_names = []

i = 1
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[fi], "r+")
s_string = "u"
slice_zonal = grab_state(s_string, jl_file)
colorrange, contour_levels = plot_helper(s_string, slice_zonal)
colorrange = (-40,40) # override
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax1 = fig[jj,ii] = Axis(fig)
mega_viz!(ax1, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label, title_string = state_names[i])

i = 2
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[fi], "r+")
s_string = "T"
slice_zonal = grab_state(s_string, jl_file)
colorrange, contour_levels = plot_helper(s_string, slice_zonal)
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax2 = fig[jj,ii] = Axis(fig)
mega_viz!(ax2, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label, title_string = state_names[i], colormap = :thermometer)
hideydecorations!(ax2, grid = false)

i = 3
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[fi], "r+")

s_string = "TT"
slice_zonal, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels = plot_helper(s_string, slice_zonal)
colorrange = (-8,48) # override
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax3 = fig[jj,ii] = Axis(fig)
mega_viz!(ax3, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label, title_string = state_names[i], colormap = :thermometer)
hideydecorations!(ax3, grid = false)

i = 4
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[fi], "r+")
s_string = "uv"
slice_zonal, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels = plot_helper(s_string, slice_zonal)
colorrange = (-60, 60)
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax4 = fig[jj,ii] = Axis(fig)
mega_viz!(ax4, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label, title_string = state_names[i])

i = 5
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[fi], "r+")
s_string = "vT"
slice_zonal, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels = plot_helper(s_string, slice_zonal)
colorrange = (-24, 24) # override
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax5 = fig[jj,ii] = Axis(fig)
mega_viz!(ax5, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label, title_string = state_names[i])
hideydecorations!(ax5, grid = false)

i = 6
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[fi], "r+")

s_string = "uu"
slice_zonal1, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels = plot_helper(s_string, slice_zonal)

s_string = "vv"
slice_zonal2, s_string = eddy_variance(s_string, jl_file)

slice_zonal = 0.5 .* (slice_zonal1 + slice_zonal2)
colorrange = extrema(slice_zonal)

s_string = "(u'u' + v'v')/2"
push!(state_names, s_string)
λ, ϕ, r, p_coord = grab_grid(jl_file)
ax6 = fig[jj,ii] = Axis(fig)
mega_viz!(ax6, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label, title_string = state_names[i], colormap = :thermometer)
hideydecorations!(ax6, grid = false)

# playin around with different plotting: 
#=
add_label = true
s_string = "T" # grab state
jl_file = jldopen(filename[1], "r+")
slice_zonal = grab_state(s_string, jl_file)
# slice_zonal, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels = plot_helper(s_string, slice_zonal)

# Colorbar(fig[1,2], hm, label = s_string, ticks = contour_levels)

# fixed state looking over files: 
fig = Figure(resolution = (1700, 1000))

i = 1
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[i], "r+")
slice_zonal = grab_state(s_string, jl_file)
colorrange, contour_levels = plot_helper(s_string, slice_zonal)
λ, ϕ, r, p_coord = grab_grid(jl_file)
ax1 = fig[jj,ii] = Axis(fig)
mega_viz!(ax1, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label, title_string = title_names[i])

i = 2
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[i], "r+")
slice_zonal = grab_state(s_string, jl_file)
colorrange, contour_levels = plot_helper(s_string, slice_zonal)
λ, ϕ, r, p_coord = grab_grid(jl_file)
ax2 = fig[jj,ii] = Axis(fig)
mega_viz!(ax2, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label, title_string = title_names[i])
hideydecorations!(ax2, grid = false)

i = 3
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[i], "r+")
slice_zonal = grab_state(s_string, jl_file)
λ, ϕ, r, p_coord = grab_grid(jl_file)
ax3 = fig[jj,ii] = Axis(fig)
mega_viz!(ax3, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label, title_string = title_names[i])
hideydecorations!(ax3, grid = false)

i = 4
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[i], "r+")
slice_zonal = grab_state(s_string, jl_file)
λ, ϕ, r, p_coord = grab_grid(jl_file)
ax4 = fig[jj,ii] = Axis(fig)
mega_viz!(ax4, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label, title_string = title_names[i])

i = 5
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[i], "r+")
slice_zonal = grab_state(s_string, jl_file)
λ, ϕ, r, p_coord = grab_grid(jl_file)
ax5 = fig[jj,ii] = Axis(fig)
mega_viz!(ax5, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label, title_string = title_names[i])
hideydecorations!(ax5, grid = false)

i = 6
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[i], "r+")
slice_zonal = grab_state(s_string, jl_file)
λ, ϕ, r, p_coord = grab_grid(jl_file)
ax6 = fig[jj,ii] = Axis(fig)
mega_viz!(ax6, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label, title_string = title_names[i])
hideydecorations!(ax6, grid = false)


tmp1, tmp2, hm = heatmap(ϕ, p_coord, slice_zonal, colorrange = colorrange, colormap = :balance, interpolate = true)
Colorbar(fig[:,4], hm, label = s_string, ticks = contour_levels)
=#
