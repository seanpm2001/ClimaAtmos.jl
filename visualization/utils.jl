using GLMakie
using LaTeXStrings
using JLD2
using Random

"""
function mega_viz!(ax, ϕ, p_coord, slice_zonal, contour_levels, colorrange; 
    add_labels = false, title_string = "", colormap = :balance,
    heuristic = 1, random_seed = 12345)

# Description 
Combines a contour plot with a heatmap and adds numbers to the contours.

"""
function contour_heatmap!(ax, ϕ, p_coord, slice_zonal, contour_levels, colorrange; 
    add_labels = false, colormap = :balance,
    heuristic = 1, random_seed = 12345)

    cplot = contour!(ax, ϕ, p_coord, slice_zonal, levels = contour_levels, color = :black, interpolate = true,  show_axis = false)
    hm = heatmap!(ax, ϕ, p_coord, slice_zonal, colorrange = colorrange, colormap = colormap, interpolate = true)

    ax.limits = (extrema(ϕ)..., extrema(p_coord)...)

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

function grab_state(s_string, jl_file)
    m_string = "moment_" * "$(length(s_string))"
    m_names = jl_file[m_string * "_names"]
    m_ind = findall(x->x==s_string, m_names)[1]
    state = jl_file[m_string][1:end-1, :, :, m_ind] 
    return mean(state, dims = 1)[1,:,:] ./ jl_file["times"]
end

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
        latex_string = L"\langle u \rangle"
    elseif s_string == "T"
        contour_levels = collect(190:10:300) # T
        latex_string = L"$\langle T$ $\rangle$"
    elseif s_string == "T'T'"
        contour_levels = collect(0:8:40) # TT
        latex_string = L"\langle T' T' \rangle"
    elseif s_string == "u'u'"
        # contour_levels = collect(0:40:160) # uu
        contour_levels = [0, 40, 80, 160, 240, 320]
        latex_string = L"\langle u' u' \rangle"
    elseif s_string =="v'T'"
        contour_levels = collect(-21:3:21) # vT
        latex_string = L"\langle u' T' \rangle"
    elseif s_string == "u'v'"
        contour_levels = collect(-30:10:30) # uv
        latex_string = L"\langle u' v' \rangle"
    else
        contour_levels = range(colorrange[1], colorrange[2], length = 10) # default
        latex_string = " "
    end
    return colorrange, contour_levels, latex_string
end