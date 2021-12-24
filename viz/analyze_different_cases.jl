
using Random
using JLD2
using Statistics

small_earth = false

filename = []
title_names = []
# Truth 
if small_earth
    push!(filename, "averages_small_earth_long_hs_he_12_hp_4_ve_12_vp_4_roefanov.jld2") # 21
else
    push!(filename, "avg_earth_hs_he_12_hp_4_ve_12_vp_4_roefanov.jld2")
end
push!(title_names, "truth")

# 2
if small_earth
    push!(filename, "averages_small_earth_long_hs_he_23_hp_1_ve_20_vp_1_roefanov.jld2")
    push!(title_names, "pₕ = 1 , pᵥ = 1")
else
    push!(filename, "avg_earth_hs_he_15_hp_1_ve_15_vp_1_roefanov.jld2")
    push!(title_names, "pₕ = 1 , pᵥ = 1")
end

# 3
if small_earth
    push!(filename, "averages_small_earth_long_hs_he_23_hp_1_ve_8_vp_4_roefanov.jld2")
    push!(title_names, "pₕ = 1 , pᵥ = 4")
else
    push!(filename, "avg_earth_hs_he_10_hp_2_ve_10_vp_2_roefanov.jld2")
    push!(title_names, "pₕ = 2 , pᵥ = 2")
end

# 4
if small_earth
    push!(filename, "averages_small_earth_long_hs_he_9_hp_4_ve_20_vp_1_roefanov.jld2")
    push!(title_names, "pₕ = 4 , pᵥ = 1")
else
    push!(filename, "avg_earth_hs_he_6_hp_4_ve_15_vp_1_roefanov.jld2")
    push!(title_names, "pₕ = 4 , pᵥ = 1")
end


# 5 
if small_earth
    push!(filename, "averages_small_earth_long_hs_he_9_hp_4_ve_8_vp_4_roefanov.jld2")
    push!(title_names, "pₕ = 4 , pᵥ = 4")
else
    push!(filename, "avg_earth_hs_he_6_hp_4_ve_6_vp_4_roefanov.jld2")
    push!(title_names, "pₕ = 4 , pᵥ = 4")
end

# 6 
if small_earth
    push!(filename, "averages_small_earth_long_hs_he_5_hp_8_ve_20_vp_1_roefanov.jld2")
    push!(title_names, "pₕ = 8, pᵥ = 1")
else
    push!(filename, "avg_earth_hs_he_5_hp_5_ve_5_vp_5_roefanov.jld2")
    push!(title_names, "pₕ = 5, pᵥ = 5")
end

# 7 
if small_earth
    push!(filename, "averages_small_earth_long_hs_he_5_hp_8_ve_20_vp_1_roefanov.jld2")
    push!(title_names, "pₕ = 8, pᵥ = 1")
else
    push!(filename, "avg_earth_hs_he_5_hp_5_ve_5_vp_5_roefanov.jld2")
    push!(title_names, "pₕ = 5, pᵥ = 5")
end

# jl_file = jldopen(filename[1], "r+")
# truth_u = grab_state("u", jl_file)

function grab_state(s_string, jl_file)
    m_string = "moment_" * "$(length(s_string))"
    m_names = jl_file[m_string * "_names"]
    m_ind = findall(x->x==s_string, m_names)[1]
    state = jl_file[m_string][1:end-1, :, :, m_ind] 
    return mean(state, dims = 1)[1,:,:] ./ jl_file["times"]
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

# comparison_field = "p"
# grab = grab_state 
field_function_list = (
                        (;comparison_field = "u", func = grab_state),
                        (;comparison_field = "T",  func = grab_state),
                        (;comparison_field = "TT", func = eddy_variance),
                        (;comparison_field = "uv", func = eddy_variance),
                        (;comparison_field = "vT", func = eddy_variance),
                        (;comparison_field = "uu", func = eddy_variance),
                        (;comparison_field = "v", func = grab_state), 
                        (;comparison_field = "vv", func = eddy_variance),
                        (;comparison_field = "wT", func = eddy_variance),
                        (;comparison_field = "ρ", func = grab_state),
                        (;comparison_field = "ρe", func = grab_state),
                        (;comparison_field = "w",  func = grab_state),
                       )
field_function_list = field_function_list[1:6]
    
error_matrix = zeros(5, length(field_function_list))
for (i, field_function) in enumerate(field_function_list)
comparison_field = field_function.comparison_field
grab = field_function.func

jl_file = jldopen(filename[1], "r+")
u_truth = grab(comparison_field , jl_file)
u_truth isa Tuple ? u_truth = u_truth[1] : nothing

jl_file = jldopen(filename[2], "r+")
u_2 = grab(comparison_field , jl_file)
u_2 isa Tuple ? u_2 = u_2[1] : nothing

jl_file = jldopen(filename[3], "r+")
u_3 = grab(comparison_field , jl_file)
u_3 isa Tuple ? u_3 = u_3[1] : nothing

jl_file = jldopen(filename[4], "r+")
u_4 = grab(comparison_field , jl_file)
u_4 isa Tuple ? u_4 = u_4[1] : nothing

jl_file = jldopen(filename[5], "r+")
u_5 = grab(comparison_field , jl_file)
u_5 isa Tuple ? u_5 = u_5[1] : nothing

jl_file = jldopen(filename[6], "r+")
u_6 = grab(comparison_field , jl_file)
u_6 isa Tuple ? u_6 = u_6[1] : nothing

jl_file = jldopen(filename[7], "r+")
u_7 = grab(comparison_field , jl_file)
u_7 isa Tuple ? u_7 = u_7[1] : nothing

error_list = Float64[]
er1 = abs.(u_truth - u_2) #  /  maximum(abs.(u_truth))
er2 = abs.(u_truth - u_3) #  / maximum(abs.(u_truth))
er3 = abs.(u_truth - u_4) #  / maximum(abs.(u_truth))
er4 = abs.(u_truth - u_5) #  / maximum(abs.(u_truth))
er5 = abs.(u_truth - u_6) #  / maximum(abs.(u_truth))
qn = 1.0 # quantile number
push!(error_list, quantile(er1[:], qn))
push!(error_list, quantile(er2[:], qn))
push!(error_list, quantile(er3[:], qn))
push!(error_list, quantile(er4[:], qn))
push!(error_list, quantile(er5[:], qn))
# push!(error_list, maximum(abs.(u_truth - u_7)))

error_matrix[:,i] .= error_list

minerror = title_names[argmin(error_list) + 1] # since u_truth is 1
println("The minimum error for " * comparison_field, " is ", minerror)
println("All the errors are ", error_list)
println("-------------------------")

end

#=
maximum(abs.(u_2 - u_3))
maximum(abs.(u_2 - u_4))
maximum(abs.(u_2 - u_5))
maximum(abs.(u_2 - u_6))

maximum(abs.(u_3 - u_4))
maximum(abs.(u_3 - u_5))
maximum(abs.(u_3 - u_6))

maximum(abs.(u_4 - u_5))
maximum(abs.(u_4 - u_6))

maximum(abs.(u_5 - u_6))
=#

#=
vizit = true
if vizit
    using GLMakie

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


    colorrange = (0, 320)
    contour_levels = [0, 40, 80, 160, 240, 300]
    fig = Figure(resolution = (1700+600, 1000+400))
    add_label = true

    for i in 1:6
        ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
        jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
        jl_file = jldopen(filename[i], "r+")

        s_string = "uu"
        slice_zonal1, s_string = eddy_variance(s_string, jl_file)
        # colorrange, contour_levels = plot_helper(s_string, slice_zonal)

        s_string = "vv"
        slice_zonal2, s_string = eddy_variance(s_string, jl_file)

        slice_zonal = 0.5 .* (slice_zonal1 + slice_zonal2)
        # colorrange = extrema(slice_zonal)

        # s_string = "(u'u' + v'v')/2"
        s_string = title_names[i]
        λ, ϕ, r, p_coord = grab_grid(jl_file)
        ax = fig[jj,ii] = Axis(fig)
        mega_viz!(ax, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label, title_string = s_string, colormap = :thermometer)
        if ii !== 1
            hideydecorations!(ax, grid = false)
        end
    end
end
=#