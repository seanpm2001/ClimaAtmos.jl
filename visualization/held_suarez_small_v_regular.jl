include("utils.jl")

filename = []
title_names = []
# 1
push!(filename, "avg_earth_hs_he_12_hp_4_ve_12_vp_4_roefanov.jld2")
push!(title_names, "Earth")

# 2
push!(filename, "averages_small_earth_long_hs_he_12_hp_4_ve_12_vp_4_roefanov.jld2")
push!(title_names, "Small Earth")

s_strings = ["wT", "ww", "TT"]
sL_strings = [L"\langle w'T' \rangle", L"\langle w'w' \rangle", L"\langle T'T' \rangle"]


colorrange = (0, 360) # modify to extrema of "truth"
contour_levels = [0, 40, 80, 160, 240, 320]
fig = Figure(resolution = (1700+600, 1000+400))
add_label = false

for i in 1:6
    ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
    jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong

    jl_file = jldopen(filename[jj], "r+")

    s_string = s_strings[ii]
    if length(s_string) == 2
        slice_zonal, s_string = eddy_variance(s_string, jl_file)
        s_string = sL_strings[ii]
    else
        slice_zonal = grab_state(s_string, jl_file)
        s_string = sL_strings[ii]
    end
    colorrange, contour_levels, s_string = plot_helper(s_string, slice_zonal)
    println(colorrange)
    
    s_string = title_names[jj]
    λ, ϕ, r, p_coord = grab_grid(jl_file)
    ax = fig[jj,ii] = Axis(fig, title = s_string, titlesize = 40)
    if ii !==1
    contour_heatmap!(ax, ϕ, p_coord, slice_zonal, contour_levels, colorrange, 
                     add_labels = add_label, colormap = :thermometer)
    else
        contour_heatmap!(ax, ϕ, p_coord, slice_zonal, contour_levels, colorrange, 
                     add_labels = add_label, colormap = :thermometer,
                     random_seed = 10)
    end
    if ii !== 1
        hideydecorations!(ax, grid = false)
    end
end