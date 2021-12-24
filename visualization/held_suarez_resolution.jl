include("utils.jl")

filename = []
title_names = []
# 1
push!(filename, "avg_earth_hs_he_12_hp_4_ve_12_vp_4_roefanov.jld2")
push!(title_names, "truth")

# 2
push!(filename, "avg_earth_hs_he_15_hp_1_ve_15_vp_1_roefanov.jld2")
push!(title_names, "pₕ = 1 , pᵥ = 1")

# 3
push!(filename, "avg_earth_hs_he_10_hp_2_ve_10_vp_2_roefanov.jld2")
push!(title_names, "pₕ = 2 , pᵥ = 2")

# 4
push!(filename, "avg_earth_hs_he_6_hp_4_ve_15_vp_1_roefanov.jld2")
push!(title_names, "pₕ = 4 , pᵥ = 1")

# 5 
push!(filename, "avg_earth_hs_he_6_hp_4_ve_6_vp_4_roefanov.jld2")
push!(title_names, "pₕ = 4 , pᵥ = 4")

# 6 
push!(filename, "avg_earth_hs_he_5_hp_5_ve_5_vp_5_roefanov.jld2")
push!(title_names, "pₕ = 5, pᵥ = 5")

# 7 
push!(filename, "avg_earth_hs_he_5_hp_5_ve_5_vp_5_roefanov.jld2")
push!(title_names, "pₕ = 5, pᵥ = 5")


colorrange = (0, 360) # modify to extrema of "truth"
contour_levels = [0, 40, 80, 160, 240, 320]
fig = Figure(resolution = (1700+600, 1000+400))
add_label = true

for i in 1:6
    ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
    jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
    jl_file = jldopen(filename[i], "r+")

    s_string = "uu"
    slice_zonal1, s_string = eddy_variance(s_string, jl_file)

    s_string = "vv"
    slice_zonal2, s_string = eddy_variance(s_string, jl_file)

    slice_zonal = 0.5 .* (slice_zonal1 + slice_zonal2)

    s_string = title_names[i]
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