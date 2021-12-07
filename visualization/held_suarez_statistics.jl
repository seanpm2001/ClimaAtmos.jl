include("utils.jl")

filename = "avg_earth_hs_he_12_hp_4_ve_12_vp_4_roefanov.jld2"

fig = Figure(resolution = (1700+600, 1000+400))
add_label = true
println("looking at ", filename)
state_names = []

i = 1
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename, "r+")
s_string = "u"
slice_zonal = grab_state(s_string, jl_file)
colorrange, contour_levels, s_string = plot_helper(s_string, slice_zonal)
colorrange = (-40,40) # override
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax1 = fig[jj,ii] = Axis(fig,title = state_names[i], titlesize = 40)
contour_heatmap!(ax1, ϕ, p_coord, slice_zonal, 
                 contour_levels, colorrange, 
                 add_labels = add_label, random_seed = 1)

i = 2
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename, "r+")
s_string = "T"
slice_zonal = grab_state(s_string, jl_file)
colorrange, contour_levels, s_string = plot_helper(s_string, slice_zonal)
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax2 = fig[jj,ii] = Axis(fig,title = state_names[i], titlesize = 40)
contour_heatmap!(ax2, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label, colormap = :thermometer)
hideydecorations!(ax2, grid = false)

i = 3
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename, "r+")

s_string = "TT"
slice_zonal, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels, s_string = plot_helper(s_string, slice_zonal)
colorrange = (-8,48) # override
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax3 = fig[jj,ii] = Axis(fig,title = state_names[i], titlesize = 40)
contour_heatmap!(ax3, ϕ, p_coord, slice_zonal, contour_levels, 
                 colorrange, add_labels = add_label,
                 colormap = :thermometer, random_seed = 10)
hideydecorations!(ax3, grid = false)

i = 4
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename, "r+")
s_string = "uv"
slice_zonal, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels, s_string = plot_helper(s_string, slice_zonal)
colorrange = (-60, 60)
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax4 = fig[jj,ii] = Axis(fig,title = state_names[i], titlesize = 40)
contour_heatmap!(ax4, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label)

i = 5
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename, "r+")
s_string = "vT"
slice_zonal, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels, s_string = plot_helper(s_string, slice_zonal)
colorrange = (-24, 24) # override
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax5 = fig[jj,ii] = Axis(fig,title = state_names[i], titlesize = 40)
contour_heatmap!(ax5, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label)
hideydecorations!(ax5, grid = false)

i = 6
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename, "r+")

s_string = "uu"
slice_zonal1, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels = plot_helper(s_string, slice_zonal)

s_string = "vv"
slice_zonal2, s_string = eddy_variance(s_string, jl_file)

slice_zonal = 0.5 .* (slice_zonal1 + slice_zonal2)
colorrange = extrema(slice_zonal)

s_string = L"\langle (u' u' + v' v')/2 \rangle"
push!(state_names, s_string)
λ, ϕ, r, p_coord = grab_grid(jl_file)
ax6 = fig[jj,ii] = Axis(fig,title = state_names[i], titlesize = 40)
contour_heatmap!(ax6, ϕ, p_coord, slice_zonal,
                 contour_levels, colorrange, add_labels = add_label, 
                 colormap = :thermometer, random_seed = 10)
hideydecorations!(ax6, grid = false)