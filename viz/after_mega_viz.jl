fig = Figure(resolution = (1700+600, 1000+400))
add_label = false
fi = 19
# fi = 21
println("looking at ", title_names[fi])
state_names = []

i = 1
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[fi], "r+")
s_string = "w"
slice_zonal = grab_state(s_string, jl_file)
colorrange, contour_levels = plot_helper(s_string, slice_zonal)
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax1 = fig[jj,ii] = Axis(fig,title = state_names[i], titlesize = 40))
mega_viz!(ax1, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label, title_string = state_names[i])
println("The extrema of ", s_string, " are ", colorrange)

i = 2
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[fi], "r+")
s_string = "wT"
slice_zonal, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels = plot_helper(s_string, slice_zonal)
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax2 = fig[jj,ii] = Axis(fig,title = state_names[i], titlesize = 40))
mega_viz!(ax2, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label,  colormap = :balance)
hideydecorations!(ax2, grid = false)
println("The extrema of ", s_string, " are ", colorrange)


i = 3
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[fi], "r+")
s_string = "ww"
slice_zonal, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels = plot_helper(s_string, slice_zonal)
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax2 = fig[jj,ii] = Axis(fig,title = state_names[i], titlesize = 40))
mega_viz!(ax2, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label,  colormap = :balance)
hideydecorations!(ax2, grid = false)
println("The extrema of ", s_string, " are ", colorrange)


i = 4
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[fi], "r+")
s_string = "uw"
slice_zonal, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels = plot_helper(s_string, slice_zonal)
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax2 = fig[jj,ii] = Axis(fig,title = state_names[i], titlesize = 40))
mega_viz!(ax2, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label,  colormap = :balance)
hideydecorations!(ax2, grid = false)
println("The extrema of ", s_string, " are ", colorrange)

i = 5
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[fi], "r+")
s_string = "vw"
slice_zonal, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels = plot_helper(s_string, slice_zonal)
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax2 = fig[jj,ii] = Axis(fig,title = state_names[i], titlesize = 40))

mega_viz!(ax2, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label,  colormap = :balance)
hideydecorations!(ax2, grid = false)
println("The extrema of ", s_string, " are ", colorrange)

i = 6
ii = (i-1)%3 + 1 # +1 on why 1 based indexing is wrong
jj = (i-1)÷3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename[fi], "r+")
s_string = "ρρ"
slice_zonal, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels = plot_helper(s_string, slice_zonal)
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax2 = fig[jj,ii] = Axis(fig,title = state_names[i], titlesize = 40))
mega_viz!(ax2, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels = add_label,  colormap = :balance)
hideydecorations!(ax2, grid = false)
println("The extrema of ", s_string, " are ", colorrange)


##
fi = 21
jl_file = jldopen(filename[fi], "r+")


p_coord = grab_state("p", jl_file)
Δp = p_coord[:,1:end-1] - p_coord[:,2:end]
λ, ϕ, r, p_coord = grab_grid(jl_file)

v = grab_state("v", jl_file)
v̅ = (v[:,1:end-1] + v[:,2:end]) .* 0.5
ṽ = r[1] / 9.81 * 2 * π * (reshape(cosd.(ϕ), (length(ϕ),1)) .* cumsum(v̅ .* Δp, dims = 2))
colorrange = extrema(ṽ)
fig, ax, cplot = contour(ϕ, p_coord[1:end-1], ṽ, color = :black, levels =21, show_axis = false)
hm = heatmap!(ax, ϕ, p_coord[1:end-1], ṽ, colorrange = colorrange, colormap = :balance, interpolate = true)
Colorbar(fig[1,2], hm, ticks = -7.2e10/20:2.4e10/20:7.2e10/20)

ax.limits = (extrema(ϕ)..., extrema(p_coord[1:end-1])...)
ax.title = "meriodonal streamfunction"
ax.titlesize = 40
ax.xlabel = "Latitude [ᵒ]"
ax.ylabel = "Stretched Height"
ax.xlabelsize = 25
ax.ylabelsize = 25 
ax.xticks = ([-80, -60, -30, 0, 30, 60, 80], ["80S", "60S", "30S", "0", "30N", "60N", "80N"])
pressure_levels = [1000, 850, 700, 550, 400, 250, 100, 10]
ax.yticks = (pressure_levels .* 1e2, string.(pressure_levels))
ax.yreversed = true


