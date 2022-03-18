using GLMakie, JLD2, Statistics
filepath = "/home/sandre/Repositories/ClimaAtmos.jl/"
filename = "interpolated_high_resolution_layer.jld2"
filename = "hr_interpolated_paper_rez_convection.jld2"
filename = "hr_interpolated_paper_rez_convection_p4.jld2"
jl_file = jldopen(filepath * filename)
println(" looking at ", filename)
tkeys = keys(jl_file["state"])

state₀ = jl_file["state"][tkeys[1]]

height_index = 94 # index 96 is where the buoyancy flux becomes zero again # first positive entry after argmin
if haskey(jl_file, "parameters")
    parameters = jl_file["parameters"]
    R_d = parameters.R_d
    cp_d = parameters.cp_d
    pₒ = parameters.pₒ
    Tₛ = parameters.Tₛ
else
    R_d = 287.0
    cp_d = 1004
    pₒ = 1e5
    Tₛ = 300.0
end

ρ₀ = state₀[:, :, :, 1]
w₀ = state₀[:, :, :, 4]
T₀ = state₀[:, :, :, 7]
p₀ = state₀[:, :, :, 6]
θ₀ = T₀ .* (pₒ ./ p₀) .^ (R_d / cp_d)

x = jl_file["grid"]["x"]
y = jl_file["grid"]["y"]
z = jl_file["grid"]["z"]

tindex = Observable(11)
state = @lift(jl_file["state"][tkeys[$tindex]])

# states
ρ = @lift($state[:, :, :, 1])
w = @lift($state[:, :, :, 4])
T = @lift($state[:, :, :, 7])
p = @lift($state[:, :, :, 6])
θ = @lift($T .* (pₒ ./ $p) .^ (R_d / cp_d))

# statistics
w̅ = @lift(mean($w, dims = (1, 2))[1, 1, :])
θ̅ = @lift(mean($θ, dims = (1, 2))[1, 1, :])
T̅ = @lift(mean($T, dims = (1, 2))[1, 1, :])

θw = @lift(mean($θ .* $w, dims = (1, 2))[1, 1, :])
ww = @lift(mean($w .* $w, dims = (1, 2))[1, 1, :])
θθ = @lift(mean($θ .* $θ, dims = (1, 2))[1, 1, :])
wT = @lift(mean($w .* $T, dims = (1, 2))[1, 1, :])

θpwp = @lift($θw - $w̅ .* $θ̅)
θpθp = @lift($θθ - $θ̅ .* $θ̅)

#=
options = (; titlesize = 30, ylabelsize = 32,
    xlabelsize = 32, xgridstyle = :dash, ygridstyle = :dash, xtickalign = 1,
    xticksize = 10, ytickalign = 1, yticksize = 10,
    xticklabelsize = 20, yticklabelsize = 20, ylabel = "z [km]")
    =#
options = (; titlesize = 30, ylabelsize = 32,
    xlabelsize = 32, xgridstyle = :dash, ygridstyle = :dash, xtickalign = 1,
    xticksize = 10, ytickalign = 1, yticksize = 10,
    xticklabelsize = 10, yticklabelsize = 10, ylabel = "z [km]")

fig = Figure(resolution = (2000, 1300));
fig = Figure(resolution = (1000, 650));
volume_width = 3

ax = LScene(fig[1:8, 2:volume_width+1], scenekw = (camera = cam3d!, show_axis = true));
stat_ax1 = Axis(fig[2:3, volume_width+2]; title = "⟨θ⟩ [K]", options...)
stat_ax2 = Axis(fig[4:5, volume_width+2]; title = "⟨θ'w'⟩ [K m s⁻¹]", options...)
stat_ax3 = Axis(fig[6:7, volume_width+2]; title = "⟨θ'θ'⟩ [K²]", options...)

cmap = :linear_kryw_0_100_c71_n256

cmapa = reverse(RGBAf.(to_colormap(cmap), 1));
updraftquantile = 0.83
cmap = vcat(fill(RGBAf(0, 0, 0, 0), floor(Int, 40 * (1 / (1 - updraftquantile) - 1))), cmapa[1:35])

plot_state = @lift($θ[:, :, 1:height_index])

clims = @lift((quantile($plot_state[:], 0.00), mean($plot_state[:, :, 1])))

v1 = volume!(ax, 0 .. 1, 0 .. 1, 0 .. 1, plot_state,
    colorrange = clims, algorithm = :absorption, absorption = 10.0f0,
    colormap = cmap)

axis = ax.scene[OldAxis]
axis[:names, :axisnames] = ("x [km]", "y [km]", "z [km]")
tstyle = axis[:names] #  get the nested attributes and work directly with them

tstyle[:textsize] = 05
tstyle[:textcolor] = (:black, :black, :black)
tstyle[:font] = "helvetica"
tstyle[:gap] = 10
axis[:ticks][:textcolor] = :black
axis[:ticks][:textsize] = 05
cbar1 = Colorbar(fig[2:7, 1], v1, label = "θ [K]", width = 25, ticklabelsize = 20,
    labelsize = 20, ticksize = 25, tickalign = 1, height = Relative(3 / 4)
)

axis[:ticks][:ranges] = ([0.0, 0.5, 1.0], [0.0, 0.5, 1.0], [0.0, 0.5, 1.0])
axis[:ticks][:labels] = (["-1.5", "0", "1.5"], ["-1.5", "0", "1.5"], ["0", "0.7", "1.4"])

line_options = (; linewidth = 3, color = cmapa[end-4])
lines!(stat_ax1, θ̅, z ./ 1e3; line_options...)
xlims!(stat_ax1, (303.1, 307.9))
ylims!(stat_ax1, (0.0, 2.5))
lines!(stat_ax2, θpwp, z ./ 1e3; line_options...) # memo to self, imposed flux is (1-exp(-L/ℓ)) this 77% factor )
xlims!(stat_ax2, (-0.03, 0.1))
lines!(stat_ax3, θpθp, z ./ 1e3; line_options...)
rotate_cam!(fig.scene.children[1], (π / 16, 0, 0))
display(fig)

#=
framerate = 30
iterations = 1:30*6
record(fig, "rotating_plot.mp4", iterations, framerate = framerate) do θ
    rotate_cam!(fig.scene.children[1], (0, 2π / length(iterations), 0))
    println("rotation θ = ", θ)
end
=#
#=
# much faster to do this 
for θ in iterations
    rotate_cam!(fig.scene.children[1], (0, 2π / length(iterations), 0))
    println("rotation θ = ", θ)
end
=#

