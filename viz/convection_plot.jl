using GLMakie, JLD2, Statistics
filepath = "/home/sandre/Repositories/ClimaAtmos.jl/"
# jl_file = jldopen("/home/sandre/Repositories/ClimaAtmos.jl/interpolated_convection_3.jld2", "a+")
# jl_file = jldopen("/home/sandre/Repositories/ClimaAtmos.jl/new_interpolated_convection_1.jld2", "a+")
# jl_file = jldopen("/home/sandre/Repositories/ClimaAtmos.jl/new_interpolated_convection_more_gridpoints.jld2", "a+")
jl_file = jldopen(filepath * "new_interpolated_convection_fixed_higher_resolution.jld2")
# jl_file = jldopen(filepath * "new_interpolated_linear_potential_temperature.jld2")
# jl_file = jldopen("/home/sandre/Repositories/ClimaAtmos.jl/new_interpolated_convection_just_tryin.jld2", "a+")
filename = "new_interpolated_linear_potential_temperature_stronger_diff.jld2"
jl_file = jldopen(filepath * filename)
println(" looking at " , filename)
tkeys = keys(jl_file["state"])
# state = jl_file["state"][tkeys[end]]
# temperature drops by about 70 Kelvin in 10km (linearish)
# pressure drops by about 28% in 8.85km (exponentialish), 7km scale height
# potential temperature increases by about 10 degrees over 2km , 15 + 273 at ground
state₀ = jl_file["state"][tkeys[1]]

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

fig = Figure(resolution = (600, 600));
ax = LScene(fig[1, 1], scenekw = (camera = cam3d!, show_axis = true));


tindex = Observable(1)
state = @lift(jl_file["state"][tkeys[$tindex]])

ρ = @lift($state[:, :, :, 1])
w = @lift($state[:, :, :, 4])
T = @lift($state[:, :, :, 7])
p = @lift($state[:, :, :, 6])
θ = @lift($T .* (pₒ ./ $p) .^ (R_d / cp_d))
θ̅ = @lift($θ .- mean($θ, dims = (1, 2)))

cmap = :balance # :Blues_9
cmapa = RGBAf.(to_colormap(cmap), 1);
cmap = vcat(cmapa[5:10], fill(RGBAf(0, 0, 0, 0), 20), cmapa[25:end])

# cmap = vcat(cmapa[1:5], fill(RGBAf(0, 0, 0, 0), 20), cmapa[25:end])
# cmap = vcat(fill(RGBAf(0, 0, 0, 0), 25), cmapa[25:end])
#=
cmap = vcat(fill(RGBAf(0, 0, 0, 0), 15), cmapa[25:end])

cmap = vcat(cmapa[1:5], fill(RGBAf(0, 0, 0, 0), 5))
for i in 1:7
    cmap = vcat(cmap..., cmapa[5*i+1:5*i+5], fill(RGBAf(0, 0, 0, 0), 5))
end
=#
# wmax = maximum(abs.(w))
clims = (-7, 7) # w
plot_state = @lift($θ[:, :, 1:div(length(z), 2)])
# clims = (351.71910286663194, 424.3662951232513)# @lift(extrema($plot_state))
# clims = (385.1117930486359, 388.04880802838943)
clims = @lift((quantile($plot_state[:], 0.01), quantile($plot_state[:], 0.95)))

v1 = volume!(ax, plot_state,
    colorrange = clims, algorithm = :absorption, absorption = 10.0f0,
    colormap = cmap)
tindex[] = length(tkeys)

# scatter(mean(θ[], dims = (1,2))[1,1,:], collect(1:128))

options = (; ylabel = "height", ylabelsize = 32,
    xlabelsize = 32, xgridstyle = :dash, ygridstyle = :dash, xtickalign = 1,
    xticksize = 30, ytickalign = 1, yticksize = 30,
    xticklabelsize = 30, yticklabelsize = 30)

# tindex[] = 2
statfig = Figure(resolution = (1800, 800));
stat_ax1 = Axis(statfig[1, 1]; title = "⟨θ⟩", options...)
stat_ax2 = Axis(statfig[1, 2]; title = "⟨θ'w'⟩", options...)
stat_ax3 = Axis(statfig[1, 3]; title = "⟨w'w'⟩", options...)
stat_ax4 = Axis(statfig[1, 4]; title = "⟨θ'θ'⟩", options...)
w̅ = mean(w[], dims = (1, 2))[1, 1, :]
θ̅ = mean(θ[], dims = (1, 2))[1, 1, :]
θw = mean(θ[] .* w[], dims = (1, 2))[1, 1, :]
ww = mean(w[] .* w[], dims = (1, 2))[1, 1, :]
θθ = mean(θ[] .* θ[], dims = (1, 2))[1, 1, :]
lines!(stat_ax1, θ̅, z)
lines!(stat_ax1, mean(θ₀, dims = (1, 2))[1, 1, :], z, color = :red)
# xlims!(stat_ax1, (305, 308))
# ylims!(stat_ax1, (0, 2000))
scatter!(stat_ax2, θw .- (w̅ .* θ̅), z)
scatter!(stat_ax3, ww .- (w̅ .* w̅), z)
scatter!(stat_ax4, θθ .- (θ̅ .* θ̅), z)


dstatfig = Figure(resolution = (900, 400));
d_ax1 = Axis(dstatfig[1, 1]; title = "∂ᶻ⟨θ⟩", options...)
dθ = (θ̅[2:end] .- θ̅[1:end-1]) ./ (z[2] - z[1])
scatter!(d_ax1, dθ, z[1:end-1])
println("The maximum potential temperature gradient is ", maximum(dθ))
println("The minimum heat flux is ", minimum(θw .- (w̅ .* θ̅)))
println("The maximum heat flux is ", maximum(θw .- (w̅ .* θ̅)))

#=
R = 287
cp = 1004
p0 = 1e5
g = 9.8
Γ = 10 / 2e3
θ0 = 300
θnew = @. θ0 + Γ * z
pnew = @. p0 * (g / (-Γ * cp) * log(θnew / θ0) + 1)^(cp / R)
scatter(pnew, z)
Tnew = @. (pnew / p0)^(R / cp) * θnew
ρnew = @. pnew / (R * Tnew)
ρnew_check = @. 1 / (-g) * p0 * (cp / R) * (g / (-Γ * cp) * log(θnew / θ0) + 1)^(cp / R - 1) * (g / (-Γ * cp) * Γ / θnew)
=#
#=
R = 287
cp = 1004
p0 = 1e5
g = 9.8
Γ = 10 / 2e3
θ0 = 300
p = p0 * (g / (Γ * cp ) * ln(1 - Γ/θ0 * z) + 1)^(cp/R)
=#
#=
framerate = 30
record(fig, "theta_convection.mp4", 1:length(tkeys), framerate = framerate) do t
    if t == 1
        # zoom!(ax.scene, 0.8)
        # zoom!(ax.scene, 0.8)
    end
    tindex[] = t
    println("at time t= ", t)
end
=#