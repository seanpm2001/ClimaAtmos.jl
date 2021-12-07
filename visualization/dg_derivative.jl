include("utils_dg.jl")

K = 6  # Number of elements

n = 4
nn = 1
r = jacobiGL(0, 0, n)
V = vandermonde(r, 0, 0, n)
fr = collect(range(-1,1,length=100))
fV = vandermonde(fr, 0, 0, n)
refine = fV * inv(V)
# cell average
avg = Diagonal(zeros(n+1))
avg[1:nn] .= 1
p0 = V * avg * inv(V)
# polyone
nn = 2
tmp = zeros(n+1)
tmp[1:nn] .= 1
avg = Diagonal(tmp)
p1 = V * avg * inv(V)
# poly two
nn = 2+1
tmp = zeros(n+1)
tmp[1:nn] .= 1
avg = Diagonal(tmp)
p2 = V * avg * inv(V)

# poly 6
Ω = (; a=0, b=2π, periodic = true)
a = Ω.a
b = Ω.b
mesh = Mesh(K, n, Ω.a, Ω.b, periodic = Ω.periodic)
x = mesh.x
weights = sum(mesh.M, dims = 1)
Δx = x[2] - x[1]

# inexact integration
DM = Diagonal(sum(mesh.M, dims = 1)[:])
mesh.M .= DM
mesh.Mi .= inv(DM)
mesh.lift[:,1] .= mesh.Mi[1,:]
mesh.lift[:,end] .= mesh.Mi[end,:]

b = Ω.b
a = Ω.a

u_first  = mesh.D *  @.sin( 2 * 2π/(b-a) * x) 
using Random
Random.seed!(123)
u_first = compute_volume_terms(u_first, mesh) + 0.1 * compute_surface_terms(randn(size(u_first)), mesh)
uo = refine * u_first
x = refine * x

DVu = refine * compute_volume_terms(u_first, mesh)
DAu = refine * compute_surface_terms(u_first, mesh)
Du = DVu + DAu

fig = Figure(resolution = (1700+600, 1000+400)) 

lims = (extrema(x)..., extrema(uo)...)
axo = Axis(fig[1,1], xlabel = L"x", ylabel = L"y", ylabelsize = 22, 
    xlabelsize= 22, xgridstyle=:dash, ygridstyle=:dash, xtickalign = 1, 
    xticksize=10, ytickalign=1, yticksize=10,  xlabelpadding = -10)

axo.limits = (extrema(x)..., (extrema(uo) .* 1.1)...)
axo.title = "Original"
axo.titlesize = 32
for i in 1:K
    lines!(axo, x[:,i], uo[:,i], linewidth = 5)
end

ax0 = Axis(fig[1,2], xlabel = L"x", ylabel = L"y", ylabelsize = 22, 
    xlabelsize= 22, xgridstyle=:dash, ygridstyle=:dash, xtickalign = 1, 
    xticksize=10, ytickalign=1, yticksize=10,  xlabelpadding = -10, ylims = (-1.1, 1.1))
lims2 = (extrema(x)..., (extrema(DVu) .* 1.1)...)
ax0.limits = lims2
ax0.title = "Volume Contribution"
ax0.titlesize = 32
for i in 1:K
    lines!(ax0, x[:,i], DVu[:,i], linewidth = 5)
end

ax1 = Axis(fig[2,1], xlabel = L"x", ylabel = L"y", ylabelsize = 22, 
    xlabelsize= 22, xgridstyle=:dash, ygridstyle=:dash, xtickalign = 1, 
    xticksize=10, ytickalign=1, yticksize=10,  xlabelpadding = -10, ylims = (-1.1, 1.1))
ax1.limits = (extrema(x)..., (extrema(DAu) .* 1.1)...)
ax1.title = "Surface Contribution"
ax1.titlesize = 32
for i in 1:K
    lines!(ax1, x[:,i], DAu[:,i], linewidth = 5)
end

ax3 = Axis(fig[2,2], xlabel = L"x", ylabel = L"y", ylabelsize = 22, 
    xlabelsize= 22, xgridstyle=:dash, ygridstyle=:dash, xtickalign = 1, 
    xticksize=10, ytickalign=1, yticksize=10,  xlabelpadding = -10, ylims = (-1.1, 1.1))
ax3.limits = (extrema(x)..., (extrema(Du) .* 1.1)...)
ax3.title = "Derivative"
ax3.titlesize = 32
for i in 1:K
    lines!(ax3, x[:,i], Du[:,i], linewidth = 5)
end