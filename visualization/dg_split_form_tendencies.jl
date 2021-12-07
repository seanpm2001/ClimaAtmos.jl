include("utils_dg.jl")

K = 7  # Number of elements

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

u_first  = @. sin(2π/(b-a) * x) *( 1+ 0.9 * sin(8 * 2π/(b-a) * x) )
using Random
Random.seed!(123)
uo = refine * u_first
x = refine * x

Duu  = refine * compute_volume_terms(u_first .* u_first, mesh)
uDu  = 2.0 .* refine * ( u_first .* compute_volume_terms(u_first, mesh) )
ESDu = 2/3 .* Duu + 1/3 .* uDu

## Plotting
fig = Figure(resolution = (1700+600, 1000+400)) 

lims = (extrema(x)..., extrema(uo)...)
axo = Axis(fig[1,1], xlabel = L"x", ylabel = L"y", ylabelsize = 22, 
    xlabelsize= 22, xgridstyle=:dash, ygridstyle=:dash, xtickalign = 1, 
    xticksize=10, ytickalign=1, yticksize=10,  xlabelpadding = -10)

axo.limits = (extrema(x)..., (extrema(uo) .* 1.1)...)
axo.title = "Function"
axo.titlesize = 32
for i in 1:K
    lines!(axo, x[:,i], uo[:,i], linewidth = 5)
end

ax0 = Axis(fig[1,2], xlabel = L"x", ylabel = L"y", ylabelsize = 22, 
    xlabelsize= 22, xgridstyle=:dash, ygridstyle=:dash, xtickalign = 1, 
    xticksize=10, ytickalign=1, yticksize=10,  xlabelpadding = -10, ylims = (-1.1, 1.1))
lims2 = (extrema(x)..., (extrema(Duu) .* 1.1)...)
ax0.limits = lims2
ax0.title = "Conservative "
ax0.titlesize = 32
for i in 1:K
    lines!(ax0, x[:,i], Duu[:,i], linewidth = 5)
end

ax1 = Axis(fig[2,1], xlabel = L"x", ylabel = L"y", ylabelsize = 22, 
    xlabelsize= 22, xgridstyle=:dash, ygridstyle=:dash, xtickalign = 1, 
    xticksize=10, ytickalign=1, yticksize=10,  xlabelpadding = -10, ylims = (-1.1, 1.1))
ax1.limits = lims2
ax1.title = "Advective"
ax1.titlesize = 32
for i in 1:K
    lines!(ax1, x[:,i], uDu[:,i], linewidth = 5)
end

ax3 = Axis(fig[2,2], xlabel = L"x", ylabel = L"y", ylabelsize = 22, 
    xlabelsize= 22, xgridstyle=:dash, ygridstyle=:dash, xtickalign = 1, 
    xticksize=10, ytickalign=1, yticksize=10,  xlabelpadding = -10, ylims = (-1.1, 1.1))
ax3.limits = lims2
ax3.title = "Entropy Stable"
ax3.titlesize = 32
for i in 1:K
    lines!(ax3, x[:,i], ESDu[:,i], linewidth = 5)
end