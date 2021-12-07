include("utils_dg.jl")

u⁰(x, a, b) = exp(-2 * (b-a) / 3 * (x - (b-a)/2)^2);
u⁰(x,a,b) = sin(2π/(b-a) * x) + 0.5
cfl = 0.3
K = 5  # Number of elements

n = 16
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
ν = 1e-8

u_first  = @.sin( 2 * 2π/(b-a) * x) 

uo = refine * u_first
x = refine * x

u0 = refine * p0 * u_first
u1 = refine * p1 * u_first
u2 = refine * p2 * u_first

## Plotting
fig = Figure(resolution = resolution = (1700+600, 1000+400)) 

lims = (extrema(x)..., (-1.4,1.4)...)
axo = Axis(fig[1,1], xlabel = L"x", ylabel = L"y", ylabelsize = 22, 
    xlabelsize= 22, xgridstyle=:dash, ygridstyle=:dash, xtickalign = 1, 
    xticksize=10, ytickalign=1, yticksize=10,  xlabelpadding = -10)
axo.limits = lims
axo.title = "Original"
axo.titlesize = 32
for i in 1:K
    lines!(axo, x[:,i], uo[:,i], linewidth = 5)
end

ax0 = Axis(fig[1,2], xlabel = L"x", ylabel = L"y", ylabelsize = 22, 
    xlabelsize= 22, xgridstyle=:dash, ygridstyle=:dash, xtickalign = 1, 
    xticksize=10, ytickalign=1, yticksize=10,  xlabelpadding = -10, ylims = (-1.1, 1.1))
ax0.limits = lims
ax0.title = "Finite Volume"
ax0.titlesize = 32
for i in 1:K
    lines!(ax0, x[:,i], u0[:,i], linewidth = 5)
end

ax1 = Axis(fig[2,1], xlabel = L"x", ylabel = L"y", ylabelsize = 22, 
    xlabelsize= 22, xgridstyle=:dash, ygridstyle=:dash, xtickalign = 1, 
    xticksize=10, ytickalign=1, yticksize=10,  xlabelpadding = -10, ylims = (-1.1, 1.1))
ax1.limits = lims
ax1.title = "Polynomial Order 1"
ax1.titlesize = 32
for i in 1:K
    lines!(ax1, x[:,i], u1[:,i], linewidth = 5)
end

ax3 = Axis(fig[2,2], xlabel = L"x", ylabel = L"y", ylabelsize = 22, 
    xlabelsize= 22, xgridstyle=:dash, ygridstyle=:dash, xtickalign = 1, 
    xticksize=10, ytickalign=1, yticksize=10,  xlabelpadding = -10, ylims = (-1.1, 1.1))
ax3.limits = lims
ax3.title = "Polynomial Order 2"
ax3.titlesize = 32
for i in 1:K
    lines!(ax3, x[:,i], u2[:,i], linewidth = 5)
end

display(fig)
