include("smallearth.jl")

# IMEX
fluxes = [:refanov, :roefanov]
fluxes = [:roefanov]
fhdofs = [60] # in a face horizontal degrees of freedom in one direction, 90 is one degreeish
vdofs = [60] # vertical degrees of freedom
hps = [4] # horizontal polynomial order
vps = [4] # vertical polynomial order

for flux in fluxes, fhdof in fhdofs, vdof in vdofs, hp in hps, vp in vps
    he = ceil(Int, fhdof / (hp+1))
    ve = ceil(Int, vdof / (vp+1))
    descriptor = (he = he, hp = hp, ve = ve, vp = vp, jld_name = "small_earth_long_hs", flux = flux)
    println("currently doing ", descriptor)
    held_suarez(model, hcfl = 0.25, sim_days = 1200.0 / small_earth_γ, he = he, hp = hp, ve = ve,
        vp = vp, jld_name = "small_earth_long_hs", flux = flux)
end


# EXPLICIT
#=
fluxes = [:rusanov]
fhdofs = [60] # in a face horizontal degrees of freedom in one direction
vdofs = [30] # vertical degrees of freedom
hps = [4] # horizontal polynomial order
vps = [3] # vertical polynomial order

for flux in fluxes, fhdof in fhdofs, vdof in vdofs, hp in hps, vp in vps
    he = ceil(Int, fhdof / (hp+1))
    ve = ceil(Int, vdof / (vp+1))
    clusterpath = "/central/scratch/jiahe/smallEarth/"
    descriptor = (he = he, hp = hp, ve = ve, vp = vp, jld_name = "small_earth_long_hs", flux = flux)
    println("currently doing ", descriptor)
    explicit_held_suarez(model, sim_days = 1200.0 / small_earth_γ, he = he, hp = hp, ve = ve,
        vp = vp, jld_name = "long_hs", flux = flux, dt = 0.1, clusterpath = clusterpath)
end
=#