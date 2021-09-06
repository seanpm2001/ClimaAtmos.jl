include("held_suarez.jl")

# held_suarez(model, he = 23, hp = 3, ve = 7, vp = 2, jld_name = "long_hs")

# held_suarez(model, he = 15, hp = 5, ve = 7, vp = 2, jld_name = "long_hs")

# held_suarez(model, he = 15, hp = 2, ve = 10, vp = 2, jld_name = "long_hs", flux = :refanov)

#=
# loop for later
fluxes = [:refanov, :roefanov]
fhdofs = [45, 90] # in a face horizontal degrees of freedom in one direction
vdofs = [20, 30] # vertical degrees of freedom
hp = [1,2,3,4,5,6] # horizontal polynomial order
vp = [1,2,3,4] # vertical polynomial order
=#

fluxes = [:refanov, :roefanov]
fhdofs = [45] # in a face horizontal degrees of freedom in one direction
vdofs = [30] # vertical degrees of freedom
hps = [3] # horizontal polynomial order
vps = [2] # vertical polynomial order

for flux in fluxes, fhdof in fhdofs, vdof in vdofs, hp in hps, vp in vps
    he = ceil(Int, fhdof / (hp+1))
    ve = ceil(Int, vdof / (vp+1))
    descriptor = (he = he, hp = hp, ve = ve, vp = vp, jld_name = "long_hs", flux = flux)
    println("currently doing ", descriptor)
    held_suarez(model, he = he, hp = hp, ve = ve, vp = vp, jld_name = "long_hs", flux = flux)
end