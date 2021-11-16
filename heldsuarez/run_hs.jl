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

fluxes = [:roefanov]
fhdofs = [30] # in a face horizontal degrees of freedom in one direction, 4x is the degree resolution
vdofs = [30] # vertical degrees of freedom
hps = [4] # horizontal polynomial order
vps = [8] # vertical polynomial order
# CHANGE DEFAULT SAVING STUFF: LAT LON NOT BEING CALLED
# CHANGE DEFAULT SAVING STUFF: JLD2 SAVE NOT BEING CALLED
for flux in fluxes, fhdof in fhdofs, vdof in vdofs, hp in hps, vp in vps
    jld_name = "quick_check"
    he = ceil(Int, fhdof / (hp+1))
    ve = ceil(Int, vdof / (vp+1))
    descriptor = (he = he, hp = hp, ve = ve, vp = vp, jld_name = jld_name, flux = flux)
    println("currently doing ", descriptor)
    local tic = time()
    held_suarez(model, sim_days = 300, dt = 60, recompute_minutes = 20, he = he, hp = hp, ve = ve, vp = vp, jld_name = jld_name, flux = flux)
    local toc = time()
    println("It took ", (toc - tic ) / (60 * 60), " hours to run ", descriptor)
end