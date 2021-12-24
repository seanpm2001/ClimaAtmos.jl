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

# paper configuration for "truth"
fluxes = [:roefanov]
fhdofs = [60] # in a face horizontal degrees of freedom in one direction, 90 is one degreeish
vdofs = [60] # vertical degrees of freedom
hps = [4] # horizontal polynomial order
vps = [4] # vertical polynomial order

# DONT FORGET TO CHANGE BACK TIMES and CALLBACKS
for flux in fluxes, fhdof in fhdofs, vdof in vdofs, hp in hps, vp in vps
    jld_name = "earth_hs"
    he = ceil(Int, fhdof / (hp+1))
    ve = ceil(Int, vdof / (vp+1))
    descriptor = (he = he, hp = hp, ve = ve, vp = vp, jld_name = jld_name, flux = flux)
    println("currently doing ", descriptor)
    local tic = time()
    held_suarez(model, sim_days = 301, hcfl = 0.2, vcfl = 15, recompute_minutes = 20, he = he, hp = hp, ve = ve, vp = vp, jld_name = jld_name, flux = flux)
    local toc = time()
    println("It took ", (toc - tic ) / (60 * 60), " hours to run ", descriptor)
end

# tests = [(; hp = 1, vp = 1), (; hp = 2, vp = 2),(; hp = 4, vp = 1),(; hp = 4, vp = 1), (; hp = 9, vp = 3)]
# tests = [(; hp = 4, vp = 4), (; hp = 5, vp = 5)]
#=
tests = [(; hp = 3, vp = 3), (; hp = 9, vp = 1)]
for test in tests
    hp = test.hp 
    vp = test.vp
    flux = :roefanov
    fhdof = 30
    vdof = 30

    jld_name = "earth_hs"
    he = ceil(Int, fhdof / (hp+1))
    ve = ceil(Int, vdof / (vp+1))
    descriptor = (he = he, hp = hp, ve = ve, vp = vp, jld_name = jld_name, flux = flux)
    println("currently doing ", descriptor)
    local tic = time()
    held_suarez(model, sim_days = 4*300, hcfl = 0.2, vcfl = 15, recompute_minutes = 20, he = he, hp = hp, ve = ve, vp = vp, jld_name = jld_name, flux = flux)
    local toc = time()
    println("It took ", (toc - tic ) / (60 * 60), " hours to run ", descriptor)

end
=#