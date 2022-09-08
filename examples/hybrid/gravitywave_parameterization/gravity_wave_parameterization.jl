function gravity_wave_cache(::Type{FT};) where {FT}
    # parameters
    gw_source_height = FT(15000) # TODO: compute source levels as f(lat)
    Bm = FT(1.2) # TODO: compute source amplitude as f(lat)
    F_S0 = FT(4e-3)
    dc = FT(0.6)
    cmax = FT(99.6)
    # create the spectrum of phase speed and horizontal wavenumber
    nc = Int(floor(2 * cmax / dc + 1))
    c0 = [(n - 1) * dc - cmax for n in 1:nc]
    nk = Int(1)
    kwv = [2π / (30 * (10^n)) * 1e3 for n in 1:nk]
    k2 = kwv .^ 2
    return (;
        gw_source_height = gw_source_height,
        gw_Bm = Bm,
        gw_F_S0 = F_S0,
        gw_c0 = c0,
        gw_nk = nk,
        gw_k = kwv,
        gw_k2 = k2,
    )
end

function gravity_wave_tendency!(Yₜ, Y, p, t)
    #unpack
    (; gw_source_height, gw_Bm, gw_F_S0, gw_c0, gw_nk, gw_k, gw_k2) = p
    (; ᶜts, ᶜT, params) = p
    ᶜρ = Y.c.ρ
    # parameters
    thermo_params = CAP.thermodynamics_params(params)
    grav = FT(CAP.grav(params))
    cp_m = @. TD.cp_m(thermo_params, ᶜts)

    # compute buoyancy frequency
    @. ᶜT = TD.air_temperature(thermo_params, ᶜts)

    grad_scalar = similar(ᶜT)  
    parent(grad_scalar) .= parent(Geometry.WVector.(ᶜgradᵥ.(ᶠinterp.(ᶜT)))) 
 
    ᶜbf = @. (grav / ᶜT) * (grad_scalar + grav / cp_m)
    ᶜbf = @. ifelse(ᶜbf < FT(2.5e-5), FT(sqrt(2.5e-5)), sqrt(ᶜbf)) # to avoid small numbers
    # alternative
    # ᶠbf = [i > 2.5e-5 ? sqrt(i) : sqrt(2.5e-5)  for i in ᶠbf]
    # TODO: create an extra layer at model top so that the gravity wave forcing
    # .     occurring between the topmost model level and the upper boundary
    # .     may be calculated
    # 

    # source level: get the index of the level that is closest to the source height (GFDL uses the fist level below instead)
    source_level = argmin(
        abs.(
            unique(parent(Fields.coordinate_field(Y.c).z) .- gw_source_height)
        ),
    )
    ᶠz = Fields.coordinate_field(Y.f).z

    # TODO: prepare input variables for gravity_wave_forcing()
    #       1. grab column wise data
    #       2. convert coviant uv to physical uv
    # need to make everything an array befere
    u_phy = Geometry.UVVector.(Y.c.uₕ).components.data.:1
    v_phy = Geometry.UVVector.(Y.c.uₕ).components.data.:2

    uforcing = ones(axes(u_phy))
    vforcing = similar(v_phy)
    # bycolume
    Fields.bycolumn(axes(ᶜρ)) do colidx
        parent(uforcing[colidx]) .= gravity_wave_forcing(
            parent(u_phy[colidx])[:, 1],
            source_level,
            gw_Bm,
            gw_c0,
            gw_nk,
            gw_k,
            gw_k2,
            parent(ᶜbf[colidx])[:, 1],
            parent(ᶜρ[colidx])[:, 1],
            parent(ᶠz[colidx])[:, 1],
        )
        parent(vforcing[colidx]) .= gravity_wave_forcing(
            parent(v_phy[colidx])[:, 1],
            source_level,
            gw_Bm,
            gw_c0,
            gw_nk,
            gw_k,
            gw_k2,
            parent(ᶜbf[colidx])[:, 1],
            parent(ᶜρ[colidx])[:, 1],
            parent(ᶠz[colidx])[:, 1],
        )
        
    end

    @. Yₜ.c.uₕ +=
        Geometry.Covariant12Vector.(Geometry.UVVector.(uforcing, vforcing))

end



function gravity_wave_forcing(
    ᶜu,
    source_level,
    gw_Bm,
    gw_c0,
    nk,
    kwv,
    k2,
    ᶜbf,
    ᶜρ,
    ᶠz,
)

    flag = Int(1)
    Bw = FT(0.4)
    Bn = FT(0.0)
    cw = FT(40.0)

    iz0 = source_level

    ampl = gw_Bm

    c0 = gw_c0
    nc = length(c0)

    # !--------------------------------------------------------------------
    # !    define wave momentum flux (B0) at source level for each phase 
    # !    speed n, and the sum over all phase speeds (Bsum), which is needed 
    # !    to calculate the intermittency. 
    # !-------------------------------------------------------------------
    # !---------------------------------------------------------------------
    # !    define wave momentum flux at source level for phase speed n. Add
    # !    the contribution from this phase speed to the previous sum.
    # !---------------------------------------------------------------------

    c0mu0 = c0 .- ᶜu[iz0]
    c = c0 * flag + c0mu0 * (1 - flag)
    Bexp = @. exp(-log(2.0) * (c / cw)^2)
    B0 = @. sign(c0mu0) * (Bw * Bexp + Bn * Bexp)
    Bsum = sum(abs.(B0))
    if (Bsum == 0.0)
        error("zero flux input at source level")
    end
    # 	!---------------------------------------------------------------------
    #   !    define the intermittency factor eps. the factor of 1.5 is currently
    #   !    unexplained.
    #   !---------------------------------------------------------------------
    eps = (ampl * 1.5 / nk) / Bsum

    ᶜdz = ᶠz[2:end] - ᶠz[1:(end - 1)]
    # append!(ᶜdz, (ᶜz[end] - ᶠz[end]))

    wv_frcng = zeros(nc)
    gwf = zeros(length(ᶜu))
    for ink in 1:nk # loop over all wave lengths

        mask = ones(nc)  # mask to determine waves that propagate upward
        for k in source_level:length(ᶜu)
            fac = 0.5 * (ᶜρ[k] / ᶜρ[source_level]) * kwv[ink] / ᶜbf[k]

            ᶜHb = - ᶜdz[k] / log(ᶜρ[k] / ᶜρ[k - 1])  # density scale height
            alp2 = 0.25 / (ᶜHb * ᶜHb)
            omc = sqrt((ᶜbf[k] * ᶜbf[k] * k2[ink]) / (k2[ink] + alp2)) # ω_r (critical frequency that marks total internal reflection)

            fm = FT(0)
            fe = FT(0)
            for n in 1:nc
                # !----------------------------------------------------------------------
                # !    check only those waves which are still propagating, i.e., msk = 1.0
                # !----------------------------------------------------------------------
                if (mask[n]) == 1.0
                    c0mu = c0[n] - ᶜu[k]
                    # !----------------------------------------------------------------------
                    # !    if phase speed matches the wind speed, remove c0(n) from the 
                    # !    set of propagating waves.
                    # !----------------------------------------------------------------------
                    if c0mu == 0.0
                        mask[n] = 0.0
                    else
                        # !---------------------------------------------------------------------
                        # !    define the criterion which determines if wave is reflected at this 
                        # !    level (test).
                        # !---------------------------------------------------------------------
                        test = abs(c0mu) * kwv[ink] - omc
                        if test >= 0.0
                            # !---------------------------------------------------------------------
                            # !    wave has undergone total internal reflection. remove it from the
                            # !    propagating set.
                            # !---------------------------------------------------------------------
                            mask[n] = 0.0
                        else
                            # !---------------------------------------------------------------------
                            # !    if wave is  not reflected at this level, determine if it is 
                            # !    breaking at this level (Foc >= 0),  or if wave speed relative to 
                            # !    windspeed has changed sign from its value at the source level 
                            # !    (c0mu0(n)*c0mu <= 0). if it is  above the source level and is
                            # !    breaking, then add its momentum flux to the accumulated sum at 
                            # !    this level, and increase the effective diffusivity accordingly. 
                            # !    set flag to remove phase speed c0(n) from the set of active waves
                            # !    moving upwards to the next level.
                            # !---------------------------------------------------------------------
                            if c0mu0[n] * c0mu <= 0.0
                                mask[n] = 0.0
                                if k < source_level
                                    fm = fm + B0[n]
                                    fe = fe + c0mu * B0[n]
                                end
                            else
                                Foc = B0[n] / (c0mu)^3 - fac
                                if Foc >= 0.0
                                    mask[n] = 0.0
                                    if k < source_level
                                        fm = fm + B0[n]
                                        fe = fe + c0mu * B0[n]
                                    end
                                end
                            end
                        end # (test >= 0.0)
                    end #(c0mu == 0.0)
                end # mask = 0

            end # phase speed loop
            # TODO: option to dump remaining flux at the top of the model 

            # !----------------------------------------------------------------------
            # !    compute the gravity wave momentum flux forcing and eddy 
            # !    diffusion coefficient obtained across the entire wave spectrum
            # !    at this level.
            # !----------------------------------------------------------------------
            if k > source_level
                rbh = sqrt(ᶜρ[k] * ᶜρ[k - 1])
                wv_frcng[k] = (ᶜρ[source_level] / rbh) * fm * eps / ᶜdz[k]
                wv_frcng[k - 1] = 0.5 * (wv_frcng[k - 1] + wv_frcng[k])
            else
                wv_frcng[k] = 0.0
            end

        end # loop over k

        for k in source_level:length(ᶜu)
            gwf[k] = gwf[k] + wv_frcng[k]
        end

    end

    # u_forcing = similar(ᶜu_field)
    # parent(u_forcing) .= reshape(gwf, (length(ᶜu), 1))
    # return u_forcing
    return reshape(gwf, (length(ᶜu), 1))
end
