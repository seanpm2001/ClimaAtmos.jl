#####
##### Precomputed quantities
#####
import CloudMicrophysics.Microphysics1M as CM1
import Cloudy as CL

# helper function to safely get precipitation from state
function qₚ(ρqₚ::FT, ρ::FT) where {FT}
    return max(FT(0), ρqₚ / ρ)
end

"""
    set_precipitation_precomputed_quantities!(Y, p, t)

Updates the precipitation terminal velocity stored in `p`
for the 1-moment or N-moment microphysics scheme
"""
function set_precipitation_precomputed_quantities!(Y, p, t, _)
    return nothing
end

function set_precipitation_precomputed_quantities!(Y, p, t, precip_model::Microphysics1Moment)
    (; ᶜwᵣ, ᶜwₛ, ᶜqᵣ, ᶜqₛ) = p.precomputed

    cmp = CAP.microphysics_precipitation_params(p.params)

    # compute the precipitation specific humidities
    @. ᶜqᵣ = qₚ(Y.c.ρq_rai, Y.c.ρ)
    @. ᶜqₛ = qₚ(Y.c.ρq_sno, Y.c.ρ)

    # compute the precipitation terminal velocity [m/s]
    @. ᶜwᵣ = CM1.terminal_velocity(
        cmp.pr,
        cmp.tv.rain,
        Y.c.ρ,
        abs(Y.c.ρq_rai / Y.c.ρ),
    )
    @. ᶜwₛ = CM1.terminal_velocity(
        cmp.ps,
        cmp.tv.snow,
        Y.c.ρ,
        abs(Y.c.ρq_sno / Y.c.ρ),
    )
    return nothing
end

function set_precipitation_precomputed_quantities!(Y, p, t, precip_model::MicrophysicsCloudy)
    (; ᶜqᵣ, pdists, weighted_vt) = p.precomputed

    clp = CAP.cloudy_params(p.params)

    # compute the precipitation specific humidities
    @. ᶜqᵣ = qₚ(Y.c.ρq_rai, Y.c.ρ)

    # update the pdists and weighted_vt
    @. pdists = get_updated_pdists(Y.c.moments, pdists, clp)
    @. weighted_vt = get_weighted_vt(FT, Y.c.moments, pdists, clp)
    return nothing
end

"""
   get_weighted_vt(moments, pdists, cloudy_params)

 - moments - current distribution moments
 - pdists - current particledistributions
 - cloudy_params - parameters for Cloudy specific stuff

Returns the moment-weighted terminal velocities corresponding to the input momnents,
where the terminal velocity parameters are specific in cloudy_params
"""
function get_weighted_vt(FT, moments, pdists, cloudy_params)
    sed_flux = CL.Sedimentation.get_sedimentation_flux(pdists, cloudy_params.vel)
    ntuple(length(moments)) do i
        ifelse(moments[i] > FT(0), -1 * sed_flux[i] * cloudy_params.mom_norms[i] / moments[i], FT(0))
    end
end

"""
    get_updated_pdists!(moments, old_pdists, cloudy_params)

 - moments - current distribution moments
 - old_pdists - previous particledistributions
 - cloudy_params - parameters for Cloudy specific stuff

Returns a new tuple of pdists with updated parameters based on the current moments
"""
function get_updated_pdists(moments, old_pdists, cloudy_params)
    mom_normed = moments ./ cloudy_params.mom_norms
    mom_i = get_dists_moments(mom_normed, cloudy_params.NProgMoms)
    ntuple(length(old_pdists)) do i
        if old_pdists[i] isa CL.ParticleDistributions.GammaPrimitiveParticleDistribution
            CL.ParticleDistributions.update_dist_from_moments(
                old_pdists[i],
                mom_i[i],
                param_range = (; :k => (1.0, 10.0)),
            )
        elseif old_pdists[i] isa CL.ParticleDistributions.LognormalPrimitiveParticleDistribution
            CL.ParticleDistributions.update_dist_from_moments(old_pdists[i], mom_i[i])
        else # Exponential or monodisperse
            CL.ParticleDistributions.update_dist_from_moments(old_pdists[i], mom_i[i][1:2])
        end
    end
end