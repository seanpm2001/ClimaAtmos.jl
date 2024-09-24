#####
##### Precomputed quantities
#####
import CloudMicrophysics.MicrophysicsNonEq as CMNe

# helper function to safely get precipitation from state
function q_cc(ρq::FT, ρ::FT) where {FT}
    return max(FT(0), ρq / ρ) # TODO - duplicated with precip
end

"""
    set_sedimentation_precomputed_quantities!(Y, p, t)

Updates the sedimentation terminal velocity stored in `p`
for the non-equilibrium microphysics scheme
"""
function set_sedimentation_precomputed_quantities!(Y, p, t)
    @assert (p.atmos.moisture_model isa MicrophysicsNonEq)

    (; ᶜwₗ, ᶜwᵢ, ᶜqₗ, ᶜqᵢ) = p.precomputed

    cmp = CAP.microphysics_sedimentation_params(p.params) # TODO - maybe don't need separate?

    # compute the precipitation specific humidities
    @. ᶜqₗ = q_cc(Y.c.ρq_liq, Y.c.ρ)
    @. ᶜqᵢ = q_cc(Y.c.ρq_ice, Y.c.ρ)

    # compute the precipitation terminal velocity [m/s]
    @. ᶜwₗ = CM1.terminal_velocity( #TODO
        cmp.pr,
        cmp.tv.rain,
        Y.c.ρ,
        abs(Y.c.ρq_rai / Y.c.ρ),
    )
    @. ᶜwᵢ = CM1.terminal_velocity( #TODO
        cmp.ps,
        cmp.tv.snow,
        Y.c.ρ,
        abs(Y.c.ρq_sno / Y.c.ρ),
    )
    return nothing
end
