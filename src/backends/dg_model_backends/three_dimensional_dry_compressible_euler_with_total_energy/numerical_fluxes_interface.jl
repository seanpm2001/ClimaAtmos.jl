struct RefanovFlux <: NumericalFluxFirstOrder end
struct RoefanovFlux <: NumericalFluxFirstOrder end
struct RusanovFlux <: NumericalFluxFirstOrder end
struct CentralVolumeFlux <: NumericalFluxFirstOrder end
struct KGVolumeFlux <: NumericalFluxFirstOrder end
struct LinearKGVolumeFlux <: NumericalFluxFirstOrder end
struct VeryLinearKGVolumeFlux <: NumericalFluxFirstOrder end
const super_hack = 0.0
# add entropy stable ones?
# https://github.com/CliMA/ClimateMachine.jl/blob/mw/entropy_stable_splittest/test/Numerics/ESDGMethods/DryAtmos/DryAtmos.jl
function numerical_volume_fluctuation_flux_first_order!(
    ::NumericalFluxFirstOrder,
    balance_law::ThreeDimensionalDryCompressibleEulerWithTotalEnergy,
    source::Grad,
    state_1::Vars,
    aux_1::Vars,
    state_2::Vars,
    aux_2::Vars,
)
    if haskey(balance_law.sources, :gravity)
        ρ_1, ρ_2 = state_1.ρ, state_2.ρ
        Φ_1, Φ_2 = aux_1.Φ, aux_2.Φ

        α = ave(ρ_1, ρ_2) * 0.5
        source.ρu -= α * (Φ_1 - Φ_2) * I

        # represents - ∇⋅({ρ}[ϕ]I) = ρ∇ϕ , analytically but not discretely

        # coriolis? 
        # ρu_1, ρu_2 = state_1.ρu, state_2.ρu
        # T_1 v = v⃗ × u⃗
        # - 2Ω⃗ × u⃗ = - 2*Ω (∇z × u⃗) = 2 Ω ϵᵢⱼₖ ∂ⱼz uₖ
        # 2Ω ϵᵢⱼₖuₖ = Tᵢⱼ => Tᵢⱼ ∂ⱼz, a non-conservative term
        # ∇⋅(T⃡ z) - z ∇⋅(T⃡) = T⃡ ∇z 
        # ρu = ave(ρu_1, ρu_2) * 2 * balance_law.parameters.Ω * 0.5
        # x̂ = @SVector{3}[0 0]
        # ρu₁ = ρ
        # T = @SMatrix [0.0 ρu[3]  -ρu[2] ;  -ρu[3]  0.0  ρu[1] ; ρu[2]  -ρu[1] 0.0]
        # factor of 1/2 cancels with factor of 2 from coriolis
        # - sign cancels out with the 
        # source.ρu += T * (aux_1.z - aux_2.z) 

        # source.ρu -= ρ_2 * Φ_1 * I * 0.5
        # source.ρu -= (ρ_1 * Φ_1 - ρ_1 * Φ_2 + ρ_2 * Φ_1 - 0 * ρ_2 * Φ_2) * I * 0.25 
    end
end

function numerical_volume_fluctuation_flux_first_order!(
    ::NumericalFluxFirstOrder,
    balance_law::DryLinearBalanceLaw,
    source::Grad,
    state_1::Vars,
    aux_1::Vars,
    state_2::Vars,
    aux_2::Vars,
)
    if haskey(balance_law.sources, :gravity)
        ρ_1, ρ_2 = state_1.ρ, state_2.ρ
        Φ_1, Φ_2 = aux_1.Φ, aux_2.Φ

        α = ave(ρ_1, ρ_2) * 0.5

        source.ρu -= α * (Φ_1 - Φ_2) * I # - sign because we take the negative divergence in code
        # source.ρu -= ρ_2 * Φ_1 * I * 0.5
        # source.ρu -= (ρ_1 * Φ_1 - ρ_1 * Φ_2 + ρ_2 * Φ_1 - 0 * ρ_2 * Φ_2) * I * 0.25 
        # coriolis? 
        # ρu_1, ρu_2 = state_1.ρu, state_2.ρu
        # T_1 v = v⃗ × u⃗
        # - 2Ω⃗ × u⃗ = - 2*Ω (∇z × u⃗) = 2 Ω ϵᵢⱼₖ ∂ⱼz uₖ
        # 2Ω ϵᵢⱼₖuₖ = Tᵢⱼ => Tᵢⱼ ∂ⱼz, a non-conservative term
        # factor of 1/2 cancels with factor of 2 from coriolis
        # - sign cancels out with the 

        # ρu = ave(ρu_1, ρu_2) * 2 * balance_law.parameters.Ω * 0.5
        # T = @SMatrix [0.0 ρu[3]  -ρu[2] ;  -ρu[3]  0.0  ρu[1] ; ρu[2]  -ρu[1] 0.0]
        # source.ρu += T * (aux_1.z - aux_2.z) 
    end
end

function numerical_volume_conservative_flux_first_order!(
    ::CentralVolumeFlux,
    m::ThreeDimensionalDryCompressibleEulerWithTotalEnergy,
    F::Grad,
    state_1::Vars,
    aux_1::Vars,
    state_2::Vars,
    aux_2::Vars,
)
    FT = eltype(F)
    F_1 = similar(F)
    flux_first_order!(m, F_1, state_1, aux_1, FT(0), EveryDirection())

    F_2 = similar(F)
    flux_first_order!(m, F_2, state_2, aux_2, FT(0), EveryDirection())

    parent(F) .= (parent(F_1) .+ parent(F_2)) ./ 2
end

function numerical_volume_conservative_flux_first_order!(
    ::KGVolumeFlux,
    balance_law::ThreeDimensionalDryCompressibleEulerWithTotalEnergy,
    F::Grad,
    state_1::Vars,
    aux_1::Vars,
    state_2::Vars,
    aux_2::Vars,
)
    eos = balance_law.equation_of_state
    parameters = balance_law.parameters

    ρ_1 = state_1.ρ
    ρu_1 = state_1.ρu
    ρe_1 = state_1.ρe
    u_1 = ρu_1 / ρ_1
    e_1 = ρe_1 / ρ_1
    p_1 = calc_pressure(eos, state_1, aux_1, parameters)

    ρ_2 = state_2.ρ
    ρu_2 = state_2.ρu
    ρe_2 = state_2.ρe
    u_2 = ρu_2 / ρ_2
    e_2 = ρe_2 / ρ_2
    p_2 = calc_pressure(eos, state_2, aux_2, parameters)

    ρ_avg = ave(ρ_1, ρ_2)
    u_avg = ave(u_1, u_2)
    e_avg = ave(e_1, e_2)
    p_avg = ave(p_1, p_2)

    # stable change to pressure calc
    Φ_avg = ave(aux_1.Φ, aux_2.Φ)
    Φ_1 = aux_1.Φ
    Φ_2 = aux_2.Φ
    α = ave(ρ_1, ρ_2) * 0.5
    p_avg -= 0.4 * super_hack * (ρ_avg * Φ_avg - ave(ρ_1 * Φ_1, ρ_2 * Φ_2))

    # fluxes
    F.ρ = ρ_avg * u_avg
    F.ρu = p_avg * I + ρ_avg * u_avg .* u_avg'
    # Kennedy-Gruber
    F.ρe = ρ_avg * u_avg * e_avg + p_avg * u_avg
    # Shima
    # F.ρe = ave(ρe_1, ρe_2) * u_avg + u_avg * (ρ_avg * (u_avg' * u_avg) - 0.5 * ave(ρu_1'*u_1, ρu_2'*u_2) - 0.5 * ave(ρ_1, ρ_2) * ave(u_1'*u_1, u_2'*u_2)) + (p_1 * u_2 + u_1 * p_2) * 0.5

    # regular instead of log averaging
    #=
    ρ_1, ρu_1, ρe_1 = state_1.ρ, state_1.ρu, state_1.ρe
    ρ_2, ρu_2, ρe_2 = state_2.ρ, state_2.ρu, state_2.ρe
    Φ_1, Φ_2 = aux_1.Φ, aux_2.Φ
    u_1 = ρu_1 / ρ_1
    u_2 = ρu_2 / ρ_2
    p_1 = calc_pressure(eos, state_1, aux_1, parameters)
    p_2 = calc_pressure(eos, state_2, aux_2, parameters)
    b_1 = ρ_1 / 2p_1
    b_2 = ρ_2 / 2p_2

    ρ_avg = ave(ρ_1, ρ_2)
    u_avg = ave(u_1, u_2)
    b_avg = ave(b_1, b_2)
    Φ_avg = ave(Φ_1, Φ_2)

    usq_avg = ave(dot(u_1, u_1), dot(u_2, u_2))

    ρ_log = ave(ρ_1, ρ_2)
    b_log = ave(b_1, b_2)

    γ = 0.4

    Fρ = u_avg * ρ_log
    Fρu = u_avg * Fρ' + ρ_avg / 2b_avg * I

    Fρe = (1 / (2 * (γ - 1) * b_log) - usq_avg / 2 + Φ_avg) * Fρ + Fρu * u_avg

    F.ρ += Fρ
    F.ρu += Fρu
    F.ρe += Fρe
    =#
end

# Shima {v⃗}( 1/2{ρ}⟨v⃗⋅v⃗⟩ + {ρe - 1/2 ρv⃗⋅ρv⃗/ρ} ) + ⟨p v⃗⟩
# ⟨a ⋅ b⟩ = 2{a}⋅{b} - {a⋅b}
# 1/2{ρ}⟨v⃗⋅v⃗⟩ + {ρe - 1/2 ρv⃗⋅ρv⃗/ρ}
# {ρe}{v⃗} + 1/2(2{ρ}{v⃗}⋅{v⃗} - {ρv⃗⋅v⃗} - {ρ} {v⃗⋅v⃗}    )
# {ρe}{v⃗} + {v⃗}({ρ}{v⃗}⋅{v⃗} - {ρv⃗⋅v⃗}) + 2{p}{v⃗} - {p v⃗} # Shima Term 
# {ρ}{e}{v⃗} + {p}{v⃗} # Kennedy-Gruber
# ave(ρe_1, ρe_2) * u_avg + u_avg * (ρ_avg * (u_avg' * u_avg) - ave(ρu_1'*ρu_1/ρ_1,ρu_2'*ρu_2/ρ_2) ) + (p_1 * v_2 + v_1 * p_2) * 0.5
# a¹b¹ + a²b² -> 2({a¹}{b¹}+ {a²}{b²}) - 2{a¹b¹+a²b²})

function numerical_volume_conservative_flux_first_order!(
    ::LinearKGVolumeFlux,
    balance_law::LinearThreeDimensionalDryCompressibleEulerWithTotalEnergy,
    F::Grad,
    state_1::Vars,
    aux_1::Vars,
    state_2::Vars,
    aux_2::Vars,
)

    eos = balance_law.equation_of_state
    parameters = balance_law.parameters

    ρu_1 = state_1.ρu
    ρuᵣ = ρu_1 * 0
    p_1 = calc_linear_pressure(eos, state_1, aux_1, parameters)

    # grab reference state
    ρᵣ_1 = aux_1.ref_state.ρ
    pᵣ_1 = aux_1.ref_state.p
    ρeᵣ_1 = aux_1.ref_state.ρe

    # only ρu fluctuates in the non-pressure terms
    u_1 = ρu_1 / ρᵣ_1
    eᵣ_1 = ρeᵣ_1 / ρᵣ_1
    ρu_2 = state_2.ρu

    ρuᵣ = ρu_2 * 0
    p_2 = calc_linear_pressure(eos, state_2, aux_2, parameters)

    # grab reference state
    ρᵣ_2 = aux_2.ref_state.ρ
    pᵣ_2 = aux_2.ref_state.p
    ρeᵣ_2 = aux_2.ref_state.ρe

    # only ρu fluctuates in the non-pressure terms
    u_2 = ρu_2 / ρᵣ_2
    eᵣ_2 = ρeᵣ_2 / ρᵣ_2

    # construct averages
    ρᵣ_avg = ave(ρᵣ_1, ρᵣ_2)
    eᵣ_avg = ave(eᵣ_1, eᵣ_2)
    pᵣ_avg = ave(pᵣ_1, pᵣ_2)

    u_avg = ave(u_1, u_2)
    p_avg = ave(p_1, p_2)

    F.ρ = ρᵣ_avg * u_avg
    F.ρu = p_avg * I + ρuᵣ .* ρuᵣ' # the latter term is needed to determine size of I
    F.ρe = (ρᵣ_avg * eᵣ_avg + pᵣ_avg) * u_avg
end

function numerical_volume_conservative_flux_first_order!(
    ::VeryLinearKGVolumeFlux,
    balance_law::VeryLinearThreeDimensionalDryCompressibleEulerWithTotalEnergy,
    F::Grad,
    state_1::Vars,
    aux_1::Vars,
    state_2::Vars,
    aux_2::Vars,
)
    eos = balance_law.equation_of_state
    parameters = balance_law.parameters

    ## State 1 Stuff 
    # unpack the perturbation state
    ρ_1 = state_1.ρ
    ρu_1 = state_1.ρu
    ρe_1 = state_1.ρe

    # grab reference state
    ρᵣ_1 = aux_1.ref_state.ρ
    ρuᵣ_1 = aux_1.ref_state.ρu
    ρeᵣ_1 = aux_1.ref_state.ρe
    pᵣ_1 = aux_1.ref_state.p

    # calculate pressure perturbation
    p_1 = calc_very_linear_pressure(eos, state_1, aux_1, parameters)

    # calculate u_1, e_1, and reference states
    u_1 = ρu_1 / ρᵣ_1 - ρ_1 * ρuᵣ_1 / (ρᵣ_1^2)
    e_1 = ρe_1 / ρᵣ_1 - ρ_1 * ρeᵣ_1 / (ρᵣ_1^2)

    uᵣ_1 = ρuᵣ_1 / ρᵣ_1
    eᵣ_1 = ρeᵣ_1 / ρᵣ_1

    ## State 2 Stuff 
    # unpack the state perubation
    ρ_2 = state_2.ρ
    ρu_2 = state_2.ρu
    ρe_2 = state_2.ρe

    # grab reference state
    ρᵣ_2 = aux_2.ref_state.ρ
    ρuᵣ_2 = aux_2.ref_state.ρu
    ρeᵣ_2 = aux_2.ref_state.ρe
    pᵣ_2 = aux_2.ref_state.p

    # calculate pressure perturbation
    p_2 = calc_very_linear_pressure(eos, state_2, aux_2, parameters)

    # calculate u_2, e_2, and reference states
    u_2 = ρu_2 / ρᵣ_2 - ρ_2 * ρuᵣ_2 / (ρᵣ_2^2)
    e_2 = ρe_2 / ρᵣ_2 - ρ_2 * ρeᵣ_2 / (ρᵣ_2^2)

    uᵣ_2 = ρuᵣ_2 / ρᵣ_2
    eᵣ_2 = ρeᵣ_2 / ρᵣ_2

    # construct averages for perturbation variables
    ρ_avg = ave(ρ_1, ρ_2)
    u_avg = ave(u_1, u_2)
    e_avg = ave(e_1, e_2)
    p_avg = ave(p_1, p_2)

    # construct averages for reference variables
    ρᵣ_avg = ave(ρᵣ_1, ρᵣ_2)
    uᵣ_avg = ave(uᵣ_1, uᵣ_2)
    eᵣ_avg = ave(eᵣ_1, eᵣ_2)
    pᵣ_avg = ave(pᵣ_1, pᵣ_2)

    # stable change to pressure calc
    Φ_avg = ave(aux_1.Φ, aux_2.Φ)
    p_avg -= 0.4 * super_hack * (ρ_avg * Φ_avg - ave(ρ_1 * aux_1.Φ, ρ_2 * aux_2.Φ))
    pᵣ_avg -= 0.4 * super_hack * (ρᵣ_avg * Φ_avg - ave(ρᵣ_1 * aux_1.Φ, ρᵣ_2 * aux_2.Φ))

    F.ρ = ρᵣ_avg * u_avg + ρ_avg * uᵣ_avg
    F.ρu = p_avg * I + ρᵣ_avg .* (uᵣ_avg .* u_avg' + u_avg .* uᵣ_avg')
    F.ρu += (ρ_avg .* uᵣ_avg) .* uᵣ_avg'
    F.ρe = (ρᵣ_avg * eᵣ_avg + pᵣ_avg) * u_avg
    F.ρe += (ρᵣ_avg * e_avg + ρ_avg * eᵣ_avg + p_avg) * uᵣ_avg

    # just use central flux
    # ρu_avg = ave(ρu_1, ρu_2)
    # F.ρ  = ρu_avg
    # F.ρu = p_avg * I + ρᵣ_avg .* (uᵣ_avg .* u_avg' + u_avg .* uᵣ_avg')
    # F.ρe = ave(ρu_1*eᵣ_1, ρu_2*eᵣ_2) + ave(ρu_1*pᵣ_1/ ρᵣ_1, ρu_2*pᵣ_2/ ρᵣ_2) 

    # Shima flux
    # F.ρe  = ave(ρe_1, ρe_2) * uᵣ_avg + ave(ρeᵣ_1, ρeᵣ_2) * u_avg 
    # F.ρe += u_avg * (ρᵣ_avg * (uᵣ_avg' * uᵣ_avg) - ave(ρuᵣ_1'* uᵣ_1, ρuᵣ_2'* uᵣ_2) )
    # F.ρe += uᵣ_avg * ( ρ_avg * (uᵣ_avg' * uᵣ_avg) + ρᵣ_avg * (u_avg' * uᵣ_avg)+ ρᵣ_avg * (uᵣ_avg' * u_avg) )
    # F.ρe += -uᵣ_avg * ( ave(ρu_1'* uᵣ_1, ρu_2'* uᵣ_2) + ave(ρuᵣ_1'* u_1, ρuᵣ_2'* u_2) )
    # F.ρe += (pᵣ_1 * u_2 + uᵣ_1 * p_2 + p_1 * uᵣ_2 + u_1 * pᵣ_2) * 0.5

end

# ave(ρe_1, ρe_2) * u_avg + u_avg * (ρ_avg * (u_avg' * u_avg) - ave(ρu_1'*ρu_1/ρ_1,ρu_2'*ρu_2/ρ_2) ) + (p_1 * v_2 + v_1 * p_2) * 0.5
# ave(ρe_1, ρe_2) * uᵣ_avg + ave(ρeᵣ_1, ρeᵣ_2) * uᵣ_avg 
# u_avg * (ρᵣ_avg * (uᵣ_avg' * uᵣ_avg) - ave(ρuᵣ_1'* uᵣ_1, ρuᵣ_2'* uᵣ_2) )
# uᵣ_avg * ( ρ_avg * (uᵣ_avg' * uᵣ_avg) + ρᵣ_avg * (u_avg' * uᵣ_avg)+ ρᵣ_avg * (uᵣ_avg' * u_avg) )
# -uᵣ_avg * ( ave(ρu_1'* uᵣ_1, ρu_2'* uᵣ_2) + ave(ρuᵣ_1'* u_1, ρuᵣ_2'* u_2) )
# (pᵣ_1 * v_2 + vᵣ_1 * p_2 + p_1 * vᵣ_2 + v_1 * pᵣ_2) * 0.5

function numerical_flux_first_order!(
    ::Nothing,
    ::ThreeDimensionalDryCompressibleEulerWithTotalEnergy,
    _...,
)
    return nothing
end

function numerical_flux_first_order!(
    ::RusanovNumericalFlux,
    balance_law::Union{ThreeDimensionalDryCompressibleEulerWithTotalEnergy,DryLinearBalanceLaw},
    fluxᵀn::Vars{S},
    normal_vector::SVector,
    state_prognostic⁻::Vars{S},
    state_auxiliary⁻::Vars{A},
    state_prognostic⁺::Vars{S},
    state_auxiliary⁺::Vars{A},
    t,
    direction,
) where {S,A}

    eos = balance_law.equation_of_state
    parameters = balance_law.parameters
    n⁻ = normal_vector
    state⁻ = state_prognostic⁻
    aux⁻ = state_auxiliary⁻
    state⁺ = state_prognostic⁺
    aux⁺ = state_auxiliary⁺

    ρ_1 = state⁻.ρ
    ρu_1 = state⁻.ρu
    ρe_1 = state⁻.ρe
    u_1 = ρu_1 / ρ_1
    e_1 = ρe_1 / ρ_1
    p_1 = calc_pressure(eos, state⁻, aux⁻, parameters)

    ρ_2 = state⁺.ρ
    ρu_2 = state⁺.ρu
    ρe_2 = state⁺.ρe
    u_2 = ρu_2 / ρ_2
    e_2 = ρe_2 / ρ_2
    p_2 = calc_pressure(eos, state⁺, aux⁺, parameters)

    ρ_avg = ave(ρ_1, ρ_2)
    u_avg = ave(u_1, u_2)
    e_avg = ave(e_1, e_2)
    p_avg = ave(p_1, p_2)
#=
    # fluxes
    fluxᵀn.ρ = (ρ_avg * u_avg)' * n⁻
    fluxᵀn.ρu = (p_avg * I + ρ_avg * u_avg .* u_avg')' * n⁻
    # Kennedy-Gruber
    fluxᵀn.ρe = (ρ_avg * u_avg * e_avg + p_avg * u_avg)' * n⁻
=#
    
        numerical_flux_first_order!(
            CentralNumericalFluxFirstOrder(),
            balance_law,
            fluxᵀn,
            normal_vector,
            state_prognostic⁻,
            state_auxiliary⁻,
            state_prognostic⁺,
            state_auxiliary⁺,
            t,
            direction,
        )
    
    #=
        numerical_volume_conservative_flux_first_order!(
            get_volume_flux(balance_law),
            balance_law,
            fluxᵀn,
            state⁻,
            aux⁻,
            state⁺,
            aux⁺,
        )
    =#
    eos = balance_law.equation_of_state
    parameters = balance_law.parameters

    cv_d = parameters.cv_d
    T_0 = parameters.T_0

    Φ = state_auxiliary⁻.Φ #Φ⁻ and Φ⁺ have the same value
    ρ⁻ = state_prognostic⁻.ρ
    ρu⁻ = state_prognostic⁻.ρu
    ρe⁻ = state_prognostic⁻.ρe
    u⁻ = ρu⁻ / ρ⁻

    p⁻ = calc_pressure(eos, state_prognostic⁻, state_auxiliary⁻, parameters)
    c⁻ = calc_sound_speed(eos, state_prognostic⁻, state_auxiliary⁻, parameters)
    h⁻ = calc_total_specific_enthalpy(eos, state_prognostic⁻, state_auxiliary⁻, parameters)

    ρ⁺ = state_prognostic⁺.ρ
    ρu⁺ = state_prognostic⁺.ρu
    ρe⁺ = state_prognostic⁺.ρe
    u⁺ = ρu⁺ / ρ⁺

    p⁺ = calc_pressure(eos, state_prognostic⁺, state_auxiliary⁺, parameters)
    c⁺ = calc_sound_speed(eos, state_prognostic⁺, state_auxiliary⁺, parameters)
    h⁺ = calc_total_specific_enthalpy(eos, state_prognostic⁺, state_auxiliary⁺, parameters)

    ρ̃ = sqrt(ρ⁻ * ρ⁺)
    ũ = roe_average(ρ⁻, ρ⁺, u⁻, u⁺)
    h̃ = roe_average(ρ⁻, ρ⁺, h⁻, h⁺)
    c̃ = sqrt(roe_average(ρ⁻, ρ⁺, c⁻^2, c⁺^2))

    ũᵀn = ũ' * normal_vector

    Δρ = ρ⁺ - ρ⁻
    Δp = p⁺ - p⁻
    Δρu = ρu⁺ - ρu⁻
    Δρe = ρe⁺ - ρe⁻
    Δu = u⁺ - u⁻
    Δuᵀn = Δu' * normal_vector

    w1 = abs(ũᵀn - c̃) * (Δp - ρ̃ * c̃ * Δuᵀn) / (2 * c̃^2)
    w2 = abs(ũᵀn + c̃) * (Δp + ρ̃ * c̃ * Δuᵀn) / (2 * c̃^2)
    w3 = abs(ũᵀn) * (Δρ - Δp / c̃^2)
    w4 = abs(ũᵀn) * ρ̃

    α = 0.0
    fluxᵀn.ρ -= c̃ * Δρ    * ( 1.0 - α ) / 2  * 0.9
    fluxᵀn.ρu -= c̃ * Δρu  * ( 1.0 - α ) / 2  * 0.9 
    fluxᵀn.ρe -= c̃ * Δρe  * ( 1.0 - α ) / 2  * 0.9 

    
    fluxᵀn.ρ -= (w1 + w2 + w3) / 2 * α
    fluxᵀn.ρu -=
        (
            w1 * (ũ - c̃ * normal_vector) +
            w2 * (ũ + c̃ * normal_vector) +
            w3 * ũ +
            w4 * (Δu - Δuᵀn * normal_vector)
        ) / 2 * α
    fluxᵀn.ρe -=
        (
            w1 * (h̃ - c̃ * ũᵀn) +
            w2 * (h̃ + c̃ * ũᵀn) +
            w3 * (ũ' * ũ / 2 + Φ - T_0 * cv_d) +
            w4 * (ũ' * Δu - ũᵀn * Δuᵀn)
        ) / 2 * α
    
end


function numerical_flux_first_order!(
    ::RoeNumericalFlux,
    balance_law::Union{ThreeDimensionalDryCompressibleEulerWithTotalEnergy,DryLinearBalanceLaw},
    fluxᵀn::Vars{S},
    normal_vector::SVector,
    state_prognostic⁻::Vars{S},
    state_auxiliary⁻::Vars{A},
    state_prognostic⁺::Vars{S},
    state_auxiliary⁺::Vars{A},
    t,
    direction,
) where {S,A}

    eos = balance_law.equation_of_state
    parameters = balance_law.parameters
    n⁻ = normal_vector
    state⁻ = state_prognostic⁻
    aux⁻ = state_auxiliary⁻
    state⁺ = state_prognostic⁺
    aux⁺ = state_auxiliary⁺

    ρ_1 = state⁻.ρ
    ρu_1 = state⁻.ρu
    ρe_1 = state⁻.ρe
    u_1 = ρu_1 / ρ_1
    e_1 = ρe_1 / ρ_1
    p_1 = calc_pressure(eos, state⁻, aux⁻, parameters)

    ρ_2 = state⁺.ρ
    ρu_2 = state⁺.ρu
    ρe_2 = state⁺.ρe
    u_2 = ρu_2 / ρ_2
    e_2 = ρe_2 / ρ_2
    p_2 = calc_pressure(eos, state⁺, aux⁺, parameters)

    ρ_avg = ave(ρ_1, ρ_2)
    u_avg = ave(u_1, u_2)
    e_avg = ave(e_1, e_2)
    p_avg = ave(p_1, p_2)

    # fluxes
    fluxᵀn.ρ = (ρ_avg * u_avg)' * n⁻
    fluxᵀn.ρu = (p_avg * I + ρ_avg * u_avg .* u_avg')' * n⁻
    # Kennedy-Gruber
    fluxᵀn.ρe = (ρ_avg * u_avg * e_avg + p_avg * u_avg)' * n⁻

    #=
        numerical_flux_first_order!(
            CentralNumericalFluxFirstOrder(),
            balance_law,
            fluxᵀn,
            normal_vector,
            state_prognostic⁻,
            state_auxiliary⁻,
            state_prognostic⁺,
            state_auxiliary⁺,
            t,
            direction,
        )
    =#
    #=
        numerical_volume_conservative_flux_first_order!(
            get_volume_flux(balance_law),
            balance_law,
            fluxᵀn,
            state⁻,
            aux⁻,
            state⁺,
            aux⁺,
        )
    =#
    eos = balance_law.equation_of_state
    parameters = balance_law.parameters

    cv_d = parameters.cv_d
    T_0 = parameters.T_0

    Φ = state_auxiliary⁻.Φ #Φ⁻ and Φ⁺ have the same value
    ρ⁻ = state_prognostic⁻.ρ
    ρu⁻ = state_prognostic⁻.ρu
    u⁻ = ρu⁻ / ρ⁻

    p⁻ = calc_pressure(eos, state_prognostic⁻, state_auxiliary⁻, parameters)
    c⁻ = calc_sound_speed(eos, state_prognostic⁻, state_auxiliary⁻, parameters)
    h⁻ = calc_total_specific_enthalpy(eos, state_prognostic⁻, state_auxiliary⁻, parameters)

    ρ⁺ = state_prognostic⁺.ρ
    ρu⁺ = state_prognostic⁺.ρu
    u⁺ = ρu⁺ / ρ⁺

    p⁺ = calc_pressure(eos, state_prognostic⁺, state_auxiliary⁺, parameters)
    c⁺ = calc_sound_speed(eos, state_prognostic⁺, state_auxiliary⁺, parameters)
    h⁺ = calc_total_specific_enthalpy(eos, state_prognostic⁺, state_auxiliary⁺, parameters)

    ρ̃ = sqrt(ρ⁻ * ρ⁺)
    ũ = roe_average(ρ⁻, ρ⁺, u⁻, u⁺)
    h̃ = roe_average(ρ⁻, ρ⁺, h⁻, h⁺)
    c̃ = sqrt(roe_average(ρ⁻, ρ⁺, c⁻^2, c⁺^2))

    ũᵀn = ũ' * normal_vector

    Δρ = ρ⁺ - ρ⁻
    Δp = p⁺ - p⁻
    Δu = u⁺ - u⁻
    Δuᵀn = Δu' * normal_vector

    w1 = abs(ũᵀn - c̃) * (Δp - ρ̃ * c̃ * Δuᵀn) / (2 * c̃^2)
    w2 = abs(ũᵀn + c̃) * (Δp + ρ̃ * c̃ * Δuᵀn) / (2 * c̃^2)
    w3 = abs(ũᵀn) * (Δρ - Δp / c̃^2)
    w4 = abs(ũᵀn) * ρ̃

    fluxᵀn.ρ -= (w1 + w2 + w3) / 2
    fluxᵀn.ρu -=
        (
            w1 * (ũ - c̃ * normal_vector) +
            w2 * (ũ + c̃ * normal_vector) +
            w3 * ũ +
            w4 * (Δu - Δuᵀn * normal_vector)
        ) / 2
    fluxᵀn.ρe -=
        (
            w1 * (h̃ - c̃ * ũᵀn) +
            w2 * (h̃ + c̃ * ũᵀn) +
            w3 * (ũ' * ũ / 2 + Φ - T_0 * cv_d) +
            w4 * (ũ' * Δu - ũᵀn * Δuᵀn)
        ) / 2
end

function numerical_flux_first_order!(
    ::LMARSNumericalFlux,
    balance_law::ThreeDimensionalDryCompressibleEulerWithTotalEnergy,
    fluxᵀn::Vars{S},
    normal_vector::SVector,
    state_prognostic⁻::Vars{S},
    state_auxiliary⁻::Vars{A},
    state_prognostic⁺::Vars{S},
    state_auxiliary⁺::Vars{A},
    t,
    direction,
) where {S,A}
    FT = eltype(fluxᵀn)
    eos = balance_law.equation_of_state
    parameters = balance_law.parameters

    ρ⁻ = state_prognostic⁻.ρ
    ρu⁻ = state_prognostic⁻.ρu
    u⁻ = ρu⁻ / ρ⁻
    uᵀn⁻ = u⁻' * normal_vector

    p⁻ = calc_pressure(eos, state_prognostic⁻, state_auxiliary⁻, parameters)
    # if the reference state is removed in the momentum equations (meaning p-p_ref is used for pressure gradient force) then we should remove the reference pressure
    #     p⁻ -= state_auxiliary⁻.ref_state.p
    # end
    c⁻ = calc_sound_speed(eos, state_prognostic⁻, state_auxiliary⁻, parameters)
    h⁻ = calc_total_specific_enthalpy(eos, state_prognostic⁻, state_auxiliary⁻, parameters)

    ρ⁺ = state_prognostic⁺.ρ
    ρu⁺ = state_prognostic⁺.ρu
    u⁺ = ρu⁺ / ρ⁺
    uᵀn⁺ = u⁺' * normal_vector

    p⁺ = calc_pressure(eos, state_prognostic⁺, state_auxiliary⁺, parameters)
    # if the reference state is removed in the momentum equations (meaning p-p_ref is used for pressure gradient force) then we should remove the reference pressure
    #     p⁺ -= state_auxiliary⁺.ref_state.p
    # end
    # c⁺ = calc_sound_speed(eos, state_prognostic⁺, state_auxiliary⁺, parameters)
    h⁺ = calc_total_specific_enthalpy(eos, state_prognostic⁺, state_auxiliary⁺, parameters)

    # Eqn (49), (50), β the tuning parameter
    β = FT(1)
    u_half = 1 / 2 * (uᵀn⁺ + uᵀn⁻) - β * 1 / (ρ⁻ + ρ⁺) / (c⁻+c⁺) * (p⁺ - p⁻)
    p_half = 1 / 2 * (p⁺ + p⁻) - β * ((ρ⁻ + ρ⁺) * c⁻) / 4 * (uᵀn⁺ - uᵀn⁻)

    # Eqn (46), (47)
    ρ_b = u_half > FT(0) ? ρ⁻ : ρ⁺
    ρu_b = u_half > FT(0) ? ρu⁻ : ρu⁺
    ρh_b = u_half > FT(0) ? ρ⁻ * h⁻ : ρ⁺ * h⁺

    # Update fluxes Eqn (18)
    fluxᵀn.ρ = ρ_b * u_half
    fluxᵀn.ρu = ρu_b * u_half .+ p_half * normal_vector
    fluxᵀn.ρe = ρh_b * u_half
end

function numerical_flux_first_order!(
    ::RefanovFlux,
    balance_law::Union{ThreeDimensionalDryCompressibleEulerWithTotalEnergy,DryLinearBalanceLaw},
    fluxᵀn::Vars{S},
    normal_vector::SVector,
    state⁻::Vars{S},
    aux⁻::Vars{A},
    state⁺::Vars{S},
    aux⁺::Vars{A},
    t,
    direction,
) where {S,A}

    numerical_flux_first_order!(
        CentralNumericalFluxFirstOrder(),
        balance_law,
        fluxᵀn,
        normal_vector,
        state⁻,
        aux⁻,
        state⁺,
        aux⁺,
        t,
        direction,
    )
    eos = balance_law.equation_of_state
    parameters = balance_law.parameters

    c⁻ = calc_ref_sound_speed(eos, state⁻, aux⁻, parameters)
    c⁺ = calc_ref_sound_speed(eos, state⁺, aux⁺, parameters)
    c = max(c⁻, c⁺)

    # - states
    ρ⁻ = state⁻.ρ
    ρu⁻ = state⁻.ρu
    ρe⁻ = state⁻.ρe

    # + states
    ρ⁺ = state⁺.ρ
    ρu⁺ = state⁺.ρu
    ρe⁺ = state⁺.ρe

    Δρ = ρ⁺ - ρ⁻
    Δρu = ρu⁺ - ρu⁻
    Δρe = ρe⁺ - ρe⁻

    fluxᵀn.ρ -= c * Δρ   
    fluxᵀn.ρu -= c * Δρu 
    fluxᵀn.ρe -= c * Δρe 

end


function numerical_flux_first_order!(
    ::RoefanovFlux,
    balance_law::Union{ThreeDimensionalDryCompressibleEulerWithTotalEnergy,DryLinearBalanceLaw},
    fluxᵀn::Vars{S},
    n⁻::SVector,
    state⁻::Vars{S},
    aux⁻::Vars{A},
    state⁺::Vars{S},
    aux⁺::Vars{A},
    t,
    direction::Tuple{EveryDirection,VerticalDirection},
) where {S,A}
    #=
        numerical_flux_first_order!(
            CentralNumericalFluxFirstOrder(),
            balance_law,
            fluxᵀn,
            n⁻,
            state⁻,
            aux⁻,
            state⁺,
            aux⁺,
            t,
            direction,
        )
    =#


    if balance_law isa ThreeDimensionalDryCompressibleEulerWithTotalEnergy
        eos = balance_law.equation_of_state
        parameters = balance_law.parameters

        ρ_1 = state⁻.ρ
        ρu_1 = state⁻.ρu
        ρe_1 = state⁻.ρe
        u_1 = ρu_1 / ρ_1
        e_1 = ρe_1 / ρ_1
        p_1 = calc_pressure(eos, state⁻, aux⁻, parameters)

        ρ_2 = state⁺.ρ
        ρu_2 = state⁺.ρu
        ρe_2 = state⁺.ρe
        u_2 = ρu_2 / ρ_2
        e_2 = ρe_2 / ρ_2
        p_2 = calc_pressure(eos, state⁺, aux⁺, parameters)

        ρ_avg = ave(ρ_1, ρ_2)
        u_avg = ave(u_1, u_2)
        e_avg = ave(e_1, e_2)
        p_avg = ave(p_1, p_2)

        # fluxes
        fluxᵀn.ρ = (ρ_avg * u_avg)' * n⁻
        fluxᵀn.ρu = (p_avg * I + ρ_avg * u_avg .* u_avg')' * n⁻
        # Kennedy-Gruber
        fluxᵀn.ρe = (ρ_avg * u_avg * e_avg + p_avg * u_avg)' * n⁻

        # fluxes
        #=
        fluxᵀn.ρ  = (ave(ρu_1, ρu_2))' * n⁻
        fluxᵀn.ρu = (p_avg * I + ave(ρu_1 *ρu_1' / ρ_1, ρu_2 * ρu_2' / ρ_2 ))' * n⁻ 
        # Kennedy-Gruber
        fluxᵀn.ρe = (ave(ρu_1 /ρ_1 * (ρe_1 + p_1), ρu_2 /ρ_2 * (ρe_2 + p_2) ))' * n⁻
        =#
    else

        eos = balance_law.equation_of_state
        parameters = balance_law.parameters

        ## State 1 Stuff 
        # unpack the perturbation state
        ρ_1 = state⁻.ρ
        ρu_1 = state⁻.ρu
        ρe_1 = state⁻.ρe

        # grab reference state
        ρᵣ_1 = aux⁻.ref_state.ρ
        ρuᵣ_1 = aux⁻.ref_state.ρu
        ρeᵣ_1 = aux⁻.ref_state.ρe
        pᵣ_1 = aux⁻.ref_state.p

        # calculate pressure perturbation
        p_1 = calc_very_linear_pressure(eos, state⁻, aux⁻, parameters)

        # calculate u_1, e_1, and reference states
        u_1 = ρu_1 / ρᵣ_1 - ρ_1 * ρuᵣ_1 / (ρᵣ_1^2)
        e_1 = ρe_1 / ρᵣ_1 - ρ_1 * ρeᵣ_1 / (ρᵣ_1^2)

        uᵣ_1 = ρuᵣ_1 / ρᵣ_1
        eᵣ_1 = ρeᵣ_1 / ρᵣ_1

        ## State 2 Stuff 
        # unpack the state perubation
        ρ_2 = state⁺.ρ
        ρu_2 = state⁺.ρu
        ρe_2 = state⁺.ρe

        # grab reference state
        ρᵣ_2 = aux⁺.ref_state.ρ
        ρuᵣ_2 = aux⁺.ref_state.ρu
        ρeᵣ_2 = aux⁺.ref_state.ρe
        pᵣ_2 = aux⁺.ref_state.p

        # calculate pressure perturbation
        p_2 = calc_very_linear_pressure(eos, state⁺, aux⁺, parameters)

        # calculate u_2, e_2, and reference states
        u_2 = ρu_2 / ρᵣ_2 - ρ_2 * ρuᵣ_2 / (ρᵣ_2^2)
        e_2 = ρe_2 / ρᵣ_2 - ρ_2 * ρeᵣ_2 / (ρᵣ_2^2)

        uᵣ_2 = ρuᵣ_2 / ρᵣ_2
        eᵣ_2 = ρeᵣ_2 / ρᵣ_2

        # construct averages for perturbation variables
        ρ_avg = ave(ρ_1, ρ_2)
        u_avg = ave(u_1, u_2)
        e_avg = ave(e_1, e_2)
        p_avg = ave(p_1, p_2)

        # construct averages for reference variables
        ρᵣ_avg = ave(ρᵣ_1, ρᵣ_2)
        uᵣ_avg = ave(uᵣ_1, uᵣ_2)
        eᵣ_avg = ave(eᵣ_1, eᵣ_2)
        pᵣ_avg = ave(pᵣ_1, pᵣ_2)

        fluxᵀn.ρ = (ρᵣ_avg * u_avg + ρ_avg * uᵣ_avg)' * n⁻
        fluxᵀn.ρu = (p_avg * I + ρᵣ_avg .* (uᵣ_avg .* u_avg' + u_avg .* uᵣ_avg'))' * n⁻
        fluxᵀn.ρu += ((ρ_avg .* uᵣ_avg) .* uᵣ_avg')' * n⁻
        fluxᵀn.ρe = ((ρᵣ_avg * eᵣ_avg + pᵣ_avg) * u_avg)' * n⁻
        fluxᵀn.ρe += ((ρᵣ_avg * e_avg + ρ_avg * eᵣ_avg + p_avg) * uᵣ_avg)' * n⁻

    end

    extra_fac = 1.0 # between 0 and 1
    eos = balance_law.equation_of_state
    parameters = balance_law.parameters

    c⁻ = calc_ref_sound_speed(eos, state⁻, aux⁻, parameters)
    c⁺ = calc_ref_sound_speed(eos, state⁺, aux⁺, parameters)
    c = max(c⁻, c⁺)

    # - states
    ρ⁻ = state⁻.ρ
    ρu⁻ = state⁻.ρu
    ρe⁻ = state⁻.ρe

    # + states
    ρ⁺ = state⁺.ρ
    ρu⁺ = state⁺.ρu
    ρe⁺ = state⁺.ρe

    Δρ = ρ⁺ - ρ⁻
    Δρu = ρu⁺ - ρu⁻
    Δρe = ρe⁺ - ρe⁻

    fluxᵀn.ρ -= c * Δρ * 0.5 * extra_fac
    fluxᵀn.ρu -= c * Δρu * 0.5 * extra_fac
    fluxᵀn.ρe -= c * Δρe * 0.5 * extra_fac

    # try it

    eos = balance_law.equation_of_state
    parameters = balance_law.parameters

    cv_d = parameters.cv_d
    T_0 = parameters.T_0

    Φ = aux⁻.Φ # Φ⁻ and Φ⁺ have the same value

    # - states
    ρ⁻ = state⁻.ρ
    ρu⁻ = state⁻.ρu
    ρe⁻ = state⁻.ρe

    ρ⁻ᵣ = aux⁻.ref_state.ρ
    ρu⁻ᵣ = aux⁻.ref_state.ρu

    # constructed states
    u⁻ = ρu⁻ / ρ⁻ᵣ - ρ⁻ * ρu⁻ᵣ / (ρ⁻ᵣ * ρ⁻ᵣ)
    u⁻ᵣ = ρu⁻ᵣ / ρ⁻ᵣ

    # in general thermodynamics
    p⁻ = calc_very_linear_pressure(eos, state⁻, aux⁻, parameters)
    c⁻ = calc_ref_sound_speed(eos, state⁻, aux⁻, parameters)

    h⁻ = calc_ref_enthalpy(eos, state⁻, aux⁻, parameters)

    # + states
    ρ⁺ = state⁺.ρ
    ρu⁺ = state⁺.ρu

    ρ⁺ᵣ = aux⁺.ref_state.ρ
    ρu⁺ᵣ = aux⁺.ref_state.ρu

    # constructed states
    u⁺ = ρu⁺ / ρ⁺ᵣ - ρ⁺ * ρu⁺ᵣ / (ρ⁺ᵣ * ρ⁺ᵣ)
    u⁺ᵣ = ρu⁺ᵣ / ρ⁺ᵣ

    # in general thermodynamics
    p⁺ = calc_very_linear_pressure(eos, state⁺, aux⁺, parameters)
    c⁺ = calc_ref_sound_speed(eos, state⁺, aux⁺, parameters)
    h⁺ = calc_ref_enthalpy(eos, state⁺, aux⁺, parameters)

    # construct roe averges
    ρ = sqrt(ρ⁻ᵣ * ρ⁺ᵣ)
    u = roe_average(ρ⁻ᵣ, ρ⁺ᵣ, u⁻ᵣ, u⁺ᵣ)
    h = roe_average(ρ⁻ᵣ, ρ⁺ᵣ, h⁻, h⁺)
    c = roe_average(ρ⁻ᵣ, ρ⁺ᵣ, c⁻, c⁺)

    # construct normal velocity
    uₙ = u' * n⁻

    # differences
    Δρ = ρ⁺ - ρ⁻
    Δp = p⁺ - p⁻
    Δu = u⁺ - u⁻
    Δuₙ = Δu' * n⁻

    # idea: keep linearize Δ's, use reference values for everything else
    # uₙ , c, u, h, ρ, as reference values 
    # w1, w2, w3, can be linear in Δp, Δu, Δρ
    # Δp, Δu, Δuₙ can use the linearized variants
    # w4 can use reference values
    # Δp, ρᵣ * c * Δu 

    # constructed values
    c⁻² = 1 / (c * c)
    w1 = abs(uₙ - c) * (Δp - ρ * c * Δuₙ) * 0.5 * c⁻²
    w2 = abs(uₙ + c) * (Δp + ρ * c * Δuₙ) * 0.5 * c⁻²
    w3 = abs(uₙ) * (Δρ - Δp * c⁻²)
    w4 = abs(uₙ) * ρ

    # fluxes

    fluxᵀn.ρ -= (w1 + w2 + w3) * 0.5 * (1.0 - extra_fac)
    fluxᵀn.ρu -=
        (
            w1 * (u - c * n⁻) +
            w2 * (u + c * n⁻) +
            w3 * u +
            w4 * (Δu - Δuₙ * n⁻)
        ) * 0.5 * (1.0 - extra_fac)
    fluxᵀn.ρe -=
        (
            w1 * (h - c * uₙ) +
            w2 * (h + c * uₙ) +
            w3 * (u' * u / 2 + Φ - T_0 * cv_d) +
            w4 * (u' * Δu - uₙ * Δuₙ)
        ) * 0.5 * (1.0 - extra_fac)

end


# Refanov Hardcoding hack
function numerical_flux_first_order!(
    ::RoefanovFlux,
    balance_law::Union{ThreeDimensionalDryCompressibleEulerWithTotalEnergy,DryLinearBalanceLaw},
    fluxᵀn::Vars{S},
    normal_vector::SVector,
    state⁻::Vars{S},
    aux⁻::Vars{A},
    state⁺::Vars{S},
    aux⁺::Vars{A},
    t,
    direction::Tuple{EveryDirection,HorizontalDirection},
) where {S,A}

    numerical_flux_first_order!(
        RoeNumericalFlux(),
        balance_law,
        fluxᵀn,
        normal_vector,
        state⁻,
        aux⁻,
        state⁺,
        aux⁺,
        t,
        direction,
    )

end

#=
function numerical_flux_first_order!(
    ::HackyRoeNumericalFlux,
    balance_law::ThreeDimensionalCompressibleEulerWithBarotropicFluid,
    fluxᵀn::Vars{S},
    n⁻::SVector,
    state⁻::Vars{S},
    aux⁻::Vars{A},
    state⁺::Vars{S},
    aux⁺::Vars{A},
    t,
    direction,
) where {S, A}

    numerical_flux_first_order!(
        CentralNumericalFluxFirstOrder(),
        balance_law,
        fluxᵀn,
        n⁻,
        state⁻,
        aux⁻,
        state⁺,
        aux⁺,
        t,
        direction,
    )
    eos = balance_law.equation_of_state
    parameters = balance_law.parameters

    cv_d = parameters.cv_d
    T_0  = parameters.T_0

    Φ = state_auxiliary⁻.Φ # Φ⁻ and Φ⁺ have the same value

    # - states
    ρ⁻ = state⁻.ρ
    ρu⁻ = state⁻.ρu
    ρe⁻ = state⁻.ρe

    ρ⁻ᵣ = aux⁻.ref_state.ρ
    ρu⁻ᵣ = aux⁻.ref_state.ρu
    ρe⁻ᵣ = aux⁻.ref_state.ρe

    # constructed states
    u⁻  = ρu⁻ / ρ⁻ᵣ - ρ⁻ * ρu⁻ᵣ / ρ⁻ᵣ^2
    u⁻ᵣ = ρu⁻ᵣ / ρ⁻ᵣ 

    # in general thermodynamics
    p⁻ = calc_very_linear_pressure(eos, state⁻, aux⁻, parameters)
    c⁻ = calc_ref_sound_speed(eos, state⁻, aux⁻, parameters)
    h⁻ = calc_ref_enthalpy(eos, state⁻, aux⁻, parameters)

    # + states
    ρ⁺ = state⁺.ρ
    ρu⁺ = state⁺.ρu
    ρe⁺ = state⁺.ρe

    # constructed states
    u⁺ = ρu⁺ / ρ⁺ᵣ - ρ⁺ * ρu⁺ᵣ / ρ⁺ᵣ^2
    u⁺ᵣ = ρu⁺ᵣ / ρ⁺ᵣ 

    # in general thermodynamics
    p⁺ = calc_very_linear_pressure(eos, state⁺, aux⁺, parameters)
    c⁺ = calc_ref_sound_speed(eos, state⁺, aux⁺, parameters)
    h⁺ = calc_ref_enthalpy(eos, state⁺, aux⁺, parameters)

    # construct roe averges
    ρ = sqrt(ρ⁻ᵣ * ρ⁺ᵣ)
    u = roe_average(ρ⁻ᵣ, ρ⁺ᵣ, u⁻ᵣ, u⁺ᵣ)
    h = roe_average(ρ⁻ᵣ, ρ⁺ᵣ, h⁻, h⁺)
    c = roe_average(ρ⁻ᵣ, ρ⁺ᵣ, c⁻, c⁺)

    # construct normal velocity
    uₙ = u' * n⁻

    # differences
    Δρ = ρ⁺ - ρ⁻
    Δp = p⁺ - p⁻
    Δu = u⁺ - u⁻
    Δuₙ = Δu' * n⁻

    # idea: keep linearize Δ's, use reference values for everything else
    # uₙ , c, u, h, ρ, as reference values 
    # w1, w2, w3, can be linear in Δp, Δu, Δρ
    # Δp, Δu, Δuₙ can use the linearized variants
    # w4 can use reference values
    # Δp, ρᵣ * c * Δu 

    # constructed values
    c⁻² = 1 / c^2
    w1 = abs(uₙ - c) * (Δp - ρ * c * Δuₙ) * 0.5 * c⁻²
    w2 = abs(uₙ + c) * (Δp + ρ * c * Δuₙ) * 0.5 * c⁻²
    w3 = abs(uₙ) * (Δρ - Δp * c⁻²)
    w4 = abs(uₙ) * ρ

    # fluxes
    fluxᵀn.ρ -= (w1 + w2 + w3) * 0.5
    fluxᵀn.ρu -=
        (
            w1 * (u - c * n⁻) +
            w2 * (u + c * n⁻) +
            w3 * u +
            w4 * (Δu - Δuₙ * n⁻)
        ) * 0.5
    fluxᵀn.ρe -=
        (
            w1 * (h - c * uₙ) +
            w2 * (h + c * uₙ) +
            w3 * (u' * u / 2 + Φ - T_0 * cv_d) +
            w4 * (u' * Δu - uₙ * Δuₙ)
        ) / 2
end
=#
function numerical_flux_second_order!(
    ::Nothing,
    ::Union{ThreeDimensionalDryCompressibleEulerWithTotalEnergy,DryLinearBalanceLaw},
    _...,
)
    return nothing
end

# utils
roe_average(ρ⁻, ρ⁺, var⁻, var⁺) =
    (sqrt(ρ⁻) * var⁻ + sqrt(ρ⁺) * var⁺) / (sqrt(ρ⁻) + sqrt(ρ⁺))

function wavespeed(
    model::ThreeDimensionalDryCompressibleEulerWithTotalEnergy,
    n⁻,
    state::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    eos = model.equation_of_state
    parameters = model.parameters
    ρ = state.ρ
    ρu = state.ρu

    u = ρu / ρ
    u_norm = abs(dot(n⁻, u))
    return u_norm + calc_sound_speed(eos, state, aux, parameters)
end