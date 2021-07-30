"""
    Declaration of state variables

    vars_state returns a NamedTuple of data types.
"""
function vars_state(balance_law::VeryLinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy, aux::Auxiliary, FT)
    @vars begin
        x::FT
        y::FT
        z::FT
        Φ::FT
        ref_state::vars_state(balance_law, balance_law.ref_state, aux, FT)
    end
end

vars_state(::VeryLinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy, ::DryReferenceState, ::Auxiliary, FT) =
    @vars(T::FT, p::FT, ρ::FT, ρu::SVector{3, FT}, ρe::FT, ρq::FT)

vars_state(::VeryLinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy, ::NoReferenceState, ::Auxiliary, FT) =
    @vars()

function vars_state(::VeryLinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy, ::Prognostic, FT)
    @vars begin
        ρ::FT
        ρu::SVector{3, FT}
        ρe::FT
        ρq::FT
    end
end

"""
    Initialization of state variables

    init_state_xyz! sets up the initial fields within our state variables
    (e.g., prognostic, auxiliary, etc.), however it seems to not initialized
    the gradient flux variables by default.
"""
function init_state_prognostic!(
        balance_law::VeryLinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy,
        state::Vars,
        aux::Vars,
        localgeo,
        t
    )
    x = aux.x
    y = aux.y
    z = aux.z

    parameters = balance_law.parameters
    ic = balance_law.initial_conditions

    # TODO!: Set to 0 by default or assign IC
    if !isnothing(ic)
        state.ρ  = ic.ρ(parameters, x, y, z)
        state.ρu = ic.ρu(parameters, x, y, z)
        state.ρe = ic.ρe(parameters, x, y, z)
        state.ρq = ic.ρq(parameters, x, y, z)
    end

    return nothing
end

function nodal_init_state_auxiliary!(
    balance_law::VeryLinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy,
    state_auxiliary,
    tmp,
    geom,
)
    init_state_auxiliary!(balance_law, balance_law.orientation, state_auxiliary, geom)
    init_state_auxiliary!(balance_law, balance_law.ref_state, state_auxiliary, geom)
end

function init_state_auxiliary!(
    balance_law::VeryLinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy,
    ::SphericalOrientation,
    state_auxiliary,
    geom,
)
    g = balance_law.parameters.g

    r = norm(geom.coord)
    state_auxiliary.x = geom.coord[1]
    state_auxiliary.y = geom.coord[2]
    state_auxiliary.z = geom.coord[3]
    state_auxiliary.Φ = g * r
end

function init_state_auxiliary!(
    balance_law::VeryLinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy,
    ::FlatOrientation,
    state_auxiliary,
    geom,
)
    g = balance_law.parameters.g

    r = geom.coord[3]
    state_auxiliary.x = geom.coord[1]
    state_auxiliary.y = geom.coord[2]
    state_auxiliary.z = geom.coord[3]
    state_auxiliary.Φ = g * r
end

function init_state_auxiliary!(
    ::VeryLinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy,
    ::NoReferenceState,
    state_auxiliary,
    geom,
) end

function init_state_auxiliary!(
    balance_law::VeryLinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy,
    ref_state::DryReferenceState,
    state_auxiliary,
    geom,
)
    FT = eltype(state_auxiliary)

    orientation = balance_law.orientation
    R_d         = balance_law.parameters.R_d
    γ           = balance_law.parameters.γ
    Φ           = state_auxiliary.Φ

    # Calculation of a dry reference state
    z = altitude(balance_law, orientation, geom)
    T, p = ref_state.temperature_profile(balance_law.parameters, z)
    ρ  = p / R_d / T
    ρu = SVector{3, FT}(0, 0, 0)
    ρe = p / (γ - 1) + dot(ρu, ρu) / 2ρ + ρ * Φ
    ρq = FT(0)

    state_auxiliary.ref_state.T  = T
    state_auxiliary.ref_state.p  = p
    state_auxiliary.ref_state.ρ  = ρ
    state_auxiliary.ref_state.ρu = ρu
    state_auxiliary.ref_state.ρe = ρe
    state_auxiliary.ref_state.ρq = ρq
end

"""
    Main model computations
"""
@inline function flux_first_order!(
    balance_law::VeryLinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    # states
    ρ   = state.ρ
    ρu  = state.ρu
    ρe  = state.ρe
    ρq  = state.ρq

    # thermodynamics
    eos = balance_law.equation_of_state
    parameters = balance_law.parameters
    p = calc_very_linear_pressure(eos, state, aux, parameters)

    # Reference states
    ρᵣ  = aux.ref_state.ρ
    ρuᵣ = aux.ref_state.ρu
    ρeᵣ = aux.ref_state.ρe
    ρqᵣ = aux.ref_state.ρq
    pᵣ  = aux.ref_state.p

    # derived states
    u = ρu / ρᵣ - ρ * ρuᵣ / (ρᵣ^2)
    q = ρq / ρᵣ - ρ * ρqᵣ / (ρᵣ^2)
    e = ρe / ρᵣ - ρ * ρeᵣ / (ρᵣ^2)

    # derived reference states
    uᵣ = ρuᵣ / ρᵣ
    qᵣ = ρqᵣ / ρᵣ
    eᵣ = ρeᵣ / ρᵣ

    # can be simplified, but written this way to look like the VeryLinearKGVolumeFlux
    flux.ρ   += ρᵣ * u + ρ * uᵣ # this is just ρu
    flux.ρu  += p * I + ρᵣ .* (uᵣ .* u' + u .* uᵣ')
    flux.ρu  += (ρ .* uᵣ) .* uᵣ'
    flux.ρe  += (ρᵣ * eᵣ + pᵣ) * u
    flux.ρe  += (ρᵣ * e + ρ * eᵣ + p) * uᵣ
    flux.ρq  += ρᵣ * qᵣ * u + (ρᵣ * q + ρ * qᵣ) * uᵣ

    nothing
end

"""
    Source computations
"""
function source!(
    balance_law::VeryLinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy,
    source,
    state_prognostic,
    state_auxiliary,
    _...
)
    sources = balance_law.sources
    #=
    ntuple(Val(length(sources))) do s
        Base.@_inline_meta
        calc_source!(source, balance_law, sources[s], state_prognostic, state_auxiliary)
    end
    =#
end

"""
    Utils
"""
function altitude(balance_law::VeryLinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy, ::SphericalOrientation, geom)
    return norm(geom.coord) - balance_law.parameters.a
end

function altitude(::VeryLinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy, ::FlatOrientation, geom)
    @inbounds geom.coord[3]
end
