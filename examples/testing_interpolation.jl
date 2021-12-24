using StaticArrays
include("../src/interface/domains.jl")
include("../src/interface/models.jl")
include("../src/interface/physics.jl")
include("../src/interface/boundary_conditions.jl")
include("../src/interface/grids.jl")
include("../src/backends/backends.jl")
include("../src/interface/timestepper_abstractions.jl")
include("../src/backends/dg_model_backends/backend_hook.jl")
include("../src/interface/simulations.jl")
include("../src/interface/callbacks.jl")
include("../src/backends/dg_model_backends/boilerplate.jl")
include("../src/utils/sphere_utils.jl")

# to be removed
using CLIMAParameters#: AbstractEarthParameterSet
struct PlanetParameterSet <: AbstractEarthParameterSet end
get_planet_parameter(p::Symbol) = getproperty(CLIMAParameters.Planet, p)(PlanetParameterSet())

# set up backend
backend = DiscontinuousGalerkinBackend(numerics = (flux = :refanov,),)

parameters = (
    a    = get_planet_parameter(:planet_radius),
    Ω    = get_planet_parameter(:Omega),
    g    = get_planet_parameter(:grav),
    κ    = get_planet_parameter(:kappa_d),
    R_d  = get_planet_parameter(:R_d),
    cv_d = get_planet_parameter(:cv_d),
    cp_d = get_planet_parameter(:cp_d),
    γ    = get_planet_parameter(:cp_d)/get_planet_parameter(:cv_d),
    H    = 30e3,
    pₒ   = 1.0e5,
    k    = 3.0,
    Γ    = 0.005,
    T_0  = 0.0,
    T_E  = 310.0,
    T_P  = 240.0,
    b    = 2.0,
    z_t  = 15e3,
    λ_c  = π / 9,
    ϕ_c  = 2 * π / 9,
    V_p  = 1.0,
    day = 86400,
    p0 = 1e5,
    T_ref = 255,
)

# Set up grid
domain = SphericalShell(
    radius = parameters.a,
    height = parameters.H,
)
discretized_domain = DiscretizedDomain(
    domain = domain,
    discretization = (
	    horizontal = SpectralElementGrid(elements = 15, polynomial_order = 2),
	    vertical = SpectralElementGrid(elements = 7, polynomial_order = 2)
	),
)

# For testing interpolation
uˡᵒⁿ(𝒫,λ,ϕ,r)   = 1.0
uˡᵃᵗ(𝒫,λ,ϕ,r)   = 1.0
uʳᵃᵈ(𝒫,λ,ϕ,r)   = 1.0

ρ₀(𝒫,λ,ϕ,r)    = r

ρuˡᵒⁿ(𝒫,λ,ϕ,r) = ρ₀(𝒫,λ,ϕ,r) * uˡᵒⁿ(𝒫,λ,ϕ,r)
ρuˡᵃᵗ(𝒫,λ,ϕ,r) = ρ₀(𝒫,λ,ϕ,r) * uˡᵃᵗ(𝒫,λ,ϕ,r)
ρuʳᵃᵈ(𝒫,λ,ϕ,r) = ρ₀(𝒫,λ,ϕ,r) * uʳᵃᵈ(𝒫,λ,ϕ,r)

ρe(𝒫,λ,ϕ,r) = ϕ /π * 180# ρ₀(𝒫,λ,ϕ,r) * (e_int(𝒫,λ,ϕ,r) + e_kin(𝒫,λ,ϕ,r) + e_pot(𝒫,λ,ϕ,r))

# Cartesian Representation (boiler plate really)
ρ₀ᶜᵃʳᵗ(𝒫, x...)  = ρ₀(𝒫, lon(x...), lat(x...), rad(x...))
ρu₀ᶜᵃʳᵗ(𝒫, x...) = (   ρuʳᵃᵈ(𝒫, lon(x...), lat(x...), rad(x...)) * r̂(x...)
                     + ρuˡᵃᵗ(𝒫, lon(x...), lat(x...), rad(x...)) * ϕ̂(x...)
                     + ρuˡᵒⁿ(𝒫, lon(x...), lat(x...), rad(x...)) * λ̂(x...) )
ρe₀ᶜᵃʳᵗ(𝒫, x...) = ρe(𝒫, lon(x...), lat(x...), rad(x...))

# set up reference state
ref_state = DryReferenceState(DecayingTemperatureProfile{FT}(parameters, FT(290), FT(220), FT(8e3)))

# Set up model
model = ModelSetup(
    equations = ThreeDimensionalEuler(
        thermodynamic_variable = TotalEnergy(),
        equation_of_state = DryIdealGas(),
        pressure_convention = Compressible(),
        sources = (
            coriolis = DeepShellCoriolis(),
            gravity = Gravity(),
        ),
        ref_state = ref_state,
    ),
    boundary_conditions = (DefaultBC(), DefaultBC()),
    initial_conditions = (
        ρ = ρ₀ᶜᵃʳᵗ, ρu = ρu₀ᶜᵃʳᵗ, ρe = ρe₀ᶜᵃʳᵗ,
    ),
    parameters = parameters,
)

# set up shadyCFL
function shady_timestep(discretized_domain::DiscretizedDomain; vcfl = 16, hcfl = 1.0, sound_speed = 330)
    # vertical cfl
    height = domain.height
    ne = discretized_domain.discretization.vertical.elements
    np = discretized_domain.discretization.vertical.polynomial_order
    vdt = height / ne / (np^2 + 1) / sound_speed * vcfl
    @info "vertical cfl implies dt=$vdt"
    # horizontal cfl
    circumference = domain.radius * 2π
    ne = discretized_domain.discretization.horizontal.elements * 4 # since 4 faces on cubed sphere equator
    np = discretized_domain.discretization.horizontal.polynomial_order
    hdt = circumference / ne / (np^2 + 1) / sound_speed * hcfl
    @info "horizontal cfl implies dt=$hdt"

    if vdt < hdt 
        dt = vdt
        @info "limited by vertical acoustic modes dt=$dt seconds"
    else
        dt = hdt
        @info "limited by horizontal acoustic modes dt=$dt seconds"

    end

    return dt
end

dt = shady_timestep(discretized_domain)

# set up simulation
simulation = Simulation(
    backend = backend,
    discretized_domain = discretized_domain,
    model = model,
    splitting = IMEXSplitting( linear_model = :linear, ),
    timestepper = (
        method = IMEX(),
        start = 0.0,
        finish = 1200 * 24 * 3600,
        timestep = dt,
    ),
    callbacks = (
        Info(),
    ),
)