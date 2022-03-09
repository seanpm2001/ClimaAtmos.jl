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

# based on https://github.com/ali-ramadhan/Atmosfoolery.jl/blob/master/sandbox/dry_convection.jl

# set up backend
backend = DiscontinuousGalerkinBackend(numerics = (flux = :roe,),)

# to be removed
using CLIMAParameters#: AbstractEarthParameterSet
struct PlanetParameterSet <: AbstractEarthParameterSet end
get_planet_parameter(p::Symbol) = getproperty(CLIMAParameters.Planet, p)(PlanetParameterSet())

parameters = (
    R_d = get_planet_parameter(:R_d),
    pₒ = get_planet_parameter(:MSLP),
    κ = get_planet_parameter(:kappa_d),
    g = get_planet_parameter(:grav),
    cp_d = get_planet_parameter(:cp_d),
    cv_d = get_planet_parameter(:cv_d),
    γ = get_planet_parameter(:cp_d) / get_planet_parameter(:cv_d),
    T_0 = 0 * get_planet_parameter(:T_0),
    xc = 5000,
    yc = 1000,
    zc = 2000,
    rc = 2000,
    xmax = 2.56e3,
    ymax = 2.56e3,
    zmax = 2.56e3,
    θₐ = 2.0,
    cₛ = 340,
    pₛ = get_planet_parameter(:MSLP),
    Tₛ = 300.0,
    ϵ = 1e-0,
)

# set up grid
x_domain = IntervalDomain(min = 0.0, max = parameters.xmax, periodic = true)
y_domain = IntervalDomain(min = 0.0, max = parameters.ymax, periodic = true)
z_domain = IntervalDomain(min = 0.0, max = parameters.zmax, periodic = false)
discretized_domain = DiscretizedDomain(
    domain = x_domain × y_domain × z_domain,
    discretization = (
        elements = (16, 16, 16),
        polynomial_order = (3, 3, 3),
        grid_stretching = nothing,),
)

## Γ = g / c_p, Tₛ = 300
## T(z) = -g/cₚ z + Tₛ
## p = R_d ρ T and dp/dz = -g ρ => ∂ᶻ( R_d ρ (-g/cₚ z + Tₛ) ) = - g ρ
## -g R_d / cₚ ρ + (-g/cₚ z  + Tₛ) R_d ∂ᶻ(ρ) = - g ρ => (-g/cₚ z  + Tₛ) ∂ᶻ(ρ) = -g( 1/R_d + 1/cₚ) ρ
## ∂ᶻ ln(ρ) = -g( 1/R_d + 1/cₚ) / (-g/cₚ z  + Tₛ) => ln(ρ) = -( cₚ /R_d + 1) ln(-g/cₚ z  + Tₛ) + C
## ρ = ρ₀*(-g/cₚ z  + Tₛ)^(-( cₚ /R_d + 1))
# set up initial condition 

# T₀(p, x, y, z) = p.Tₛ - p.g / p.cp_d * z
# ρ₀(p, x, y, z) = (-p.g / p.cp_d * z / p.Tₛ + 1.0)^(-(p.cp_d / p.R_d + 1))
# p₀(p, x, y, z) = p.R_d * ρ₀(p, x, y, z) * T₀(p, x, y, z)

ρ₀(p, x, y, z) = 1.0 - 0.3 * z / p.zmax
p₀(p, x, y, z) = -p.g * (z - 0.3 / 2 * z^2 / p.zmax) + p.pₒ
T₀(p, x, y, z) = p₀(p, x, y, z) / (ρ₀(p, x, y, z) * p.R_d)

ρu₀(p, x, y, z) = 0.01 * @SVector [randn(), randn(), randn()]

e_pot(p, x, y, z) = p.g * z
e_int(p, x, y, z) = p.cv_d * (T₀(p, x, y, z) - p.T_0)
e_kin(p, x, y, z) = 0.0

ρe₀(p, x, y, z) = ρ₀(p, x, y, z) * (e_kin(p, x, y, z) + e_int(p, x, y, z) + e_pot(p, x, y, z)) + rand()

#=
π_exn(p, x, y, z) = 1.0 - p.g / (p.cp_d * θ₀(p, x, y, z)) * z

e_pot(p, x, y, z) = p.g * z
e_int(p, x, y, z) = p.cv_d * (θ₀(p, x, y, z) * π_exn(p, x, y, z) - p.T_0)
e_kin(p, x, y, z) = 0.0

ρ₀(p, x, y, z) = p.pₒ / (p.R_d * θ₀(p, x, y, z)) * (π_exn(p, x, y, z))^(p.cv_d / p.R_d)
ρu₀(p, x, y, z) = @SVector [0.0, 0.0, 0.0]
ρe₀(p, x, y, z) = ρ₀(p, x, y, z) * (e_kin(p, x, y, z) + e_int(p, x, y, z) + e_pot(p, x, y, z))
=#

# Radiative Forcing 
struct ConvectiveForcing{S} <: AbstractForcing
    parameters::S
end

convective_parameters = (;
    nunm = 1.0
)

function calc_source!(
    source,
    balance_law::ThreeDimensionalDryCompressibleEulerWithTotalEnergy,
    hsf::ConvectiveForcing,
    state,
    aux,
)
    FT = eltype(state)

    # Extract the state
    ρ = state.ρ
    ρu = state.ρu
    ρe = state.ρe
    Φ = aux.Φ

    x = aux.x
    y = aux.y
    z = aux.z
    coord = @SVector [x, y, z]

    Q₀ = 100.0
    ℓ = 100.0

    radiation_profile = Q₀ / ℓ * exp(-z / ℓ) * 1e1

    L = 2.56e3
    γ = 0.2 * L
    damping_profile = -exp(-(L - z) / γ)
    λ = 1 / 10.0

    # Apply convective forcing
    source.ρu += λ * damping_profile * ρu 
    source.ρe += ρ * radiation_profile 

    return nothing
end

# set up model

model = ModelSetup(
    equations = ThreeDimensionalEuler(
        thermodynamic_variable = TotalEnergy(),
        equation_of_state = DryIdealGas(),
        pressure_convention = Compressible(),
        sources = (
            # coriolis = DeepShellCoriolis(),
            gravity = Gravity(),
            forcing = ConvectiveForcing(convective_parameters),
        ),
        ref_state = NoReferenceState(),
    ),
    boundary_conditions = (DefaultBC(), DefaultBC(), DefaultBC(), DefaultBC(), DefaultBC(), DefaultBC()),
    initial_conditions = (
        ρ = ρ₀, ρu = ρu₀, ρe = ρe₀,
    ), parameters = parameters,
)


end_time = 60 * 60 * 8
Δt = 0.02
iteration_partition = 100
iterations = floor(Int, end_time / Δt / iteration_partition)


numerical_grid = create_grid(backend, discretized_domain);
cₛ = 330
Δxᵥ = min_node_distance(numerical_grid, VerticalDirection())
Δxₕ = min_node_distance(numerical_grid, HorizontalDirection())
Δt = (Δxᵥ / cₛ) * 0.1
println("The timestep is ", Δt)
vCFL = Δt / (Δxᵥ / cₛ)
hCFL = Δt / (Δxₕ / cₛ)

iteration_partition = 1000
iterations = floor(Int, end_time / Δt / iteration_partition)


println("vcfl is ", vCFL)
println("hcfl is ", hCFL)

filename = "convection_1.jld2"
jl_cb = JLD2State(iteration = iterations, filepath = filename)

# set up simulation
simulation = Simulation(
    backend = backend,
    discretized_domain = discretized_domain,
    model = model,
    timestepper = (
        method = SSPRK22Heuns,
        start = 0.0,
        finish = end_time,
        timestep = Δt,
    ),
    callbacks = (
        Info(),
        jl_cb,
        CFL(),
    ),
)

# run the simulation
tic = Base.time()
initialize!(simulation)
evolve!(simulation)
toc = Base.time()
println("The amount of time for the simulation is ", (toc - tic)/60, " minutes")

include("convection_to_uniform.jl")
