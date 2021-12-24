abstract type AbstractTimestepper end
abstract type AbstractRate end
abstract type AbstractSplitting end 
abstract type AbstractCFL end

# struct CFL{𝒜} <: AbstractCFL
#     cfl::𝒜
# end

# Base.@kwdef struct Timestepper{𝒜,ℬ,𝒞,𝒟} <: AbstractTimestepper
#     method::𝒜
#     start::ℬ
#     finish::𝒞
#     timestep::𝒟
# end

# function Timestepper(; method = SSPRK22Heuns, start = 0.0, finish, timestep = CFL(1.0))

# end

struct NoSplitting <: AbstractSplitting end

Base.@kwdef struct IMEXSplitting{𝒜,ℬ,𝒞} <: AbstractSplitting
    linear_model::𝒞 = :linear
    implicit_method::𝒜 = LinearBackwardEulerSolver(ManyColumnLU(); isadjustable = false)
    split_explicit_implicit::ℬ = false
end

# TODO: Add more methods here such as MultiRate, Explicit [can't reuse word]
Base.@kwdef struct IMEX{ℱ}
    method::ℱ
end
# ARK1ForwardBackwardEuler
# ARK2ImplicitExplicitMidpoint
# ARK2GiraldoKellyConstantinescu # by far the best
# ARK548L2SA2KennedyCarpenter
# ARK437L2SA1KennedyCarpenter

function IMEX()
    return IMEX(ARK2GiraldoKellyConstantinescu )
end

function construct_odesolver(::NoSplitting, simulation)
    method        = simulation.timestepper.method
    start         = simulation.timestepper.start
    timestep      = simulation.timestepper.timestep
    rhs           = simulation.rhs
    state         = simulation.state

    ode_solver = method(
        rhs,
        state;
        dt = timestep,
        t0 = start,
    )

    return ode_solver
end

function construct_odesolver(splitting::IMEXSplitting, simulation; t0 = 0, split_explicit_implicit = false)
    method       = simulation.timestepper.method.method
    start        = simulation.timestepper.start
    timestep     = simulation.timestepper.timestep
    state        = simulation.state 

    explicit_rhs    = simulation.rhs[1]
    implicit_rhs    = simulation.rhs[2]

    implicit_method         = splitting.implicit_method
    split_explicit_implicit = splitting.split_explicit_implicit

    odesolver = method(
        explicit_rhs,
        implicit_rhs,
        implicit_method,
        state;
        dt = timestep,
        t0 = start,
        split_explicit_implicit = split_explicit_implicit,
    )

    return odesolver
end
