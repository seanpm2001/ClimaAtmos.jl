abstract type DryLinearBalanceLaw <: BalanceLaw end
abstract type MoistLinearBalanceLaw <: BalanceLaw end

Base.@kwdef struct ThreeDimensionalCompressibleEulerWithBarotropicFluid{𝒜,ℬ,𝒞,𝒟,ℰ,ℱ} <: BalanceLaw 
    equation_of_state::𝒜
    sources::ℬ
    boundary_conditions::𝒞
    initial_conditions::𝒟
    ref_state::ℰ
    parameters::ℱ
end

Base.@kwdef struct ThreeDimensionalDryCompressibleEulerWithTotalEnergy{𝒜,ℬ,𝒞,𝒟,ℰ,ℱ,𝒢} <: BalanceLaw
    orientation::𝒜
    equation_of_state::ℬ
    sources::𝒞
    boundary_conditions::𝒟
    initial_conditions::ℰ
    ref_state::ℱ 
    parameters::𝒢
end

Base.@kwdef struct LinearThreeDimensionalDryCompressibleEulerWithTotalEnergy{𝒜,ℬ,𝒞,𝒟,ℰ,ℱ,𝒢} <: DryLinearBalanceLaw
    orientation::𝒜
    equation_of_state::ℬ
    sources::𝒞 # may not need
    boundary_conditions::𝒟
    initial_conditions::ℰ
    ref_state::ℱ 
    parameters::𝒢
end

Base.@kwdef struct VeryLinearThreeDimensionalDryCompressibleEulerWithTotalEnergy{𝒜,ℬ,𝒞,𝒟,ℰ,ℱ,𝒢} <: DryLinearBalanceLaw
    orientation::𝒜
    equation_of_state::ℬ
    sources::𝒞 # may not need
    boundary_conditions::𝒟
    initial_conditions::ℰ
    ref_state::ℱ 
    parameters::𝒢
end

Base.@kwdef struct ThreeDimensionalMoistCompressibleEulerWithTotalEnergy{𝒜,ℬ,𝒞,𝒟,ℰ,ℱ,𝒢} <: BalanceLaw
    orientation::𝒜
    equation_of_state::ℬ
    sources::𝒞
    boundary_conditions::𝒟
    initial_conditions::ℰ
    ref_state::ℱ 
    parameters::𝒢
end

Base.@kwdef struct LinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy{𝒜,ℬ,𝒞,𝒟,ℰ,ℱ,𝒢} <: MoistLinearBalanceLaw
    orientation::𝒜
    equation_of_state::ℬ
    sources::𝒞 # may not need
    boundary_conditions::𝒟
    initial_conditions::ℰ
    ref_state::ℱ 
    parameters::𝒢
end

Base.@kwdef struct VeryLinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy{𝒜,ℬ,𝒞,𝒟,ℰ,ℱ,𝒢} <: MoistLinearBalanceLaw
    orientation::𝒜
    equation_of_state::ℬ
    sources::𝒞 # may not need
    boundary_conditions::𝒟
    initial_conditions::ℰ
    ref_state::ℱ 
    parameters::𝒢
end

function create_balance_law(model::ModelSetup{𝒜}, _...) where 
    {𝒜 <: ThreeDimensionalEuler{Density, BarotropicFluid, Compressible}}

    return ThreeDimensionalCompressibleEulerWithBarotropicFluid(
        equation_of_state = model.equations.equation_of_state,
        sources = model.equations.sources,
        boundary_conditions = model.boundary_conditions,
        initial_conditions = model.initial_conditions,
        ref_state = model.equations.ref_state,
        parameters = model.parameters, 
    )
end

function create_balance_law(model::ModelSetup{𝒜}, domain) where 
    {𝒜 <: ThreeDimensionalEuler{TotalEnergy, DryIdealGas, Compressible}}

    return ThreeDimensionalDryCompressibleEulerWithTotalEnergy(
        orientation = create_orientation(domain),
        equation_of_state = model.equations.equation_of_state,
        sources = model.equations.sources,
        boundary_conditions = model.boundary_conditions,
        initial_conditions = model.initial_conditions,
        ref_state = model.equations.ref_state,
        parameters = model.parameters, 
    )
end

function create_balance_law(model::ModelSetup{𝒜}, domain) where 
    {𝒜 <: ThreeDimensionalEuler{TotalEnergy, MoistIdealGas, Compressible}}

    return ThreeDimensionalMoistCompressibleEulerWithTotalEnergy(
        orientation = create_orientation(domain),
        equation_of_state = model.equations.equation_of_state,
        sources = model.equations.sources,
        boundary_conditions = model.boundary_conditions,
        initial_conditions = model.initial_conditions,
        ref_state = model.equations.ref_state,
        parameters = model.parameters, 
    )
end

function linearize_balance_law(balance_law::ThreeDimensionalDryCompressibleEulerWithTotalEnergy) 

    return LinearThreeDimensionalDryCompressibleEulerWithTotalEnergy(
        orientation = balance_law.orientation,
        equation_of_state = balance_law.equation_of_state,
        sources = (gravity = Gravity(),),
        boundary_conditions = (DefaultBC(), DefaultBC()),
        initial_conditions = balance_law.initial_conditions,
        ref_state = balance_law.ref_state,
        parameters = balance_law.parameters, 
    )
end

function linearize_balance_law(balance_law::ThreeDimensionalMoistCompressibleEulerWithTotalEnergy) 

    return LinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy(
        orientation = balance_law.orientation,
        equation_of_state = balance_law.equation_of_state,
        sources = (gravity = Gravity(),),
        boundary_conditions = (DefaultBC(), DefaultBC()),
        initial_conditions = balance_law.initial_conditions,
        ref_state = balance_law.ref_state,
        parameters = balance_law.parameters, 
    )
end

function verylinearize_balance_law(balance_law::ThreeDimensionalDryCompressibleEulerWithTotalEnergy) 

    return VeryLinearThreeDimensionalDryCompressibleEulerWithTotalEnergy(
        orientation = balance_law.orientation,
        equation_of_state = balance_law.equation_of_state,
        sources = (gravity = Gravity(),),
        boundary_conditions = (DefaultBC(), DefaultBC()),
        initial_conditions = balance_law.initial_conditions,
        ref_state = balance_law.ref_state,
        parameters = balance_law.parameters, 
    )
end

function verylinearize_balance_law(balance_law::ThreeDimensionalMoistCompressibleEulerWithTotalEnergy) 

    return VeryLinearThreeDimensionalMoistCompressibleEulerWithTotalEnergy(
        orientation = balance_law.orientation,
        equation_of_state = balance_law.equation_of_state,
        sources = (gravity = Gravity(),),
        boundary_conditions = (DefaultBC(), DefaultBC()),
        initial_conditions = balance_law.initial_conditions,
        ref_state = balance_law.ref_state,
        parameters = balance_law.parameters, 
    )
end

function create_numerical_flux(surface_flux)
    if surface_flux == :lmars
        return LMARSNumericalFlux()
    elseif surface_flux == :roe
        return RoeNumericalFlux()
    elseif surface_flux == :refanov
        return RefanovFlux()
    elseif surface_flux == :roefanov
        return RoefanovFlux()
    elseif surface_flux == :rusanov
        return RusanovNumericalFlux()
    else
        return nothing
    end
end

create_orientation(::ProductDomain) = FlatOrientation()
create_orientation(::SphericalShell) = SphericalOrientation()