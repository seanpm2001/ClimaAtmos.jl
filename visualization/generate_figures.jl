# Discontinuous Galerkin Intro visualizations 
include("dg_projection.jl")
save("DG_projections.png", fig)
include("dg_derivative.jl")
save("DG_derivatives.png", fig)
include("dg_split_form_tendencies.jl")
save("DG_split_form_tendencies.png", fig)

# Held-Suarez visualizations
include("held_suarez_statistics.jl")
save("Held_Suarez_Statistics.png", fig)
include("held_suarez_resolution.jl")
save("Held_Suarez_Resolutions.png", fig)
include("small_held_suarez_statistics.jl")
save("Small_Held_Suarez_Resolutions.png", fig)