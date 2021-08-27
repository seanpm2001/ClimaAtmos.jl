# use after running held_suarez
using ClimateMachine.Mesh.Interpolation

grid = simulation.rhs[1].grid

vert_range = grid1d(
        domain.radius, 
        domain.radius + domain.height, 
        nothing,
        nelem = discretized_domain.discretization.vertical.elements,
)

lat_grd = [-80, 0, 80] .* 1.0
long_grd = [-170, -10, 10, 170] .* 1.0
rad_grd = [domain.radius + 1e3, domain.radius + 2e3] .* 1.0
tic = time()
println("creating interpolation object")
interpol = InterpolationCubedSphere(
        grid,
        vert_range,
        discretized_domain.discretization.horizontal.elements,
        lat_grd,
        long_grd,
        rad_grd
) 
toc = time()
println("done creating interpolation object: took ", toc - tic, " seconds")

bl = simulation.rhs[1]
n_states = size(simulation.state)[2]
istate = similar(simulation.state.data, interpol.Npl, n_states)
interpolate_local!(interpol, simulation.state.data, istate) # istate has linear indexing

#=
# project to usual velocity components if wanted
_ρu, _ρv, _ρw = 2, 3, 4
project_cubed_sphere!(interpol, istate, (_ρu, _ρv, _ρw))
=#

# longitude latitude radial
all_state_data = accumulate_interpolated_data(MPI.COMM_WORLD, interpol, istate)

