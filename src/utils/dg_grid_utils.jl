
"""
function coordinates(grid::DiscontinuousSpectralElementGrid)
# Description
Gets the (x,y,z) coordinates corresponding to the grid
# Arguments
- `grid`: DiscontinuousSpectralElementGrid
# Return
- `x, y, z`: views of x, y, z coordinates
"""
function coordinates(grid::DiscontinuousSpectralElementGrid)
    x = view(grid.vgeo, :, grid.x1id, :)   # x-direction	
    y = view(grid.vgeo, :, grid.x2id, :)   # y-direction	
    z = view(grid.vgeo, :, grid.x3id, :)   # z-direction
    return x, y, z
end

"""
function massmatrix(grid; M = nothing)
# Description
Get the mass matrix of the grid
# Arguments
- `grid`: DiscontinuousSpectralElementGrid
# Return
- Tuple of cell-centers
"""
function massmatrix(grid)
    return view(grid.vgeo, :, grid.Mid, :)
end

# for computing fluxes and checking budgets, the horizontal mass matrix is useful
# Mᴴ = reshape(grid.numerical.vgeo[:, grid.numerical.MHid, :], (n_ijk..., n_e[3], 6*n_e[1]*n_e[2]))

using ClimateMachine.Mesh.Interpolation
import ClimateMachine.Mesh.Interpolation: InterpolationBrick

function InterpolationBrick(
    simulation::Simulation;
    xlength = 64,
    ylength = 64,
    zlength = 64
)
    @info "Recreating Grid"
    grid = create_grid(simulation.backend, simulation.discretized_domain)
    @info "Creating Interpolation Object"
    simulation.discretized_domain

    xmin = discretized_domain.domain[1].min .* 1.0
    xmax = discretized_domain.domain[1].max .* 1.0
    xgrid = collect(0:xlength-1) ./ (xlength - 1) * (xmax - xmin) .+ xmin

    ymin = discretized_domain.domain[2].min .* 1.0
    ymax = discretized_domain.domain[2].max .* 1.0
    ygrid = collect(0:ylength-1) ./ (ylength - 1) * (ymax - ymin) .+ ymin

    zmin = discretized_domain.domain[3].min .* 1.0
    zmax = discretized_domain.domain[3].max .* 1.0
    zgrid = collect(0:zlength-1) ./ (zlength - 1) * (zmax - zmin) .+ zmin

    boundary = [xmin xmin zmin; xmax ymax zmax] .* 1.0
    tic = time()
    ib = InterpolationBrick(
        grid,
        boundary,
        xgrid,
        ygrid,
        zgrid,
    )
    toc = time()
    @info "took $(toc-tic) seconds"
    return ib
end