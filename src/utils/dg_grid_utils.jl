
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