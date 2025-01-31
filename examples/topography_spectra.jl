using LazyArtifacts
using ClimaAtmos.AtmosArtifacts
import ClimaAtmos as CA
import ClimaCore as CC
import ClimaComms
import ClimaParams as CP
import ClimaCore.Quadratures
using ClimaUtilities: SpaceVaryingInputs
using ClimaCore.Hypsography: diffuse_surface_elevation!
using ClimaCore: Remapping, Geometry, Operators, Spaces
using CairoMakie

import ClimaCoreSpectra: power_spectrum_2d

const AA = AtmosArtifacts


function mask(x::FT) where {FT}
    return x * FT(x > 0)
end

"""
    generate_spaces(; h_elem=16)
For a given number of elements `h_elem`, this function generates a 
horizontal spectral space (ClimaCore.SpectralElementSpace2D) and
generates an interpolated (using SpaceVaryingInputs) orography field
on this space. The `diffuse_surface_elevation` function is called 
to apply Laplacian diffusion (canonical diffusion equation) passes,
with a timescale defined by the attenuation required at the smallest 
resolved lengthscale.
Returns a `Field` containing surface orography. 
"""
function generate_spaces(;
    h_elem = 16,
    n_attenuation = 1000,
    planet_radius = 6.378e6,
)
    FT = Float32
    cubed_sphere_mesh =
        CA.cubed_sphere_mesh(; radius = FT(planet_radius), h_elem)
    quad = Quadratures.GLL{4}()
    comms_ctx = ClimaComms.context()
    h_space = CA.make_horizontal_space(cubed_sphere_mesh, quad, comms_ctx, true)
    Δh_scale = Spaces.node_horizontal_length_scale(h_space)
    @assert h_space isa CC.Spaces.SpectralElementSpace2D
    coords = CC.Fields.coordinate_field(h_space)
    target_field = CC.Fields.zeros(h_space)
    elev_from_file = SpaceVaryingInputs.SpaceVaryingInput(
        AA.earth_orography_file_path(; context = comms_ctx),
        "z",
        h_space,
    )
    elev_from_file = @. mask(elev_from_file)
    # Fractional damping of smallest resolved scale
    # Approximated as k₀ ≈ 1/Δx, with n_attenuation 
    # the factor by which we wish to damp wavenumber
    # k₀. Assume dt == 1 for this example.
    diff_courant = FT(0.05)
    κ = diff_courant * Δh_scale^2
    maxiter = Int(round(log(n_attenuation) / diff_courant))
    diffuse_surface_elevation!(elev_from_file; κ, dt = FT(1), maxiter)
    elev_from_file = @. mask(elev_from_file)
    return elev_from_file
end

"""
    remap_to_array(test_var, hcoords)
Given an Array of `hcoords`, this function remaps a ClimaCore
Field `test_var` onto `hcoords`. 
Returns an `Array` containing the values which make up `test_var`.
"""
function remap_to_array(test_var, hcoords)
    remapper = Remapping.Remapper(axes(test_var), hcoords)
    orog = Array(Remapping.interpolate(remapper, test_var))
end

"""
    gen_spectra(test_var)
Given a Field `test_var`, this function first calls `remap_to_array`
to remap this spectral element field onto a latitude-longitude
grid, and computes the spectra.
Returns (`w_numbers`, `power_spectrum`) which contain the spherical
wavenumbers and the spectral energy corresponding to this array
of `w_numbers`. 
"""
function gen_spectra(test_var)
    FT = eltype(test_var)
    Δh_scale = Spaces.node_horizontal_length_scale(Spaces.axes(test_var))
    planet_radius = FT(6.378e6)
    npts = Int(round(2π * planet_radius / Δh_scale))
    npts = ifelse(rem(npts, 2) == 0, npts, npts + 1)
    longpts = range(FT(-180), FT(180.0), Int(npts))
    latpts = range(FT(-90), FT(90), Int(npts // 2))
    hcoords =
        [Geometry.LatLongPoint(lat, long) for long in longpts, lat in latpts]
    test_var = remap_to_array(test_var, hcoords)
    # Use ClimaCoreSpectra to generate power spectra for orography data. 
    len1 = size(test_var)[1]
    len2 = size(test_var)[2]
    @assert len1 == 2len2
    mass_weight = FT(1) # No weighting applied 
    spectrum_data, wave_numbers, _spherical, mesh_info =
        power_spectrum_2d(FT, test_var, mass_weight)
    power_spectrum =
        dropdims(sum(spectrum_data, dims = 1), dims = 1)[begin:(end - 1), :]
    w_numbers = collect(0:1:(mesh_info.num_spherical - 1))
    @info "Returning wavenumber array (spherical) and orography power spectrum"
    return w_numbers, power_spectrum
end

"""
    generate_all_spectra(;h_elem=16)
Uses the spectral calculator and space generation tools in this example file
to generate a series of surface orography fields with different extents of 
Laplacian smoothing. 
Returns a `CairoMakie.Figure`. 
"""
function generate_all_spectra(; h_elem = 16)
    fig = Figure(; size = (1200, 900))
    ax1 = Axis(
        fig[1, 1],
        xlabel = "spherical wavenumber",
        ylabel = "elevation spectra",
        title = "Surface elevation spectra",
        xticks = [2, 4, 8, 16, 32, 64, 128, 256],
        yticks = [10^-1, 10^0, 10^1, 10^2, 10^3, 10^4],
        xscale = log10,
        yscale = log10,
        xgridvisible = true,
        ygridvisible = false,
    )
    for ii in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512)
        n_attenuation = ii
        test_var = generate_spaces(; h_elem, n_attenuation)
        Δh_scale = Spaces.node_horizontal_length_scale(Spaces.axes(test_var))
        diff_courant = 0.05 # Arbitrary example value.
        κ = diff_courant * Δh_scale^2
        maxiter = Int(round(log(n_attenuation) / diff_courant))
        sph_wn, psd = gen_spectra(test_var)
        scatterlines!(sph_wn[2:end], psd[2:end], label = "iter = $(maxiter)")
    end
    CairoMakie.Legend(fig[1, 2], ax1)
    return fig
end
