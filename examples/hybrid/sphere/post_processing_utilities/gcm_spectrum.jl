"""
    SpectralSphericalMesh
Struct with mesh information.
"""
mutable struct SpectralSphericalMesh
    # grid info
    num_fourier::Int64
    num_spherical::Int64
    nλ::Int64
    nθ::Int64
    nd::Int64
    Δλ::Float64
    qwg::Array{Float64, 3}
    qnm::Array{Float64, 3}   # n,m coordinates
    wave_numbers::Array{Int64, 2}

    # variables
    var_grid::Array{Float64, 3}
    var_fourier::Array{ComplexF64, 3}
    var_spherical::Array{ComplexF64, 4}
    var_spectrum::Array{Float64, 3}

end
function SpectralSphericalMesh(nθ::Int64, nd::Int64)
    nλ = 2nθ
    Δλ = 2π / nλ

    num_fourier = Int64(floor((2 * nθ - 1) / 3)) # number of truncated zonal wavenumbers (m): minimum truncation given nθ - e.g.: nlat = 32 -> T21 (can change manually for more a severe truncation)
    num_spherical = Int64(num_fourier + 1) # number of total wavenumbers (n)

    radius = Float64(6371000)
    wave_numbers = compute_wave_numbers(num_fourier, num_spherical)

    qwg = zeros(Float64, num_fourier + 1, num_spherical + 1, nθ)
    qnm = zeros(Float64, num_fourier + 1, num_spherical + 2, nθ)

    var_fourier = zeros(Complex{Float64}, nλ, nθ, nd)
    var_grid = zeros(Float64, nλ, nθ, nd)
    nθ_half = div(nθ, 2)
    var_spherical =
        zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd, nθ_half)
    var_spectrum = zeros(Float64, num_fourier + 1, num_spherical + 1, nd)

    SpectralSphericalMesh(
        num_fourier,
        num_spherical,
        nλ,
        nθ,
        nd,
        Δλ,
        qwg,
        qnm,
        wave_numbers,
        var_grid,
        var_fourier,
        var_spherical,
        var_spectrum,
    )
end


"""
    power_spectrum_1d(var_grid, z, lat, lon, weight)
For a variable `var_grid` on a (lon,lat,z) grid, given an array of 
`weight`s, compute the zonal (1D) power spectrum using a Fourier 
traansform at each latitude, from a 3D velocity field. The input velocities
must be intepolated to a Gaussian grid.
"""
function power_spectrum_1d(var_grid, z, lat, lon, weight)
    num_lev = length(z)
    num_lat = length(lat)
    num_lon = length(lon)
    num_fourier = Int64(num_lon)

    # get number of positive Fourier coefficients incl. 0
    if mod(num_lon, 2) == 0 # even
        num_pfourier = div(num_lon, 2)
    else # odd
        num_pfourier = div(num_lon, 2) + 1
    end

    zon_spectrum = zeros(Float64, num_pfourier, num_lat, num_lev)
    freqs = zeros(Float64, num_pfourier, num_lat, num_lev)

    for k in 1:num_lev
        for j in 1:num_lat
            # compute fft frequencies for each latitude
            x = lon ./ 180.0 .* π
            dx = (lon[2] - lon[1]) ./ 180.0 .* π

            freqs_ = fftfreq(num_fourier, 1.0 / dx) # 0,+ve freq,-ve freqs (lowest to highest)
            freqs[:, j, k] = freqs_[1:num_pfourier] .* 2.0 .* π

            # compute the fourier coefficients for all latitudes
            fourier = fft(var_grid[:, j, k]) # e.g. vcos_grid, ucos_grid
            fourier = (fourier / num_fourier)

            # convert to energy spectra
            zon_spectrum[1, j, k] =
                zon_spectrum[1, j, k] +
                weight[k] * fourier[1] .* conj(fourier[1])

            for m in 2:num_pfourier
                zon_spectrum[m, j, k] =
                    zon_spectrum[m, j, k] +
                    2.0 * weight[k] * fourier[m] * conj(fourier[m]) # factor 2 for neg freq contribution
            end
        end
    end
    return zon_spectrum, freqs
end

"""
    power_spectrum_2d(var_grid, mass_weight)
- transform variable on grid to the 2d spectral space using fft on latitude circles
(as for the 1D spectrum) and Legendre polynomials for meridians, and calculate spectra
# Arguments
- var_grid: variable (typically u or v) on a Gausian (lon, lat, z) grid to be transformed
- mass_weight: weight for mass-weighted calculations
# References
- [Baer1972](@cite)
"""
function power_spectrum_2d(var_grid, mass_weight)
    #  initialize spherical mesh variables
    nθ, nd = (size(var_grid, 2), size(var_grid, 3))
    mesh = SpectralSphericalMesh(nθ, nd)
    var_spectrum = mesh.var_spectrum
    var_spherical = mesh.var_spherical

    sinθ, wts = compute_gaussian!(mesh.nθ) # latitude weights using Gaussian quadrature, to orthogonalize Legendre polynomials upon summation
    mesh.qnm =
        compute_legendre!(mesh.num_fourier, mesh.num_spherical, sinθ, mesh.nθ) #  normalized associated Legendre polynomials

    for k in 1:(mesh.nd)
        # apply Gaussian quadrature weights
        for i in 1:(mesh.nθ)
            mesh.qwg[:, :, i] .= mesh.qnm[:, :, i] * wts[i] * mass_weight[k]
        end

        # Transform variable using spherical harmonics
        var_spherical[:, :, k, :] =
            trans_grid_to_spherical!(mesh, var_grid[:, :, k]) # var_spherical[m,n,k,sinθ]

        # Calculate energy spectra
        var_spectrum[:, :, k] =
            2.0 .* sum(var_spherical[:, :, k, :], dims = 3) .*
            conj(sum(var_spherical[:, :, k, :], dims = 3))  # var_spectrum[m,n,k] # factor 2 to account for negative Fourier frequencies
        var_spectrum[1, :, k] = var_spectrum[1, :, k] ./ 2.0 # m=0
    end
    return var_spectrum, mesh.wave_numbers, var_spherical, mesh
end

"""
    compute_legendre!(num_fourier, num_spherical, sinθ, nθ)
Normalized associated Legendre polynomials, P_{m,l} = qnm
# Arguments:
- num_fourier
- num_spherical
- sinθ
- nθ
# References:
- Ehrendorfer, M. (2011) Spectral Numerical Weather Prediction Models, Appendix B, Society for Industrial and Applied Mathematics
- Winch, D. (2007) Spherical harmonics, in Encyclopedia of Geomagnetism and Paleomagnetism, Eds Gubbins D. and Herrero-Bervera, E., Springer
# Details (using notation and Eq. references from Ehrendorfer, 2011):
    l=0,1...∞    and m = -l, -l+1, ... l-1, l
    P_{0,0} = 1, such that 1/4π ∫∫YYdS = δ (where Y = spherical harmonics, S = domain surface area)
    P_{m,m} = sqrt((2m+1)/2m) cosθ P_{m-1m-1}
    P_{m+1,m} = sqrt(2m+3) sinθ P_{m m}
    sqrt((l^2-m^2)/(4l^2-1))P_{l,m} = P_{l-1, m} -  sqrt(((l-1)^2-m^2)/(4(l-1)^2 - 1))P_{l-2,m}
    THe normalization assures that 1/2 ∫_{-1}^1 P_{l,m}(sinθ) P_{n,m}(sinθ) dsinθ = δ_{n,l}
    Julia index starts with 1, so qnm[m+1,l+1] = P_l^m
"""
function compute_legendre!(num_fourier, num_spherical, sinθ, nθ)
    qnm = zeros(Float64, num_fourier + 1, num_spherical + 2, nθ)

    cosθ = sqrt.(1 .- sinθ .^ 2)
    ε = zeros(Float64, num_fourier + 1, num_spherical + 2)

    qnm[1, 1, :] .= 1
    for m in 1:num_fourier
        qnm[m + 1, m + 1, :] = -sqrt((2m + 1) / (2m)) .* cosθ .* qnm[m, m, :] # Eq. B.20
        qnm[m, m + 1, :] = sqrt(2m + 1) * sinθ .* qnm[m, m, :] # Eq. B.22
    end
    qnm[num_fourier + 1, num_fourier + 2, :] =
        sqrt(2 * (num_fourier + 2)) * sinθ .*
        qnm[num_fourier + 1, num_fourier + 1, :]

    for m in 0:num_fourier
        for l in (m + 2):(num_spherical + 1)
            ε1 = sqrt(((l - 1)^2 - m^2) ./ (4 * (l - 1)^2 - 1))
            ε2 = sqrt((l^2 - m^2) ./ (4 * l^2 - 1))
            qnm[m + 1, l + 1, :] =
                (sinθ .* qnm[m + 1, l, :] - ε1 * qnm[m + 1, l - 1, :]) / ε2 # Eq. B.18
        end
    end

    return qnm[:, 1:(num_spherical + 1), :]
end

"""
    compute_gaussian!(n)
Compute sin(latitude) and the weight factors for Gaussian integration
# Arguments
- n: number of latitudes
# References
- Ehrendorfer, M., Spectral Numerical Weather Prediction Models, Appendix B, Society for Industrial and Applied Mathematics, 2011
# Details (following notation from Ehrendorfer, 2011):
    Pn(x) is an odd function
    solve half of the n roots and weightes of Pn(x) # n = 2n_half
    P_{-1}(x) = 0
    P_0(x) = 1
    P_1(x) = x
    nP_n(x) = (2n-1)xP_{n-1}(x) - (n-1)P_{n-2}(x)
    P'_n(x) = n/(x^2-1)(xP_{n}(x) - P_{n-1}(x))
    x -= P_n(x)/P'_{n}()
    Initial guess xi^{0} = cos(π(i-0.25)/(n+0.5))
    wi = 2/(1-xi^2)/P_n'(xi)^2
"""
function compute_gaussian!(n)
    itermax = 10000
    tol = 1.0e-15

    sinθ = zeros(Float64, n)
    wts = zeros(Float64, n)

    n_half = Int64(n / 2)
    for i in 1:n_half
        dp = 0.0
        z = cos(pi * (i - 0.25) / (n + 0.5))
        for iter in 1:itermax
            p2 = 0.0
            p1 = 1.0

            for j in 1:n
                p3 = p2 # Pj-2
                p2 = p1 # Pj-1
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j  #Pj
            end
            # P'_n
            dp = n * (z * p1 - p2) / (z * z - 1.0)
            z1 = z
            z = z1 - p1 / dp
            if (abs(z - z1) <= tol)
                break
            end
            if iter == itermax
                @error("Compute_Gaussian! does not converge!")
            end
        end

        sinθ[i], sinθ[n - i + 1], = -z, z
        wts[i] = wts[n - i + 1] = 2.0 / ((1.0 - z * z) * dp * dp)
    end

    return sinθ, wts
end

"""
    trans_grid_to_spherical!(mesh::SpectralSphericalMesh, pfield::Array{Float64,2})
Transforms a variable on a Gaussian grid (pfield[nλ, nθ]) into the spherical harmonics domain (var_spherical2d[num_fourier+1, num_spherical+1])
Here λ = longitude, θ = latitude, η = sinθ, m = zonal wavenumber, n = total wavenumber:
var_spherical2d = F_{m,n}    # Output variable in spectral space (Complex{Float64}[num_fourier+1, num_spherical+1])
qwg = P_{m,n}(η)w(η)         # Weighted Legendre polynomials (Float64[num_fourier+1, num_spherical+1, nθ])
var_fourier2d = g_{m, θ}     # Untruncated Fourier transformation (Complex{Float64} [nλ, nθ])
pfield = F(λ, η)             # Input variable on Gaussian grid Float64[nλ, nθ]
# Arguments
- mesh: struct with mesh information
- pfield: variable on Gaussian grid to be transformed
# References
- Ehrendorfer, M., Spectral Numerical Weather Prediction Models, Appendix B, Society for Industrial and Applied Mathematics, 2011
- [Wiin1967](@cite)
"""
function trans_grid_to_spherical!(
    mesh::SpectralSphericalMesh,
    pfield::Array{Float64, 2},
)

    num_fourier, num_spherical = mesh.num_fourier, mesh.num_spherical
    var_fourier2d, var_spherical2d =
        mesh.var_fourier[:, :, 1] * 0.0, mesh.var_spherical[:, :, 1, :] * 0.0
    nλ, nθ, nd = mesh.nλ, mesh.nθ, mesh.nd

    # Retrieve weighted Legendre polynomials
    qwg = mesh.qwg # qwg[m,n,nθ]

    # Fourier transformation
    for j in 1:nθ
        var_fourier2d[:, j] = fft(pfield[:, j], 1) / nλ
    end

    # Complete spherical harmonic transformation
    @assert(nθ % 2 == 0)
    nθ_half = div(nθ, 2)
    for m in 1:(num_fourier + 1)
        for n in m:num_spherical
            var_fourier2d_t = transpose(var_fourier2d[m, :])  # truncates var_fourier(nlon, nhlat) to (nfourier,nlat)
            if (n - m) % 2 == 0
                var_spherical2d[m, n, :] .=
                    (
                        var_fourier2d_t[1:nθ_half] .+
                        var_fourier2d_t[nθ:-1:(nθ_half + 1)]
                    ) .* qwg[m, n, 1:nθ_half] ./ 2.0
            else
                var_spherical2d[m, n, :] .=
                    (
                        var_fourier2d_t[1:nθ_half] .-
                        var_fourier2d_t[nθ:-1:(nθ_half + 1)]
                    ) .* qwg[m, n, 1:nθ_half] ./ 2.0
            end
        end
    end

    return var_spherical2d
end

function compute_wave_numbers(num_fourier::Int64, num_spherical::Int64)
    """
    See wave_numers[i,j] saves the wave number of this basis
    """
    wave_numbers = zeros(Int64, num_fourier + 1, num_spherical + 1)

    for m in 0:num_fourier
        for n in m:num_spherical
            wave_numbers[m + 1, n + 1] = n
        end
    end

    return wave_numbers

end



"""
    compute_legendre!(num_fourier, num_spherical, sinθ, nθ)
Normalized associated Legendre polynomials, P_{m,l} = qnm
# Arguments:
- num_fourier
- num_spherical
- sinθ
- nθ
# References:
- Ehrendorfer, M. (2011) Spectral Numerical Weather Prediction Models, Appendix B, Society for Industrial and Applied Mathematics
- Winch, D. (2007) Spherical harmonics, in Encyclopedia of Geomagnetism and Paleomagnetism, Eds Gubbins D. and Herrero-Bervera, E., Springer
# Details (using notation and Eq. references from Ehrendorfer, 2011):
    l=0,1...∞    and m = -l, -l+1, ... l-1, l
    P_{0,0} = 1, such that 1/4π ∫∫YYdS = δ (where Y = spherical harmonics, S = domain surface area)
    P_{m,m} = sqrt((2m+1)/2m) cosθ P_{m-1m-1}
    P_{m+1,m} = sqrt(2m+3) sinθ P_{m m}
    sqrt((l^2-m^2)/(4l^2-1))P_{l,m} = P_{l-1, m} -  sqrt(((l-1)^2-m^2)/(4(l-1)^2 - 1))P_{l-2,m}
    THe normalization assures that 1/2 ∫_{-1}^1 P_{l,m}(sinθ) P_{n,m}(sinθ) dsinθ = δ_{n,l}
    Julia index starts with 1, so qnm[m+1,l+1] = P_l^m
"""
function compute_legendre!(num_fourier, num_spherical, sinθ, nθ)
    qnm = zeros(Float64, num_fourier + 1, num_spherical + 2, nθ)

    cosθ = sqrt.(1 .- sinθ .^ 2)
    ε = zeros(Float64, num_fourier + 1, num_spherical + 2)

    qnm[1, 1, :] .= 1
    for m in 1:num_fourier
        qnm[m + 1, m + 1, :] = -sqrt((2m + 1) / (2m)) .* cosθ .* qnm[m, m, :] # Eq. B.20
        qnm[m, m + 1, :] = sqrt(2m + 1) * sinθ .* qnm[m, m, :] # Eq. B.22
    end
    qnm[num_fourier + 1, num_fourier + 2, :] =
        sqrt(2 * (num_fourier + 2)) * sinθ .*
        qnm[num_fourier + 1, num_fourier + 1, :]

    for m in 0:num_fourier
        for l in (m + 2):(num_spherical + 1)
            ε1 = sqrt(((l - 1)^2 - m^2) ./ (4 * (l - 1)^2 - 1))
            ε2 = sqrt((l^2 - m^2) ./ (4 * l^2 - 1))
            qnm[m + 1, l + 1, :] =
                (sinθ .* qnm[m + 1, l, :] - ε1 * qnm[m + 1, l - 1, :]) / ε2 # Eq. B.18
        end
    end

    return qnm[:, 1:(num_spherical + 1), :]
end

"""
    compute_gaussian!(n)
Compute sin(latitude) and the weight factors for Gaussian integration
# Arguments
- n: number of latitudes
# References
- Ehrendorfer, M., Spectral Numerical Weather Prediction Models, Appendix B, Society for Industrial and Applied Mathematics, 2011
# Details (following notation from Ehrendorfer, 2011):
    Pn(x) is an odd function
    solve half of the n roots and weightes of Pn(x) # n = 2n_half
    P_{-1}(x) = 0
    P_0(x) = 1
    P_1(x) = x
    nP_n(x) = (2n-1)xP_{n-1}(x) - (n-1)P_{n-2}(x)
    P'_n(x) = n/(x^2-1)(xP_{n}(x) - P_{n-1}(x))
    x -= P_n(x)/P'_{n}()
    Initial guess xi^{0} = cos(π(i-0.25)/(n+0.5))
    wi = 2/(1-xi^2)/P_n'(xi)^2
"""
function compute_gaussian!(n)
    itermax = 10000
    tol = 1.0e-15

    sinθ = zeros(Float64, n)
    wts = zeros(Float64, n)

    n_half = Int64(n / 2)
    for i in 1:n_half
        dp = 0.0
        z = cos(pi * (i - 0.25) / (n + 0.5))
        for iter in 1:itermax
            p2 = 0.0
            p1 = 1.0

            for j in 1:n
                p3 = p2 # Pj-2
                p2 = p1 # Pj-1
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j  #Pj
            end
            # P'_n
            dp = n * (z * p1 - p2) / (z * z - 1.0)
            z1 = z
            z = z1 - p1 / dp
            if (abs(z - z1) <= tol)
                break
            end
            if iter == itermax
                @error("Compute_Gaussian! does not converge!")
            end
        end

        sinθ[i], sinθ[n - i + 1], = -z, z
        wts[i] = wts[n - i + 1] = 2.0 / ((1.0 - z * z) * dp * dp)
    end

    return sinθ, wts
end

"""
    trans_grid_to_spherical!(mesh::SpectralSphericalMesh, pfield::Array{Float64,2})
Transforms a variable on a Gaussian grid (pfield[nλ, nθ]) into the spherical harmonics domain (var_spherical2d[num_fourier+1, num_spherical+1])
Here λ = longitude, θ = latitude, η = sinθ, m = zonal wavenumber, n = total wavenumber:
var_spherical2d = F_{m,n}    # Output variable in spectral space (Complex{Float64}[num_fourier+1, num_spherical+1])
qwg = P_{m,n}(η)w(η)         # Weighted Legendre polynomials (Float64[num_fourier+1, num_spherical+1, nθ])
var_fourier2d = g_{m, θ}     # Untruncated Fourier transformation (Complex{Float64} [nλ, nθ])
pfield = F(λ, η)             # Input variable on Gaussian grid Float64[nλ, nθ]
# Arguments
- mesh: struct with mesh information
- pfield: variable on Gaussian grid to be transformed
# References
- Ehrendorfer, M., Spectral Numerical Weather Prediction Models, Appendix B, Society for Industrial and Applied Mathematics, 2011
- [Wiin1967](@cite)
"""
function trans_grid_to_spherical!(
    mesh::SpectralSphericalMesh,
    pfield::Array{Float64, 2},
)

    num_fourier, num_spherical = mesh.num_fourier, mesh.num_spherical
    var_fourier2d, var_spherical2d =
        mesh.var_fourier[:, :, 1] * 0.0, mesh.var_spherical[:, :, 1, :] * 0.0
    nλ, nθ, nd = mesh.nλ, mesh.nθ, mesh.nd

    # Retrieve weighted Legendre polynomials
    qwg = mesh.qwg # qwg[m,n,nθ]

    # Fourier transformation
    for j in 1:nθ
        var_fourier2d[:, j] = fft(pfield[:, j], 1) / nλ
    end

    # Complete spherical harmonic transformation
    @assert(nθ % 2 == 0)
    nθ_half = div(nθ, 2)
    for m in 1:(num_fourier + 1)
        for n in m:num_spherical
            var_fourier2d_t = transpose(var_fourier2d[m, :])  # truncates var_fourier(nlon, nhlat) to (nfourier,nlat)
            if (n - m) % 2 == 0
                var_spherical2d[m, n, :] .=
                    (
                        var_fourier2d_t[1:nθ_half] .+
                        var_fourier2d_t[nθ:-1:(nθ_half + 1)]
                    ) .* qwg[m, n, 1:nθ_half] ./ 2.0
            else
                var_spherical2d[m, n, :] .=
                    (
                        var_fourier2d_t[1:nθ_half] .-
                        var_fourier2d_t[nθ:-1:(nθ_half + 1)]
                    ) .* qwg[m, n, 1:nθ_half] ./ 2.0
            end
        end
    end

    return var_spherical2d
end

function compute_wave_numbers(num_fourier::Int64, num_spherical::Int64)
    """
    See wave_numbers[i,j] saves the wave number of this basis
    """
    wave_numbers = zeros(Int64, num_fourier + 1, num_spherical + 1)

    for m in 0:num_fourier
        for n in m:num_spherical
            wave_numbers[m + 1, n + 1] = n
        end
    end
    return wave_numbers
end