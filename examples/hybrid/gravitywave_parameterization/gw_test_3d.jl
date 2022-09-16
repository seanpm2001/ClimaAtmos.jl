using NCDatasets
using Dates
using Interpolations
using Statistics
using Plots
include("gravity_wave_parameterization.jl")
const FT = Float64

# test Figure 8 of the Alexander and Dunkerton (1999) paper:
# https://journals.ametsoc.org/view/journals/atsc/56/24/1520-0469_1999_056_4167_aspomf_2.0.co_2.xml?tab_body=pdf

face_z = 0:1e3:0.47e5
center_z = 0.5.*(face_z[1:end-1].+face_z[2:end])


# compute the source parameters
params = gravity_wave_cache(FT; Bm = 0.4, cmax = 150, kwv = 2π/100e3 )
source_level = argmin( abs.(center_z .- params.gw_source_height) )

# ERA5 data 1973 Jan
ds = NCDataset("./era5-monthly.nc")

lon = ds["longitude"][:] 
lat = ds["latitude"][:] 
lev = ds["level"][:] .* 100
time = ds["time"][:]

gZ = ds["z"][:] 
T = ds["t"][:] 
u = ds["u"][:]

# compute density and buoyancy frequency
R_d = 287.0
grav = 9.8
cp_d = 1004.0

Z = gZ./grav
ρ = ones(size(T)) .* reshape(lev,(1,1,length(lev),1)) ./ T / R_d

dTdz = zeros(size(T))
@. dTdz[:, :, 1, :] = (T[:,:,2,:] - T[:,:,1,:]) / (Z[:,:,2,:] - Z[:,:,1,:])
@. dTdz[:, :, end, :] = (T[:,:,end,:] - T[:,:,end-1,:]) / (Z[:,:,end,:] - Z[:,:,end-1,:])
@. dTdz[:, :, 2:end-1, :] = (T[:,:, 3:end, :] - T[:,:, 1:end-2, :])/(Z[:,:, 3:end, :] - Z[:,:, 1:end-2, :])
bf = @. (grav / T) * (dTdz + grav / cp_d)
bf = @. ifelse(bf < 2.5e-5, sqrt(2.5e-5), sqrt(abs(bf)))

# interpolation to center_z grid
center_u = zeros(length(lon), length(lat), length(center_z), length(time))
center_bf = zeros(length(lon), length(lat), length(center_z), length(time))
center_ρ = zeros(length(lon), length(lat), length(center_z), length(time))
for i in 1:length(lon)
  for j in 1:length(lat)
	for it in 1:length(time)
		interp_linear = LinearInterpolation(Z[i,j,:,it][end:-1:1], u[i,j,:,it][end:-1:1], extrapolation_bc=Line())
		center_u[i,j,:,it] = interp_linear.(center_z)

		interp_linear = LinearInterpolation(Z[i,j,:,it][end:-1:1], bf[i,j,:,it][end:-1:1], extrapolation_bc=Line())
		center_bf[i,j,:,it] = interp_linear.(center_z)

		interp_linear = LinearInterpolation(Z[i,j,:,it][end:-1:1], ρ[i,j,:,it][end:-1:1], extrapolation_bc=Line())
		center_ρ[i,j,:,it] = interp_linear.(center_z)
	end
  end
end

# compute zonal mean profile first and apply parameterization 
center_u_zonalave = mean(center_u, dims = 1)[1,:,:,:]
center_bf_zonalave = mean(center_bf, dims = 1)[1,:,:,:]
center_ρ_zonalave = mean(center_ρ, dims = 1)[1,:,:,:]

# Jan
month = Dates.month.(time)

Jan_u = mean(center_u_zonalave[:, :, month.==1], dims=3)[:,:,1]
Jan_bf = mean(center_bf_zonalave[:, :, month.==1], dims=3)[:,:,1]
Jan_ρ = mean(center_ρ_zonalave[:, :, month.==1], dims=3)[:,:,1]
Jan_uforcing = zeros(length(lat), length(center_z))
for j in 1:length(lat)
	Jan_uforcing[j,:] = gravity_wave_forcing(
	    Jan_u[j,:],
	    source_level,
		params.gw_F_S0,
	    params.gw_Bm,
	    params.gw_c,
		params.gw_cw,
		params.gw_c0,
	    params.gw_nk,
	    params.gw_k,
	    params.gw_k2,
	    Jan_bf[j,:],
	    Jan_ρ[j,:],
	    face_z,
	)
end

png(
	contourf( lat[end:-1:1], center_z[source_level:end], 86400*Jan_uforcing[end:-1:1, source_level:end]', color=:balance, clim=(-1, 1) ),
	"./test-fig8.png"
)
