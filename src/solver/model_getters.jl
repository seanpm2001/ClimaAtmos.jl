function get_moisture_model(parsed_args)
    moisture_name = parsed_args["moist"]
    @assert moisture_name in ("dry", "equil", "nonequil")
    return if moisture_name == "dry"
        DryModel()
    elseif moisture_name == "equil"
        EquilMoistModel()
    elseif moisture_name == "nonequil"
        NonEquilMoistModel()
    elseif moisture_name == "cloudy"
        CloudyMoistModel()
    end
end

function get_model_config(parsed_args)
    config = parsed_args["config"]

    valid_configurations = ("sphere", "column", "box", "plane")

    if !(config ∈ valid_configurations)
        error_message = string(
            "config = $config is not one of the ",
            "valid configurations $valid_configurations",
        )
        throw(ArgumentError(error_message))
    end

    return if config == "sphere"
        SphericalModel()
    elseif config == "column"
        SingleColumnModel()
    elseif config == "box"
        BoxModel()
    elseif config == "plane"
        PlaneModel()
    end
end

function get_sfc_temperature_form(parsed_args)
    surface_temperature = parsed_args["surface_temperature"]
    @assert surface_temperature in ("ZonallyAsymmetric", "ZonallySymmetric")
    return if surface_temperature == "ZonallyAsymmetric"
        ZonallyAsymmetricSST()
    elseif surface_temperature == "ZonallySymmetric"
        ZonallySymmetricSST()
    end
end

function get_hyperdiffusion_model(parsed_args, ::Type{FT}) where {FT}
    hyperdiff_name = parsed_args["hyperdiff"]
    if hyperdiff_name in ("ClimaHyperdiffusion", "true", true)
        ν₄_vorticity_coeff =
            FT(parsed_args["vorticity_hyperdiffusion_coefficient"])
        ν₄_scalar_coeff = FT(parsed_args["scalar_hyperdiffusion_coefficient"])
        divergence_damping_factor = FT(parsed_args["divergence_damping_factor"])
        return ClimaHyperdiffusion(;
            ν₄_vorticity_coeff,
            ν₄_scalar_coeff,
            divergence_damping_factor,
        )
    elseif hyperdiff_name in ("CAM_SE",)
        # To match hyperviscosity coefficients in:
        #    https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2017MS001257
        #    for equation A18 and A19
        # Need to scale by (1.1e5 / (sqrt(4 * pi / 6) * 6.371e6 / (3*30)) )^3  ≈ 1.238
        ν₄_vorticity_coeff = FT(0.150 * 1.238)
        ν₄_scalar_coeff = FT(0.751 * 1.238)
        divergence_damping_factor = FT(5)
        return ClimaHyperdiffusion(;
            ν₄_vorticity_coeff,
            ν₄_scalar_coeff,
            divergence_damping_factor,
        )
    elseif hyperdiff_name in ("none", "false", false)
        return nothing
    else
        error("Uncaught hyperdiffusion model type.")
    end
end

function get_vertical_diffusion_model(
    diffuse_momentum,
    parsed_args,
    params,
    ::Type{FT},
) where {FT}
    vert_diff_name = parsed_args["vert_diff"]
    return if vert_diff_name in ("false", false, "none")
        nothing
    elseif vert_diff_name in ("true", true, "VerticalDiffusion")
        VerticalDiffusion{diffuse_momentum, FT}(; C_E = params.C_E)
    elseif vert_diff_name in ("FriersonDiffusion",)
        FriersonDiffusion{diffuse_momentum, FT}()
    else
        error("Uncaught diffusion model `$vert_diff_name`.")
    end
end

function get_surface_model(parsed_args)
    prognostic_surface_name = parsed_args["prognostic_surface"]
    return if prognostic_surface_name in
              ("false", false, "PrescribedSurfaceTemperature")
        PrescribedSurfaceTemperature()
    elseif prognostic_surface_name in
           ("true", true, "PrognosticSurfaceTemperature")
        PrognosticSurfaceTemperature()
    else
        error("Uncaught surface model `$prognostic_surface_name`.")
    end
end

function get_surface_albedo_model(parsed_args, params, ::Type{FT}) where {FT}
    albedo_name = parsed_args["albedo_model"]
    return if albedo_name in ("ConstantAlbedo",)
        ConstantAlbedo{FT}(; α = params.idealized_ocean_albedo)
    elseif albedo_name in ("RegressionFunctionAlbedo",)
        isnothing(parsed_args["rad"]) && error(
            "Radiation model not specified, so cannot use RegressionFunctionAlbedo",
        )
        RegressionFunctionAlbedo{FT}(; n = params.water_refractive_index)
    elseif albedo_name in ("CouplerAlbedo",)
        CouplerAlbedo()
    else
        error("Uncaught surface albedo model `$albedo_name`.")
    end
end

function get_viscous_sponge_model(parsed_args, params, ::Type{FT}) where {FT}
    vs_name = parsed_args["viscous_sponge"]
    return if vs_name in ("false", false, "none")
        nothing
    elseif vs_name in ("true", true, "ViscousSponge")
        zd = params.zd_viscous
        κ₂ = params.kappa_2_sponge
        ViscousSponge{FT}(; zd, κ₂)
    else
        error("Uncaught viscous sponge model `$vs_name`.")
    end
end

function get_smagorinsky_lilly_model(parsed_args, params, ::Type{FT}) where {FT}
    is_model_on = parsed_args["smagorinsky_lilly"]
    Cs = parsed_args["c_smag"]
    @assert is_model_on in (true, false)
    return if is_model_on == true
        SmagorinskyLilly{FT}(; Cs)
    else
        nothing
    end
end

function get_rayleigh_sponge_model(parsed_args, params, ::Type{FT}) where {FT}
    rs_name = parsed_args["rayleigh_sponge"]
    return if rs_name in ("false", false)
        nothing
    elseif rs_name in ("true", true, "RayleighSponge")
        zd = params.zd_rayleigh
        α_uₕ = params.alpha_rayleigh_uh
        α_w = params.alpha_rayleigh_w
        RayleighSponge{FT}(; zd, α_uₕ, α_w)
    else
        error("Uncaught rayleigh sponge model `$rs_name`.")
    end
end

function get_non_orographic_gravity_wave_model(
    parsed_args,
    model_config,
    ::Type{FT},
) where {FT}
    nogw_name = parsed_args["non_orographic_gravity_wave"]
    @assert nogw_name in (true, false)
    return if nogw_name == true
        if model_config isa SingleColumnModel
            NonOrographyGravityWave{FT}(; Bw = 1.2, Bn = 0.0, Bt_0 = 4e-3)
        elseif model_config isa SphericalModel
            NonOrographyGravityWave{FT}(;
                Bw = 0.4,
                Bn = 0.0,
                cw = 35.0,
                cw_tropics = 35.0,
                cn = 2.0,
                Bt_0 = 0.0043,
                Bt_n = 0.0,
                Bt_eq = 0.0043,
                Bt_s = 0.0,
                ϕ0_n = 15,
                ϕ0_s = -15,
                dϕ_n = 10,
                dϕ_s = -10,
            )
        else
            error("Uncaught case")
        end
    else
        nothing
    end
end

function get_orographic_gravity_wave_model(parsed_args, ::Type{FT}) where {FT}
    ogw_name = parsed_args["orographic_gravity_wave"]
    @assert ogw_name in (nothing, "gfdl_restart", "raw_topo")
    return if ogw_name == "gfdl_restart"
        OrographicGravityWave{FT, String}()
    elseif ogw_name == "raw_topo"
        OrographicGravityWave{FT, String}(topo_info = "raw_topo")
    else
        nothing
    end
end

function get_perf_mode(parsed_args)
    return if parsed_args["perf_mode"] == "PerfExperimental"
        PerfExperimental()
    else
        PerfStandard()
    end
end

function get_radiation_mode(parsed_args, ::Type{FT}) where {FT}
    idealized_h2o = parsed_args["idealized_h2o"]
    @assert idealized_h2o in (true, false)
    idealized_insolation = parsed_args["idealized_insolation"]
    @assert idealized_insolation in (true, false)
    idealized_clouds = parsed_args["idealized_clouds"]
    @assert idealized_clouds in (true, false)
    radiation_name = parsed_args["rad"]
    @assert radiation_name in (
        nothing,
        "nothing",
        "clearsky",
        "gray",
        "allsky",
        "allskywithclear",
        "DYCOMS_RF01",
        "TRMM_LBA",
    )
    return if radiation_name == "clearsky"
        RRTMGPI.ClearSkyRadiation(
            idealized_h2o,
            idealized_insolation,
            idealized_clouds,
        )
    elseif radiation_name == "gray"
        RRTMGPI.GrayRadiation(
            idealized_h2o,
            idealized_insolation,
            idealized_clouds,
        )
    elseif radiation_name == "allsky"
        RRTMGPI.AllSkyRadiation(
            idealized_h2o,
            idealized_insolation,
            idealized_clouds,
        )
    elseif radiation_name == "allskywithclear"
        RRTMGPI.AllSkyRadiationWithClearSkyDiagnostics(
            idealized_h2o,
            idealized_insolation,
            idealized_clouds,
        )
    elseif radiation_name == "DYCOMS_RF01"
        RadiationDYCOMS_RF01{FT}()
    elseif radiation_name == "TRMM_LBA"
        RadiationTRMM_LBA(FT)
    else
        nothing
    end
end

function get_precipitation_model(parsed_args)
    precip_model = parsed_args["precip_model"]
    return if precip_model == nothing || precip_model == "nothing"
        NoPrecipitation()
    elseif precip_model == "0M"
        Microphysics0Moment()
    elseif precip_model == "1M"
        Microphysics1Moment()
    elseif precip_model == "Cloudy"
        MicrophysicsCloudy()
    else
        error("Invalid precip_model $(precip_model)")
    end
end

function get_cloud_model(parsed_args)
    cloud_model = parsed_args["cloud_model"]
    return if cloud_model == "grid_scale"
        GridScaleCloud()
    elseif cloud_model == "quadrature"
        QuadratureCloud()
    else
        error("Invalid cloud_model $(cloud_model)")
    end
end

function get_forcing_type(parsed_args)
    forcing = parsed_args["forcing"]
    @assert forcing in (nothing, "held_suarez")
    return if forcing == nothing
        nothing
    elseif forcing == "held_suarez"
        HeldSuarezForcing()
    end
end

struct CallCloudDiagnosticsPerStage end
function get_call_cloud_diagnostics_per_stage(parsed_args)
    ccdps = parsed_args["call_cloud_diagnostics_per_stage"]
    @assert ccdps in (nothing, true, false)
    return if ccdps in (nothing, false)
        nothing
    elseif ccdps == true
        CallCloudDiagnosticsPerStage()
    end
end

function get_subsidence_model(parsed_args, radiation_mode, FT)
    subsidence = parsed_args["subsidence"]
    subsidence == nothing && return nothing

    prof = if subsidence == "Bomex"
        APL.Bomex_subsidence(FT)
    elseif subsidence == "LifeCycleTan2018"
        APL.LifeCycleTan2018_subsidence(FT)
    elseif subsidence == "Rico"
        APL.Rico_subsidence(FT)
    elseif subsidence == "DYCOMS"
        @assert radiation_mode isa RadiationDYCOMS_RF01
        z -> -z * radiation_mode.divergence
    else
        error("Uncaught case")
    end
    return Subsidence(prof)
end

function get_large_scale_advection_model(parsed_args, ::Type{FT}) where {FT}
    ls_adv = parsed_args["ls_adv"]
    ls_adv == nothing && return nothing

    (prof_dTdt₀, prof_dqtdt₀) = if ls_adv == "Bomex"
        (APL.Bomex_dTdt(FT), APL.Bomex_dqtdt(FT))
    elseif ls_adv == "LifeCycleTan2018"
        (APL.LifeCycleTan2018_dTdt(FT), APL.LifeCycleTan2018_dqtdt(FT))
    elseif ls_adv == "Rico"
        (APL.Rico_dTdt(FT), APL.Rico_dqtdt(FT))
    elseif ls_adv == "ARM_SGP"
        (APL.ARM_SGP_dTdt(FT), APL.ARM_SGP_dqtdt(FT))
    elseif ls_adv == "GATE_III"
        (APL.GATE_III_dTdt(FT), APL.GATE_III_dqtdt(FT))
    else
        error("Uncaught case")
    end
    # See https://clima.github.io/AtmosphericProfilesLibrary.jl/dev/
    # for which functions accept which arguments.
    prof_dqtdt = if ls_adv in ("Bomex", "LifeCycleTan2018", "Rico", "GATE_III")
        (thermo_params, ᶜts, t, z) -> prof_dqtdt₀(z)
    elseif ls_adv == "ARM_SGP"
        (thermo_params, ᶜts, t, z) ->
            prof_dqtdt₀(TD.exner(thermo_params, ᶜts), t, z)
    end
    prof_dTdt = if ls_adv in ("Bomex", "LifeCycleTan2018", "Rico")
        (thermo_params, ᶜts, t, z) ->
            prof_dTdt₀(TD.exner(thermo_params, ᶜts), z)
    elseif ls_adv == "ARM_SGP"
        (thermo_params, ᶜts, t, z) -> prof_dTdt₀(t, z)
    elseif ls_adv == "GATE_III"
        (thermo_params, ᶜts, t, z) -> prof_dTdt₀(z)
    end

    return LargeScaleAdvection(prof_dTdt, prof_dqtdt)
end

function get_external_forcing_model(parsed_args)
    external_forcing = parsed_args["external_forcing"]
    @assert external_forcing in (nothing, "GCM")
    return if isnothing(external_forcing)
        nothing
    elseif external_forcing == "GCM"
        DType = Float64  # TODO: Read from `parsed_args`
        GCMForcing{DType}(parsed_args["external_forcing_file"])
    end
end

function get_edmf_coriolis(parsed_args, ::Type{FT}) where {FT}
    edmf_coriolis = parsed_args["edmf_coriolis"]
    edmf_coriolis == nothing && return nothing
    (prof_u, prof_v) = if edmf_coriolis == "Bomex"
        (APL.Bomex_geostrophic_u(FT), z -> FT(0))
    elseif edmf_coriolis == "LifeCycleTan2018"
        (APL.LifeCycleTan2018_geostrophic_u(FT), z -> FT(0))
    elseif edmf_coriolis == "Rico"
        (APL.Rico_geostrophic_ug(FT), APL.Rico_geostrophic_vg(FT))
    elseif edmf_coriolis == "ARM_SGP"
        (z -> FT(10), z -> FT(0))
    elseif edmf_coriolis == "DYCOMS_RF01"
        (z -> FT(7), z -> FT(-5.5))
    elseif edmf_coriolis == "DYCOMS_RF02"
        (z -> FT(5), z -> FT(-5.5))
    elseif edmf_coriolis == "GABLS"
        (APL.GABLS_geostrophic_ug(FT), APL.GABLS_geostrophic_vg(FT))
    else
        error("Uncaught case")
    end

    coriolis_params = Dict()
    coriolis_params["Bomex"] = FT(0.376e-4)
    coriolis_params["LifeCycleTan2018"] = FT(0.376e-4)
    coriolis_params["Rico"] = FT(4.5e-5)
    coriolis_params["ARM_SGP"] = FT(8.5e-5)
    coriolis_params["DYCOMS_RF01"] = FT(0) # TODO: check this
    coriolis_params["DYCOMS_RF02"] = FT(0) # TODO: check this
    coriolis_params["GABLS"] = FT(1.39e-4)
    coriolis_param = coriolis_params[edmf_coriolis]
    return EDMFCoriolis(prof_u, prof_v, coriolis_param)
end

function get_turbconv_model(FT, parsed_args, turbconv_params)
    turbconv = parsed_args["turbconv"]
    @assert turbconv in
            (nothing, "edmfx", "prognostic_edmfx", "diagnostic_edmfx")

    return if turbconv == "prognostic_edmfx"
        N = parsed_args["updraft_number"]
        TKE = parsed_args["prognostic_tke"]
        PrognosticEDMFX{N, TKE}(turbconv_params.min_area)
    elseif turbconv == "diagnostic_edmfx"
        N = parsed_args["updraft_number"]
        TKE = parsed_args["prognostic_tke"]
        DiagnosticEDMFX{N, TKE}(FT(0.1), turbconv_params.min_area)
    else
        nothing
    end
end

function get_entrainment_model(parsed_args)
    entr_model = parsed_args["edmfx_entr_model"]
    return if entr_model == nothing || entr_model == "nothing"
        NoEntrainment()
    elseif entr_model == "PiGroups"
        PiGroupsEntrainment()
    elseif entr_model == "Generalized"
        GeneralizedEntrainment()
    elseif entr_model == "GeneralizedHarmonics"
        GeneralizedHarmonicsEntrainment()
    else
        error("Invalid entr_model $(entr_model)")
    end
end

function get_detrainment_model(parsed_args)
    detr_model = parsed_args["edmfx_detr_model"]
    return if detr_model == nothing || detr_model == "nothing"
        NoDetrainment()
    elseif detr_model == "PiGroups"
        PiGroupsDetrainment()
    elseif detr_model == "Generalized"
        GeneralizedDetrainment()
    elseif detr_model == "GeneralizedHarmonics"
        GeneralizedHarmonicsDetrainment()
    elseif detr_model == "ConstantArea"
        ConstantAreaDetrainment()
    else
        error("Invalid detr_model $(detr_model)")
    end
end

function get_surface_thermo_state_type(parsed_args)
    dict = Dict()
    dict["GCMSurfaceThermoState"] = GCMSurfaceThermoState()
    return dict[parsed_args["surface_thermo_state_type"]]
end

function get_tracers(parsed_args)
    aerosol_names = Tuple(parsed_args["prescribed_aerosols"])
    return (; aerosol_names)
end

function get_tendency_model(parsed_args)
    zero_tendency_name = parsed_args["zero_tendency"]
    @assert zero_tendency_name in (nothing, "grid_scale", "subgrid_scale")
    return if zero_tendency_name == "grid_scale"
        NoGridScaleTendency()
    elseif zero_tendency_name == "subgrid_scale"
        NoSubgridScaleTendency()
    elseif isnothing(zero_tendency_name)
        UseAllTendency()
    end
end
