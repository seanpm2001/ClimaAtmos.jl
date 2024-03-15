using Adapt
using Dates: DateTime, @dateformat_str
using Dierckx
using Interpolations
import NCDatasets
import ClimaCore: InputOutput, Meshes, Spaces, Quadratures
import ClimaAtmos.RRTMGPInterface as RRTMGPI
import ClimaAtmos as CA
import LinearAlgebra
import ClimaCore.Fields
import ClimaTimeSteppers as CTS
import DiffEqCallbacks as DECB
import DiffEqBase: KeywordArgSilent

function get_atmos(config::AtmosConfig, params)
    (; turbconv_params) = params
    (; parsed_args) = config
    FT = eltype(config)
    moisture_model = get_moisture_model(parsed_args)
    precip_model = get_precipitation_model(parsed_args)
    cloud_model = get_cloud_model(parsed_args)
    radiation_mode = get_radiation_mode(parsed_args, FT)
    forcing_type = get_forcing_type(parsed_args)

    diffuse_momentum = !(forcing_type isa HeldSuarezForcing)

    advection_test = parsed_args["advection_test"]
    @assert advection_test in (false, true)

    gs_tendency = parsed_args["gs_tendency"]
    @assert gs_tendency in (false, true)

    edmfx_entr_model = get_entrainment_model(parsed_args)
    edmfx_detr_model = get_detrainment_model(parsed_args)

    edmfx_sgs_mass_flux = parsed_args["edmfx_sgs_mass_flux"]
    @assert edmfx_sgs_mass_flux in (false, true)

    edmfx_sgs_diffusive_flux = parsed_args["edmfx_sgs_diffusive_flux"]
    @assert edmfx_sgs_diffusive_flux in (false, true)

    edmfx_nh_pressure = parsed_args["edmfx_nh_pressure"]
    @assert edmfx_nh_pressure in (false, true)

    edmfx_velocity_relaxation = parsed_args["edmfx_velocity_relaxation"]
    @assert edmfx_velocity_relaxation in (false, true)

    implicit_diffusion = parsed_args["implicit_diffusion"]
    @assert implicit_diffusion in (true, false)

    model_config = get_model_config(parsed_args)
    vert_diff =
        get_vertical_diffusion_model(diffuse_momentum, parsed_args, params, FT)
    atmos = AtmosModel(;
        moisture_model,
        model_config,
        perf_mode = get_perf_mode(parsed_args),
        radiation_mode,
        subsidence = get_subsidence_model(parsed_args, radiation_mode, FT),
        ls_adv = get_large_scale_advection_model(parsed_args, FT),
        edmf_coriolis = get_edmf_coriolis(parsed_args, FT),
        advection_test,
        gs_tendency,
        edmfx_entr_model,
        edmfx_detr_model,
        edmfx_sgs_mass_flux,
        edmfx_sgs_diffusive_flux,
        edmfx_nh_pressure,
        edmfx_velocity_relaxation,
        precip_model,
        cloud_model,
        forcing_type,
        turbconv_model = get_turbconv_model(FT, parsed_args, turbconv_params),
        non_orographic_gravity_wave = get_non_orographic_gravity_wave_model(
            parsed_args,
            model_config,
            FT,
        ),
        orographic_gravity_wave = get_orographic_gravity_wave_model(
            parsed_args,
            FT,
        ),
        hyperdiff = get_hyperdiffusion_model(parsed_args, FT),
        vert_diff,
        diff_mode = implicit_diffusion ? Implicit() : Explicit(),
        viscous_sponge = get_viscous_sponge_model(parsed_args, params, FT),
        rayleigh_sponge = get_rayleigh_sponge_model(parsed_args, params, FT),
        sfc_temperature = get_sfc_temperature_form(parsed_args),
        surface_model = get_surface_model(parsed_args),
        surface_albedo = get_surface_albedo_model(parsed_args, FT),
        numerics = get_numerics(parsed_args),
    )
    @assert !@any_reltype(atmos, (UnionAll, DataType))

    @info "AtmosModel: \n$(summary(atmos))"
    return atmos
end

function get_numerics(parsed_args)
    test_dycore =
        parsed_args["test_dycore_consistency"] ? TestDycoreConsistency() :
        nothing

    energy_upwinding = Val(Symbol(parsed_args["energy_upwinding"]))
    tracer_upwinding = Val(Symbol(parsed_args["tracer_upwinding"]))
    edmfx_upwinding = Val(Symbol(parsed_args["edmfx_upwinding"]))
    edmfx_sgsflux_upwinding =
        Val(Symbol(parsed_args["edmfx_sgsflux_upwinding"]))

    limiter = parsed_args["apply_limiter"] ? CA.QuasiMonotoneLimiter() : nothing

    # wrap each upwinding mode in a Val for dispatch
    numerics = AtmosNumerics(;
        energy_upwinding,
        tracer_upwinding,
        edmfx_upwinding,
        edmfx_sgsflux_upwinding,
        limiter,
        test_dycore_consistency = test_dycore,
        use_reference_state = parsed_args["use_reference_state"],
    )
    @info "numerics $(summary(numerics))"

    return numerics
end

function get_spaces(parsed_args, params, comms_ctx)

    FT = eltype(params)
    z_elem = Int(parsed_args["z_elem"])
    z_max = FT(parsed_args["z_max"])
    dz_bottom = FT(parsed_args["dz_bottom"])
    dz_top = FT(parsed_args["dz_top"])
    topography = parsed_args["topography"]
    bubble = parsed_args["bubble"]
    deep = parsed_args["deep_atmosphere"]

    @assert topography in ("NoWarp", "DCMIP200", "Earth", "Agnesi", "Schar")
    if topography == "DCMIP200"
        warp_function = topography_dcmip200
    elseif topography == "Agnesi"
        warp_function = topography_agnesi
    elseif topography == "Schar"
        warp_function = topography_schar
    elseif topography == "NoWarp"
        warp_function = nothing
    elseif topography == "Earth"
        data_path = joinpath(topo_elev_dataset_path(), "ETOPO1_coarse.nc")
        array_type = ClimaComms.array_type(comms_ctx.device)
        earth_spline = NCDatasets.NCDataset(data_path) do data
            zlevels = Array(data["elevation"])
            lon = Array(data["longitude"])
            lat = Array(data["latitude"])
            # Apply Smoothing
            smooth_degree = Int(parsed_args["smoothing_order"])
            esmth = CA.gaussian_smooth(zlevels, smooth_degree)
            Adapt.adapt(
                array_type,
                linear_interpolation(
                    (lon, lat),
                    esmth,
                    extrapolation_bc = (Periodic(), Flat()),
                ),
            )
        end
        @info "Generated interpolation stencil"
        warp_function = generate_topography_warp(earth_spline)
    end
    @info "Topography" topography


    h_elem = parsed_args["h_elem"]
    radius = CAP.planet_radius(params)
    center_space, face_space = if parsed_args["config"] == "sphere"
        nh_poly = parsed_args["nh_poly"]
        quad = Quadratures.GLL{nh_poly + 1}()
        horizontal_mesh = cubed_sphere_mesh(; radius, h_elem)
        h_space =
            make_horizontal_space(horizontal_mesh, quad, comms_ctx, bubble)
        z_stretch = if parsed_args["z_stretch"]
            Meshes.GeneralizedExponentialStretching(dz_bottom, dz_top)
        else
            Meshes.Uniform()
        end
        if warp_function == nothing
            make_hybrid_spaces(h_space, z_max, z_elem, z_stretch; deep)
        else
            make_hybrid_spaces(
                h_space,
                z_max,
                z_elem,
                z_stretch;
                parsed_args = parsed_args,
                surface_warp = warp_function,
                deep,
            )
        end
    elseif parsed_args["config"] == "column" # single column
        @warn "perturb_initstate flag is ignored for single column configuration"
        FT = eltype(params)
        Δx = FT(1) # Note: This value shouldn't matter, since we only have 1 column.
        quad = Quadratures.GL{1}()
        horizontal_mesh = periodic_rectangle_mesh(;
            x_max = Δx,
            y_max = Δx,
            x_elem = 1,
            y_elem = 1,
        )
        if bubble
            @warn "Bubble correction not compatible with single column configuration. It will be switched off."
            bubble = false
        end
        h_space =
            make_horizontal_space(horizontal_mesh, quad, comms_ctx, bubble)
        z_stretch = if parsed_args["z_stretch"]
            Meshes.GeneralizedExponentialStretching(dz_bottom, dz_top)
        else
            Meshes.Uniform()
        end
        make_hybrid_spaces(h_space, z_max, z_elem, z_stretch; parsed_args)
    elseif parsed_args["config"] == "box"
        FT = eltype(params)
        nh_poly = parsed_args["nh_poly"]
        quad = Quadratures.GLL{nh_poly + 1}()
        x_elem = Int(parsed_args["x_elem"])
        x_max = FT(parsed_args["x_max"])
        y_elem = Int(parsed_args["y_elem"])
        y_max = FT(parsed_args["y_max"])
        horizontal_mesh = periodic_rectangle_mesh(;
            x_max = x_max,
            y_max = y_max,
            x_elem = x_elem,
            y_elem = y_elem,
        )
        h_space =
            make_horizontal_space(horizontal_mesh, quad, comms_ctx, bubble)
        z_stretch = if parsed_args["z_stretch"]
            Meshes.GeneralizedExponentialStretching(dz_bottom, dz_top)
        else
            Meshes.Uniform()
        end
        make_hybrid_spaces(
            h_space,
            z_max,
            z_elem,
            z_stretch;
            parsed_args,
            surface_warp = warp_function,
            deep,
        )
    elseif parsed_args["config"] == "plane"
        FT = eltype(params)
        nh_poly = parsed_args["nh_poly"]
        quad = Quadratures.GLL{nh_poly + 1}()
        x_elem = Int(parsed_args["x_elem"])
        x_max = FT(parsed_args["x_max"])
        horizontal_mesh =
            periodic_line_mesh(; x_max = x_max, x_elem = x_elem)
        h_space =
            make_horizontal_space(horizontal_mesh, quad, comms_ctx, bubble)
        z_stretch = if parsed_args["z_stretch"]
            Meshes.GeneralizedExponentialStretching(dz_bottom, dz_top)
        else
            Meshes.Uniform()
        end
        make_hybrid_spaces(
            h_space,
            z_max,
            z_elem,
            z_stretch;
            parsed_args,
            surface_warp = warp_function,
            deep,
        )
    end
    ncols = Fields.ncolumns(center_space)
    ndofs_total = ncols * z_elem
    hspace = Spaces.horizontal_space(center_space)
    quad_style = Spaces.quadrature_style(hspace)
    Nq = Quadratures.degrees_of_freedom(quad_style)

    @info "Resolution stats: " Nq h_elem z_elem ncols ndofs_total
    return (;
        center_space,
        face_space,
        horizontal_mesh,
        quad,
        z_max,
        z_elem,
        z_stretch,
    )
end

function get_spaces_restart(Y)
    center_space = axes(Y.c)
    face_space = axes(Y.f)
    return (; center_space, face_space)
end

function get_state_restart(comms_ctx)
    @assert haskey(ENV, "RESTART_FILE")
    reader = InputOutput.HDF5Reader(ENV["RESTART_FILE"], comms_ctx)
    Y = InputOutput.read_field(reader, "Y")
    t_start = InputOutput.HDF5.read_attribute(reader.file, "time")
    return (Y, t_start)
end

function get_initial_condition(parsed_args)
    if parsed_args["initial_condition"] in [
        "DryBaroclinicWave",
        "MoistBaroclinicWave",
        "DecayingProfile",
        "MoistBaroclinicWaveWithEDMF",
        "MoistAdiabaticProfileEDMFX",
    ]
        return getproperty(ICs, Symbol(parsed_args["initial_condition"]))(
            parsed_args["perturb_initstate"],
        )
    elseif parsed_args["initial_condition"] in [
        "Nieuwstadt",
        "GABLS",
        "GATE_III",
        "Soares",
        "Bomex",
        "LifeCycleTan2018",
        "ARM_SGP",
        "DYCOMS_RF01",
        "DYCOMS_RF02",
        "Rico",
        "TRMM_LBA",
        "SimplePlume",
    ]
        return getproperty(ICs, Symbol(parsed_args["initial_condition"]))(
            parsed_args["prognostic_tke"],
        )
    elseif parsed_args["initial_condition"] in [
        "IsothermalProfile",
        "AgnesiHProfile",
        "DryDensityCurrentProfile",
        "RisingThermalBubbleProfile",
        "ScharProfile",
        "PrecipitatingColumn",
    ]
        return getproperty(ICs, Symbol(parsed_args["initial_condition"]))()
    else
        error(
            "Unknown `initial_condition`: $(parsed_args["initial_condition"])",
        )
    end
end

function get_surface_setup(parsed_args)
    return getproperty(SurfaceConditions, Symbol(parsed_args["surface_setup"]))()
end

jac_kwargs(ode_algo, Y, p, parsed_args) =
    if ode_algo isa CTS.IMEXAlgorithm
        use_exact_jacobian = parsed_args["use_exact_jacobian"]
        always_update_exact_jacobian =
            parsed_args["n_steps_update_exact_jacobian"] == 0
        approximate_solve_iters = parsed_args["approximate_linear_solve_iters"]
        jacobian_algorithm = if parsed_args["debug_approximate_jacobian"]
            DebugJacobian(
                ApproxJacobian(; approximate_solve_iters);
                use_exact_jacobian,
                always_update_exact_jacobian,
            )
        else
            use_exact_jacobian ?
            ExactJacobian(; always_update_exact_jacobian) :
            ApproxJacobian(; approximate_solve_iters)
        end
        A = ImplicitEquationJacobian(jacobian_algorithm, Y, p)
        @info "Jacobian algorithm: $(dump_string(jacobian_algorithm))"
        (; jac_prototype = A, Wfact = update_jacobian!)
    else
        (;)
    end

#=
    ode_configuration(Y, parsed_args)

Returns the ode algorithm
=#
function ode_configuration(::Type{FT}, parsed_args) where {FT}
    ode_name = parsed_args["ode_algo"]
    ode_algo_name = getproperty(CTS, Symbol(ode_name))
    @info "Using ODE config: `$ode_algo_name`"

    if ode_algo_name <: CTS.ERKAlgorithmName
        return CTS.ExplicitAlgorithm(ode_algo_name())
    else
        @assert ode_algo_name <: CTS.IMEXARKAlgorithmName
        newtons_method = CTS.NewtonsMethod(;
            max_iters = parsed_args["max_newton_iters_ode"],
            krylov_method = if parsed_args["use_krylov_method"]
                CTS.KrylovMethod(;
                    jacobian_free_jvp = CTS.ForwardDiffJVP(;
                        step_adjustment = FT(
                            parsed_args["jvp_step_adjustment"],
                        ),
                    ),
                    forcing_term = if parsed_args["use_dynamic_krylov_rtol"]
                        α = FT(parsed_args["eisenstat_walker_forcing_alpha"])
                        CTS.EisenstatWalkerForcing(; α)
                    else
                        CTS.ConstantForcing(FT(parsed_args["krylov_rtol"]))
                    end,
                )
            else
                nothing
            end,
            convergence_checker = if parsed_args["use_newton_rtol"]
                norm_condition = CTS.MaximumRelativeError(
                    FT(parsed_args["newton_rtol"]),
                )
                CTS.ConvergenceChecker(; norm_condition)
            else
                nothing
            end,
        )
        return CTS.IMEXAlgorithm(ode_algo_name(), newtons_method)
    end
end

thermo_state_type(::DryModel, ::Type{FT}) where {FT} = TD.PhaseDry{FT}
thermo_state_type(::EquilMoistModel, ::Type{FT}) where {FT} = TD.PhaseEquil{FT}
thermo_state_type(::NonEquilMoistModel, ::Type{FT}) where {FT} =
    TD.PhaseNonEquil{FT}

function get_sim_info(config::AtmosConfig)
    (; parsed_args) = config
    FT = eltype(config)

    job_id = if isnothing(parsed_args["job_id"])
        job_id_from_config(parsed_args)
    else
        parsed_args["job_id"]
    end
    default_output = haskey(ENV, "CI") ? job_id : joinpath("output", job_id)
    out_dir = parsed_args["output_dir"]
    output_dir = isnothing(out_dir) ? default_output : out_dir
    mkpath(output_dir)

    sim = (;
        output_dir,
        restart = haskey(ENV, "RESTART_FILE"),
        job_id,
        dt = FT(time_to_seconds(parsed_args["dt"])),
        start_date = DateTime(parsed_args["start_date"], dateformat"yyyymmdd"),
        t_end = FT(time_to_seconds(parsed_args["t_end"])),
    )
    n_steps = floor(Int, sim.t_end / sim.dt)
    @info(
        "Time info:",
        dt = parsed_args["dt"],
        t_end = parsed_args["t_end"],
        floor_n_steps = n_steps,
    )

    return sim
end

function get_diagnostics(parsed_args, atmos_model, cspace)

    # We either get the diagnostics section in the YAML file, or we return an empty list
    # (which will result in an empty list being created by the map below)
    yaml_diagnostics = get(parsed_args, "diagnostics", [])

    # ALLOWED_REDUCTIONS is the collection of reductions we support. The keys are the
    # strings that have to be provided in the YAML file. The values are tuples with the
    # function that has to be passed to reduction_time_func and the one that has to passed
    # to pre_output_hook!

    # We make "nothing" a string so that we can accept also the word "nothing", in addition
    # to the absence of the value
    #
    # NOTE: Everything has to be lowercase in ALLOWED_REDUCTIONS (so that we can match
    # "max" and "Max")
    ALLOWED_REDUCTIONS = Dict(
        "nothing" => (nothing, nothing), # nothing is: just dump the variable
        "max" => (max, nothing),
        "min" => (min, nothing),
        "average" => ((+), CAD.average_pre_output_hook!),
    )

    hdf5_writer = CAD.HDF5Writer()

    if !isnothing(parsed_args["netcdf_interpolation_num_points"])
        num_netcdf_points =
            tuple(parsed_args["netcdf_interpolation_num_points"]...)
    else
        # TODO: Once https://github.com/CliMA/ClimaCore.jl/pull/1567 is merged,
        # dispatch over the Grid type
        num_netcdf_points = (180, 90, 50)
    end

    netcdf_writer = CAD.NetCDFWriter(;
        cspace,
        num_points = num_netcdf_points,
        disable_vertical_interpolation = parsed_args["netcdf_output_at_levels"],
    )
    writers = (hdf5_writer, netcdf_writer)

    # The default writer is HDF5
    ALLOWED_WRITERS = Dict(
        "nothing" => netcdf_writer,
        "h5" => hdf5_writer,
        "hdf5" => hdf5_writer,
        "nc" => netcdf_writer,
        "netcdf" => netcdf_writer,
    )

    diagnostics_ragged = map(yaml_diagnostics) do yaml_diag
        short_names = yaml_diag["short_name"]
        output_name = get(yaml_diag, "output_name", nothing)

        if short_names isa Vector
            isnothing(output_name) || error(
                "Diagnostics: cannot have multiple short_names while specifying output_name",
            )
        else
            short_names = [short_names]
        end

        ret_value = map(short_names) do short_name
            # Return "nothing" if "reduction_time" is not in the YAML block
            #
            # We also normalize everything to lowercase, so that can accept "max" but
            # also "Max"
            reduction_time_yaml =
                lowercase(get(yaml_diag, "reduction_time", "nothing"))

            if !haskey(ALLOWED_REDUCTIONS, reduction_time_yaml)
                error("reduction $reduction_time_yaml not implemented")
            else
                reduction_time_func, pre_output_hook! =
                    ALLOWED_REDUCTIONS[reduction_time_yaml]
            end

            writer_ext = lowercase(get(yaml_diag, "writer", "nothing"))

            if !haskey(ALLOWED_WRITERS, writer_ext)
                error("writer $writer_ext not implemented")
            else
                writer = ALLOWED_WRITERS[writer_ext]
            end

            haskey(yaml_diag, "period") ||
                error("period keyword required for diagnostics")

            period_seconds = time_to_seconds(yaml_diag["period"])

            if isnothing(output_name)
                output_short_name = CAD.descriptive_short_name(
                    CAD.get_diagnostic_variable(short_name),
                    period_seconds,
                    reduction_time_func,
                    pre_output_hook!,
                )
            end

            if isnothing(reduction_time_func)
                compute_every = period_seconds
            else
                compute_every = :timestep
            end

            return CAD.ScheduledDiagnosticTime(
                variable = CAD.get_diagnostic_variable(short_name),
                output_every = period_seconds,
                compute_every = compute_every,
                reduction_time_func = reduction_time_func,
                pre_output_hook! = pre_output_hook!,
                output_writer = writer,
                output_short_name = output_short_name,
            )
        end
        return ret_value
    end

    # Flatten the array of arrays of diagnostics
    diagnostics = vcat(diagnostics_ragged...)

    if parsed_args["output_default_diagnostics"]
        t_end = time_to_seconds(parsed_args["t_end"])
        return [
            CAD.default_diagnostics(
                atmos_model,
                t_end;
                output_writer = netcdf_writer,
            )...,
            diagnostics...,
        ],
        writers
    else
        return collect(diagnostics), writers
    end
end

function args_integrator(parsed_args, Y, p, tspan, ode_algo, callback)
    (; dt) = p
    dt_save_to_sol = time_to_seconds(parsed_args["dt_save_to_sol"])

    s = @timed_str begin
        func = if parsed_args["split_ode"]
            implicit_func = SciMLBase.ODEFunction(
                implicit_tendency!;
                jac_kwargs(ode_algo, Y, p, parsed_args)...,
                tgrad = (∂Y∂t, Y, p, t) -> (∂Y∂t .= 0),
            )
            CTS.ClimaODEFunction(;
                T_exp_T_lim! = remaining_tendency!,
                T_imp! = implicit_func,
                # Can we just pass implicit_tendency! and jac_prototype etc.?
                lim! = limiters_func!,
                dss!,
                post_explicit! = set_precomputed_quantities!,
                post_implicit! = set_precomputed_quantities!,
            )
        else
            remaining_tendency! # should be total_tendency!
        end
    end
    @info "Define ode function: $s"
    problem = SciMLBase.ODEProblem(func, Y, tspan, p)
    saveat = if dt_save_to_sol == Inf
        tspan[2]
    elseif tspan[2] % dt_save_to_sol == 0
        dt_save_to_sol
    else
        [tspan[1]:dt_save_to_sol:tspan[2]..., tspan[2]]
    end # ensure that tspan[2] is always saved
    @info "dt_save_to_sol: $dt_save_to_sol, length(saveat): $(length(saveat))"
    args = (problem, ode_algo)
    kwargs = (;
        saveat,
        callback,
        dt,
        kwargshandle = KeywordArgSilent, # allow custom kwargs
        adjustfinal = true,
    )
    return (args, kwargs)
end

import ClimaComms, Logging, NVTX
function get_comms_context(parsed_args)
    device = if parsed_args["device"] == "auto"
        ClimaComms.device()
    elseif parsed_args["device"] == "CUDADevice"
        ClimaComms.CUDADevice()
    elseif parsed_args["device"] == "CPUMultiThreaded" || Threads.nthreads() > 1
        ClimaComms.CPUMultiThreaded()
    else
        ClimaComms.CPUSingleThreaded()
    end
    comms_ctx = ClimaComms.context(device)
    ClimaComms.init(comms_ctx)
    if ClimaComms.iamroot(comms_ctx)
        Logging.global_logger(Logging.ConsoleLogger(stderr, Logging.Info))
    else
        Logging.global_logger(Logging.NullLogger())
    end
    @info "Running on $(nameof(typeof(device)))."
    if comms_ctx isa ClimaComms.SingletonCommsContext
        @info "Setting up single-process ClimaAtmos run"
    else
        @info "Setting up distributed ClimaAtmos run" nprocs =
            ClimaComms.nprocs(comms_ctx)
    end
    if NVTX.isactive()
        # makes output on buildkite a bit nicer
        if ClimaComms.iamroot(comms_ctx)
            atexit() do
                println("--- Saving profiler information")
            end
        end
    end

    return comms_ctx
end

function get_simulation(config::AtmosConfig)
    params = create_parameter_set(config)
    atmos = get_atmos(config, params)

    sim_info = get_sim_info(config)
    job_id = sim_info.job_id
    output_dir = sim_info.output_dir

    CP.log_parameter_information(
        config.toml_dict,
        joinpath(output_dir, "$(job_id)_parameters.toml"),
        strict = true,
    )
    YAML.write_file(joinpath(output_dir, "$job_id.yml"), config.parsed_args)

    if sim_info.restart
        s = @timed_str begin
            (Y, t_start) = get_state_restart(config.comms_ctx)
            spaces = get_spaces_restart(Y)
            @warn "Progress estimates do not support restarted simulations"
        end
        @info "Allocating Y: $s"
    else
        spaces = get_spaces(config.parsed_args, params, config.comms_ctx)
    end

    initial_condition = get_initial_condition(config.parsed_args)
    surface_setup = get_surface_setup(config.parsed_args)

    if !sim_info.restart
        s = @timed_str begin
            Y = ICs.atmos_state(
                initial_condition(params),
                atmos,
                spaces.center_space,
                spaces.face_space,
            )
            t_start = Spaces.undertype(axes(Y.c))(0)
        end
        @info "Allocating Y: $s"
    end

    s = @timed_str begin
        p = build_cache(Y, atmos, params, surface_setup, sim_info)
    end
    @info "Allocating cache (p): $s"

    if config.parsed_args["discrete_hydrostatic_balance"]
        set_discrete_hydrostatic_balanced_state!(Y, p)
    end

    FT = Spaces.undertype(axes(Y.c))
    s = @timed_str begin
        ode_algo = ode_configuration(FT, config.parsed_args)
    end
    @info "ode_configuration: $s"

    s = @timed_str begin
        callback = get_callbacks(config, sim_info, atmos, params, Y, p, t_start)
    end
    @info "get_callbacks: $s"

    # Initialize diagnostics
    s = @timed_str begin
        diagnostics, writers =
            get_diagnostics(config.parsed_args, atmos, Spaces.axes(Y.c))
    end
    @info "initializing diagnostics: $s"

    length(diagnostics) > 0 && @info "Computing diagnostics:"

    for writer in writers
        writer_str = nameof(typeof(writer))
        diags_with_writer =
            filter((x) -> getproperty(x, :output_writer) == writer, diagnostics)
        diags_outputs = [
            getproperty(diag, :output_short_name) for diag in diags_with_writer
        ]
        @info "$writer_str: $diags_outputs"
    end

    # First, we convert all the ScheduledDiagnosticTime into ScheduledDiagnosticIteration,
    # ensuring that there is consistency in the timestep and the periods and translating
    # those periods that depended on the timestep
    diagnostics_iterations =
        [CAD.ScheduledDiagnosticIterations(d, sim_info.dt) for d in diagnostics]

    # For diagnostics that perform reductions, the storage is used for the values computed
    # at each call. Reductions also save the accumulated value in diagnostic_accumulators.
    diagnostic_storage = Dict()
    diagnostic_accumulators = Dict()
    diagnostic_counters = Dict()

    s = @timed_str begin
        diagnostics_functions = CAD.get_callbacks_from_diagnostics(
            diagnostics_iterations,
            diagnostic_storage,
            diagnostic_accumulators,
            diagnostic_counters,
            output_dir,
        )
    end
    @info "Prepared diagnostic callbacks: $s"

    # It would be nice to just pass the callbacks to the integrator. However, this leads to
    # a significant increase in compile time for reasons that are not known. For this
    # reason, we only add one callback to the integrator, and this function takes care of
    # executing the other callbacks. This single function is orchestrate_diagnostics

    function orchestrate_diagnostics(integrator)
        for d in diagnostics_functions
            if d.cbf.n > 0 && integrator.step % d.cbf.n == 0
                d.f!(integrator)
            end
        end
    end

    diagnostic_callbacks =
        call_every_n_steps(orchestrate_diagnostics, skip_first = true)

    # The generic constructor for SciMLBase.CallbackSet has to split callbacks into discrete
    # and continuous. This is not hard, but can introduce significant latency. However, all
    # the callbacks in ClimaAtmos are discrete_callbacks, so we directly pass this
    # information to the constructor
    continuous_callbacks = tuple()
    discrete_callbacks = (callback..., diagnostic_callbacks)

    s = @timed_str begin
        all_callbacks =
            SciMLBase.CallbackSet(continuous_callbacks, discrete_callbacks)
    end
    @info "Prepared SciMLBase.CallbackSet callbacks: $s"
    steps_cycle_non_diag = n_steps_per_cycle_per_cb(all_callbacks, sim_info.dt)
    steps_cycle_diag =
        n_steps_per_cycle_per_cb_diagnostic(diagnostics_functions)
    steps_cycle = lcm([steps_cycle_non_diag..., steps_cycle_diag...])
    @info "n_steps_per_cycle_per_cb (non diagnostics): $steps_cycle_non_diag"
    @info "n_steps_per_cycle_per_cb_diagnostic: $steps_cycle_diag"
    @info "n_steps_per_cycle (non diagnostics): $steps_cycle"

    tspan = (t_start, sim_info.t_end)
    s = @timed_str begin
        integrator_args, integrator_kwargs = args_integrator(
            config.parsed_args,
            Y,
            p,
            tspan,
            ode_algo,
            all_callbacks,
        )
    end

    s = @timed_str begin
        integrator = SciMLBase.init(integrator_args...; integrator_kwargs...)
    end
    @info "init integrator: $s"
    reset_graceful_exit(output_dir)

    s = @timed_str begin
        for diag in diagnostics_iterations
            variable = diag.variable
            try
                # The first time we call compute! we use its return value. All
                # the subsequent times (in the callbacks), we will write the
                # result in place
                diagnostic_storage[diag] =
                    variable.compute!(nothing, Y, p, t_start)
                diagnostic_counters[diag] = 1
                # If it is not a reduction, call the output writer as well
                if isnothing(diag.reduction_time_func)
                    CAD.write_field!(
                        diag.output_writer,
                        diagnostic_storage[diag],
                        diag,
                        Y,
                        p,
                        t_start,
                        output_dir,
                    )
                else
                    # Add to the accumulator

                    # We use similar + .= instead of copy because CUDA 5.2 does
                    # not supported nested wrappers with view(reshape(view))
                    # objects. See discussion in
                    # https://github.com/CliMA/ClimaAtmos.jl/pull/2579 and
                    # https://github.com/JuliaGPU/Adapt.jl/issues/21
                    diagnostic_accumulators[diag] =
                        similar(diagnostic_storage[diag])
                    diagnostic_accumulators[diag] .=
                        diagnostic_storage[diag]
                end
            catch e
                error("Could not compute diagnostic $(variable.long_name): $e")
            end
        end
    end
    @info "Init diagnostics: $s"

    if config.parsed_args["warn_allocations_diagnostics"]
        for diag in diagnostics_iterations
            # We write over the storage space we have already prepared (and filled) before
            allocs = @allocated diag.variable.compute!(
                diagnostic_storage[diag],
                Y,
                p,
                t_start,
            )
            if allocs > 10 * 1024
                @warn "Diagnostics $(diag.output_short_name) allocates $allocs bytes"
            end
        end
    end

    return AtmosSimulation(
        job_id,
        output_dir,
        sim_info.start_date,
        sim_info.t_end,
        writers,
        integrator,
    )
end

# Compatibility with old get_integrator
function get_integrator(config::AtmosConfig)
    Base.depwarn(
        "get_integrator is deprecated, use get_simulation instead",
        :get_integrator,
    )
    return get_simulation(config).integrator
end
