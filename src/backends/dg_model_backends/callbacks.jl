function create_callback(::Info, simulation::Simulation{<:DiscontinuousGalerkinBackend}, odesolver)
    Q = simulation.state
    timeend = simulation.timestepper.finish
    mpicomm = MPI.COMM_WORLD

    starttime = Ref(now())
    cbinfo = ClimateMachine.GenericCallbacks.EveryXWallTimeSeconds(
        60,
        mpicomm,
    ) do (s = false)
        if s
            starttime[] = now()
        else
            energy = norm(Q)
            @info @sprintf(
                """Update
                simtime = %8.4f / %8.4f
                runtime = %s
                norm(Q) = %.16e""",
                ClimateMachine.ODESolvers.gettime(odesolver),
                timeend,
                Dates.format(
                    convert(Dates.DateTime, Dates.now() - starttime[]),
                    Dates.dateformat"HH:MM:SS",
                ),
                energy
            )

            if isnan(energy)
                error("NaNs")
            end
        end
    end

    return cbinfo
end

function create_callback(::CFL, simulation::Simulation{<:DiscontinuousGalerkinBackend}, odesolver)
    Q = simulation.state
    # timeend = simulation.time.finish
    # mpicomm = MPI.COMM_WORLD
    # starttime = Ref(now())
    cbcfl = EveryXSimulationSteps(100) do
            simtime = gettime(odesolver)

            @views begin
                ρ = Array(Q.data[:, 1, :])
                ρu = Array(Q.data[:, 2, :])
                ρv = Array(Q.data[:, 3, :])
                ρw = Array(Q.data[:, 4, :])
            end

            u = ρu ./ ρ
            v = ρv ./ ρ
            w = ρw ./ ρ

            # TODO! transform onto sphere

            ue = extrema(u)
            ve = extrema(v)
            we = extrema(w)

            @info @sprintf """CFL
                    simtime = %.16e
                    u = (%.4e, %.4e)
                    v = (%.4e, %.4e)
                    w = (%.4e, %.4e)
                    """ simtime ue... ve... we...
        end

    return cbcfl
end

function create_callback(callback::StateCheck, simulation::Simulation{<:DiscontinuousGalerkinBackend}, _...)
    sim_length = simulation.time.finish - simulation.time.start
    timestep = simulation.timestepper.timestep
    nChecks = callback.number_of_checks

    nt_freq = floor(Int, sim_length / timestep / nChecks)

    cbcs_dg = ClimateMachine.StateCheck.sccreate(
        [(simulation.state, "state")],
        nt_freq,
    )

    return cbcs_dg
end

function create_callback(output::JLD2State, simulation::Simulation{<:DiscontinuousGalerkinBackend}, odesolver)
    # Initialize output
    output.overwrite &&
        isfile(output.filepath) &&
        rm(output.filepath; force = output.overwrite)

    Q = simulation.state
    mpicomm = MPI.COMM_WORLD
    iteration = output.iteration

    steps = ClimateMachine.ODESolvers.getsteps(odesolver)
    time = ClimateMachine.ODESolvers.gettime(odesolver)

    file = jldopen(output.filepath, "a+")
    JLD2.Group(file, "state")
    JLD2.Group(file, "time")
    file["state"][string(steps)] = Array(Q)
    file["time"][string(steps)] = time
    close(file)


    jldcallback = ClimateMachine.GenericCallbacks.EveryXSimulationSteps(
        iteration,
    ) do (s = false)
        steps = ClimateMachine.ODESolvers.getsteps(odesolver)
        time = ClimateMachine.ODESolvers.gettime(odesolver)
        @info steps, time
        file = jldopen(output.filepath, "a+")
        file["state"][string(steps)] = Array(Q)
        file["time"][string(steps)] = time
        close(file)
        return nothing
    end

    return jldcallback
end

function create_callback(output::AveragedState, simulation::Simulation{<:DiscontinuousGalerkinBackend}, odesolver)
    # Initialize output
    output.overwrite &&
        isfile(output.filepath) &&
        rm(output.filepath; force = output.overwrite)

    Q = simulation.state
    mpicomm = MPI.COMM_WORLD
    iteration = output.iteration

    steps = ClimateMachine.ODESolvers.getsteps(odesolver)
    time = ClimateMachine.ODESolvers.gettime(odesolver)

    file = jldopen(output.filepath, "a+")
    # JLD2.Group(file, "state")
    # JLD2.Group(file, "time")
    
    if output.start_iteration <= 0
        file["state"] = Array(Q)
        # file["times"] = 0
    else
        file["state"] = Array(Q) .* 0.0

        # file["times"] = []
    end
    file["times"] = 1
    close(file)


    jldcallback = ClimateMachine.GenericCallbacks.EveryXSimulationSteps(
        iteration,
    ) do (s = false)
        steps = ClimateMachine.ODESolvers.getsteps(odesolver)
        time = ClimateMachine.ODESolvers.gettime(odesolver)
        @info steps, time
        # a bit hacky but gets the job done. removes old file and creates new one
        if steps > output.start_iteration
            @info "accumulating average"
            # open old file and grab data
            file = jldopen(output.filepath, "a+")
            oldQ = copy(file["state"])
            oldt = file["times"]
            close(file)
            rm(output.filepath)
            # put data in new file as a part of running average
            new_file = jldopen(output.filepath, "a+")
            new_file["state"] = Array(Q) + oldQ
            new_file["times"] = oldt + 1
            close(new_file)
        end
        return nothing
    end

    return jldcallback
end

function create_callback(output::VTKState, simulation::Simulation{<:DiscontinuousGalerkinBackend}, odesolver)
    # Initialize output
    output.overwrite &&
        isfile(output.filepath) &&
        rm(output.filepath; force = output.overwrite)
    mkpath(output.filepath)

    state = simulation.state
    if simulation.rhs isa Tuple
        if simulation.rhs[1] isa AbstractRate 
            model = simulation.rhs[1].model
        else
            model = simulation.rhs[1]
        end
    else
        model = simulation.rhs
    end
    # model = (simulation.rhs isa Tuple) ? simulation.rhs[1] : simulation.rhs 

    function do_output(counter, model, state)
        mpicomm = MPI.COMM_WORLD
        balance_law = model.balance_law
        aux_state = model.state_auxiliary

        outprefix = @sprintf(
            "%s/mpirank%04d_step%04d",
            output.filepath,
            MPI.Comm_rank(mpicomm),
            counter[1],
        )

        @info "doing VTK output" outprefix

        state_names =
            flattenednames(vars_state(balance_law, Prognostic(), eltype(state)))
        aux_names =
            flattenednames(vars_state(balance_law, Auxiliary(), eltype(state)))

        writevtk(outprefix, state, model, state_names, aux_state, aux_names)

        counter[1] += 1

        return nothing
    end

    do_output(output.counter, model, state)
    cbvtk =
        ClimateMachine.GenericCallbacks.EveryXSimulationSteps(output.iteration) do (
            init = false
        )
            do_output(output.counter, model, state)
            return nothing
        end

    return cbvtk
end

function create_callback(filter::PositivityPreservingCallback, simulation::Simulation{<:DiscontinuousGalerkinBackend}, odesolver)
    Q = simulation.state
    rhs = simulation.rhs
    if rhs isa SpaceDiscretization
        grid = simulation.rhs.grid
    elseif rhs isa Tuple 
        grid = simulation.rhs[1].grid
    else
        println("rhs error => fail to initialize PositivityPreservingCallback")
    end
    grid = simulation.rhs[1].grid
    tmar_filter = EveryXSimulationSteps(1) do
        ClimateMachine.Mesh.Filters.apply!(Q, filter.filterstates, grid, TMARFilter())
        end
    return tmar_filter
end

# helper function 
function update_ref_state!(
    balance_law::BalanceLaw,
    state::Vars,
    aux::Vars,
    t::Real,
)
    eos = balance_law.equation_of_state
    parameters = balance_law.parameters
    ρ = state.ρ
    ρu = state.ρu
    ρe = state.ρe

    aux.ref_state.ρ = ρ
    aux.ref_state.ρu = ρu # @SVector[0.0,0.0,0.0]
    aux.ref_state.ρe = ρe
    aux.ref_state.p = calc_pressure(eos, state, aux, parameters)
end

function create_callback(update_ref::ReferenceStateUpdate, simulation::Simulation{<:DiscontinuousGalerkinBackend}, odesolver)
    Q = simulation.state
    step =  update_ref.recompute
    dg = simulation.rhs[2]
    balance_law = dg.balance_law

    relinearize = EveryXSimulationSteps(step) do       
        t = gettime(odesolver)
        
        update_auxiliary_state!(update_ref_state!, dg, balance_law, Q, t)

        α = odesolver.dt * odesolver.RKA_implicit[2, 2]
        # hack
        be_solver = odesolver.implicit_solvers[odesolver.RKA_implicit[2, 2]][1]
        update_backward_Euler_solver!(be_solver, Q, α)
        nothing
    end
    return relinearize
end

# lat lon callback need some helper functions

function get_state(ρQ, geopotential, parameters)
    A_ρQ = Array(ρQ)
    size_ρQ = size(A_ρQ)
    Q = zeros(size_ρQ[1], 7, size_ρQ[3])

    # define indices
    _ρ  = 1 #density
    _ρu = 2 #x-velocity
    _ρv = 3 #y-velocity
    _ρw = 4 #z-velocity
    _ρe = 5 # total energy density
    _p  = 6 # pressure
    _T  = 7 # temperature

    # grab density and total energy
    Q[:,_ρ,:]  .= A_ρQ[:,_ρ,:] 
    Q[:,_ρe,:] .= A_ρQ[:,_ρe,:] 

    # grab velocities 
    Q[:,_ρu,:]  .= A_ρQ[:,_ρu,:] ./ A_ρQ[:, _ρ,:]
    Q[:,_ρv,:]  .= A_ρQ[:,_ρv,:] ./ A_ρQ[:, _ρ,:]
    Q[:,_ρw,:]  .= A_ρQ[:,_ρw,:] ./ A_ρQ[:, _ρ,:]

    # calculate pressure 
    Q[:,_p,:] .= pressure(A_ρQ, Array(geopotential)[:,1,:], parameters)

    # calculate temperature (kelvin)
    Q[:,_T,:] .= Q[:,_p,:] ./ parameters.R_d ./ A_ρQ[:, _ρ,:]

    # string for name 
    state_names = ["ρ", "u", "v", "w", "e", "p", "T"]

    return Q, state_names
end

function get_second_moments(Q, state_names)

    s_Q = size(Q)
    s_states = size(Q)[2]
    s_M = floor(Int, s_states * (s_states+1) / 2) # second momements
    QQ = zeros(s_Q[1], s_M, s_Q[3])
    
    clk = [1]
    moment_names = []
    for i in 1:s_states 
        for j in i:s_states
            QQ[:,clk[1],:] .= Q[:,i,:] .* Q[:,j,:]
            push!(moment_names, state_names[i] * state_names[j])
            clk .= clk .+ 1
        end
    end

    return QQ, moment_names
end

"""
get second moments of the lat lon version of things
"""
function get_second_moments_ll(Q, state_names)

    s_Q = size(Q)
    s_states = s_Q[end]
    s_M = floor(Int, s_states * (s_states+1) / 2) # second moments
    QQ = zeros(s_Q[1:end-1]..., s_M)
    
    clk = [1]
    moment_names = []
    for i in 1:s_states 
        for j in i:s_states
            QQ[:,:,:, clk[1]] .= Q[:,:,:,i] .* Q[:,:,:,j]
            push!(moment_names, state_names[i] * state_names[j])
            clk .= clk .+ 1
        end
    end

    return QQ, moment_names
end

function pressure(ρQ, geopotential, parameters)
    γ = 1 / (parameters.cv_d / parameters.R_d) + 1
    ρ  = ρQ[:,1,:]
    ρu = ρQ[:,2,:]
    ρv = ρQ[:,3,:]
    ρw = ρQ[:,4,:]
    ρe = ρQ[:,5,:]

    ϕ  = geopotential
    p = (γ-1) .* (ρe - 0.5 * (ρu .^2 + ρv .^2 + ρw .^2) ./ ρ .- ρ .* ϕ)
    return p
end

function create_callback(output::LatLonDiagnostics, simulation::Simulation{<:DiscontinuousGalerkinBackend}, odesolver)
    # Initialize output
    output.overwrite &&
        isfile(output.filepath) &&
        rm(output.filepath; force = output.overwrite)

    Q = simulation.state
    geopotential = simulation.rhs[1].state_auxiliary.Φ # capital \Phi

    # immediately grab state 
    moment_1, moment_1_names = get_state(Q, geopotential, parameters)

    latitude  = output.latitude
    longitude = output.longitude
    raditude  = output.radius

    interpol = InterpolationCubedSphere(simulation, latitude = latitude, longitude = longitude, raditude = raditude)

    mpicomm = MPI.COMM_WORLD
    iteration = output.iteration

    steps = ClimateMachine.ODESolvers.getsteps(odesolver)
    time = ClimateMachine.ODESolvers.gettime(odesolver)

    file = jldopen(output.filepath, "a+")
    
    _ρu, _ρv, _ρw = 2, 3, 4

    istate = ClimateMachine.CUDA.CuArray(similar(Q, interpol.Npl, 7)) # 7 because, ρ, ρu, ρv, ρw, ρe, p, T
    
    # get moment and second moments in spherical coordinates
    interpolate_local!(interpol, ClimateMachine.CUDA.CuArray(moment_1), istate) 
    project_cubed_sphere!(interpol, istate, (_ρu, _ρv, _ρw))
    moment_1_ll = Array(accumulate_interpolated_data(MPI.COMM_WORLD, interpol, istate))
    moment_2_ll, moment_2_names = get_second_moments_ll(moment_1_ll, moment_1_names)

    # save
    if output.start_iteration <= 0
        file["moment_1"] = moment_1_ll
        file["moment_2"] = moment_2_ll  
        file["times"] = 1
    else
        file["moment_1"] = moment_1_ll .* 0.0
        file["moment_2"] = moment_2_ll .* 0.0
        file["times"] = 0
    end
    file["moment_1_names"] = moment_1_names
    file["moment_2_names"] = moment_2_names

    JLD2.Group(file, "grid")
    file["grid"]["latitude"] = latitude
    file["grid"]["longitude"] = longitude
    file["grid"]["radius"] = raditude

    close(file)


    jldcallback = ClimateMachine.GenericCallbacks.EveryXSimulationSteps(
        iteration,
    ) do (s = false)

        steps = ClimateMachine.ODESolvers.getsteps(odesolver)
        time = ClimateMachine.ODESolvers.gettime(odesolver)
        @info steps, time/86400
        # a bit hacky but gets the job done. removes old file and creates new one
        if steps > output.start_iteration
            @info "accumulating average"
            # open old file and grab data
            file = jldopen(output.filepath, "a+")
            old_moment_1 = copy(file["moment_1"])
            old_moment_2 = copy(file["moment_2"])
            oldt = file["times"]
            close(file)
            rm(output.filepath)
            moment_1, moment_1_names = get_state(Q, geopotential, parameters)
            interpolate_local!(interpol, ClimateMachine.CUDA.CuArray(moment_1), istate) 
            project_cubed_sphere!(interpol, istate, (_ρu, _ρv, _ρw))
            moment_1_ll = Array(accumulate_interpolated_data(MPI.COMM_WORLD, interpol, istate))
            moment_2_ll, moment_2_names = get_second_moments_ll(moment_1_ll, moment_1_names)
            # put data in new file as a part of running average
            new_file = jldopen(output.filepath, "a+")
            new_file["moment_1"] = moment_1_ll + old_moment_1
            new_file["moment_2"] = moment_2_ll + old_moment_2
            
            new_file["moment_1_names"] = moment_1_names
            new_file["moment_2_names"] = moment_2_names
            new_file["times"] = oldt + 1

            JLD2.Group(new_file, "grid")
            new_file["grid"]["latitude"] = latitude
            new_file["grid"]["longitude"] = longitude
            new_file["grid"]["radius"] = raditude

            close(new_file)
        end

        return nothing
    end

    return jldcallback
end
