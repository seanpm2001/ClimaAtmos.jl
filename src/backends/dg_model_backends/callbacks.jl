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

#=
Base.@kwdef struct AveragedState{𝒜, ℬ, 𝒞, 𝒟} <: AbstractCallback
    iteration::𝒜
    filepath::ℬ
    overwrite::𝒞 = true
    start_iteration::𝒟
end
=#

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
            new_file["times"] = time + 1
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