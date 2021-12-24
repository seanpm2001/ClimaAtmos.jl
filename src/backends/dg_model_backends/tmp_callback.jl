# goal: write a function that computes all second order statistics
# u, v, w, p, T
# uu, uv, uw, uT, vv, vw, vT, ww, wT, ww, TT
# check: how to extract geopotential from dg object
# step 1, just output velocities to start
# step 2, output Temperature and Pressure 
# step 3, output second momemnts 

# design decisions: store everything in one big state?
# design decisions: netcdf output instead?

Base.@kwdef struct DefaultDiagnostics{𝒜, ℬ, 𝒞, 𝒟} <: AbstractCallback
    iteration::𝒜
    filepath::ℬ
    overwrite::𝒞 = true
    start_iteration::𝒟 = 0
end

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

function create_callback(output::DefaultDiagnostics, simulation::Simulation{<:DiscontinuousGalerkinBackend}, odesolver)
    # Initialize output
    output.overwrite &&
        isfile(output.filepath) &&
        rm(output.filepath; force = output.overwrite)

    Q = simulation.state
    geopotential = simulation.rhs[1].state_auxiliary.Φ # capital \Phi

    mpicomm = MPI.COMM_WORLD
    iteration = output.iteration

    steps = ClimateMachine.ODESolvers.getsteps(odesolver)
    time = ClimateMachine.ODESolvers.gettime(odesolver)

    file = jldopen(output.filepath, "a+")
    # JLD2.Group(file, "state")
    # JLD2.Group(file, "time")
    
    moment_1, moment_1_names = get_state(Q, geopotential, parameters)
    moment_2, moment_2_names = get_second_moments(moment_1, moment_1_names)
    if output.start_iteration <= 0
        file["moment_1"] = moment_1
        file["moment_2"] = moment_2  
        file["times"] = 1
    else
        file["moment_1"] = moment_1 .* 0.0
        file["moment_2"] = moment_2 .* 0.0
        file["times"] = 0
    end
    file["moment_1_names"] = moment_1_names
    file["moment_2_names"] = moment_2_names

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
            old_moment_1 = copy(file["moment_1"])
            old_moment_2 = copy(file["moment_2"])
            oldt = file["times"]
            close(file)
            rm(output.filepath)
            moment_1, moment_1_names = get_state(Q, geopotential, parameters)
            moment_2, moment_2_names = get_second_moments(moment_1, moment_1_names)
            # put data in new file as a part of running average
            new_file = jldopen(output.filepath, "a+")
            new_file["moment_1"] = moment_1 + old_moment_1
            new_file["moment_2"] = moment_2 + old_moment_2
            
            new_file["moment_1_names"] = moment_1_names
            new_file["moment_2_names"] = moment_2_names
            new_file["times"] = oldt + 1
            close(new_file)
        end
        return nothing
    end

    return jldcallback
end

# Lat lon callback 

Base.@kwdef struct LatLonDiagnostics{𝒜, ℬ, 𝒞, 𝒟, ℰ} <: AbstractCallback
    iteration::𝒜
    filepath::ℬ
    overwrite::𝒞 = true
    start_iteration::𝒟 = 0
    latitude::ℰ
    longitude::ℰ
    radius::ℰ
end

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





