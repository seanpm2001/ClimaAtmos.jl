

dg_fs = simulation.rhs[1]
dg_sd = simulation.rhs[2]
test_state = simulation.state
old_state = copy(test_state)
state = copy(test_state)

odesolver = construct_odesolver(IMEXSplitting(), simulation);
be_solver! = odesolver.implicit_solvers[odesolver.RKA_implicit[2, 2]][1];
be_solver!.isadjustable = true

println("----")
println("for Nq⃗=", (hp + 1, hp + 1, vp + 1), ", Kh = ", he, ", Kv =", ve)

#=
tic = Base.time()
for i in 1:10
    α = odesolver.dt * odesolver.RKA_implicit[2, 2]
    # hack
    be_solver = odesolver.implicit_solvers[odesolver.RKA_implicit[2, 2]][1]
    update_backward_Euler_solver!(be_solver, test_state, α)
end
toc = Base.time()
createlu = [(toc - tic) * 1e8]
=#
createlu = @benchmark begin
    α = odesolver.dt * odesolver.RKA_implicit[2, 2]
    # hack
    be_solver = odesolver.implicit_solvers[odesolver.RKA_implicit[2, 2]][1]
    update_backward_Euler_solver!(be_solver, test_state, α)
end

solvetimes = @benchmark begin
    ClimateMachine.ODESolvers.dostep!(test_state, odesolver, nothing, 0.0, nothing, nothing, nothing)
    test_state .= old_state
    odesolver.t = 0.0
end
#okay a bit mislabeled
meanlu = median(createlu.times)
meansolvetimes = median(solvetimes.times)
println("the time in nanoseconds to do the lufactorization of linear system is ", meanlu)
println("the time in nanoseconds to take one timestep is ", meansolvetimes)
println("the ratio of lu construction to solving one step is ", meanlu / meansolvetimes)

rhstimes = @benchmark begin
    dg_fs(state, test_state, nothing, false, false)
    dg_sd(state, test_state, nothing, false, false)
end

ldivtimes = @benchmark begin
    odesolver.rhs_implicit!(test_state, state, nothing, dt, increment=false)
    be_solver!(test_state, state, false, nothing, odesolver.dt * odesolver.RKA_implicit[2, 2])
end
meanldiv = median(ldivtimes.times)
meanrhs = median(rhstimes.times)
println("the time in nanoseconds to solve the factorized linear system is ", meanldiv)
println("the time in nanoseconds to evaluate the explicit rhs is ", meanrhs)
println("the ratio of implicit to explicit is ", meanldiv / meanrhs)
current_timings = [meanlu, meansolvetimes, meanldiv, meanrhs]
println("----")

##



