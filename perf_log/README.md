# Long-term performace log
We need to be able to answerthe following:
- what is the performance improvement over the last few months, which PRs made the most difference (maybe last week is enough - pre #620 and perf tables thereafter)
- how this compares to Tempest
- what is the log of long-term change in walltime (graph?)
- what speedup we're hoping to get for the planned items

# Tempest Questions
- unthreaded not possible? (only Cartesian) 
- MPI ntasks = 6 ??; Maximum number of patches currently equals communicator size

# Log
- Tempest:
    - walltime
    - mem alloc
- ClimaAtmos: 
    - walltime
    - mem alloc

# TODO
- integrate automated Buildkite logging (with rough estimates pre #620)
- plan with each item including a minimal p.o.c. () 





