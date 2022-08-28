#! /bin/bash
set -euxo pipefail
set +x

nargs=$#

if (( nargs != 3 ))
then
    echo "please provide arguments for \"nprocs\", \"resolution\", (high/low) and \"profiling\" (enable/disable)"
fi

nprocs=$1
resolution="$2"
profiling="$3"

job_id="sphere_held_suarez_${resolution}_res_rhoe_$nprocs"

profiling_params="nsys profile --trace=nvtx,mpi,osrt --mpi-impl=openmpi --output=${job_id}/report.%q{NPROCS}.%q{OMPI_COMM_WORLD_RANK}"

if [[ "$resolution" == "low" ]]
then
    sim_params="julia --color=yes --project=examples examples/hybrid/driver.jl --job_id $job_id --forcing held_suarez --enable_threading false --upwinding none --t_end 10days --dt 400secs --z_elem 10 --h_elem 4 --kappa_4 2e17"
else
    sim_params="julia --color=yes --project=examples examples/hybrid/driver.jl --job_id $job_id --forcing held_suarez --enable_threading false --upwinding none --t_end 1days --dt 50secs --z_elem 45 --h_elem 24 --kappa_4 5e14"
fi

if [[ "$profiling" == "enable" ]]
then
    module load cuda/11.3 nsight-systems/2022.2.1
    mpiexec $profiling_params $sim_params
else
    mpiexec $sim_params
fi
