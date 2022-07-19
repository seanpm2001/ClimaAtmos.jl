# Tempest benchmark

# setup dir location
HOME = "/home/elencz/"
TEMPEST_HOME="$HOME/tempestmodel"
OUTPUT_DIR="/central/groups/esm/lenka/tempest/test/"

# clone Tempest
run(`git clone https://github.com/paullric/tempestmodel`)
run(`cd $TEMPEST_HOME`)

# delete unnecessary breaking line
run(`sed '10d' ${TEMPEST_HOME}/src/base/netcdf.cpp`)

# load necessary modules
run(`module load julia/1.7.2 hdf5/1.10.1 netcdf-c/4.6.1 openmpi/4.0.1`)

# configure and build Tempest using the Makefile
run(`make`)

# run Tempest

# 1. TODO: unthreaded (nDistributedPatches = nPatchCount = 1 not possible > ("nMinPatchCount must be >= 6"))
# [Paul] - should compare with ntasks6?

# 2. threaded (min = 6 so that 1 per panel); TODO: time the time loop
run(`mpirun -np 1 $TEMPEST_HOME/test/nonhydro_sphere/HeldSuarezTest
        --resolution 12 --order 4 --levels 45 --ztop 30000 --endtime 1d --dt 400s --outputtime 1200s
        --nu 7e15 --nud 7e15 --nuv 7e15 --output_x 90 --output_y 45 
        --output_dir $OUTPUT_DIR`)



# 3. TODO: MPI (np = 8) threaded (no patches = no comms)

#=
# Appendix A: CLI args
Parameters:
  --output_none <bool> [false] 
  --output_dir <string> ["/central/groups/esm/lenka/tempest/test/"] 
  --output_prefix <string> ["out"] 
  --restart_file <string> [""] 
  --perturb_restart <bool> [false] 
  --output_perfile <integer> [-1] 
  --output_restart_dt <dtime> [] 
  --output_x <integer> [90] 
  --output_y <integer> [45] 
  --output_z <integer> [0] 
  --output_vort <bool> [false] 
  --output_div <bool> [false] 
  --output_temp <bool> [false] 
  --output_ps <bool> [false] 
  --output_Ri <bool> [false] 
  --norefstate <bool> [false] 
  --notracers <bool> [false] 
  --nohypervis <bool> [false] 
  --hypervisorder <integer> [4] 
  --nu <double> [3.150000e+14] 
  --nud <double> [3.150000e+14] 
  --nuv <double> [3.150000e+14] 
  --inud <double> [0.000000] 
  --explicitvertical <bool> [false] 
  --vstagger <string> ["LOR"] (LEV | INT | LOR | CPH)
  --vdisc <string> ["FE"] (FE | FV)
  --vmassfluxlevels <bool> [false] 
  --vstretch <string> ["uniform"] 
  --vhypervisorder <integer> [0] 
  --timescheme <string> ["strang"] 
  --hmethod <string> ["V1"] (V1 | V2 | SPEX)
  --vmethod <string> ["V1"] (V1 | V2 | SCHUR | NONE)
  --resolution <integer> [4] 
  --levels <integer> [10] 
  --outputtime <dtime> [400s] 
  --dt <dtime> [400s] 
  --endtime <time> [10d] 
  --order <integer> [4] 
  --vertorder <integer> [1] 
  --ztop <double> [30000.000000] 
  --zh <double> [30000.000000] 
  --tau0 <double> [25.000000] 
  =#