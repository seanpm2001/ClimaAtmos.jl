import ClimaCore
using ClimaCore:
    level,
    Geometry,
    Operators,
    InputOutput
import Thermodynamics as TD

import ClimaCore.Utilities: half
using ClimaCorePlots, Plots

filein = "data/day50.0.hdf5"
reader = InputOutput.HDF5Reader(filein)
Y = InputOutput.read_field(reader, "Y")
 
# atmos config to generate parameters
config = ...
params = ...
thermo_params = CAP.thermodynamics_params(params)

# extract all supersaturated grid cells

# check number of sa iterations needed
ts = TD.thermostate()
