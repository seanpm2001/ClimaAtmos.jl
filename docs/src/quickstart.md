# Quickstart 
I want to run a model setup of my own! 
Atmos has several types of models including single columns and sphere (global) models. Why study SCM cases vs global cases?

General structure of what goes into simulation within the CliMA architecture
 - Physical core (ClimaCore.jl)
 - Numerical methods to handle timestepping (ClimaCore.jl)
 - parameterizations (ClimaAtmos.jl, Microphysics...)
   - parameters (ClimaParams.jl)
... other things

## Structure of ClimaAtmos
Everything is run from a driver file (examples/hybrid/driver.jl) with the details of the configuration specified by the config (.yml) file, e.g. column or sphere (global) run...

# SCM Tutorial 
At minimum we have to tell `julia` which environment to run, the driver, and config file (which typically calls a toml file). 
```bash
julia --project=examples examples/hybrid/driver.jl --config_file config/default_configs/default_config.yml --job_id first_experiment
```
`default_config.yml` sets all the parameters unless overridden by the user. If overridding you need to make sure that the arguments don't conflict otherwise unexpected behavior may occur. For example we could specify a dry convection case and then ask for the 1 moment microphysics scheme. However, if the simulation is dry we won't use microphysics! 

## REPL 
For developing a model it is helpful to debug inside the repl and also allows us to avoid repeated, annoying `julia` precompilation. From command line with your environment of choice (e.g. copy examples and add your own stuff) run:
```bash
julia --project=. --color=yes
```
Then in REPL we can do: 
```julia
using Revise # allows file changes to be incorporated real-time into our REPL session
using YAML
import ClimaAtmos as CA 

config_dict = YAML.load_file("...")
config = CA.AtmosConfig(config_dict; job_id = ...)
simulation = CA.get_simulation(config)

# run simulation
sol_res = CA.solve_atmos!(simulation)
```

To examine what objects are created or avaiable (e.g. examining the cache of what's stored) I find it helpful to unpack the simulation object to examine the data this way 
```julia
(; integrator) = simulation;
(; atmos, params) = integrator.p;
(; p) = integrator;
```
For example we might be interested in looking at the density field. We can access this with `integrator.u.c.ρ`.

## Basic arguments 
time step, resolution, etc - tradeoffs etc 

## Diagnostics 
What are they and how are they useful - point to the documentation. Specify what output we want from the simulation. How do we control reductions etc. Point to where we can add diagnostics

## TOML
Describe what types of things go in the toml (e.g. we want to test different parameters) and is also how we interface with EKP.jl / calibration frameworks.

# SCM "Case study": State of the art TC parameterization - EDMF 
After intro we could run users through an example which could be EDMF or something of your choice. We would include some plots and maybe some general plotting code that uses simple `ClimaAnalysis.jl` tools (e.g. time, lat/lon reductions and windowing) 

Could talk about prognostic vs diagnostic EDMF. What are the differences, specify with the `turbconv` argument in yml.

## Changing external forcing 
GCM driven, ERA5 driven (currently being developed). Typically it is useful to add new structs for each case that we try and add (e.g. GCMDrivenInsolation) and link them to an abstract type (e.g. AbstractInsolation) as we want to reuse some of the utility across structs.

### cfsites and LES data 
Not sure if we want to go this detailed for a tutorial.

## Calibration 
??