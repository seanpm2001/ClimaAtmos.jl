module Aleph

using Pkg
export add_climate_machine

"""
add_climate_machine()

# Description
Adds the right branch of ClimateMachine

# Return 
- nothing

"""
function add_climate_machine()
    Pkg.add(url = "https://github.com/CliMA/ClimateMachine.jl.git#tb/refactoring_ans_sphere")
    return nothing
end

end # module