using LinearAlgebra

# positions
x⃗ = [-1, -0.25, 0, 0.25, 1] .* 1.0 

# position operator
x̂ = Diagonal(x⃗)

# construct differentiation matrix / position operator
X⃗ = [x⃗[i]^(j-1) for i in eachindex(x⃗), j in eachindex(x⃗)]
dX⃗ = [(j-1) * x⃗[i]^(j-2) for i in eachindex(x⃗), j in eachindex(x⃗)]
dX⃗[:,1] .= 0.0 # get rid of NaNs
# p̂ * X⃗ = dX⃗, so just invert to find p̂
p̂ = dX⃗ / X⃗

# fake identity
# won't be identiy for polynomial orders ≥ length(x⃗)-1
I_s =  p̂ * x̂ - x̂ * p̂

# and check that its correct to machine precision
check = 3 .+ x⃗ .+ x⃗ .^2 .+ x⃗ .^3 
I_s * check - check
