using SpecialFunctions # for gamma function reasons
using LinearAlgebra    # for Gauss quadrature
using OffsetArrays     # for notational clarity
using SparseArrays     # for data structures

using GLMakie          # plotting
using LaTeXStrings     # plotting

"""
unimesh1D(xmin, xmax, K)
# Description
    Generates a uniform 1D mesh
# Arguments
    xmin: smallest value of array
    xmax: largest values of array
    K: number of elements in an array
# Return Values: VX, EtoV
    VX: vertex values | an Array of size K+1
    EtoV: element to node connectivity | a Matrix of size Kx2
# Example
xmin = -1
xmax =  1
K    =  4
VX, EtoV = unimesh1D(xmin, xmax, K)
"""
function unimesh1D(xmin, xmax, K)
    VX = collect(0:K) ./ K .* (xmax - xmin) .+ xmin
    EtoV = Int.(ones(K, 2))
    for i = 1:K
        EtoV[i,1] = Int(i)
        EtoV[i,2] = Int(i+1)
    end
    return VX, EtoV
end

# Mathy aliases
const Γ = gamma

# Coefficients in the Jacobi polynomial recurrence relations.
aᴾ(α, β, n) = 2/(2n+α+β) * √(n * (n+α+β) * (n+α) * (n+β) / (2n+α+β-1) / (2n+α+β+1))
bᴾ(α, β, n) = -(α^2 - β^2) / (2n+α+β) / (2n+α+β+2)

# code checked against the matlab code
"""
jacobi(x, α, β, n)
# Description
- Evaluates the jacobi polynomial at the point x
# Arguments
- `x`: point at which you will evaluate the jacobi polynomial
- `α`: first parameter for Jacobi polynomials
- `β`: second parameter for Jacobi polynomials
- `n` : order
# Return
-  `y`: the value of the of the Jacobi polynomial
"""
function jacobi(x, α, β, n::Int)
    Pᵅᵝ = n <= 1 ? OffsetArray(zeros(2), 0:1) : OffsetArray(zeros(n+1), 0:n)
    Pᵅᵝ[0] = √(2.0^-(α+β+1) * Γ(α+β+2) / Γ(α+1) / Γ(β+1))
    Pᵅᵝ[1] = Pᵅᵝ[0]/2 * √((α+β+3) / (α+1) / (β+1)) * ((α+β+2)*x + α - β)
    for n′ in 1:n-1
        Pᵅᵝ[n′+1] = ((x - bᴾ(α,β,n′)) * Pᵅᵝ[n′] - aᴾ(α,β,n′) * Pᵅᵝ[n′-1]) / aᴾ(α, β, n′+1)
    end
    return Pᵅᵝ[n]
end

"""
djacobi(x, α, β, n)
# Description
- Evaluates the derivative of the jacobi polynomial at the point x
# Arguments
- `x`: point at which you will evaluate the derivative of the jacobi polynomial
- `α`: first parameter for Jacobi polynomials
- `β`: second parameter for Jacobi polynomials
- `n` : order
# Return
-  `y`: the derivative of the of the Jacobi polynomial
"""
djacobi(x, α, β, n::Int) = √(n * (n+α+β+1)) * jacobi(x, α+1, β+1, n-1)

"""
vandermonde(x, α, β, N)
# Description
    Return vandermonde matrix of order N at the values x
    Allocates a little bit of memory
# Arguments
-   `x`: points at which to evaluate the Jacobi polynomials
-   `α`: first parameter for Jacobi polynomials
-   `β`: second parameter for Jacobi polynomials
-   `N`: maximum order of Jacobi polynomial to include
# Return Values
-   `v`: vandermonde matrix
# Example
    See LegendreTests.jl
"""
function vandermonde(x, α, β, N)
    # compute first two coefficients
    γ0 = 2^(α + β + 1) * factorial(α) * factorial(β) / ((α + β + 1) * factorial(α + β))
    γ1 = (α + 1) * (β + 1) / (α + β + 3) * γ0

    # create view to assign values
    v = zeros(length(x), N+1)
    v1 = view(v, :, 1)
    @. v1 = 1 / sqrt(γ0)

    # explicitly compute second coefficient
    if N == 0
        return v
    end

    v2 = view(v, :, 2)
    @. v2 = ( (α + β + 2) * x/2 + (α - β)/2) / sqrt(γ1)

    if N == 1
        return v
    end

    aʲ = 2 / (2 + α + β) * sqrt((α+1) * (β+1) / (α + β + 3))

    for i in 3:(N+1)
        # get views for ith, i-1th, and i-2th columns
        vi = view(v, :, i)
        vM1 = view(v, :, i-1)
        vM2 = view(v, :, i-2)

        # compute new a and b values
        h1 = 2 * (i-2) + α + β
        aⁱ = 2 / (h1 + 2) * sqrt((i-1) * (i-1 + α + β) * (i-1 + α) * (i-1 + β) / ((h1 + 1) * (h1 + 3)))
        bⁱ = - (α^2 - β^2) / (h1 * (h1 + 2))

        # compute coefficients for ith column
        @. vi = 1 / aⁱ * (-aʲ * vM2 + (x - bⁱ) * vM1)

        # save a coefficient for next iteration
        aʲ = aⁱ
    end

    return v
end

"""
dvandermonde(x, α, β, N)
# Description
    Return the gradient of the vandermonde matrix of order N at the values x
    Allocates a little bit of memory
# Arguments
-   `x`: points at which to evaluate the Jacobi polynomials
-   `α`: first parameter for Jacobi polynomials
-   `β`: second paramater for Jacobi polynomials
-   `N`: maximum order of Jacobi polynomial to include
# Return Values
-   `vr`: gradient of vandermonde matrix
# Example
    See LegendreTests.jl
"""
function dvandermonde(x, α, β, N)
    # create empty matrix (also handles first set of derivatives)
    vr = zeros(length(x), N+1)

    if N == 0
        return vr
    end

    # set values using vandermonde matrix
    v = vandermonde(x, α+1, β+1, N)
    for i in 1:N
        vi = view(v, :, i)
        vrP1 = view(vr, :, i+1)
        @. vrP1 = sqrt(i * (α + β + i+1)) * vi
    end

    return vr
end

"""
dmatrix(x, α, β, N)
# Description
    Return the differentiation matrix of order N at the values x
    Allocates too much memory
# Arguments
-   `x`: points at which to evaluate the Jacobi polynomials
-   `α`: first parameter for Jacobi polynomials
-   `β`: second paramater for Jacobi polynomials
-   `N`: maximum order of Jacobi polynomial to include
# Return Values
-   `D`: the differentiation matrix
# Example
    See LegendreTests.jl
"""
function dmatrix(x, α, β, N)
    # calculate vandermonde matrix and grad of vandermonde matrix
    vr = dvandermonde(x, α, β, N)
    v  =  vandermonde(x, α, β, N)

    # calculate values using D = vr * v^-1
    d = vr / v

    return d
end

"""
lift1D(V, y)
for computing fluxes
helps compute a surface integral of a quantity
note that the parentheses are necessary to prevent too much multiplcation
the E function takes the surface integrals are presents it
with respect to the full space inside an element
the entire operator represents how fluxes flow
into the interior of an element
"""
function lift1D(V)
    m,n = size(V)

    E = zeros(m , 2)
    E[1,1] = 1.0
    E[m,2] = 1.0

    return V * (transpose(V) * E)
end

"""
lift1D_v2(V, y)
for computing fluxes
nodal form
helps compute a surface integral of a quantity
note that the parantheses are necessary to prevent too much multiplcation
the E function takes the surface integrals are presents it
with respect to the full space inside an element
the entire operator represents how fluxes flow
into the interior of an element
"""
function lift1D_v2(V)
    m,n = size(V)

    E = zeros(m , 2)
    E[1,1] = 1.0
    E[m,2] = 1.0

    return E
end


"""
jacobiGQ(α, β, N)
# Description
    Guass Quadrature points and weights for the Jacobi Polynomial (α,β)
# Input
α, β: Jacobi polynomial descriptors
N:    order of quadrature points
# Return: x,w
x: quadrature points | array of size N+1
w: quadrature weights | array of size N+1
#Example
α = 0
β = 0
N = 4
x, w = jacobiGQ(α, β, N)
"""
function jacobiGQ(α, β, N)
    N == 0 && return [(α-β) / (α+β+2)], [2]

    # Form symmetric matrix from recurrence.
    dv = OffsetArray(zeros(N+1), 0:N)  # diagonal vector
    ev = OffsetArray(zeros(N+1), 0:N)  # sub/super-diagonal vector

    for n in 0:N
        dv[n] = bᴾ(α, β, n)
        ev[n] = aᴾ(α, β, n)
    end

    # Create full matrix combining the two.
    # Need to pass arrays that are not offset.
    J = SymTridiagonal(dv[0:N], ev[1:N])
    (α + β) ≈ 0 && (J[1, 1] = 0)

    # Compute quadrature points and weights by eigenvalue solve.
    x, V = eigen(J)
    w = @. V[1, :]^2 * 2^(α+β+1) / (α+β+1)
    @. w *= factorial(α) * factorial(β) / factorial(α+β)

    return x, w
end

"""
jacobiGL(α, β, N)
# Description
    Guass Labatto quadrature points for the Jacobi Polynomial (α,β)
    The quadrature weights are computed as well (but not returned)
# Arguments
- `α, β`: Jacobi polynomial descriptors
- `N`:    order of quadrature
# Return: x
- `x`: quadrature points  | array of size N+1
# Examples
```julia-repl
julia> x = jacobiGL(0, 0, 4)
5-element Array{Float64,1}:
 -1.0
 -0.6546536707079759
  4.440892098500626e-16
  0.6546536707079771
  1.0
```
"""
function jacobiGL(α, β, N)
    N == 0 && error("What are you doing? Go back to finite volume land.")
    N == 1 && return [-1, 1]

    x = zeros(N+1)
    x[1], x[N+1] = -1, 1

    x_GQ, _ = jacobiGQ(α+1, β+1, N-2)
    x[2:N] .= x_GQ

    return x
end

abstract type AbstractMesh end

"""
gridvalues1D(xmin, xmax, K)

# Description

    Generates physical gridpoints with each element

# Arguments

    VX: vertex values | an Array of size K+1

    EtoV: element to node connectivity | a Matrix of size Kx2

    r: LGL nodes in reference element | an array

# Return Values: x

    x: physical coordinates of solution

# Example (uses ../utils.jl as well)

xmin = 0
xmax = 2π
K = 4
# call functions
VX, EtoV = unimesh1D(xmin, xmax, K)
r = jacobiGL(0, 0, 4)
x = gridvalues1D(VX, EtoV, r)
# x[:,1] is the physical coordinates within the first element
# for plotting
f(x) = sin(x)
plot(x, f.(x))
# scatter(x,f.(x)) tends to work better
"""
function gridvalues1D(VX, EtoV, r)
    # get low and high edges
    va = view(EtoV, :, 1)
    vb = view(EtoV, :, 2)

    # compute physical coordinates of the grid points
    x = ones(length(r),1) * (VX[va]') .+ 0.5 .* (r .+ 1 ) * ((VX[vb] - VX[va])')
    return x
end

"""
facemask1D(r)

# Description

    creates face mask

# Arguments

    r: GL points

# Return Values: x

    fmask1: standard facemask
    fmask2: alternate form

# Example | ../utils.jl

r = jacobiGL(0, 0, 4)
fmask = fmask1D(r)

"""
function fmask1D(r)
    # check if index is left or right edge
    fm1 = @. abs(r+1) < eps(1.0);
    fm2 = @. abs(r-1) < eps(1.0);
    fmask1 = (fm1,fm2)

    # alternate form
    tmp = collect(1:length(r))
    fmask2  = [tmp[fm1]; tmp[fm2]]

    fmask = (fmask1, fmask2)
    return fmask
end

"""
edgevalues1D(fmask, x)

# Description

    calculates edge values

# Arguments

    fmask: face mask for GL edges

    x:  physical coordinates of solution on each element

# Return Values: x

    fx: face values of x

# Example | ../utils.jl

r = jacobiGL(0, 0, 4)
x = gridvalues1D(VX, EtoV, r)
fmask = fmask1D(r)[1]
fx = edgevalues1D(fmask,x)

# the locations of the edges in element 1 is fx[:, 1]


"""
function edgevalues1D(fmask, x)
    # compute x values at selected indices
    fx1 = x[fmask[1],:]
    fx2 = x[fmask[2],:]

    # return list of physical edge positions
    fx = [fx1; fx2]
    return fx
end

"""
normals1D(K)

# Description

    calculates face normals

# Arguments

    K: number of elements

# Return Values: normals

    normals: face normals along each grid

# Example

"""
function normals1D(K)
    normals  = ones(2,K)
    @. normals[1,:] *= -1
    return normals
end

"""
geometric_factors(x, Dʳ)

# Description

    computes the geometric factors for local mappings of 1D elements

# Arguments

    x: physical coordinates of solution for each element

    Dʳ:

# Return Values: rx, J

    rx: inverse jacobian

    J: jacobian (in 1D a scalar)

# Example

"""
function geometric_factors(x, Dʳ)
    J = Dʳ * x
    rx = 1 ./ J # for 1D
    return rx, J
end

"""
connect1D(EtoV)

# Description

    builds global connectivity arrays for 1D

# Arguments

    EtoV: element to node connectivity | a Matrix of size Kx2

# Return Values: EtoE, EtoF

    EtoE: element to element connectivity
    EtoF: element to face connectivity

# Example

"""
function connect1D(EtoV)
    nFaces = 2 # for 1d elements

    # Find number of elements and vertices
    K = size(EtoV,1)
    total_faces = nFaces * K
    Nv = K+1

    # list of local face to local vertex connections
    vn = [1, 2]

    # build global face to vertex array
    FtoV = Int.(spzeros(total_faces, Nv))
    let sk = 1
        for k = 1:K
            for faces = 1:nFaces
                FtoV[sk, EtoV[k, vn[faces]]] = 1;
                sk += 1
            end
        end
    end

    # build global face to face array
    FtoF = FtoV * (FtoV') - sparse(I, total_faces, total_faces)

    # find all face to face connections
    # check
    #(faces1, faces2) = findnz(FtoF)
    faces1, faces2 = findnz(FtoF .== 1)

    # convert global face number to element and face numbers
    element1 = @. floor(Int, (faces1 - 1) / nFaces) + 1
    face1    = @. Int( mod( (faces1 - 1),  nFaces) + 1)
    element2 = @. floor(Int, (faces2 - 1) / nFaces) + 1
    face2    = @. Int( mod( (faces2 - 1),  nFaces) + 1)

    # Rearrange into Nelement x Nfaces sized arrays
    ind = diag( LinearIndices(ones(K, nFaces))[element1,face1] ) # this line is a terrible idea.
    EtoE = collect(1:K) * ones(1, nFaces)
    EtoF = ones(K, 1) * (collect(1:nFaces)')
    EtoE[ind] = copy(element2);
    EtoF[ind] = face2;
    return EtoE, EtoF
end

"""
buildmaps1D(K, nGL, nFP, nFaces, fmask, EtoE, EtoF, x)
# Description

    connectivity matrices for element to elements and elements to face

# Arguments

-   `K`: number of elements
-   `nGL`: number of points within an element (polynomial degree + 1)
-   `nFP`: 1
-   `nFaces`: 2
-   `fmask`: an element by element mask to extract edge values
-   `EtoE`: element to element connectivity
-   `EtoF`: element to face connectivity
-   `x`: Guass lobatto points

# Return Values: vmapM, vmapP, vmapB, mapB, mapI, mapO, vmapI, vmapO

-   `vmapM`: vertex indices, (used for interior u values)
-   `vmapP`: vertex indices, (used for exterior u values)
-   `vmapB`: vertex indices, corresponding to boundaries
-   `mapB`: use to extract vmapB from vmapM
-   `mapI`: Index of left boundary condition
-   `mapO`: Index of right boundary condition

# Example | uses ../utils.jl

K = 3
n = 3; α = 0; β = 0; xmin = 0; xmax = 2π;
nGL = n + 1
nFP = 1
nFaces = 2

r = jacobiGL(α, β, n)

VX, EtoV = unimesh1D(xmin, xmax, K)
EtoE, EtoF = connect1D(EtoV)
x = gridvalues1D(VX, EtoV, r)
fx = edgevalues1D(r,x)

vmapM, vmapP, vmapB, mapB, mapI, mapO, vmapI, vmapO = buildmaps1D(K, nGL, nFP, nFaces, fmask, EtoE, EtoF, x)
"""
function buildmaps1D(K, nGL, nFP, nFaces, fmask, EtoE, EtoF, x)
    # number volume nodes consecutively
    nodeids = reshape(collect(1:(K*nGL)), nGL, K)
    vmapM = zeros(nFP, nFaces, K)
    vmapP = zeros(nFP, nFaces, K)
    # find index of face nodes wrt volume node ordering
    for k1 in 1:K
        for f1 in 1:nFaces
            vmapM[:, f1, k1] = nodeids[fmask[f1], k1]
        end
    end

    for k1 = 1:K
        for f1 = 1:nFaces
            # find neighbor
            k2 = Int.( EtoE[k1, f1])
            f2 = Int.( EtoF[k1, f1])

            # find volume node numbers of left and right nodes
            vidM = Int.( vmapM[:, f1, k1])
            vidP = Int.( vmapM[:, f2, k2])

            x1 = x[vidM]
            x2 = x[vidP]

            # compute distance matrix
            # need to figure out this part
            D = @. (x1 - x2)^2
            m = length(x1)
            for j = 1:m
                if D[j] < eps(1.0)*10^5
                    vmapP[j, f1, k1] = vidP[j]
                end
            end
        end
    end

    # reshape arrays
    vmapP = Int.( reshape(vmapP, length(vmapP)) )
    vmapM = Int.( reshape(vmapM, length(vmapM)) )

    # Create list of boundary nodes
    mapB = Int.( collect(1:length(vmapP))[vmapP .== vmapM] )
    vmapB = Int.( vmapM[mapB] )

    # inflow and outflow maps
    mapI = 1
    mapO = K * nFaces
    vmapI = 1
    vmapO = K*nGL
    return vmapM, vmapP, vmapB, mapB, mapI, mapO, vmapI, vmapO
end

"""
make_periodic1D!(vmapP, u)

# Description

    makes the grid periodic by modifying vmapP.
    Assumes that the first node is connected to the last.

# Arguments

    vmapP: exterior vertex map
    u: vertex vector

# Return Values: none

# Example

"""
function make_periodic1D!(vmapP, u)
    vmapP[1] = length(u)
    vmapP[end] = 1

    return nothing
end

struct Mesh{T,S,U,W} <: AbstractMesh
    # inputs
    K::S
    n::S

    # face stuff
    nFP::S
    nFaces::S

    # GL points
    r::U
    x::T

    # vertex maps
    vmapM::W
    vmapP::W
    vmapB::W
    mapB::W

    # inflow/outflow maps
    mapI::S
    mapO::S
    vmapI::S
    vmapO::S

    # structures for computation
    D::T
    M::T
    Mi::T
    lift::T
    rx::T
    normals::T
    fscale::T
end

"""
mesh(K, n, xmin, xmax)

# Description

    outer_constructor mesh struct

# Arguments

- `K`: number of elements
- `n`: polynomial order
- `xmin`: lower bound
- `xmax`: upper bound


# Return Values: x
    return grid values

"""
function Mesh(KK, nn, xmin, xmax; periodic = false)
    # initialize parameters
    K = KK
    α = 0; β = 0;
    n = nn

    # number of vertices
    nGL = n+1
    nFP = 1
    nFaces = 2

    # compute Gauss Lobatto grid
    r = jacobiGL(α, β, n)

    # build grid
    VX, EtoV = unimesh1D(xmin, xmax, K)

    # build coordinates of all the nodes
    x = gridvalues1D(VX, EtoV, r)

    # build connectivity matrix
    EtoE, EtoF = connect1D(EtoV)

    # build face masks
    fmask1,fmask2 = fmask1D(r)

    # build connectivity maps
    vmapM,vmapP,vmapB,mapB, mapI,mapO,vmapI,vmapO = buildmaps1D(K, nGL,nFP,nFaces, fmask1, EtoE,EtoF, x)

    # build differentiation matrix
    D = dmatrix(r, α, β, n)

    # build surface integral terms
    V = vandermonde(r, α, β, n)
    lift = lift1D(V)

    # build mass matrix and inverse of mass matrix
    Mi = V * V'
    M = inv(Mi)

    # calculate geometric factors
    rx,J = geometric_factors(x, D)

    # build surface normals
    normals = normals1D(K)

    # build inverse metric at the surface
    fscale = 1 ./ J[fmask2,:]

    # hack for periodicity
    if periodic
        vmapP[1] = vmapM[end]
        vmapP[end] = vmapM[1]
    end

    return Mesh{typeof(x),typeof(K),typeof(r),typeof(vmapP)}(K, n, nFP, nFaces, r, x, vmapM, vmapP, vmapB, mapB, mapI, mapO, vmapI, vmapO, D, M, Mi, lift, rx, normals, fscale)
end

function compute_volume_terms(data::AbstractArray, mesh::Mesh)
    q = mesh.D * data
    @. q *= mesh.rx
    return q
end

function compute_surface_terms(data::AbstractArray, 𝒢::Mesh)
    diffs = reshape( (data[𝒢.vmapM] - data[𝒢.vmapP]), (𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
    # Include factor of 2 for the weak-strong form
    @. diffs *= 1.0 / 2.0
    # Compute Lift Operator
    lifted = - 𝒢.lift * (𝒢.fscale .* 𝒢.normals .* diffs)
    return lifted
end