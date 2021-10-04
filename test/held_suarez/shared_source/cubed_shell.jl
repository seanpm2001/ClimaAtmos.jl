using CubedSphere, Rotations

# This function was taken from https://github.com/CliMA/ClimateMachine.jl/blob/master/src/Numerics/Mesh/Topologies.jl
# Author: Valeria Barra <valeriabarra21@gmail.com>

function cubedshellwarp(a, b, c, R = max(abs(a), abs(b), abs(c)))

    fdim = argmax(abs.((a, b, c)))
    M = max(abs.((a, b, c))...)
    if fdim == 1 && a < 0
        # left face
        x1, x2, x3 = conformal_cubed_sphere_mapping(-b / M, c / M)
        x1, x2, x3 = RotX(π / 2) * RotY(-π / 2) * [x1, x2, x3]
    elseif fdim == 2 && b < 0
        # front face
        x1, x2, x3 = conformal_cubed_sphere_mapping(a / M, c / M)
        x1, x2, x3 = RotX(π / 2) * [x1, x2, x3]
    elseif fdim == 1 && a > 0
        # right face
        x1, x2, x3 = conformal_cubed_sphere_mapping(b / M, c / M)
        x1, x2, x3 = RotX(π / 2) * RotY(π / 2) * [x1, x2, x3]
    elseif fdim == 2 && b > 0
        # back face
        x1, x2, x3 = conformal_cubed_sphere_mapping(a / M, -c / M)
        x1, x2, x3 = RotX(-π / 2) * [x1, x2, x3]
    elseif fdim == 3 && c > 0
        # top face
        x1, x2, x3 = conformal_cubed_sphere_mapping(a / M, b / M)
    elseif fdim == 3 && c < 0
        # bottom face
        x1, x2, x3 = conformal_cubed_sphere_mapping(a / M, -b / M)
        x1, x2, x3 = RotX(π) * [x1, x2, x3]
    else
        error("invalid case for cubed_sphere_warp(::ConformalCubedSphere): $a, $b, $c")
    end

    return x1 * R, x2 * R, x3 * R


end
