"""
    make_function_space(domain::Column)
"""
function make_function_space(domain::Column{FT}) where {FT}
    column = ClimaCore.Domains.IntervalDomain(
        Geometry.ZPoint{FT}(domain.zlim[1]),
        Geometry.ZPoint{FT}(domain.zlim[2]);
        boundary_tags = (:bottom, :top),
    )
    mesh = Meshes.IntervalMesh(column; nelems = domain.nelements)
    center_space = Spaces.CenterFiniteDifferenceSpace(mesh)
    face_space = Spaces.FaceFiniteDifferenceSpace(center_space)

    return center_space, face_space
end

"""
    make_function_space(domain::Plane)
"""
function make_function_space(domain::Plane{FT}) where {FT}
    rectangle = ClimaCore.Domains.RectangleDomain(
        Interval(
            Geometry.XPoint(domain.xlim[1]),
            Geometry.XPoint(domain.xlim[2]),
        ),
        Interval(
            Geometry.YPoint(domain.ylim[1]),
            Geometry.YPoint(domain.ylim[2]),
        ),
        x1periodic = domain.periodic[1],
        x2periodic = domain.periodic[2],
    )
    mesh = Meshes.EquispacedRectangleMesh(
        rectangle,
        domain.nelements[1],
        domain.nelements[2],
    )
    grid_topology = Topologies.GridTopology(mesh)
    quad = Spaces.Quadratures.GLL{domain.npolynomial + 1}()
    space = Spaces.SpectralElementSpace2D(grid_topology, quad)

    return space
end

"""
    make_function_space(domain::HybridPlane)
"""
function make_function_space(domain::HybridPlane{FT}) where {FT}
    vertdomain = ClimaCore.Domains.IntervalDomain(
        Geometry.ZPoint{FT}(domain.zlim[1]),
        Geometry.ZPoint{FT}(domain.zlim[2]);
        boundary_tags = (:bottom, :top),
    )

    vertmesh = Meshes.IntervalMesh(vertdomain, nelems = domain.nelements[2])
    vert_center_space = Spaces.CenterFiniteDifferenceSpace(vertmesh)

    horzdomain = ClimaCore.Domains.RectangleDomain(
        Interval(
            Geometry.XPoint(domain.xlim[1]),
            Geometry.XPoint(domain.xlim[2]),
        ),
        Interval(Geometry.YPoint(-0), Geometry.YPoint(0)),
        x1periodic = true,
        x2boundary = (:a, :b),
    )
    horzmesh =
        Meshes.EquispacedRectangleMesh(horzdomain, domain.nelements[1], 1)
    horztopology = Topologies.GridTopology(horzmesh)

    quad = Spaces.Quadratures.GLL{domain.npolynomial + 1}()
    horzspace = Spaces.SpectralElementSpace1D(horztopology, quad)

    hv_center_space =
        Spaces.ExtrudedFiniteDifferenceSpace(horzspace, vert_center_space)
    hv_face_space = Spaces.FaceExtrudedFiniteDifferenceSpace(hv_center_space)

    return hv_center_space, hv_face_space
end

"""
    make_function_space(domain::Sphere)
"""
function make_function_space(domain::Sphere{FT}) where {FT}
    sphere = ClimaCore.Domains.SphereDomain(domain.radius)
    mesh =
        Meshes.Mesh2D(sphere, Meshes.EquiangularSphereWarp(), domain.nelements)
    grid_topology = Topologies.Grid2DTopology(mesh)
    quad = Spaces.Quadratures.GLL{domain.npolynomial + 1}()
    space = Spaces.SpectralElementSpace2D(grid_topology, quad)

    return space
end
