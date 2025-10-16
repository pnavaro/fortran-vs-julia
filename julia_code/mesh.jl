struct Mesh

    dimx::Float64
    nx::Int
    dx::Float64
    dimy::Float64
    ny::Int
    dy::Float64

    function Mesh(dimx::Float64, nx::Int, dimy::Float64, ny::Int)

        dx = dimx / nx
        dy = dimy / ny

        return new(dimx, nx, dx, dimy, ny, dy)

    end

end

struct MeshFields

    mesh::Mesh
    ex::Array{Float64, 2}
    ey::Array{Float64, 2}
    bz::Array{Float64, 2}

    function MeshFields(mesh::Mesh)

        nx, ny = mesh.nx, mesh.ny
        ex = zeros(Float64, (nx, ny + 1))
        ey = zeros(Float64, (nx + 1, ny))
        bz = zeros(Float64, (nx + 1, ny + 1))
        return new(mesh, ex, ey, bz)

    end

    function MeshFields(mesh::Mesh, nx::Int, ny::Int)

        ex = zeros(Float64, (nx + 1, ny + 1))
        ey = zeros(Float64, (nx + 1, ny + 1))
        bz = zeros(Float64, (nx + 1, ny + 1))
        return new(mesh, ex, ey, bz)

    end

end
