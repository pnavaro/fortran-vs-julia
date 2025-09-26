struct Mesh

    dimx::Float64
    nx::Int64
    dx::Float64
    dimy::Float64
    ny::Int64
    dy::Float64
    kx::Any
    ky::Any

    function Mesh(dimx, nx, dimy, ny)

        dx = dimx / nx
        dy = dimy / ny
        kx = 2π ./ dimx .* vcat(0:(nx ÷ 2 - 1), (-nx ÷ 2):-1)
        ky = 2π ./ dimy .* vcat(0:(ny ÷ 2 - 1), (-ny ÷ 2):-1) |> transpose

        return new(dimx, nx, dx, dimy, ny, dy, kx, ky)

    end

end

struct MeshFields

    mesh::Mesh
    ex::Array{Float64, 2}
    ey::Array{Float64, 2}
    bz::Array{Float64, 2}

    function MeshFields(mesh)

        nx, ny = mesh.nx, mesh.ny
        ex = zeros(Float64, (nx, ny + 1))
        ey = zeros(Float64, (nx + 1, ny))
        bz = zeros(Float64, (nx + 1, ny + 1))
        return new(mesh, ex, ey, bz)

    end

    function MeshFields(mesh, nx, ny)

        ex = zeros(Float64, (nx + 1, ny + 1))
        ey = zeros(Float64, (nx + 1, ny + 1))
        bz = zeros(Float64, (nx + 1, ny + 1))
        return new(mesh, ex, ey, bz)

    end

end
