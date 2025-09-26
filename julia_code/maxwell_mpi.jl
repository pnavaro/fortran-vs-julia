using MPI
using Printf

const c = 1.0 # speed of light
const csq = c * c

struct Mesh

    nx::Int64
    dx::Float64
    ny::Int64
    dy::Float64

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

end


function faraday!(fields, dt)

    dx, dy = fields.mesh.dx, fields.mesh.dy
    nx, ny = fields.mesh.nx, fields.mesh.ny

    for j in 1:ny, i in 1:nx
        dex_dy = (fields.ex[i, j + 1] - fields.ex[i, j]) / dy
        dey_dx = (fields.ey[i + 1, j] - fields.ey[i, j]) / dx
        fields.bz[i, j] = fields.bz[i, j] + dt * (dex_dy - dey_dx)
    end

    return
end

function ampere_maxwell!(fields, dt)

    dx, dy = fields.mesh.dx, fields.mesh.dy
    nx, ny = fields.mesh.nx, fields.mesh.ny

    for j in 2:(ny + 1), i in 1:nx
        dbz_dy = (fields.bz[i, j] - fields.bz[i, j - 1]) / dy
        fields.ex[i, j] = fields.ex[i, j] + dt * csq * dbz_dy
    end

    for j in 1:ny, i in 2:(nx + 1)
        dbz_dx = (fields.bz[i, j] - fields.bz[i - 1, j]) / dx
        fields.ey[i, j] = fields.ey[i, j] - dt * csq * dbz_dx
    end

    return
end

include("plot_mpi.jl")

function main(nstep)

    cfl = 0.4    # Courant-Friedrich-Levy
    tfinal = 1.0    # final time
    nstepmax = 1000 # max steps
    md = 2          # md : wave number x (initial condition)
    nd = 2          # nd : wave number y (initial condition)
    nx = 1200       # x number of points
    ny = 1200       # y number of points
    dx = 0.01       # width
    dy = 0.01       # height
    dimx = nx * dx
    dimy = ny * dy

    comm = MPI.COMM_WORLD
    proc = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)

    MPI.Barrier(comm)

    tcpu1 = MPI.Wtime()

    dims = [0 0]
    ndims = length(dims)
    periods = [1 1]
    reorder = 1

    MPI.Dims_create!(proc, dims)
    comm2d = MPI.Cart_create(comm, dims, periods, reorder)

    @assert MPI.Comm_size(comm2d) == proc

    north, south = MPI.Cart_shift(comm2d, 0, 1)
    west, east = MPI.Cart_shift(comm2d, 1, 1)

    coords = MPI.Cart_coords(comm2d)

    nxp, nyp = dims

    mx = nx ÷ nxp
    my = ny ÷ nyp

    dt = cfl / sqrt(1 / dx^2 + 1 / dy^2) / c

    nstep = min(nstepmax, nstep)

    MPI.Barrier(comm)

    # Origin of local mesh
    xp = coords[1] * dimx / nyp
    yp = coords[2] * dimy / nyp

    mesh = Mesh(mx, dx, my, dy)

    fields = MeshFields(mesh)

    omega = c * sqrt((md * pi / dimx)^2 + (nd * pi / dimy)^2)
    for j in 1:my, i in 1:mx
        x = xp + (i - 0.5) * dx
        y = yp + (j - 0.5) * dy
        fields.bz[i, j] =
            (-cos(md * pi * x / dimx) * cos(nd * pi * y / dimy) * cos(omega * 0.5 * dt))
    end

    tag = 1111

    for istep in 1:nstep # Loop over time

        # E(n) [1:mx]*[1:my] --> B(n+1/2) [1:mx-1]*[1:my-1]

        faraday!(fields, dt)


        # Send to North  and receive from South
        MPI.Sendrecv!(
            view(fields.bz, 1, 1:my),
            north,
            tag,
            view(fields.bz, mx + 1, 1:my),
            south,
            tag,
            comm2d,
        )

        # Send to West and receive from East
        MPI.Sendrecv!(
            view(fields.bz, 1:mx, 1),
            west,
            tag,
            view(fields.bz, 1:mx, my + 1),
            east,
            tag,
            comm2d,
        )

        # Bz(n+1/2) [1:mx]*[1:my] --> Ex(n+1) [1:mx]*[2:my+1]
        # Bz(n+1/2) [1:mx]*[1:my] --> Ey(n+1) [2:mx+1]*[1:my]

        ampere_maxwell!(fields, dt)

        # Send to East and receive from West
        MPI.Sendrecv!(
            view(fields.ex, :, my + 1),
            east,
            tag,
            view(fields.ex, :, 1),
            west,
            tag,
            comm2d,
        )

        # Send to South and receive from North
        MPI.Sendrecv!(
            view(fields.ey, mx + 1, :),
            south,
            tag,
            view(fields.ey, 1, :),
            north,
            tag,
            comm2d,
        )

        # plot_fields(mesh, rank, proc, fields.bz, xp, yp, istep )


    end # next time step


    err_l2 = 0.0
    time = (nstep - 0.5) * dt
    for j in 1:my, i in 1:mx
        x = xp + (i - 0.5) * dx
        y = yp + (j - 0.5) * dy
        th_bz = (-cos(md * pi * x / dimx) * cos(nd * pi * y / dimy) * cos(omega * time))
        err_l2 += (fields.bz[i, j] - th_bz)^2
    end

    for k in 0:(proc - 1)
        if rank == k
            println("----")
            println("$rank : $mx, $my  $err_l2 ")
            println("$rank : $dx, $dy  $err_l2 ")
        end
        MPI.Barrier(comm)
    end


    sum_err_l2 = MPI.Allreduce(err_l2, +, comm2d)

    return sqrt(sum_err_l2)

end


MPI.Init()
tbegin = MPI.Wtime()

err_l2 = main(1)
@time err_l2 = main(1000)

tend = MPI.Wtime()

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    println(" error : $(err_l2) ")
    println(" time : $(tend - tbegin) ")
end

MPI.Finalize()
