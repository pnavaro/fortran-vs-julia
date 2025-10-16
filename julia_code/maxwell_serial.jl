const csq = 1.0

include("mesh.jl")
include("fdtd.jl")

function @main(ARGS)

	if length(ARGS) == 0
		println(Core.stdout, "No command line arguments provided. nstep = 1")
        nstep = 1
	else
        nstep = parse(Int, first(ARGS))
	end

	println(Core.stdout, nstep)

    c = 1.0 # speed of light
    csq = c * c

    cfl = 0.1        # Courant-Friedrich-Levy
    tfinal = 10.0    # final time
    nstepmax = 1000  # max steps
    md = 2           # md : wave number x (initial condition)
    nd = 2           # nd : wave number y (initial condition)
    nx = 1200        # x number of points
    ny = 1200        # y number of points
    dimx = 1.0       # width
    dimy = 1.0       # height

    dx = dimx / nx
    dy = dimy / ny

    mesh = Mesh(dimx, nx, dimy, ny)

    dx, dy = mesh.dx, mesh.dy

    dt = cfl / sqrt(1 / dx^2 + 1 / dy^2) / c

    nstep = min(nstepmax, nstep)

    fields = MeshFields(mesh)

    omega = c * sqrt((md * pi / dimx)^2 + (nd * pi / dimy)^2)

    # Ex and Ey are set at t = 0.0
    # Bz is set at  t = -dt/2

    for j in 1:ny, i in 1:nx
        fields.bz[i, j] = (
            -cos(md * pi * ((i - 0.5) * dx / dimx)) *
                cos(nd * pi * ((j - 0.5) * dy / dimy)) *
                cos(omega * (-0.5 * dt))
        )
    end

    time = 0.0
    err_l2 = Inf

    for istep in 1:nstep # Loop over time

        # Ex(n) [1:nx]*[1:ny+1] --> B(n+1/2) [1:nx]*[1:ny]
        # Ey(n) [1:nx+1]*[1:ny] --> B(n+1/2) [1:nx]*[1:ny]

        faraday!(fields, 1, nx, 1, ny, dt)

        time += 0.5dt

        periodic_bc!(fields, 1, nx, 1, ny, dt)

        ampere_maxwell!(fields, 1, nx, 1, ny, dt)


        err_l2 = 0.0
        time = (istep - 0.5) * dt
        for j in 1:ny, i in 1:nx
            th_bz = (
                -cos(md * pi * ((i - 0.5) * dx / dimx)) *
                 cos(nd * pi * ((j - 0.5) * dy / dimy)) *
                 cos(omega * time)
            )
            err_l2 += (fields.bz[i, j] - th_bz)^2
        end

        time += 0.5dt

        println(Core.stdout, sqrt(err_l2))

    end # next time step


    return 0
end
