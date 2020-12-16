include("mesh.jl")
include("pstd.jl")

function main( nstep )

    @show nx   = 1024	   # x number of points
    @show ny   = 1024	   # y number of points
    @show dimx = 2π	   # width
    @show dimy = 2π	   # height

    dx = dimx / nx
    dy = dimy / ny

    mesh = Mesh( dimx, nx, dimy, ny)
    pstd = PSTD( mesh )

    dx, dy = mesh.dx, mesh.dy

    x = LinRange( 0, dimx, nx+1)[1:end-1] |> collect 
    y = LinRange( 0, dimy, ny+1)[1:end-1] |> collect |> transpose
    
    dt = (dx + dy) / ( π * sqrt(2))
    
    ex = zeros( nx, ny)
    ey = zeros( nx, ny)
    bz = zeros( nx, ny)
    
    omega = sqrt(2)

    # Ex and Ey are set at t = 0.0
    # Bz is set at  t = -dt/2

    bz .= - cos.(x) .* cos.(y) .* cos.(omega*(-0.5*dt))
    
    for istep in 1:nstep

        # Bz(n-1/2) --> B(n+1/2) 
        faraday!(bz, pstd, ex, ey, dt)   

        # Ex(n) --> B(n+1/2) 
        # Ey(n) --> B(n+1/2)
        ampere_maxwell!(ex, ey, pstd, bz, dt) 

    end

    @show time = (nstep-0.5)*dt
    err_l2 = 0.0
    for j = 1:ny, i = 1:nx
        th_bz = - cos(x[i]) * cos(y[j]) * cos(omega*time)
        err_l2 += (bz[i,j] - th_bz)^2
    end

    return sqrt(err_l2)


end
# -

main( 1 ) # trigger building

@time println(main( 400 ))
