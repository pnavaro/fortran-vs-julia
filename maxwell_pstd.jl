using LinearAlgebra

const c = 1.0 # speed of light 
const csq = c * c

include("mesh.jl")
include("pstd.jl")


function main( nstep )

    @show cfl    = 0.1     # Courant-Friedrich-Levy
    @show tfinal = 10.     # final time
    @show nstepmax = 1000  # max steps
    @show md = 2	   # md : wave number x (initial condition)
    @show nd = 2	   # nd : wave number y (initial condition)
    @show nx = 1024	   # x number of points
    @show ny = 1024	   # y number of points
    @show dimx = 1.0	   # width
    @show dimy = 1.0	   # height

    dx = dimx / nx
    dy = dimy / ny

    mesh = Mesh( dimx, nx, dimy, ny)

    dx, dy = mesh.dx, mesh.dy

    x = LinRange( 0, dimx, nx+1)[1:end-1]
    y = LinRange( 0, dimy, ny+1)[1:end-1]
    
    dt = cfl / sqrt(1/dx^2+1/dy^2) / c
    
    @show nstep  = min( nstepmax, nstep)
    
    ex = zeros(ComplexF64, nx, ny)
    ey = zeros(ComplexF64, nx, ny)
    bz = zeros(ComplexF64, nx, ny)
    
    omega = c * sqrt((md*pi/dimx)^2+(nd*pi/dimy)^2)

    # Ex and Ey are set at t = 0.0
    # Bz is set at  t = -dt/2

    bz .= ( - cos.(md*pi*x/dimx) 
           .* cos.(nd*pi*y/dimy)
           .* cos.(omega*(-0.5*dt)) )
    
    for istep = 1:nstep # Loop over time
    
       # Ex(n) [1:nx]*[1:ny+1] --> B(n+1/2) [1:nx]*[1:ny]
       # Ey(n) [1:nx+1]*[1:ny] --> B(n+1/2) [1:nx]*[1:ny]

       ampere_maxwell!(ex, ey, mesh, bz, dt) 

       faraday!(bz, mesh, ex, ey, dt)   
    
       time = (istep-0.5)*dt
       err_l2 = 0.0
       for j = 1:ny, i = 1:nx
           th_bz = (- cos(md*pi*x[i]/dimx)
                    * cos(nd*pi*y[j]/dimy)
                    * cos(omega*time))
           err_l2 += (real(bz[i,j]) - th_bz)^2
       end

       println(sqrt(err_l2))

    end # next time step

end

main( 1 ) # trigger building

@time println(main( 1000 ))
