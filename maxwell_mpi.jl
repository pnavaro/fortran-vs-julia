using MPI
using Printf

const c = 1.0 # speed of light 
const csq = c * c

include("mesh.jl")
include("fdtd.jl")

function plot_fields(mesh, rank, proc, field, xp, yp, iplot )

    dx, dy = mesh.dx, mesh.dy
    ix, jx = 1, mesh.mx
    iy, jy = 1, mesh.my

    if iplot == 1
        mkpath("data/$rank")
    end

    io = open("data/$(rank)/$(iplot)", "w")
    for j=iy:jy
        for i=ix:jx
            @printf( io, "%f %f %f \n", xp+(i-0.5)*dx, yp+(j-1)*dy, field[i,j])
        end
        @printf( io, "\n")
    end
    close(io)
   
    # write master file

    if rank == 0

      if iplot == 1 
         io = open( "field.gnu", "w" )
         write(io, "set xr[-0.1:1.1]\n")
         write(io, "set yr[-0.1:1.1]\n")
         write(io, "set zr[-1.1:1.1]\n")
         write(io, "set surf\n")
      else
         io = open( "field.gnu", "a" )
      end
      write(io, "set title '$(iplot)' \n")
      write(io, "splot 'data/$(rank)/$(iplot)' u 1:2:3 w lines")
   
      for p = 1:proc-1
         write(io, ", 'data/$(p)/$(iplot)' u 1:2:3 w lines")
      end
      write( io, "\n")
      write( io, "set term gif \n")
      write( io, "set output 'image$(lpad(iplot,3,"0")).gif'\n")
      write( io, "replot\n")

      close(io)

    end

end 

function main( nstep )

    cfl    = 0.2    # Courant-Friedrich-Levy
    tfinal = 1.	    # final time
    nstepmax = 500  # max steps
    md = 2          # md : wave number x (initial condition)
    nd = 2          # nd : wave number y (initial condition)
    nx = 1000       # x number of points
    ny = 1000       # y number of points
    dimx = 1.0      # width
    dimy = 1.0      # height

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
    
    north, south = MPI.Cart_shift(comm2d,0,1)
    west,  east  = MPI.Cart_shift(comm2d,1,1)
    
    coords = MPI.Cart_coords(comm2d)
    
    nxp, nyp = dims

    mx = nx ÷ nxp
    my = ny ÷ nyp

    @assert nx % nxp == 0
    @assert ny % nyp == 0

    global_mesh = Mesh( dimx, nx, dimy, ny)
    dx, dy = global_mesh.dx, global_mesh.dy

    dt = cfl / sqrt(1/dx^2+1/dy^2) / c
    
    nstep  = min( nstepmax, nstep)
    
    MPI.Barrier(comm)
    
    # Origin of local mesh
    xp = coords[1]/nxp * dimx 
    yp = coords[2]/nyp * dimy 

    local_mesh = Mesh( dimx/nxp, mx, dimy/nyp, my)

    @assert global_mesh.dx == local_mesh.dx
    @assert global_mesh.dy == local_mesh.dy
    
    fields = MeshFields(local_mesh)
    
    omega = c * sqrt((md*pi/dimx)^2+(nd*pi/dimy)^2)
    for j=1:my, i=1:mx
        fields.bz[i,j] = (- cos(md*pi*(xp+(i-0.5)*dx/dimx)) 
                      * cos(nd*pi*(yp+(j-0.5)*dy/dimy))
                      * cos(omega*(-0.5*dt)) )
    end  
    
    tag = 1111
    
    for istep = 1:nstep # Loop over time
    
       # E(n) [1:mx]*[1:my] --> B(n+1/2) [1:mx-1]*[1:my-1]
       faraday!(fields, 1, mx, 1, my, dt)   
    
       # Bz(n+1/2) [1:mx]*[1:my] --> Ex(n+1) [1:mx]*[2:my]
       # Bz(n+1/2) [1:mx]*[1:my] --> Ey(n+1) [2:mx]*[1:my]
       ampere_maxwell!(fields, 1, mx, 1, my, dt) 
    
       # Send to East and receive from West
       MPI.Sendrecv!(fields.ex[1:mx,my+1], east, tag,
                     fields.ex[1:mx,   1], west, tag, comm2d)
    
       # Send to South and receive from North
       MPI.Sendrecv!(fields.ey[mx+1,1:my], south, tag,
                     fields.ey[   1,1:my], north, tag, comm2d)
    
       err_l2 = 0.0
       time = (istep-0.5)*dt
       for j = 1:my, i = 1:mx
           th_bz = (- cos(md*pi*(xp+(i-0.5)*dx)/dimx)
                    * cos(nd*pi*(yp+(j-0.5)*dy)/dimy)
                    * cos(omega*time))
           err_l2 += (fields.bz[i,j] - th_bz)^2
       end
    
       sum_err_l2 = MPI.Reduce(err_l2, +, 0, comm2d)
       if rank == 0 
           println(sqrt(sum_err_l2))
       end
    
       #plot_fields(mesh, rank, proc, fields.bz, xp, yp, istep )
    
    
    end # next time step
    
    MPI.Barrier(comm)

end

MPI.Init()

main( 1 ) # trigger building

tbegin = MPI.Wtime()

main( 500 )

tend = MPI.Wtime()

println(" time : $(tend -tbegin) ")

MPI.Finalize()
