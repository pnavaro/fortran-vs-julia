using MPI
using Printf

struct MeshFields

   ex :: Array{Float64, 2}
   ey :: Array{Float64, 2}
   bz :: Array{Float64, 2}

    function MeshFields( mx, my)
        ex = zeros(Float64, (mx+1,my+1))
        ey = zeros(Float64, (mx+1,my+1))
        bz = zeros(Float64, (mx+1,my+1))
        new( ex, ey, bz)
    end

end 

c = 1.0 # speed of light 
csq = c * c

cfl    = 0.2 	# Courant-Friedrich-Levy
tfinal = 1.	# final time
nstepmax = 500  # max steps
md = 2		# md : wave number x
nd = 2		# nd : wave number y
nx = 100	# x number of points
ny = 100	# y number of points
dimx = 1.0	# width
dimy = 1.0	# height


function faraday!( tm, ix, jx, iy, jy )

    for j=iy:jy, i=ix:jx
       dex_dy     = (tm.ex[i,j+1]-tm.ex[i,j]) / dy
       dey_dx     = (tm.ey[i+1,j]-tm.ey[i,j]) / dx
       tm.bz[i,j] = tm.bz[i,j] + dt * (dex_dy - dey_dx)
    end

end

function ampere_maxwell!( tm, ix, jx, iy, jy )

    for j=iy+1:jy, i=ix:jx
       dbz_dy = (tm.bz[i,j]-tm.bz[i,j-1]) / dy
       tm.ex[i,j] = tm.ex[i,j] + dt*csq*dbz_dy 
    end

    for j=iy:jy, i=ix+1:jx
       dbz_dx = (tm.bz[i,j]-tm.bz[i-1,j]) / dx
       tm.ey[i,j] = tm.ey[i,j] - dt*csq*dbz_dx 
    end

end 

function periodic_bc(tm, ix, jx, iy, jy)

    for i = ix:jx
       dbz_dy = (tm.bz[i,iy]-tm.bz[i,jy]) / dy
       tm.ex[i,iy] = tm.ex[i,iy] + dt*csq*dbz_dy 
       tm.ex[i,jy+1] = tm.ex[i,iy] 
    end
    
    for j = iy:jy
       dbz_dx = (tm.bz[ix,j]-tm.bz[jx,j]) / dx
       tm.ey[ix,j] = tm.ey[ix,j] - dt*csq*dbz_dx 
       tm.ey[jx+1,j] = tm.ey[ix,j]
    end

end 


function plot_fields(rank, proc, field, ix, jx, iy, jy, xp, yp, iplot )

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

MPI.Init()
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

disp = 1

north, south = MPI.Cart_shift(comm2d,0,1)
west,  east  = MPI.Cart_shift(comm2d,1,1)

coords = MPI.Cart_coords(comm2d)

nxp, nyp = dims

mx = nx ÷ nxp
dx = dimx / mx / nxp
my = ny ÷ nyp
dy = dimy / my / nyp

dt = cfl / sqrt(1/dx^2+1/dy^2) / c

nstep = floor(Int64, tfinal/dt)

nstep  = min( nstepmax, nstep)

MPI.Barrier(comm)

xp = coords[1]/nxp * dimx 
yp = coords[2]/nyp * dimy 

tm = MeshFields(mx, my)

omega = c * sqrt((md*pi/dimx)^2+(nd*pi/dimy)^2)
for j=1:my, i=1:mx
    tm.bz[i,j] = (- cos(md*pi*(xp+(i-0.5)*dx/dimx)) 
                  * cos(nd*pi*(yp+(j-0.5)*dy/dimy))
                  * cos(omega*(-0.5*dt)) )
end  

tag = 1111

for istep = 1:nstep #*** Loop over time

   # E(n) [1:mx]*[1:my] --> B(n+1/2) [1:mx-1]*[1:my-1]
   faraday!(tm, 1, mx, 1, my)   

   # Send to North  and receive from South
   MPI.Sendrecv!(tm.bz[ 1, 1:my+1], north, tag,
                 tm.bz[mx, 1:my+1], south, tag, comm2d)

   # Send to West and receive from East
   MPI.Sendrecv!(tm.bz[1:mx+1, 1], west, tag,
                 tm.bz[1:mx+1,my], east, tag, comm2d)

   # Bz(n+1/2) [1:mx]*[1:my] --> Ex(n+1) [1:mx]*[2:my]
   # Bz(n+1/2) [1:mx]*[1:my] --> Ey(n+1) [2:mx]*[1:my]
   ampere_maxwell!(tm, 1, mx, 1, my) 

   # Send to East and receive from West
   MPI.Sendrecv!(tm.ex[1:mx+1,my+1], east, tag,
                 tm.ex[1:mx+1,   1], west, tag, comm2d)

   # Send to South and receive from North
   MPI.Sendrecv!(tm.ey[mx+1,1:my+1], south, tag,
                 tm.ey[   1,1:my+1], north, tag, comm2d)

   err_l2 = 0.0
   time = (istep-0.5)*dt
   for j = 1:my, i = 1:mx
       th_bz = (- cos(md*pi*(xp+(i-0.5)*dx/dimx))
                * cos(nd*pi*(yp+(j-0.5)*dy/dimy))
                * cos(omega*time))
       err_l2 += (tm.bz[i,j] - th_bz)^2
   end

   sum_err_l2 = MPI.Reduce(err_l2, +, 0, comm2d)
   if rank == 0 
       println(sum_err_l2)
   end

   plot_fields(rank, proc, tm.bz, 1, mx, 1, my, xp, yp, istep )


end # next time step

MPI.Barrier(comm)
MPI.Finalize()
