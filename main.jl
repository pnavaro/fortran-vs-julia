using MPI

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
idiag  = 1	# output period
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

#coords = MPI.Cart_coords(comm2d, rank)
coords = MPI.Cart_coords(comm2d, ndims)

nxp, nyp = dims

mx = nx ÷ nxp
dx = dimx / mx / nxp
my = ny ÷ nyp
dy = dimy/ my / nyp

dt = cfl / sqrt(1/dx^2+1/dy^2) / c

nstep = floor(Int64, tfinal/dt)

nstep  = min( nstepmax, nstep)

MPI.Barrier(comm)

MPI_REAL8 = MPI.Datatype(MPI.MPI_Datatype(MPI.MPI_DOUBLE))

column_type = MPI.Type_Contiguous(mx+1, MPI_REAL8)
MPI.Type_Commit!(column_type)
row_type = MPI.Type_Vector(my+1,1,mx+1,MPI_REAL8)
MPI.Type_Commit(row_type)

#=

!column datatype

CALL MPI_TYPE_CONTIGUOUS(mx+1,MPI_REAL8,type_colonne,code)
CALL MPI_TYPE_COMMIT(type_colonne,code)
!type ligne
CALL MPI_TYPE_VECTOR(my+1,1,mx+1,MPI_REAL8,type_ligne,code)
CALL MPI_TYPE_COMMIT(type_ligne,code)

!Initialisation des champs 
tm.ex(:,:) = 0d0; tm.ey(:,:) = 0d0; tm.bz(:,:) = 0d0

xp = dble(coords(1))/nxp * dimx 
yp = dble(coords(2))/nyp * dimy 

omega = c * sqrt((md*pi/dimx)**2+(nd*pi/dimy)**2)
for j=1,my
   for i=1,mx
      tm.bz(i,j) = - cos(md*pi*(xp+(i-0.5)*dx/dimx))  &
                   * cos(nd*pi*(yp+(j-0.5)*dy/dimy))  &
                   * cos(omega*(-0.5*dt))
      !r2 = (xp+(i-0.5)*dx)**2 +  (yp+(j-0.5)*dy - 0.5*dimy)**2
      !tm.bz(i,j) =   exp(-r2/0.02)
   end  
end  

for istep = 1, nstep !*** Loop over time

   !E(n) [1:mx]*[1:my] --> B(n+1/2) [1:mx-1]*[1:my-1]
   call faraday(tm, 1, mx+1, 1, my+1)   

   time = time + 0.5*dt

   for j=1,my
   for i=1,mx
      th%bz(i,j) =   - cos(md*pi*(xp+(i-0.5)*dx/dimx))  &
                     * cos(nd*pi*(yp+(j-0.5)*dy/dimy))  &
                     * cos(omega*time)
   end  
   end  

   !for i = 1, mx
      !print"(11f7.2)", ( (i-0.5)*dx, j= 1,my)
   !end

   !Envoi au voisin N et reception du voisin S
   CALL MPI_SENDRECV(tm.bz(   1,   1),1,type_ligne,voisin(N),tag, 	&
                     tm.bz(mx+1,   1),1,type_ligne,voisin(S),tag, 	&
                     comm2d, statut, code)

   !Envoi au voisin W et reception du voisin E
   CALL MPI_SENDRECV(tm.bz(   1,   1),1,type_colonne,voisin(W),tag,	&
                     tm.bz(   1,my+1),1,type_colonne,voisin(E),tag,	&
                     comm2d, statut, code)

   !Envoi au voisin NW et reception du voisin SE
   !CALL MPI_SENDRECV(tm.bz(   1,   1),1,MPI_REAL8,voisin(NW),tag,	&
   !                  tm.bz(  mx,  my),1,MPI_REAL8,voisin(SE),tag,	&
   !                  comm2d, statut, code)

   !Bz(n+1/2) [1:mx]*[1:my] --> Ex(n+1) [1:mx]*[2:my]
   !Bz(n+1/2) [1:mx]*[1:my] --> Ey(n+1) [2:mx]*[1:my]
   call ampere_maxwell(tm, 1, mx+1, 1, my+1) 

   !Envoi au voisin E et reception du voisin W
   CALL MPI_SENDRECV(tm.ex(   1,my+1),1,type_colonne,voisin(E),tag,	&
                     tm.ex(   1,   1),1,type_colonne,voisin(W),tag,	&
                     comm2d, statut, code)

   !Envoi au voisin S et reception du voisin N
   CALL MPI_SENDRECV(tm.ey(mx+1,   1),1,type_ligne,voisin(S),tag, 	&
                     tm.ey(   1,   1),1,type_ligne,voisin(N),tag, 	&
                     comm2d, statut, code)

   time = time + 0.5*dt

   for j=1,my
   for i=1,mx
      th%ex(i,j) = + (csq*nd*pi)/(omega*dimy)   	&
                    * cos(md*pi*(xp+(i-0.5)*dx/dimx)) 	&
                    * sin(nd*pi*(yp+(j-1.0)*dy/dimy)) 	&
                    * sin(omega*time)
      th%ey(i,j) = - (csq*md*pi)/(omega*dimx)   	&
                    * sin(md*pi*(xp+(i-1.0)*dx/dimx)) 	&
                    * cos(nd*pi*(yp+(j-0.5)*dy/dimy)) 	&
                    * sin(omega*time)
   end  
   end  

   !-----------------------------------------------------------!
   !    Sorties graphiques (les champs sont connus au temps n) ! 
   !    pour le solveur de MAXWELL                             !
   !-----------------------------------------------------------!

   if ( istep==1 .or. mod(istep,idiag) == 0.0) then
      iplot = iplot + 1

      call plot_fields(rang, nproc, tm, th, 1, mx, 1, my, 	&
                       xp, yp, iplot, time )
      err_l2 = 0.0
      for j = 1, my
      for i = 1, mx
         err_l2 = err_l2 + (tm.bz(i,j) - th%bz(i,j))**2
      end
      end

      call MPI_REDUCE (err_l2,sum_l2,1,MPI_REAL8,MPI_SUM,0,comm2d,code)

      if (rang == 0) then
         open(17,file="thf.dat",position="append")
         if (istep==1) rewind(17)
         write(17,*) time, err_l2
         close(17)
         write(*,"(10x,' istep = ',I6)",advance="no") istep
         write(*,"(' time = ',e17.3,' ns')",advance="no") time*1e09
         write(*,"(' erreur L2 = ',g10.3)") sqrt(sum_l2)
      end

      !call write_domains(rang, xp, yp, dx, dy, mx, my, tm, iplot)
      !if (rang == 0) call write_master(nproc, iplot)

   end

end ! next time step

call FLUSH(6)
call MPI_BARRIER(MPI_COMM_WORLD, code)
tcpu2 = MPI_WTIME()
if (rang == 0) &
   write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") (tcpu2-tcpu1)*nproc

=#
MPI.Barrier(comm)
MPI.Finalize()
