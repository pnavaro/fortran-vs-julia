program maxyee_mpi

use mpi

use commun
use sorties
use solveur_yee

implicit none

type(mesh_fields) :: fields
integer :: i, j, iproc, istep, iplot
real(8) :: tcpu1, tcpu2, time, err_l2, sum_l2, xp, yp
real(8) :: x, y, dtloc
real(8) :: th_bz
integer,dimension(MPI_STATUS_SIZE) :: statut
integer                  :: rang, nproc, code, comm2d
integer,parameter        :: tag=1111
integer,dimension(8)     :: voisin
integer,parameter        :: N =1, S =2, W =3, E =4
integer,parameter        :: ndims = 2
integer,dimension(ndims) :: dims, coords
logical                  :: reorder
logical,dimension(ndims) :: periods
integer                  :: nxp, nyp, mx, my

!Initialisation de MPI
CALL MPI_INIT(code)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,code)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
CALL MPI_BARRIER(MPI_COMM_WORLD, code)

tcpu1 = MPI_WTIME()

!Nombre de processus suivant x et y
dims(:) = 0
CALL MPI_DIMS_CREATE(nproc,ndims,dims,code)
nxp = dims(1)
nyp = dims(2)

!Creation de la grille 2D periodique en x et y
periods(1) = .true.
periods(2) = .true.
reorder    = .true.

CALL MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims,periods,reorder,comm2d,code)

!Initialisation du tableau voisin 
voisin(:) = MPI_PROC_NULL

!Recherche de mes voisins Sud et Nord
CALL MPI_CART_SHIFT(comm2d,0,1,voisin(N),voisin(S),code)
!Recherche de mes voisins Ouest et Est
CALL MPI_CART_SHIFT(comm2d,1,1,voisin(W),voisin(E),code)
!Connaitre mes coordonnees dans la topologie
CALL MPI_COMM_RANK(comm2d,rang,code)
CALL MPI_CART_COORDS(comm2d,rang,ndims,coords,code)

if (rang == 0) call readin( )

call MPI_BCAST(nx,    1,MPI_INTEGER, 0,comm2d,code)
call MPI_BCAST(ny,    1,MPI_INTEGER, 0,comm2d,code)
call MPI_BCAST(nd,    1,MPI_INTEGER, 0,comm2d,code)
call MPI_BCAST(md,    1,MPI_INTEGER, 0,comm2d,code)
call MPI_BCAST(idiag, 1,MPI_INTEGER, 0,comm2d,code)
call MPI_BCAST(dimx,  1,MPI_REAL8,   0,comm2d,code)
call MPI_BCAST(dimy,  1,MPI_REAL8,   0,comm2d,code)
call MPI_BCAST(pi,    1,MPI_REAL8,   0,comm2d,code)
call MPI_BCAST(csq,   1,MPI_REAL8,   0,comm2d,code)
call MPI_BCAST(c,     1,MPI_REAL8,   0,comm2d,code)

mx = nx/nxp; dx = dimx/nx
my = ny/nyp; dy = dimy/ny

call MPI_BCAST(cfl,   1,MPI_REAL8,   0,comm2d,code)
dtloc = cfl  / sqrt (1./(dx*dx)+1./(dy*dy)) / c
call MPI_ALLREDUCE(dtloc,dt,1,MPI_REAL8,MPI_MAX,comm2d,code)

call MPI_BCAST(tfinal,1,MPI_REAL8,   0,comm2d,code)
call MPI_BCAST(nstepmax,1,MPI_INTEGER,   0,comm2d,code)
nstep = floor(tfinal/dt)
if( nstep > nstepmax ) nstep = nstepmax

allocate(fields%ex(mx,my+1)); fields%ex(:,:) = 0.0d0
allocate(fields%ey(mx+1,my)); fields%ey(:,:) = 0.0d0
allocate(fields%bz(mx+1,my+1)); fields%bz(:,:) = 0.0d0

xp = dble(coords(1))/nxp * dimx 
yp = dble(coords(2))/nyp * dimy 

print*, "proc= ",rang,": ",mx,"x",my," mes coords sont ",coords(:)
print*, "proc= ",rang,": (xp,yp) = (",sngl(xp),";", sngl(yp), ")"
print*, "proc= ",rang,": mes voisins sont ", voisin(1:4)
print*, "proc= ",rang,": Nombre d'iteration nstep = ", nstep
print*, "proc= ",rang,": dx = ", sngl(dx), " dy = ", sngl(dy), " dt = ", sngl(dt)

CALL FLUSH(6)
call MPI_BARRIER(MPI_COMM_WORLD, code)

time  = 0.
iplot = 0

!Initialisation des champs 
fields%ex(:,:) = 0d0; fields%ey(:,:) = 0d0; fields%bz(:,:) = 0d0

xp = coords(1) * dimx /nxp
yp = coords(2) * dimy /nyp

omega = c * sqrt((md*pi/dimx)**2+(nd*pi/dimy)**2)
do j=1,my
   do i=1,mx
      x = xp + (i-.5) * dx
      y = yp + (j-.5) * dy
      fields%bz(i,j) = - cos(md*pi*x/dimx)  &
                       * cos(nd*pi*y/dimy)  &
                       * cos(omega*(-0.5*dt))
   end do  
end do  

do istep = 1, nstep !*** Loop over time

   !E(n) [1:mx]*[1:my] --> B(n+1/2) [1:mx-1]*[1:my-1]
   call faraday(fields, 1, mx, 1, my)   

   time = time + 0.5*dt

   !Envoi au voisin N et reception du voisin S
   CALL MPI_SENDRECV(fields%bz(   1,  :), my+1,MPI_REAL8,voisin(N),tag, 	&
                     fields%bz(mx+1,  :), my+1,MPI_REAL8,voisin(S),tag, 	&
                     comm2d, statut, code)

   !Envoi au voisin W et reception du voisin E
   CALL MPI_SENDRECV(fields%bz(   :,   1), mx+1,MPI_REAL8,voisin(W),tag,	&
                     fields%bz(   :,my+1), mx+1,MPI_REAL8,voisin(E),tag,	&
                     comm2d, statut, code)

   !Bz(n+1/2) [1:mx]*[1:my] --> Ex(n+1) [1:mx]*[2:my]
   !Bz(n+1/2) [1:mx]*[1:my] --> Ey(n+1) [2:mx]*[1:my]
   call ampere_maxwell(fields, 1, mx, 1, my) 

   !Envoi au voisin E et reception du voisin W
   CALL MPI_SENDRECV(fields%ex(   :,my+1), mx, MPI_REAL8,voisin(E),tag,	&
                     fields%ex(   :,   1), mx, MPI_REAL8,voisin(W),tag,	&
                     comm2d, statut, code)

   !Envoi au voisin S et reception du voisin N
   CALL MPI_SENDRECV(fields%ey(mx+1,   :), my, MPI_REAL8,voisin(S),tag, 	&
                     fields%ey(   1,   :), my, MPI_REAL8,voisin(N),tag, 	&
                     comm2d, statut, code)


   !-----------------------------------------------------------!
   !    Sorties graphiques (les champs sont connus au temps n) ! 
   !    pour le solveur de MAXWELL                             !
   !-----------------------------------------------------------!
   err_l2 = 0.0
   do j = 1, my
       do i = 1, mx
          x = xp + (i-.5) * dx
          y = yp + (j-.5) * dy
          th_bz = - cos(md*pi*x/dimx)  &
                  * cos(nd*pi*y/dimy)  &
                  * cos(omega*time)
          err_l2 = err_l2 + (fields%bz(i,j) - th_bz)**2
       end do
   end do

   call MPI_REDUCE (err_l2,sum_l2,1,MPI_REAL8,MPI_SUM,0,comm2d,code)

   if (rang == 0) then
       write(*,"(10x,' istep = ',I6)",advance="no") istep
       write(*,"(' time = ',e17.3,' ns')",advance="no") time
       write(*,"(' erreur L2 = ',g10.3)") sqrt(sum_l2)
   end if

   if ( istep==1 .or. mod(istep,idiag) == 0.0) then
      iplot = iplot + 1
      call plot_fields(rang, nproc, fields, 1, mx, 1, my, xp, yp, iplot, time )
   end if

   time = time + 0.5*dt

end do ! next time step

call FLUSH(6)
call MPI_BARRIER(MPI_COMM_WORLD, code)
tcpu2 = MPI_WTIME()
if (rang == 0) &
   write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") (tcpu2-tcpu1)*nproc

!Desactivation de MPI
CALL MPI_FINALIZE(code)
stop

1000 format(10f8.4)

end program maxyee_mpi
