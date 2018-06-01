program Maxwell_Yee_2d
use commun, only: nx, ny, tm_mesh_fields, pi, csq, nstep, omega, &
		  dimx, dimy, dx, dy, dt, c, md, nd, idiag, 	 &
		  readin, tfinal, cfl, nstepmax
use sorties, only: plot_fields
use solveur_yee, only: faraday, ampere_maxwell
#ifdef HAVE_SILO
use silo_module,only: write_master, write_domains
#endif

implicit none
include 'mpif.h'

type(tm_mesh_fields) :: tm, th
integer :: i, j, iproc, istep, iplot
real(8) :: tcpu1, tcpu2, time, err_l2, sum_l2, xp, yp
real(8) :: x1, y1, r2, dtloc
integer,dimension(MPI_STATUS_SIZE) :: statut
integer                  :: rang, nproc, code,comm2d
integer,parameter        :: tag=1111
integer,dimension(8)     :: voisin
integer,parameter        :: N =1, S =2, W =3, E =4
integer,parameter        :: NW=5, SW=6, NE=7, SE=8
integer,parameter        :: ndims = 2
integer,dimension(ndims) :: dims, coords
integer,dimension(ndims) :: coordsse, coordssw, coordsne, coordsnw
logical                  :: reorder
logical,dimension(ndims) :: periods
integer                  :: nxp, nyp, mx, my
integer                  :: type_ligne, type_colonne

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

if (voisin(N) /= MPI_PROC_NULL) then
   coordsnw(1) = coords(1)-1; coordsnw(2) = coords(2)-1
   coordsne(1) = coords(1)-1; coordsne(2) = coords(2)+1
   CALL MPI_CART_RANK(comm2d,coordsnw,voisin(NW),code)
   CALL MPI_CART_RANK(comm2d,coordsne,voisin(NE),code)
end if

if (voisin(S) /= MPI_PROC_NULL) then
   coordssw(1) = coords(1)+1; coordssw(2) = coords(2)-1
   coordsse(1) = coords(1)+1; coordsse(2) = coords(2)+1
   CALL MPI_CART_RANK(comm2d,coordssw,voisin(SW),code)
   CALL MPI_CART_RANK(comm2d,coordsse,voisin(SE),code)
end if

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

mx = nx/nxp; dx = dimx/mx/nxp
my = ny/nyp; dy = dimy/my/nyp

call MPI_BCAST(cfl,   1,MPI_REAL8,   0,comm2d,code)
dtloc = cfl  / sqrt (1./(dx*dx)+1./(dy*dy)) / c
call MPI_ALLREDUCE(dtloc,dt,1,MPI_REAL8,MPI_MAX,comm2d,code)

call MPI_BCAST(tfinal,1,MPI_REAL8,   0,comm2d,code)
call MPI_BCAST(nstepmax,1,MPI_INTEGER,   0,comm2d,code)
nstep = floor(tfinal/dt)
if( nstep > nstepmax ) nstep = nstepmax

allocate(tm%ex(mx+1,my+1)); tm%ex(:,:) = 0.0d0
allocate(tm%ey(mx+1,my+1)); tm%ey(:,:) = 0.0d0
allocate(tm%bz(mx+1,my+1)); tm%bz(:,:) = 0.0d0
allocate(th%ex(mx,my),th%ey(mx,my),th%bz(mx,my))

xp = dble(coords(1))/nxp * dimx 
yp = dble(coords(2))/nyp * dimy 

if ( rang == 0) then

   open(4,file="yee.mtv")
   write(4,*)"$DATA=CURVE3D"
   write(4,*)"%equalscale=T"
   write(4,*)"%xmin=",-0.1*dimx, " xmax = ", 1.1*dimx
   write(4,*)"%ymin=",-0.1*dimy, " ymax = ", 1.1*dimy
   write(4,"(3f7.3)") 0.   , 0.   , 0.0
   write(4,"(3f7.3)") dimx , 0.   , 0.0
   write(4,"(3f7.3)") dimx , dimy , 0.0
   write(4,"(3f7.3)") 0.   , dimy , 0.0
   write(4,"(3f7.3)") 0.   , 0.   , 0.0
   write(4,*)

   do iproc = 0, nproc-1
      if (iproc > 0) then
         call MPI_RECV(xp,1,MPI_REAL8,iproc,tag,comm2d,statut,code)
         call MPI_RECV(yp,1,MPI_REAL8,iproc,tag,comm2d,statut,code)
      end if
      do i = 1, mx
      do j = 1, my
         x1 = xp + (i-1)*dx; y1 = yp + (j-1)*dy
         write(4,"('%linecolor=',i2)")iproc+1
         write(4,"(3f7.3)") x1   , y1   , tm%bz(i,j)
         write(4,"(3f7.3)") x1+dx, y1   , tm%bz(i,j)
         write(4,"(3f7.3)") x1+dx, y1+dy, tm%bz(i,j)
         write(4,"(3f7.3)") x1   , y1+dy, tm%bz(i,j)
         write(4,"(3f7.3)") x1   , y1   , tm%bz(i,j)
         write(4,*)
      end do
      end do
      !write(4,*)"@text x1=",xp+.25*dimx,' y1=',yp+.25*dimy,    &
      !          " z1=0.0 linelabel='proc:",iproc,"'"
   end do

   write(4,"('$END')")
   close(4)

else

   call MPI_SEND(xp, 1, MPI_REAL8, 0, tag, comm2d, code)
   call MPI_SEND(yp, 1, MPI_REAL8, 0, tag, comm2d, code)

end if

print*, "proc= ",rang,": ",mx,"x",my," mes coords sont ",coords(:)
print*, "proc= ",rang,": (xp,yp) = (",sngl(xp),";", sngl(yp), ")"
print*, "proc= ",rang,": mes voisins sont ", voisin(1:4)
print*, "proc= ",rang,": Nombre d'iteration nstep = ", nstep
print*, "proc= ",rang,": dx = ", sngl(dx), " dy = ", sngl(dy), " dt = ", sngl(dt)

CALL FLUSH(6)
call MPI_BARRIER(MPI_COMM_WORLD, code)

time  = 0.
iplot = 0

!type colonne
CALL MPI_TYPE_CONTIGUOUS(mx+1,MPI_REAL8,type_colonne,code)
CALL MPI_TYPE_COMMIT(type_colonne,code)
!type ligne
CALL MPI_TYPE_VECTOR(my+1,1,mx+1,MPI_REAL8,type_ligne,code)
CALL MPI_TYPE_COMMIT(type_ligne,code)

!Initialisation des champs 
tm%ex(:,:) = 0d0; tm%ey(:,:) = 0d0; tm%bz(:,:) = 0d0

xp = dble(coords(1))/nxp * dimx 
yp = dble(coords(2))/nyp * dimy 

omega = c * sqrt((md*pi/dimx)**2+(nd*pi/dimy)**2)
do j=1,my
   do i=1,mx
      tm%bz(i,j) = - cos(md*pi*(xp+(i-0.5)*dx/dimx))  &
                   * cos(nd*pi*(yp+(j-0.5)*dy/dimy))  &
                   * cos(omega*(-0.5*dt))
      !r2 = (xp+(i-0.5)*dx)**2 +  (yp+(j-0.5)*dy - 0.5*dimy)**2
      !tm%bz(i,j) =   exp(-r2/0.02)
   end do  
end do  

do istep = 1, nstep !*** Loop over time

   !E(n) [1:mx]*[1:my] --> B(n+1/2) [1:mx-1]*[1:my-1]
   call faraday(tm, 1, mx+1, 1, my+1)   

   time = time + 0.5*dt

   do j=1,my
   do i=1,mx
      th%bz(i,j) =   - cos(md*pi*(xp+(i-0.5)*dx/dimx))  &
                     * cos(nd*pi*(yp+(j-0.5)*dy/dimy))  &
                     * cos(omega*time)
   end do  
   end do  

   !do i = 1, mx
      !print"(11f7.2)", ( (i-0.5)*dx, j= 1,my)
   !end do

   !Envoi au voisin N et reception du voisin S
   CALL MPI_SENDRECV(tm%bz(   1,   1),1,type_ligne,voisin(N),tag, 	&
                     tm%bz(mx+1,   1),1,type_ligne,voisin(S),tag, 	&
                     comm2d, statut, code)

   !Envoi au voisin W et reception du voisin E
   CALL MPI_SENDRECV(tm%bz(   1,   1),1,type_colonne,voisin(W),tag,	&
                     tm%bz(   1,my+1),1,type_colonne,voisin(E),tag,	&
                     comm2d, statut, code)

   !Envoi au voisin NW et reception du voisin SE
   !CALL MPI_SENDRECV(tm%bz(   1,   1),1,MPI_REAL8,voisin(NW),tag,	&
   !                  tm%bz(  mx,  my),1,MPI_REAL8,voisin(SE),tag,	&
   !                  comm2d, statut, code)

   !Bz(n+1/2) [1:mx]*[1:my] --> Ex(n+1) [1:mx]*[2:my]
   !Bz(n+1/2) [1:mx]*[1:my] --> Ey(n+1) [2:mx]*[1:my]
   call ampere_maxwell(tm, 1, mx+1, 1, my+1) 

   !Envoi au voisin E et reception du voisin W
   CALL MPI_SENDRECV(tm%ex(   1,my+1),1,type_colonne,voisin(E),tag,	&
                     tm%ex(   1,   1),1,type_colonne,voisin(W),tag,	&
                     comm2d, statut, code)

   !Envoi au voisin S et reception du voisin N
   CALL MPI_SENDRECV(tm%ey(mx+1,   1),1,type_ligne,voisin(S),tag, 	&
                     tm%ey(   1,   1),1,type_ligne,voisin(N),tag, 	&
                     comm2d, statut, code)

   time = time + 0.5*dt

   do j=1,my
   do i=1,mx
      th%ex(i,j) = + (csq*nd*pi)/(omega*dimy)   	&
                    * cos(md*pi*(xp+(i-0.5)*dx/dimx)) 	&
                    * sin(nd*pi*(yp+(j-1.0)*dy/dimy)) 	&
                    * sin(omega*time)
      th%ey(i,j) = - (csq*md*pi)/(omega*dimx)   	&
                    * sin(md*pi*(xp+(i-1.0)*dx/dimx)) 	&
                    * cos(nd*pi*(yp+(j-0.5)*dy/dimy)) 	&
                    * sin(omega*time)
   end do  
   end do  

   !-----------------------------------------------------------!
   !    Sorties graphiques (les champs sont connus au temps n) ! 
   !    pour le solveur de MAXWELL                             !
   !-----------------------------------------------------------!

   if ( istep==1 .or. mod(istep,idiag) == 0.0) then
      iplot = iplot + 1

      call plot_fields(rang, nproc, tm, th, 1, mx, 1, my, 	&
                       xp, yp, iplot, time )
      err_l2 = 0.0
      do j = 1, my
      do i = 1, mx
         err_l2 = err_l2 + (tm%bz(i,j) - th%bz(i,j))**2
      end do
      end do

      call MPI_REDUCE (err_l2,sum_l2,1,MPI_REAL8,MPI_SUM,0,comm2d,code)

      if (rang == 0) then
         open(17,file="thf.dat",position="append")
         if (istep==1) rewind(17)
         write(17,*) time, err_l2
         close(17)
         write(*,"(10x,' istep = ',I6)",advance="no") istep
         write(*,"(' time = ',e17.3,' ns')",advance="no") time*1e09
         write(*,"(' erreur L2 = ',g10.3)") sqrt(sum_l2)
      end if

      !call write_domains(rang, xp, yp, dx, dy, mx, my, tm, iplot)
      !if (rang == 0) call write_master(nproc, iplot)

   end if

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

end program Maxwell_Yee_2d
