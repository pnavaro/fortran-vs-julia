program Maxwell_Yee_2d

use commun
use sorties
use solveur_yee

implicit none

integer :: i, j
real(8) :: tcpu, x0, y0
type(tm_mesh_fields) :: tm, th
real(8) :: time, err_l2, r2
integer :: istep, iplot

call cpu_time(tcpu)
call readin( )

allocate(th%ex(nx,ny),tm%ex(nx,ny)) 
allocate(th%ey(nx,ny),tm%ey(nx,ny)) 
allocate(th%bz(nx,ny),tm%bz(nx,ny))

dx = dimx / (nx-1d0)
dy = dimy / (ny-1d0)
dt = cfl  / sqrt (1./(dx*dx)+1./(dy*dy)) / c
nstep = floor(tfinal/dt)
write(*,*)
write(*,*) " dx = ", dx
write(*,*) " dy = ", dy
write(*,*) " dt = ", dt
write(*,*)
if( nstep > nstepmax ) nstep = nstepmax
write(*,*) " Nombre d'iteration nstep = ", nstep

time  = 0.
iplot = 0

omega = c * sqrt((md*pi/dimx)**2+(nd*pi/dimy)**2)

tm%ex(:,:) = 0.0d0
tm%ey(:,:) = 0.0d0

do i=1,nx
do j=1,ny
   tm%bz(i,j) =   - cos(md*pi*(i-0.5)*dx/dimx)    &
                  * cos(nd*pi*(j-0.5)*dy/dimy)	  &
                  * cos(omega*(-0.5*dt))
!   r2 = ((i-0.5)*dx - 0.25*dimx)**2 +  ((j-0.5)*dy - 0.5*dimy)**2
!   tm%bz(i,j) =   exp(-r2/0.02)
end do  
end do  

open(4,file="yee.mtv") 
write(4,*)"$DATA=CURVE3D"
write(4,*)"%equalscale=T"
write(4,*)"%xmin=",-0.2, " xmax = ", 1.2
write(4,*)"%ymin=",-0.2, " ymax = ", 1.2
do i = 1, nx
do j = 1, ny
   x0 = (i-1)*dx; y0 = (j-1)*dy
   write(4,"(3f7.3)") x0   , y0   , 0.0
   write(4,"(3f7.3)") x0+dx, y0   , 0.0
   write(4,"(3f7.3)") x0+dx, y0+dy, 0.0
   write(4,"(3f7.3)") x0   , y0+dy, 0.0
   write(4,"(3f7.3)") x0   , y0   , 0.0
   write(4,*) 
end do
end do
write(4,*)"$DATA=CONTOUR" 
write(4,*)"%contstyle=",3 
write(4,*)"%nx=",nx 
write(4,*)"%ny=",ny 
write(4,*)"%xmin=",0.0
write(4,*)"%ymin=",0.0
write(4,*)"%xmax=",dimx 
write(4,*)"%ymax=",dimy
write(4,*)"%nsteps=",10 
do i = 1,nx 
   write(4,*)(tm%bz(i,j),j=1,ny) 
end do
write(4,*)"$END"
close(4)


do istep = 1, nstep !*** Loop over time

   !*** Calcul de B(n+1/2) sur les pts interieurs   
   call faraday(tm, 1, nx, 1, ny)   !Calcul de B(n-1/2)--> B(n+1/2)

   time = time + 0.5*dt

   do i=1,nx
   do j=1,ny
      th%bz(i,j) =   - cos(md*pi*(i-0.5)*dx/dimx)    &
                     * cos(nd*pi*(j-0.5)*dy/dimy)    &
                     * cos(omega*time)
   end do  
   end do  

   call cl_periodiques(tm, 1, nx, 1, ny)
   !Conducteur parfait 
   !call cl_condparfait(tm, 0, nx+1, 0, ny+1,'N')
   !call cl_condparfait(tm, 0, nx+1, 0, ny+1,'S')
   !call cl_condparfait(tm, 0, nx+1, 0, ny+1,'E')
   !call cl_condparfait(tm, 0, nx+1, 0, ny+1,'W')
   !Absorbantes 
   !call silver_muller(tm, 0, nx+1, 0, ny+1,'N')
   !call silver_muller(tm, 0, nx+1, 0, ny+1,'S')
   !call silver_muller(tm, 0, nx+1, 0, ny+1,'E')
   !call silver_muller(tm, 0, nx+1, 0, ny+1,'W')

   !*** Calcul de E(t=n+1) sur les pts interieurs
   call ampere_maxwell(tm, 1, nx, 1, ny) 

   time = time + 0.5*dt

   do i=1,nx
   do j=1,ny

      th%ex(i,j) = + (csq*nd*pi)/(omega*dimy)   	&
                    * cos(md*pi*(i-0.5)*dx/dimx) 	&
                    * sin(nd*pi*(j-1.0)*dy/dimy) 	&
                    * sin(omega*time)

      th%ey(i,j) = -(csq*md*pi)/(omega*dimx)   	&
                    * sin(md*pi*(i-1.0)*dx/dimx) 	&
                    * cos(nd*pi*(j-0.5)*dy/dimy) 	&
                    * sin(omega*time)

   end do  
   end do  

   !-----------------------------------------------------------!
   !*** Sorties graphiques (les champs sont connus au temps n) ! 
   !*** pour le solveur de MAXWELL                             !
   !-----------------------------------------------------------!

   if ( istep==1 .or. mod(istep,idiag) == 0.0) then
      iplot = iplot + 1
      call plot_fields(0, 1, tm, th, 1, nx, 1, ny, 0d0, 0d0, iplot, time )
      err_l2 = 0.0
      do j = 1, ny
      do i = 1, nx
         err_l2 = err_l2 + (tm%bz(i,j) - th%bz(i,j))**2
         !err_l2 = err_l2 + tm%bz(i,j)**2
      end do
      end do
      write(*,"(10x,' istep = ',I6)",advance="no") istep
      write(*,"(' time = ',e15.3,' sec')",advance="no") time
      write(*,"(' erreur L2 = ',g10.5)") sqrt(err_l2)
   end if

end do ! next time step

call cpu_time(tcpu)
write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") tcpu

stop

1000 format(19f8.4)
1001 format(19i8)

end program Maxwell_Yee_2d
