program Maxwell_Yee_2d

use commun
use sorties
use solveur_yee

implicit none

integer :: i, j
real(8) :: tcpu, x0, y0
real(8) :: time, err_l2, th_bz
integer :: istep, iplot

type(mesh_fields) :: fields, th

call cpu_time(tcpu)
call readin( )

allocate(fields%ex(nx,ny+1)) 
allocate(fields%ey(nx+1,ny)) 
allocate(fields%bz(nx,ny))

dx = dimx / real(nx,kind=8)
dy = dimy / real(ny,kind=8)
dt = cfl  / sqrt (1./(dx*dx)+1./(dy*dy)) / c
nstep = floor(tfinal/dt)
write(*,*)
write(*,*) " dx = ", dx
write(*,*) " dy = ", dy
write(*,*) " dt = ", dt
write(*,*)
write(*,*) nx * dx
if( nstep > nstepmax ) nstep = nstepmax
write(*,*) " Nombre d'iteration nstep = ", nstep

time  = 0.
iplot = 0

omega = c * sqrt((md*pi/dimx)**2+(nd*pi/dimy)**2)

fields%ex(:,:) = 0.0d0
fields%ey(:,:) = 0.0d0

do j=1,ny
do i=1,nx
   fields%bz(i,j) =   - cos(md*pi*(i-0.5)*dx/dimx)    &
                  * cos(nd*pi*(j-0.5)*dy/dimy)	  &
                  * cos(omega*(-0.5*dt))
end do  
end do  

do istep = 1, nstep !*** Loop over time

   !*** Calcul de B(n+1/2) sur les pts interieurs   
   call faraday(fields, 1, nx, 1, ny)   !Calcul de B(n-1/2)--> B(n+1/2)
 
   time = time + 0.5*dt

   call cl_periodiques(fields, 1, nx, 1, ny)

   !*** Calcul de E(t=n+1) sur les pts interieurs
   call ampere_maxwell(fields, 1, nx, 1, ny) 

   !*** diagnostics ***
   err_l2 = 0.0
   do j = 1, ny
   do i = 1, nx
      th_bz = - cos(md*pi*(i-0.5)*dx/dimx)    &
              * cos(nd*pi*(j-0.5)*dy/dimy)    &
              * cos(omega*time)
      err_l2 = err_l2 + (fields%bz(i,j) - th_bz)**2
   end do
   end do
   write(*,"(10x,' istep = ',I6)",advance="no") istep
   write(*,"(' time = ',e15.3,' sec')",advance="no") time
   write(*,"(' erreur L2 = ',g10.5)") sqrt(err_l2)

   if (mod(istep,idiag) == 0.0) then
      iplot = iplot + 1
      call plot_fields(0, 1, fields, 1, nx, 1, ny, 0d0, 0d0, iplot, time )
   end if

   time = time + 0.5*dt

end do ! next time step

call cpu_time(tcpu)
write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") tcpu

stop

1000 format(19f8.4)
1001 format(19i8)

end program Maxwell_Yee_2d
