program Maxwell_Yee_2d
   use, intrinsic :: iso_fortran_env, only: output_unit
   use commun
   use sorties, only: plot_fields
   use solveur_yee, only: faraday, ampere_maxwell, cl_periodiques
   implicit none(external)

   integer(ilp) :: i, j
   real(dp) :: tcpu, x0, y0
   real(dp) :: time, err_l2, th_bz
   integer(ilp) :: istep, iplot
   type(mesh_fields) :: fields, th

   call cpu_time(tcpu)
   call readin()

   allocate (fields%ex(nx, ny + 1), source=0.0_dp)
   allocate (fields%ey(nx + 1, ny), source=0.0_dp)
   allocate (fields%bz(nx, ny))

   dx = dimx/real(nx, kind=dp)
   dy = dimy/real(ny, kind=dp)
   dt = cfl/sqrt(1.0_dp/(dx*dx) + 1.0_dp/(dy*dy))/c
   nstep = floor(tfinal/dt)
   write (output_unit, *)
   write (output_unit, *) " dx = ", dx
   write (output_unit, *) " dy = ", dy
   write (output_unit, *) " dt = ", dt
   write (output_unit, *)
   write (output_unit, *) nx*dx
   nstep = min(nstep, nstepmax)
   write (output_unit, *) " Nombre d'iteration nstep = ", nstep

   time = 0.0_dp
   iplot = 0

   omega = c*sqrt((md*pi/dimx)**2 + (nd*pi/dimy)**2)

   do concurrent(i=1:nx, j=1:ny)
      fields%bz(i, j) = -cos(md*pi*(i - 0.5_dp)*dx/dimx) &
                        *cos(nd*pi*(j - 0.5_dp)*dy/dimy) &
                        *cos(omega*(-0.5_dp*dt))
   end do

   do istep = 1, nstep !*** Loop over time

      !*** Calcul de B(n+1/2) sur les pts interieurs
      call faraday(fields, 1, nx, 1, ny)   !Calcul de B(n-1/2)--> B(n+1/2)

      time = time + 0.5_dp*dt

      call cl_periodiques(fields, 1, nx, 1, ny)

      !*** Calcul de E(t=n+1) sur les pts interieurs
      call ampere_maxwell(fields, 1, nx, 1, ny)

      if (mod(istep, idiag) == 0) then
         iplot = iplot + 1
         call plot_fields(0, 1, fields, 1, nx, 1, ny, 0.0_dp, 0.0_dp, iplot, time)
      end if

      time = time + 0.5_dp*dt

   end do ! next time step

   err_l2 = 0.0_dp
   do concurrent(i=1:nx, j=1:ny) reduce(+:err_l2) local(th_bz)
      th_bz = -cos(md*pi*(i - 0.5)*dx/dimx) &
              *cos(nd*pi*(j - 0.5)*dy/dimy) &
              *cos(omega*(time - 0.5*dt))
      err_l2 = err_l2 + (fields%bz(i, j) - th_bz)**2
   end do

   call cpu_time(tcpu)

   write (output_unit, "(10x,' istep = ',I6)", advance="no") istep
   write (output_unit, "(' time = ',e15.3,' sec')", advance="no") time
   write (output_unit, "(' erreur L2 = ',g10.5)") sqrt(err_l2)
   write (output_unit, "(//10x,' Temps CPU = ', G15.3, ' sec' )") tcpu

1000 format(19f8.4)
1001 format(19i8)

end program Maxwell_Yee_2d
