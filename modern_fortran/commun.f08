module commun
   use, intrinsic:: iso_fortran_env, only: dp => real64, ilp => int32
   implicit none(external)
   public

   type mesh_fields
      real(dp), dimension(:, :), pointer :: ex, ey
      real(dp), dimension(:, :), pointer :: bz
   end type mesh_fields

   real(dp), parameter :: c = 1.0_dp, csq = c**2
   real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)

   integer(ilp) :: md, nd
   integer(ilp) :: nx, ny
   integer(ilp) :: nstep, nstepmax
   integer(ilp) :: idiag

   real(dp) :: dt, dx, dy
   real(dp) :: dimx, dimy
   real(dp) :: cfl
   real(dp) :: tfinal
   real(dp) :: omega

contains

   subroutine readin()
      implicit none(external)
      integer :: io
      namelist /donnees/ cfl, & !nbre de Courant
         tfinal, &              !duree maxi
         nstepmax, &            !nbre d'iterations maxi
         idiag, &               !frequence des diagnostics
         md, nd, &              !nombre d'onde de la solution initiale
         nx, ny, &              !nombre de points dans les deux directions
         dimx, dimy             !dimensions du domaine

      open (newunit=io, file="input_data.nml", status='old')
      read (io, donnees)
      close (io)

      write (*, "(a,g12.3)") " Width dimx              = ", dimx
      write (*, "(a,g12.3)") " Length dimy             = ", dimy
      write (*, "(a,g12.3)") " Final time              = ", tfinal
      write (*, "(a,g12.3)") " Number nx               = ", nx
      write (*, "(a,g12.3)") " Number ny               = ", ny
      write (*, "(a,g12.3)") " CFL                     = ", cfl
      write (*, "(a,g12.3)") " Diagnostic interval     = ", idiag
   end subroutine readin

end module commun
