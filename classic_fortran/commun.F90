module commun

   type mesh_fields
      real(8), dimension(:, :), pointer :: ex, ey
      real(8), dimension(:, :), pointer :: bz
   end type mesh_fields

   real(8) :: c, csq
   real(8) :: pi

   integer :: md, nd
   integer :: nx, ny
   integer :: nstep, nstepmax
   integer :: idiag

   integer, private :: i, j

   real(8) :: dt, dx, dy
   real(8) :: dimx, dimy
   real(8) :: cfl
   real(8) :: tfinal
   real(8) :: omega

contains

   subroutine readin()

      implicit none

      namelist /donnees/ cfl, &            !nbre de Courant
         tfinal, &            !duree maxi
         nstepmax, &        !nbre d'iterations maxi
         idiag, &        !frequence des diagnostics
         md, nd, &        !nombre d'onde de la solution initiale
         nx, ny, &        !nombre de points dans les deux directions
         dimx, dimy                !dimensions du domaine

!***Initialisation  des valeurs pas default

      pi = 4.*atan(1.)

      c = 1d0                !celerite de la lumiere
      csq = c*c

      open (10, file="input_data.nml", status='old')
      read (10, donnees)
      close (10)

      write (*, "(a,g12.3)") " largeur dimx              = ", dimx
      write (*, "(a,g12.3)") " longueur dimy             = ", dimy
      write (*, "(a,g12.3)") " temps                     = ", tfinal
      write (*, "(a,g12.3)") " nombre nx                 = ", nx
      write (*, "(a,g12.3)") " nombre ny                 = ", ny
      write (*, "(a,g12.3)") " nombre de Courant         = ", cfl
      write (*, "(a,g12.3)") " frequence des diagnostics = ", idiag

   end subroutine readin

end module commun
