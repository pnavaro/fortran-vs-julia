module solveur_yee
   use commun, only: dp, ilp, mesh_fields, dx, dy, csq, c, dt
   implicit none(external)
   private

   interface
      pure module subroutine faraday(tm, ix, jx, iy, jy)
         implicit none(external)
         type(mesh_fields), intent(inout) :: tm
         integer(ilp), intent(in) :: ix, jx, iy, jy
      end subroutine faraday
      pure module subroutine ampere_maxwell(tm, ix, jx, iy, jy)
         implicit none(external)
         type(mesh_fields), intent(inout) :: tm
         integer(ilp), intent(in) :: ix, jx, iy, jy
      end subroutine ampere_maxwell
      pure module subroutine cl_periodiques(tm, ix, jx, iy, jy)
         implicit none(external)
         type(mesh_fields), intent(inout) :: tm
         integer(ilp), intent(in) :: ix, jx, iy, jy
      end subroutine cl_periodiques
   end interface
   public :: faraday, ampere_maxwell, cl_periodiques

contains

   module procedure faraday
   integer(ilp) :: i, j
   real(dp) :: dex_dy, dey_dx
   !*** On utilise l'equation de Faraday sur un demi pas
   !*** de temps pour le calcul du champ magnetique  Bz
   !*** a l'instant n puis n+1/2
   !Ex[1:nx,1:ny+1,]*Ey[1:nx+1,1:ny] -> Bz[1:nx,1:ny]
   do concurrent(i=ix:jx, j=iy:jy) local(dex_dy, dey_dx)
      dex_dy = (tm%ex(i, j + 1) - tm%ex(i, j))/dy
      dey_dx = (tm%ey(i + 1, j) - tm%ey(i, j))/dx
      tm%bz(i, j) = tm%bz(i, j) + dt*(dex_dy - dey_dx)
   end do
   end procedure faraday

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module procedure ampere_maxwell
   integer(ilp) :: i, j
   real(dp) :: dbz_dy, dbz_dx
   !*** Calcul du champ electrique E au temps n+1
   !*** sur les points internes du maillage
   !*** Ex aux points (i+1/2,j)
   !*** Ey aux points (i,j+1/2)
   do concurrent(i=ix:jx, j=iy + 1:jy) local(dbz_dy)
      dbz_dy = (tm%bz(i, j) - tm%bz(i, j - 1))/dy
      tm%ex(i, j) = tm%ex(i, j) + dt*csq*dbz_dy
   end do
   do concurrent(i=ix + 1:jx, j=iy:jy) local(dbz_dx)
      dbz_dx = (tm%bz(i, j) - tm%bz(i - 1, j))/dx
      tm%ey(i, j) = tm%ey(i, j) - dt*csq*dbz_dx
   end do
   end procedure ampere_maxwell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module procedure cl_periodiques
   integer(ilp) :: i, j
   real(dp) :: dbz_dy, dbz_dx
   do concurrent(i=ix:jx) local(dbz_dy)
      dbz_dy = (tm%bz(i, iy) - tm%bz(i, jy))/dy
      tm%ex(i, iy) = tm%ex(i, iy) + dt*csq*dbz_dy
      tm%ex(i, jy + 1) = tm%ex(i, iy)
   end do
   do concurrent(j=iy:jy) local(dbz_dx)
      dbz_dx = (tm%bz(ix, j) - tm%bz(jx, j))/dx
      tm%ey(ix, j) = tm%ey(ix, j) - dt*csq*dbz_dx
      tm%ey(jx + 1, j) = tm%ey(ix, j)
   end do
   end procedure cl_periodiques

end module solveur_yee
