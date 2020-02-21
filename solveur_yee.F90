module solveur_yee

use commun, only: mesh_fields, dx, dy, csq , c, dt

implicit none

integer, private :: i, j
real(8), private :: dex_dx, dey_dy
real(8), private :: dex_dy, dey_dx
real(8), private :: dbz_dx, dbz_dy

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine faraday( tm, ix, jx, iy, jy )

type( mesh_fields ) :: tm
integer, intent(in) :: ix, jx, iy, jy

!*** On utilise l'equation de Faraday sur un demi pas
!*** de temps pour le calcul du champ magnetique  Bz 
!*** a l'instant n puis n+1/2 
!Ex[1:nx,1:ny+1,]*Ey[1:nx+1,1:ny] -> Bz[1:nx,1:ny]

do j=iy,jy
do i=ix,jx
   dex_dy     = (tm%ex(i,j+1)-tm%ex(i,j)) / dy
   dey_dx     = (tm%ey(i+1,j)-tm%ey(i,j)) / dx
   tm%bz(i,j) = tm%bz(i,j) + dt * (dex_dy - dey_dx)
end do
end do

end subroutine faraday

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ampere_maxwell( tm, ix, jx, iy, jy )

type( mesh_fields ) :: tm
integer, intent(in) :: ix, jx, iy, jy

!*** Calcul du champ electrique E au temps n+1
!*** sur les points internes du maillage
!*** Ex aux points (i+1/2,j)
!*** Ey aux points (i,j+1/2)

do j=iy+1,jy
do i=ix,jx
   dbz_dy = (tm%bz(i,j)-tm%bz(i,j-1)) / dy
   tm%ex(i,j) = tm%ex(i,j) + dt*csq*dbz_dy 
end do
end do

do j=iy,jy
do i=ix+1,jx
   dbz_dx = (tm%bz(i,j)-tm%bz(i-1,j)) / dx
   tm%ey(i,j) = tm%ey(i,j) - dt*csq*dbz_dx 
end do
end do

end subroutine ampere_maxwell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cl_periodiques(tm, ix, jx, iy, jy)

type( mesh_fields ) :: tm
integer, intent(in) :: ix, jx, iy, jy

do i = ix, jx
   dbz_dy = (tm%bz(i,iy)-tm%bz(i,jy)) / dy
   tm%ex(i,iy) = tm%ex(i,iy) + dt*csq*dbz_dy 
   tm%ex(i,jy+1) = tm%ex(i,iy) 
end do

     
do j = iy, jy
   dbz_dx = (tm%bz(ix,j)-tm%bz(jx,j)) / dx
   tm%ey(ix,j) = tm%ey(ix,j) - dt*csq*dbz_dx 
   tm%ey(jx+1,j) = tm%ey(ix,j)
end do

end subroutine

end module solveur_yee
