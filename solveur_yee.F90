module solveur_yee

use commun, only: mesh_fields, dx, dy, dt, csq , c

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cl_condparfait(tm, ix, jx, iy, jy, cchar)

character, intent(in) :: cchar
type( mesh_fields ) :: tm
integer, intent(in) :: ix, jx, iy, jy

select case(cchar)
case('S')
   do i = ix, jx
      tm%ex(i,iy) = 0.d0
   end do
case('N')
   do i = ix, jx
      tm%ex(i,jy) = 0.d0
      tm%bz(i,jy) = tm%bz(i,jy-1)	
   end do
case('W')
   do j = iy, jy
      tm%ey(ix,j) = 0.d0
   end do
case('E')
   do j = iy, jy
      tm%ey(jx,j) = 0.d0
      tm%bz(jx,j) = tm%bz(jx-1,j)
   end do
end select

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine silver_muller( tm, ix, jx, iy, jy, ccall )

character, intent(in) :: ccall
integer, intent(in) :: ix, jx, iy, jy
type( mesh_fields ) :: tm
real(8) :: a11,a12,a21,a22,b1,b2,dis


!Conditions de Silver-Muller
!------------------------------------
!Ey = -c Bz sur la frontiere ouest
!Ey =  c Bz sur la frontiere est
!Ex = -c Bz sur la frontiere nord
!Ex =  c Bz sur la frontiere sud
   
!On effectue le calcul de B sur les points fictifs du maillage
!simultanement avec la prise en compte des conditions limites sur
!E. Pour cela, on resout sur les points frontieres, l'equation de la
!condition limite en moyennant en temps pour E et en espace pour B puis
!l'equation d'Ampere

select case (ccall)

case ('N')
   !Frontiere Nord : Ex = -c Bz 
   do i = ix, jx
         
      a11 = 1.;	a12 = + c
      a21 = 1./dt;	a22 = - csq / dy
      b1  = - tm%ex(i,jy) - c * tm%bz(i,jy-1)
      b2  =   tm%ex(i,jy)/dt - csq/dy*tm%bz(i,jy-1)
         
      dis = a11*a22-a21*a12 
         
      !tm%ex(i,jy) = (b1*a22-b2*a12)/dis
      tm%bz(i,jy) = (a11*b2-a21*b1)/dis
         
   end do
      
case ('S')

   !Frontiere Sud : Ex =  c Bz
   do i = ix, jx
         
      a11 = 1.;	a12 = - c
      a21 = 1./dt;	a22 = csq / dy
      b1  = - tm%ex(i,iy) + c * tm%bz(i,iy+1)
      b2  = tm%ex(i,iy)/dt + csq / dy * tm%bz(i,iy+1) 
         
      dis = a11*a22-a21*a12 
         
      tm%ex(i,iy) = (b1*a22-b2*a12)/dis
      !tm%bz(i,iy) = (a11*b2-a21*b1)/dis
         
   end do
      
case ('E')

   !Frontiere Est : Ey =  c Bz
   do j = iy, jy
         
      a11 = 1.;	a12 = - c
      a21 = 1./dt; a22 = + csq / dx
      b1  = - tm%ey(jx,j) + c * tm%bz(jx-1,j)
      b2  = tm%ey(jx,j)/dt + csq/dx*tm%bz(jx-1,j) 
         
      dis = a11*a22-a21*a12 
         
      !tm%ey(jx,j) = (b1*a22-b2*a12)/dis
      tm%bz(jx,j) = (a11*b2-a21*b1)/dis
      
   end do
      
case ('W')

   !Frontiere Ouest : Ey = -c Bz
   do j = iy, jy
      
      a11 = 1.;	a12 = + c
      a21 = 1./dt;	a22 = - csq / dx
      b1  = - tm%ey(ix,j) - c * tm%bz(ix+1,j)
      b2  =   tm%ey(ix,j)/dt - csq/dx*tm%bz(ix+1,j) 
      
      dis = a11*a22-a21*a12 
   
      tm%ey(ix,j) = (b1*a22-b2*a12)/dis
      !tm%bz(ix,j) = (a11*b2-a21*b1)/dis
      
   end do

end select

end subroutine silver_muller

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module solveur_yee
