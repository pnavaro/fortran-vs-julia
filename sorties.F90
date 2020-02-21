module sorties

use commun, only: mesh_fields, dx, dy

implicit none

integer, private :: i, j, k
character(len=10) :: outdir = 'data/'

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plot_fields(rang, nproc, f, ix, jx, iy, jy, &
                       xp, yp, iplot, time )

integer, intent(in) :: rang, ix, jx, iy, jy, nproc
type (mesh_fields), intent(in) :: f
integer :: iplot, i, j
real(8) :: time, xp, yp
integer :: kk0, kk1, kk2, kk3, kk4
character(len=4) :: fin
character(len=1) :: aa,bb,cc,dd

if (iplot == 1) then
   k = len_trim(outdir)
   outdir(k:) = "/"//char(rang+48)//"/"
   call system("mkdir -p "//outdir)
end if

kk0 = iplot
kk1 = kk0/1000
aa  = char(kk1 + 48)
kk2 = (kk0 - kk1*1000)/100
bb  = char(kk2 + 48)
kk3 = (kk0 - (kk1*1000) - (kk2*100))/10
cc  = char(kk3 + 48)
kk4 = (kk0 - (kk1*1000) - (kk2*100) - (kk3*10))/1
dd  = char(kk4 + 48)
fin = aa//bb//cc//dd

open( 80, file = trim(outdir)//fin//".dat" )
   do i=ix,jx
      do j=iy,jy
         write(80,"(4e10.2)") xp+(i-0.5)*dx, yp+(j-.5)*dy, f%bz(i,j)
      end do
      write(80,*) 
   end do
close(80)
   
if (rang == 0) then
   do k = 1, 3
      open( 90, file = 'bz.gnu', position="append" )
      if ( iplot == 1 ) then
         rewind(90)
         write(90,*)"set xr[-0.1:1.1]"
         write(90,*)"set yr[-0.1:1.1]"
         write(90,*)"set zr[-1.1:1.1]"
         !write(90,*)"set cbrange[-1:1]"
         !write(90,*)"set pm3d"
         write(90,*)"set surf"
         !write(90,*)"set term x11"
      end if
      write(90,*)"set title 'Time = ",time,"'"
      write(90,"(a)",advance='no')"splot '" 	&
      & //trim(outdir)//fin//".dat' u 1:2:3 w lines"
   
      do j = 1, nproc - 1
         write(90,"(a)",advance='no')"&
         & ,'"//outdir(1:5)//char(j+48)//"/"//fin// &
         & ".dat' u 1:2:3 w lines" 
      end do
      write(90,*)
      close(90)
   end do
end if

end subroutine plot_fields

end module sorties
