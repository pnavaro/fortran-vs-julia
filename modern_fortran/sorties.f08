module sorties
   use commun, only: dp, ilp, mesh_fields, dx, dy
   implicit none(external)
   private

   character(len=10) :: outdir = '../data/'

   interface
      module subroutine plot_fields(rang, nproc, f, ix, jx, iy, jy, &
                                    xp, yp, iplot, time)
         implicit none(external)
         integer(ilp), intent(in) :: rang, ix, jx, iy, jy, nproc, iplot
         type(mesh_fields), intent(in) :: f
         real(dp), intent(in) :: xp, yp, time
      end subroutine plot_fields
   end interface
   public :: plot_fields

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   module procedure plot_fields
   integer :: i, j, k
   integer :: kk0, kk1, kk2, kk3, kk4
   character(len=4) :: fin
   character(len=1) :: aa, bb, cc, dd
   integer :: data_unit, gnu_unit

   if (iplot == 1) then
      k = len_trim(outdir)
      outdir(k:) = "/"//char(rang + 48)//"/"
      call system("mkdir -p "//outdir)
   end if

   kk0 = iplot
   kk1 = kk0/1000
   aa = char(kk1 + 48)
   kk2 = (kk0 - kk1*1000)/100
   bb = char(kk2 + 48)
   kk3 = (kk0 - (kk1*1000) - (kk2*100))/10
   cc = char(kk3 + 48)
   kk4 = (kk0 - (kk1*1000) - (kk2*100) - (kk3*10))/1
   dd = char(kk4 + 48)
   fin = aa//bb//cc//dd

   open (newunit=data_unit, file=trim(outdir)//fin//".dat")
   do i = ix, jx
      do j = iy, jy
         write (data_unit, "(4e10.2)") xp + (i - 0.5)*dx, yp + (j - .5)*dy, f%bz(i, j)
      end do
      write (data_unit, *)
   end do
   close (data_unit)

   if (rang == 0) then
      do k = 1, 3
         open (newunit=gnu_unit, file='bz.gnu', position="append")
         if (iplot == 1) then
            rewind (gnu_unit)
            write (gnu_unit, *) "set xr[-0.1:1.1]"
            write (gnu_unit, *) "set yr[-0.1:1.1]"
            write (gnu_unit, *) "set zr[-1.1:1.1]"
            !write(gnu_unit,*)"set cbrange[-1:1]"
            !write(gnu_unit,*)"set pm3d"
            write (gnu_unit, *) "set surf"
            !write(gnu_unit,*)"set term x11"
         end if
         write (gnu_unit, *) "set title 'Time = ", time, "'"
         write (gnu_unit, "(a)", advance='no') "splot '"         &
         & //trim(outdir)//fin//".dat' u 1:2:3 w lines"

         do j = 1, nproc - 1
            write (gnu_unit, "(a)", advance='no') "&
            & ,'"//outdir(1:5)//char(j + 48)//"/"//fin// &
            & ".dat' u 1:2:3 w lines"
         end do
         write (gnu_unit, *)
         close (gnu_unit)
      end do
   end if

   end procedure plot_fields

end module sorties
