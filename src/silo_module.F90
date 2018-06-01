module silo_module

use commun, only:tm_mesh_fields
implicit none
#ifdef HAVE_SILO
include "silo.inc"
integer, private :: i, j
integer :: kk0, kk1, kk2, kk3, kk4
character(len=4) :: fin
character(len=1) :: aa,bb,cc,dd


contains


subroutine write_domains(imesh, xc, yc, dx, dy, mx, my, phi, iplot)
integer ::  err, ierr, i, dom, ndims, imesh, iplot
integer :: dims(2), mx, my, dbfile
type(tm_mesh_fields), intent(in) :: phi
real(8), dimension(mx) :: x
real(8), dimension(my) :: y
real(8) ::  xc, yc, dx, dy
character(len=20) :: filename 
integer :: lfilename

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

ndims = 2
dims(1) = mx
dims(2) = my
! Poke a number into the filename.
filename = "data/"//char(48+imesh)//"/"//fin
lfilename = len_trim(filename)
! Create a new silo file.
ierr = dbcreate(filename, lfilename, DB_CLOBBER, DB_LOCAL,	&
                "multimesh data", 14, DB_HDF5, dbfile)
if (dbfile.eq.-1) then
   write (6,*) 'Could not create Silo file!\n'
   return
end if
! Displace the coordinates
do i=1,mx
   x(i) = xc + (i-1)*dx
end do
do j=1,my
   y(j) = yc + (j-1)*dy
end do

! Write the multimesh
err = dbputqm (dbfile, "mesh", 4, "x", 1,		&	
     "y", 1, "z", 1, x, y, DB_F77NULL, dims, ndims,	&	
     DB_DOUBLE, DB_COLLINEAR, DB_F77NULL, ierr)

err = dbputqv1(dbfile, "Ex", 2, "mesh", 4, phi%ex(1:mx,1:my), dims, &
               ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
err = dbputqv1(dbfile, "Ey", 2, "mesh", 4, phi%ey(1:mx,1:my), dims, &
               ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
err = dbputqv1(dbfile, "Bz", 2, "mesh", 4, phi%bz(1:mx,1:my), dims, &
               ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
! Close the Silo file
ierr = dbclose(dbfile)

end subroutine write_domains

subroutine write_master(nmesh, iplot)

integer :: err, ierr, nmesh, oldlen, nvar, dbfile, iplot
character(len=20), dimension(nmesh) :: meshnames
character(len=20), dimension(nmesh) :: varnames1, varnames2, varnames3
integer, dimension(nmesh) :: meshtypes , lmeshnames
integer, dimension(nmesh) :: lvarnames, vartypes

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

! Create a new silo file
ierr = dbcreate("data/"//fin//".root",14, DB_CLOBBER, DB_LOCAL, &
                "multimesh root", 14, DB_HDF5, dbfile)

if(dbfile.eq.-1) then
   write (6,*) "Could not create Silo file!\n"
   return
endif

! Set the maximum string length to 20 since that's how long our
! strings are
oldlen = dbget2dstrlen()
err = dbset2dstrlen(20)

do i = 1, nmesh
   meshnames(i)  = char(i+47)//"/"//fin//":mesh"
   meshtypes(i)  = DB_QUAD_RECT
   lmeshnames(i) = len_trim(meshnames(i))
   varnames1(i)  = char(i+47)//"/"//fin//":Ex"
   varnames2(i)  = char(i+47)//"/"//fin//":Ey"
   varnames3(i)  = char(i+47)//"/"//fin//":Bz"
   lvarnames(i)  = 9
   vartypes(i)   = DB_QUADVAR
end do
! Write the multimesh object.
err = dbputmmesh(dbfile, "mesh", 4, nmesh, meshnames,	&
                 lmeshnames, meshtypes, DB_F77NULL, ierr)
nvar = nmesh
err = dbputmvar(dbfile, "Ex", 2, nvar, varnames1, lvarnames,	&
                vartypes, DB_F77NULL, ierr)
err = dbputmvar(dbfile, "Ey", 2, nvar, varnames2, lvarnames,	&
                vartypes, DB_F77NULL, ierr)
err = dbputmvar(dbfile, "Bz", 2, nvar, varnames3, lvarnames,	&
                vartypes, DB_F77NULL, ierr)

! Restore the previous value for maximum string length
err = dbset2dstrlen(oldlen)
! Close the Silo file
ierr = dbclose(dbfile)

end subroutine write_master

#endif
end module silo_module
