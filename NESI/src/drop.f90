!--------------------------------------------------------------------------------
! To handle the dropping of the spheroid to the well bottom.
! The blob takes the shape of a squashed sphere, with a piece cut off the bottom.
! The parameters of the shape are adrop, bdrop, cdrop.  The centre of the whole
! deformed sphere is at (x0,y0,z0drop), the minimum (i.e. the cut plane corresponding
! to the well bottom) is at zmin = the minimum cell position at the time of dropping.
! Subsequently R is computed for the equivalent sphere, i.e. R = (3V/4)^1/3,
! and z0drop = zmin + (bdrop-cdrop)*R
!--------------------------------------------------------------------------------

module drop

use global
use boundary

implicit none

type path_type
	integer :: x, y
	integer :: nstack
end type
real(REAL_KIND), allocatable :: bdist(:,:,:), cdist(:,:), rz(:)
integer, allocatable :: nstack(:,:)
logical, allocatable :: usable(:,:)
contains

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
logical function OutsideSquashedSphere(site) result(out)
integer :: site(3)
real(REAL_KIND) :: r, z, cos2, sin2

z = site(3) + cdrop - zmin
if (z < 0 .or. z > 2*bdrop*Radius) then
	out = .true.
	return
endif
r = sqrt((site(1) - Centre(1))**2 + (site(2) - Centre(2))**2)
cos2 = (1 - z/(bdrop*Radius))**2
sin2 = 1 - cos2
out = (r > adrop*Radius*sqrt(sin2))
end function

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
integer function CountOutside result(nout)
integer :: kcell

nout = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	if (OutsideSquashedSphere(cell_list(kcell)%site)) nout = nout+1
enddo
end function

!--------------------------------------------------------------------------------
! (1) At (x,y) with r <= a.R.sin(theta0), the column of cells is dropped until
! the lowest cell is at z = zmin.
! (2) At (x,y) with a.R.sin(theta0) < r < aR, drop the column to zbnd(x,y)
! where zbnd = the lower boundary z for (x,y).  
!--------------------------------------------------------------------------------
subroutine dropper
real(REAL_KIND) :: sintheta0, Rcontact, Rc2, Ra2, r, r2, rb, cosa, sina, theta
integer :: x, y, z, z1, z2, dz, kcell, zbmax, zmax, xv, yv, nv, npath, nvtot, nstot, newtot
integer :: zlow, zb(NX,NY), noutside, incontact
real(REAL_KIND) :: z0drop	! drop centre is at x0,y0,z0drop = zmin + (bdrop-cdrop)*R
logical :: ok
type(path_type) :: path(NX)

allocate(bdist(NX,NY,NZ))
allocate(cdist(NX,NY))
allocate(rz(NZ))
allocate(nstack(NX,NY))
allocate(usable(NX,NY))

!adrop = 1.105
!bdrop = 0.829
!cdrop = 0.111

call SetRadius(Nsites)
zmin = GetZmin()
write(*,*) 'Radius, zmin: ',Radius,zmin
z0drop = zmin + (bdrop-cdrop)*Radius
z0 = z0drop
Centre(3) = z0
sintheta0 = sqrt(1 - (1-cdrop/bdrop)**2)
Rcontact = adrop*Radius*sintheta0
Rc2 = Rcontact*Rcontact
Ra2 = (adrop*Radius)**2
write(*,*) 'z0drop,sintheta0,Rcontact: ',z0drop,sintheta0,Rcontact
! drop cells in contact disc
zmax = 0
do x = 1,NX
	do y = 1,NY
		r2 = (x-x0)*(x-x0) + (y-y0)*(y-y0)
!		if (r2 > Rc2) cycle							! outside the contact disc
		if (r2 > Ra2) cycle
		r = sqrt(r2)
		if (r2 > Rc2) then
			incontact = 0
			! we need to find the desired lower boundary z=zlow at (x,y)
			! r = a.R.sin(theta) ==> theta = asin(r/aR)
			theta = asin(r/(adrop*Radius))
			zlow = bdrop*Radius*(1-cos(theta)) + zmin - cdrop + 1
		else
			incontact = 1
			zlow = zmin
		endif
		do z = zlow+1,NZ
			if (occupancy(x,y,z)%indx(1) > 0) then
				zmax = max(z,zmax)
			endif
		enddo
		if (occupancy(x,y,zlow)%indx(1) > 0) cycle	! already in contact
!		write(*,*) x,y,r2,occupancy(x,y,zmin)%indx(1)
		z1 = 0
		z2 = 0
		do z = zlow+1,NZ
			if (occupancy(x,y,z)%indx(1) > 0) then
				if (z1 == 0) then
					z1 = z
				endif
				z2 = z
			elseif (z1 > 0) then
				exit
			endif
		enddo
		if (z1 <= 0) cycle
!		if (incontact == 0) then
!			write(*,'(6i4,f6.1)') x,y,incontact,zlow,z1,z2,r
!		endif
!		write(*,*) 'x,y,z1,z2: ',x,y,z1,z2
		! cells span z1 <= z <= z2
		dz = z1 - zlow
		do z = z1,z2
			! move the cell at (x,y,z) to (x,y,z-dz)
			kcell = occupancy(x,y,z)%indx(1)
			cell_list(kcell)%site = (/x,y,z-dz/)
			occupancy(x,y,z-dz)%indx(1) = kcell
			occupancy(x,y,z)%indx(1) = OUTSIDE_TAG
		enddo
	enddo
enddo
	
! Compute initial nstack, bdist, cdist
bdist = -1
zbmax = 0
nstot = 0
do x = 1,NX
	do y = 1,NY
		! calculate maximum z: upper surface of squashed spheroid
		zb(x,y) = GetZb(x,y,zmin)
		zbmax = max(zbmax,zb(x,y))
	enddo
enddo

! Compute squashed spheroid radius at each z
rz = 0
newtot = 0
do z = zmin,zbmax
	cosa = 1 - (z - zmin + cdrop*Radius)/(bdrop*Radius)	
	sina = sqrt(1 - cosa*cosa)
	rz(z) = adrop*Radius*sina
	newtot = newtot + 3.14159*rz(z)*rz(z)
enddo
write(*,*) 'Approximate number of sites in the squashed spheroid: ',newtot

do x = 1,NX
	do y = 1,NY
		nstack(x,y) = 0
		do z = 1,NZ
			if (occupancy(x,y,z)%indx(1) <= 0) cycle
			if (z > zb(x,y)) nstack(x,y) = nstack(x,y) + 1
		enddo
		nstot = nstot + nstack(x,y)
		r = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0))
		cdist(x,y) = r
!		if (nstack(x,y) > 0) then
!			write(*,*) x,y,nstack(x,y),r
!		endif
		do z = zmin,zmax
			if (z > zmin + Radius*(2*bdrop - cdrop)) then
				bdist(x,y,z) = -1
			else
				if (rz(z) >= r) then
					bdist(x,y,z) = rz(z) - r
!					write(*,*) x,y,z,rb,r
				else
					bdist(x,y,z) = -1
				endif
			endif
		enddo
	enddo
enddo

write(*,*) 'zbmax: ',zbmax, zmin + (2*bdrop-cdrop)*Radius
write(*,*) 'zmax: ',zmax
nvtot = 0
do z = zmin+1,zmax
	call GetNearestVacantSite(z,xv,yv,nv)
!	write(*,'(a,4i4)') 'z: xv,yv,nv: ',z,xv,yv,nv
	nvtot = nvtot + nv
enddo
!write(*,*) 'nstot,nvtot: ',nstot,nvtot

noutside = CountOutside()
write(*,*) 'noutside, Ncells: ',noutside,Ncells

if (.true.) then
! move cells to the expanded radius
do z = zmin+1,zmax
!	write(*,*) 'z: ',z
	usable = .true.		! initially set all sites as usable
	do 
		call GetNearestVacantSite(z,xv,yv,nv)
		if (nv == 0) exit
		call GetCentrePath(xv,yv,path,npath,ok)
		if (ok) then
!			write(*,'(a,2f8.1,3i5)') 'path: ',x0,y0,xv,yv,npath
!			write(*,'(10i6)') path(1:npath)
			call UseCentrePath(xv,yv,z,path,npath)
		else
			usable(xv,yv) = .false.
			cycle
		endif
	enddo
enddo
endif

deallocate(bdist)
deallocate(cdist)
deallocate(rz)
deallocate(nstack)
deallocate(usable)

noutside = CountOutside()
write(*,*) 'noutside, Ncells: ',noutside,Ncells

is_squashed = .true.
occupancy(:,:,1:zmin-1)%indx(1) = UNREACHABLE_TAG
! Presumably we need to first delete the bdrylist
call DestroyBdryList
call CreateBdryList

call CheckBdryList('after dropper')

end subroutine

!--------------------------------------------------------------------------------
! Find the shortest path from (xv,yv) to the centre of the squashed spheroid,
! in the z plane, i.e. to (x0,y0)
! If all nstack are 0, no use can be made of the path, ok = .false.
!--------------------------------------------------------------------------------
subroutine GetCentrePath(xv,yv,path,npath,ok)
integer :: xv, yv, npath
type(path_type) :: path(NX)
logical :: ok
integer :: x, y, dx, dy, xmin, ymin, xlast, ylast, nstot
real(REAL_KIND) :: r2, r2min, r2prev

ok = .true.
nstot = 0
xlast = xv
ylast = yv
r2prev = 1.0e10
npath = 0
do
	r2min = 1.0e10
	do dx = -1,1,2
		do dy = -1,1,2
			x = xlast + dx
			y = ylast + dy
			r2 = (x-x0)*(x-x0) + (y-y0)*(y-y0)
			if (r2 < r2min) then
				r2min = r2
				xmin = x
				ymin = y
			endif
		enddo
	enddo
!	write(*,*) xmin,ymin,r2min
	if (r2min < r2prev) then
		npath = npath + 1
		path(npath)%x = xmin
		path(npath)%y = ymin
		path(npath)%nstack = nstack(xmin,ymin)
		nstot = nstot + nstack(xmin,ymin)
		r2prev = r2min
		xlast = xmin
		ylast = ymin
	else
		ok = (nstot > 0)
		return
	endif
enddo		
end subroutine

!--------------------------------------------------------------------------------
! Select the path point with the max nstack, move cells along the path from this
! site, then drop cells down to fill the created vacancy.
! Decrement nstack().
!--------------------------------------------------------------------------------
subroutine UseCentrePath(xv,yv,zslice,path,npath)
integer :: xv, yv, zslice, npath
type(path_type) :: path(NX)
integer :: kpath, nsmax, np, kcell, site(3), z, z1, z2


nsmax = 0
do kpath = 1,npath
	if (path(kpath)%nstack > nsmax) then
		nsmax = path(kpath)%nstack
		np = kpath
	endif
enddo

! Move cells in the z plane
do kpath = 1,np
	kcell = occupancy(path(kpath)%x,path(kpath)%y,zslice)%indx(1)
!	write(*,*) 'loc, kcell: ',kpath,np,path(kpath)%x,path(kpath)%y,zslice,kcell
	cell_list(kcell)%site = (/xv,yv,zslice/)
	occupancy(xv,yv,zslice)%indx(1) = kcell
	xv = path(kpath)%x
	yv = path(kpath)%y
	occupancy(xv,yv,zslice)%indx(1) = OUTSIDE_TAG
enddo

! Drop cells to fill vacancy at (xv,yv)
z1 = 0
z2 = 0
do z = zslice+1,NZ
	if (occupancy(xv,yv,z)%indx(1) > 0) then
		if (z1 == 0) then
			z1 = z
		endif
		z2 = z
	elseif (z1 > 0) then
		exit
	endif
enddo
do z = z1,z2
	! move the cell at (xv,yv,z) to (xv,yv,z-1)
	kcell = occupancy(xv,yv,z)%indx(1)
	cell_list(kcell)%site = (/xv,yv,z-1/)
	occupancy(xv,yv,z-1)%indx(1) = kcell
	occupancy(xv,yv,z)%indx(1) = OUTSIDE_TAG
enddo
nstack(xv,yv) = nstack(xv,yv) - 1
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine GetNearestVacantSite(z,xv,yv,nv)
integer :: z, xv, yv, nv
integer :: x, y
real(REAL_KIND) :: cdistmin

cdistmin = 999
xv = 0
yv = 0
nv = 0
do x = 1,NX
	do y = 1,NY
		if (bdist(x,y,z) > 0) then										! inside the squashed sphere
			if (usable(x,y) .and. occupancy(x,y,z)%indx(1) <= 0) then	! usable vacant site
				nv = nv + 1
				if (cdist(x,y) < cdistmin) then
					cdistmin = cdist(x,y)
					xv = x	
					yv = y
				endif
			endif
		endif
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
! Get z position of the lowest cell.  This becomes the blob lower limit (i.e. just
! above the surface that the blob is sitting on.)
!--------------------------------------------------------------------------------
integer function GetZmin()
integer :: minval
integer :: kcell

GetZmin = 99999
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	GetZmin = min(GetZmin,cell_list(kcell)%site(3))
enddo
end function

!--------------------------------------------------------------------------------
! Get upper bound of squashed spheroid at (x,y).
! A value of -1 flags a site outside the projection of the squashed spheroid.
!--------------------------------------------------------------------------------
integer function GetZb(x, y, zmin)
integer :: x, y, zmin
real(REAL_KIND) :: r, sina, cosa

r = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0))
if (r > adrop*Radius) then
	GetZb = -1
	return
endif
sina = r/(adrop*Radius)
cosa = sqrt(1 - sina*sina)
GetZb = zmin + Radius*(bdrop*(1 + cosa) - cdrop)
end function

end module

