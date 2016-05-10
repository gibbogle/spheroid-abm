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

z = site(3) + cdrop*Radius - zmin
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
integer function CountOutside() result(nout)
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
integer :: x, y, z, z1, z2, dz, kcell, zbmax, zmax, xv, yv, nv, npath, nvtot, nstot, newtot, nbtot
integer :: zlow, zb(NX,NY), noutside, incontact
real(REAL_KIND) :: z0drop	! drop centre is at x0,y0,z0drop = zmin + (bdrop-cdrop)*R
real(REAL_KIND) :: Cext(MAX_CHEMO)
logical :: ok
type(path_type) :: path(NX)

allocate(bdist(NX,NY,NZ))
allocate(cdist(NX,NY))
allocate(rz(NZ))
allocate(nstack(NX,NY))
allocate(usable(NX,NY))

!call SetRadius(Nsites)
zmin = GetZmin()
z0drop = zmin + (bdrop-cdrop)*Radius
z0 = z0drop
Centre(3) = z0
sintheta0 = sqrt(1 - (1-cdrop/bdrop)**2)
Rcontact = adrop*Radius*sintheta0
Rc2 = Rcontact*Rcontact
Ra2 = (adrop*Radius)**2

! Stage 1
!--------
! drop cells in contact disc
zmax = 0
do x = 1,NX
	do y = 1,NY
		r2 = (x-x0)*(x-x0) + (y-y0)*(y-y0)
		if (r2 > Ra2) cycle
		r = sqrt(r2)
		if (r2 > Rc2) then
			incontact = 0
			! we need to find the desired lower boundary z=zlow at (x,y)
			! r = a.R.sin(theta) ==> theta = asin(r/aR)
			theta = asin(r/(adrop*Radius))
			zlow = bdrop*Radius*(1-cos(theta)) + zmin - cdrop*Radius + 1
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
		! cells span z1 <= z <= z2
		dz = z1 - zlow
		do z = z1,z2
			! move the site+cell at (x,y,z) to (x,y,z-dz)
			kcell = occupancy(x,y,z)%indx(1)
			Cext = occupancy(x,y,z)%C
			cell_list(kcell)%site = (/x,y,z-dz/)
			occupancy(x,y,z-dz)%indx(1) = kcell
			occupancy(x,y,z-dz)%C = Cext
			occupancy(x,y,z)%indx(1) = OUTSIDE_TAG
		enddo
	enddo
enddo
	
! Stage 2
!--------
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
!write(nflog,*) 'Approximate number of sites in the squashed spheroid: ',newtot

bdist = -1
nbtot = 0
do x = 1,NX
	do y = 1,NY
		nstack(x,y) = 0
		if (zb(x,y) <= 0) cycle
		do z = 1,NZ
			if (occupancy(x,y,z)%indx(1) <= 0) cycle
			if (z > zb(x,y)) nstack(x,y) = nstack(x,y) + 1
		enddo
		nstot = nstot + nstack(x,y)
		r = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0))
		cdist(x,y) = r
		do z = zmin,zmax
			if (z <= zmin + Radius*(2*bdrop - cdrop)) then
				if (rz(z) >= r) then
					bdist(x,y,z) = rz(z) - r
					nbtot = nbtot + 1
				endif
			endif
		enddo
	enddo
enddo

!write(nflog,*) 'Actual number of sites in the squashed spheroid: ',nbtot
nvtot = 0
do z = zmin+1,zmax
	usable = .true.		! initially set all sites in this layer as usable
	call GetNearestVacantSite(z,xv,yv,nv)
	nvtot = nvtot + nv
enddo

noutside = CountOutside()

if (.true.) then
! move cells to the expanded radius
do z = zmin+1,zmax
	usable = .true.		! initially set all sites in this layer as usable
	do 
		call GetNearestVacantSite(z,xv,yv,nv)
		if (nv == 0) exit
		call GetBestPath(xv,yv,z,path,npath,ok)
		if (ok) then
			call UsePath(xv,yv,z,path,npath)
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

is_dropped = .true.
occupancy(:,:,1:zmin-1)%indx(1) = UNREACHABLE_TAG
! Presumably we need to first delete the bdrylist
call DestroyBdryList
call CreateBdryList

call CheckBdryList('after dropper')

end subroutine


!--------------------------------------------------------------------------------
! The best path is the shortest path from (xv,yv) to the best site.  
! The best site is the site (x,y) within a circle that has the highest nstack(x,y)/dist,
! where dist is the distance from (xv,yv)  
! The circle is defined such that (xv,yv) and (x0,y0) are at ends of a diameter,
! and both lie in the circle.
!--------------------------------------------------------------------------------
subroutine GetBestPath(xv,yv,z,path,npath,ok)
integer :: xv, yv, z, npath
type(path_type) :: path(NX)
logical :: ok
integer :: ix0, iy0, x, y, xbest, ybest
real(REAL_KIND) :: cx0, cy0, cr, cr2, r2, d, val, vmax

ix0 = x0
iy0 = y0

! Define the circle: radius and centre
cr = sqrt((xv-ix0)**2 + (yv-iy0)**2 + 2.0)/2
cr2 = cr*cr
cx0 = (xv + ix0)/2.
cy0 = (yv + iy0)/2.
! Find best site inside the circle
vmax = 0
do x = 1,NX
	do y = 1,NY
		if (x == xv .and. y == yv) cycle
		if (bdist(x,y,z) > 0) then	! inside the squashed sphere
			r2 = (x - cx0)**2 + (y - cy0)**2
			if (r2 <= cr2 .and. nstack(x,y) > 0) then	! inside the circle
				d = sqrt((x-xv)**2. + (y-yv)**2.)
				val = nstack(x,y)/d
				if (val > vmax) then
					vmax = val
					xbest = x
					ybest = y
				endif
			endif
		endif
	enddo
enddo
if (vmax == 0) then
	ok = .false.
	return
endif
ok = .true.
! Now we need the path from (xv,yv) to (xbest,ybest)
call GetPath(xv,yv,z,xbest,ybest,path,npath)
end subroutine

!--------------------------------------------------------------------------------
! Get shortest path from (xv,yv) to (xe,ye)
!--------------------------------------------------------------------------------
subroutine GetPath(xv,yv,z,xe,ye,path,npath)
integer :: xv, yv, z, xe, ye, npath
type(path_type) :: path(NX)
integer :: x, y, dx, dy, xmin, ymin, xlast, ylast, nstot
real(REAL_KIND) :: r2, r2min, r2prev

nstot = 0
xlast = xv
ylast = yv
r2prev = 1.0e10
npath = 0
do
	r2min = 1.0e10
	do dx = -1,1	!,2
		do dy = -1,1	!,2
			if (dx == 0 .and. dy == 0) cycle
			x = xlast + dx
			y = ylast + dy
			r2 = (x-xe)*(x-xe) + (y-ye)*(y-ye)
			if (r2 < r2min) then
				r2min = r2
				xmin = x
				ymin = y
			endif
		enddo
	enddo
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
		return
	endif
enddo		
end subroutine

!--------------------------------------------------------------------------------
! Find the shortest path from (xv,yv) to the centre of the squashed spheroid,
! in the z plane, i.e. to (x0,y0)
! If all nstack are 0, no use can be made of the path, ok = .false.
! NOT USED
!--------------------------------------------------------------------------------
subroutine GetCentrePath(xv,yv,z,path,npath,ok)
integer :: xv, yv, z,npath
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
	do dx = -1,1	!,2
		do dy = -1,1	!,2
			if (dx == 0 .and. dy == 0) cycle
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
subroutine UsePath(xv,yv,zslice,path,npath)
integer :: xv, yv, zslice, npath
type(path_type) :: path(NX)
integer :: kpath, nsmax, np, kcell, site(3), z, z1, z2
real(REAL_KIND) :: Cext(MAX_CHEMO)

nsmax = 0
do kpath = 1,npath
	if (path(kpath)%nstack > nsmax) then
		nsmax = path(kpath)%nstack
		np = kpath
	endif
enddo

! Move site+cells in the z plane
do kpath = 1,np
	kcell = occupancy(path(kpath)%x,path(kpath)%y,zslice)%indx(1)
	if (kcell <= 0) then
		write(logmsg,*) 'Error: UsePath: kcell: ',kcell,path(kpath)%x,path(kpath)%y,zslice
		call logger(logmsg)
		stop
	endif
	Cext = occupancy(path(kpath)%x,path(kpath)%y,zslice)%C
	cell_list(kcell)%site = (/xv,yv,zslice/)
	occupancy(xv,yv,zslice)%indx(1) = kcell
	occupancy(xv,yv,zslice)%C = Cext
	xv = path(kpath)%x
	yv = path(kpath)%y
	occupancy(xv,yv,zslice)%indx(1) = 0
enddo

! Drop site+cells to fill vacancy at (xv,yv)
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
	Cext = occupancy(xv,yv,z)%C
	cell_list(kcell)%site = (/xv,yv,z-1/)
	occupancy(xv,yv,z-1)%indx(1) = kcell
	occupancy(xv,yv,z-1)%C = Cext
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
cosa = -sqrt(1 - sina*sina)
GetZb = zmin + Radius*(bdrop*(1 - cosa) - cdrop)
end function

end module

