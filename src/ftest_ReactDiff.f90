! To test cvReactDiff_bnd

module ReactDiff
use iso_c_binding
implicit none

integer, parameter :: dp = kind(1.d0)

logical, parameter :: SPHERE = .true.
real(dp), parameter :: RADIUS = 15
real(dp), parameter :: K01 = 0.001
!#define CADV  0.0			 /* coeff of advection term */

real(dp), parameter :: XMIN = 0           ! grid boundaries in x  
real(dp), parameter :: XMAX = 20.0
real(dp), parameter :: YMIN = 30.0        ! grid boundaries in y  
real(dp), parameter :: YMAX = 50.0
real(dp), parameter :: ZMIN = 0
real(dp), parameter :: ZMAX = 20.0

!static void SetIC(int ndim, int mx, int my, int mz, int NG, int nvars, int *xmap, int *ymap, int *zmap, realtype *v);
!static void	MakeMaps(int ndim, int mx, int my, int mz, int *NG, int **xmap, int **ymap, int **zmap);
!void React(realtype *C, realtype *dCdt);
!void ReactJac(realtype C[], realtype dfdC[], int nvars);

interface
	subroutine Solve(ndim, mx, my, mz, NG, nvars, vbnd, cdiff, &
		xmap, ymap, zmap, v, t0, t1, dtout, nout) bind(C)
	use iso_c_binding
	integer(c_int), VALUE :: ndim, mx, my, mz, NG, nvars, nout
	integer(c_int) :: xmap(*), ymap(*), zmap(*)
	real(c_double) :: vbnd(*), cdiff(*), v(*)
	real(c_double), VALUE :: t0, t1, dtout
	end subroutine Solve
end interface

contains

!------------------------------------------------------------------------------
! This function is to generate test geometries.
! In actual use, NG, xmap, ymap, zmap will be passed from the invoking program.
! ndim, MX, MY, MZ will also be supplied.
! Note that all "inside" gridcells must have 1 <= x <= MX etc.
!------------------------------------------------------------------------------
subroutine MakeMaps(ndim, MX, MY, MZ, NG, xmap, ymap, zmap)
integer :: ndim, MX, MY, MZ, NG
integer, allocatable :: xmap(:), ymap(:), zmap(:)
integer :: x, y, z, kg, nin, n2, n3
real(dp) :: x0, y0, z0, r2
logical, allocatable :: inside(:)
logical :: in

if (ndim == 2) then
	if (SPHERE) then
		x0 = (MX + 1)/2.0
		y0 = (MY + 1)/2.0
		allocate(inside(0:(MX+2)*(MY+2)))
		nin = 0
		do x = 0, MX+2-1
			do y = 0, MY+2-1
				r2 = (x - x0)*(x-x0) + (y-y0)*(y-y0)
				if (r2 < RADIUS*RADIUS) then
					inside(y + x*(MY+2)) = .TRUE.
					nin = nin+1
				else
					inside(y + x*(MY+2)) = .FALSE.
				endif
			enddo
		enddo
	else
		nin = MX*MY
	endif
	allocate(xmap(0:nin-1))
	allocate(ymap(0:nin-1))
	allocate(zmap(1))
	NG = nin
	kg = -1
	do x = 0, MX+2-1
		do y = 0, MY+2-1
			n2 = y + x*(MY+2)
			in = .FALSE.
			if (SPHERE) then
				if (inside(n2)) in = .TRUE.
			else
				if (x >= 1 .and. x <= MX .and. y >= 1 .and. y <= MY) in = .TRUE.
			endif
			if (in) then
				kg = kg + 1
				xmap(kg) = x
				ymap(kg) = y
			endif
		enddo
	enddo
	! kg is the grid index in the list.
	if (SPHERE) then
		deallocate(inside)
	endif
elseif (ndim == 3) then
	if (SPHERE) then
		x0 = (MX + 1)/2.0
		y0 = (MY + 1)/2.0
		z0 = (MZ + 1)/2.0
		allocate(inside(0:(MX+2)*(MY+2)*(MZ+2)))
		nin = 0
		do x = 0, MX+2-1
			do y = 0, MY+2-1
				do z = 0, MZ+2-1
					r2 = (x - x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0)
					if (r2 < RADIUS*RADIUS) then
						inside(z + y*(MZ+2) + x*(MZ+2)*(MY+2)) = .TRUE.
						nin = nin + 1
					else
						inside(z + y*(MZ+2) + x*(MZ+2)*(MY+2)) = .FALSE.
					endif
				enddo
			enddo
		enddo
	else
		nin = MX*MY*MZ
	endif
	allocate(xmap(0:nin-1))
	allocate(ymap(0:nin-1))
	allocate(zmap(0:nin-1))
	NG = nin
	kg = -1
	do x = 0, MX+2-1
		do y = 0, MY+2-1
			do z = 0, MZ+2-1
				n3 = z + y*(MZ+2) + x*(MZ+2)*(MY+2)
				in = .FALSE.
				if (SPHERE) then
					if (inside(n3)) in = .TRUE.
				else
					if (x >= 1 .and. x <= MX &
						.and. y >= 1 .and. y <= MY &
						.and. z >= 1 .and. z <= MZ) in = .TRUE.
				endif
				if (in) then
					kg = kg+1
					xmap(kg) = x
					ymap(kg) = y
					zmap(kg) = z
				endif
			enddo
		enddo
	enddo
	if (SPHERE) then
		deallocate(inside)
	endif
endif
end subroutine

!----------------------------------------------------------------------------
! Set test-case initial conditions in v
!----------------------------------------------------------------------------
subroutine SetIC(ndim, mx, my, mz, NG, nvars, xmap, ymap, zmap, v)
integer :: ndim, mx, my, mz, NG, nvars
integer :: xmap(:), ymap(:), zmap(:)
real(dp) :: v(:)
integer :: kg, ke
real(dp) :: x, y, z, dx, dy, dz
real(dp) :: xmid, ymid, zmid, ax, ay, az

dx = (XMAX-XMIN)/mx  
dy = (YMAX-YMIN)/my
xmid = 0.5*(XMIN + XMAX)
ymid = 0.5*(YMIN + YMAX)
if (ndim == 2) then
	do kg = 0, NG-1
		x = XMIN + dx*(xmap(kg) - 0.5)
		y = YMIN + dy*(ymap(kg) - 0.5)
		ke = kg*nvars
		if (x < xmid) then
			ax = (x-XMIN)/(xmid-XMIN)
		else
			ax = (XMAX-x)/(XMAX-xmid)
		endif
		if (y < ymid) then
			ay = (y-YMIN)/(ymid-YMIN)
		else
			ay = (YMAX-y)/(YMAX-ymid)
		endif
		v(ke)   = 	ax*ay
		v(ke+1) = 10*ax*ay
	enddo
else if (ndim == 3) then
	dz = (ZMAX-ZMIN)/mz
	xmid = 0.5*(XMIN + XMAX)
	ymid = 0.5*(YMIN + YMAX)
	zmid = 0.5*(ZMIN + ZMAX)
	do kg = 0, NG-1
		x = XMIN + dx*(xmap(kg) - 0.5)
		y = YMIN + dy*(ymap(kg) - 0.5)
		z = ZMIN + dz*(zmap(kg) - 0.5)
		ke = kg*nvars
		if (x < xmid) then
			ax = (x-XMIN)/(xmid-XMIN)
		else
			ax = (XMAX-x)/(XMAX-xmid)
		endif
		if (y < ymid) then
			ay = (y-YMIN)/(ymid-YMIN)
		else
			ay = (YMAX-y)/(YMAX-ymid)
		endif
		if (z < zmid) then
			az = (z-ZMIN)/(zmid-ZMIN)
		else
			az = (ZMAX-z)/(ZMAX-zmid)
		endif
		v(ke)   = 20*ax*ay*az
		v(ke+1) = 10*ax*ay*az
	enddo
endif
end subroutine

!------------------------------------------------------------------------------
! Rates of change of consituents are determined from constituent concentrations
! (later we may wish to add sources and sinks)
!------------------------------------------------------------------------------
subroutine React(C, dCdt) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: React
use, intrinsic :: iso_c_binding
real(c_double) :: C(*), dCdt(*)

dCdt(0) = K01*C(0)*C(1)
dCdt(1) = -K01*C(0)*C(1)
end subroutine

!------------------------------------------------------------------------------
! dfdC[k'][k] = df[k']/dC[k]
! df(jc)/dC(ic) = dfdC[k], index k = jc + ic*nvars
!------------------------------------------------------------------------------
subroutine ReactJac(C, dfdC, nvars)
!DEC$ ATTRIBUTES DLLEXPORT :: ReactJac
use, intrinsic :: iso_c_binding
real(c_double) :: C(*), dfdC(*)
integer(c_int) :: nvars
integer :: ic, jc, k

!dfdC[0][0] = K01*C[1];
!dfdC[0][1] = K01*C[0];
!dfdC[1][0] = -K01*C[1];
!dfdC[1][1] = -K01*C[0];
jc = 0
ic = 0
k = jc + ic*nvars
dfdC(k) = K01*C(1)
jc = 0
ic = 1
k = jc + ic*nvars
dfdC(k) = K01*C(0)
jc = 1
ic = 0
k = jc + ic*nvars
dfdC(k) = -K01*C(1)
jc = 1
ic = 1
k = jc + ic*nvars
dfdC(k) = -K01*C(0)
end subroutine

end module

program main
use ReactDiff

!clock_t begin, end;
!double time_spent;
integer :: ndim, mx, my, mz, NG, nvars, nout
integer, allocatable :: xmap(:), ymap(:), zmap(:)
real(dp), allocatable :: v(:)
real(dp) :: dx, cdiff(3), vbnd(2), t0, t1, dtout
real(dp) :: DIFFCOEF = 1.0e-1		 ! coeff of diffusion

!begin = clock();

ndim = 3
mx = 30
my = 30
mz = 30
nvars = 2
vbnd(0) = 0
vbnd(1) = 0
dx = (XMAX-XMIN)/mx  
cdiff(0) = DIFFCOEF/(dx*dx)
cdiff(1) = DIFFCOEF/(dx*dx)
cdiff(2) = DIFFCOEF/(dx*dx)

t0 = 0
t1 = 10
dtout = t1
nout = 20
call MakeMaps(ndim,mx,my,mz,NG,xmap,ymap,zmap)

allocate(v(NG*nvars))
call SetIC(ndim,mx,my,mz,NG,nvars,xmap,ymap,zmap,v)
call Solve(ndim,mx,my,mz,NG,nvars,vbnd,cdiff,xmap,ymap,zmap,v,t0,t1,dtout,nout)

!end = clock();
!time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
!printf("Solve time: %6.1f\n",time_spent);
end program
