! Global definitions

module global

use omp_lib
use par_zig_mod

implicit none

!INTEGER,  PARAMETER  ::  DP=SELECTED_REAL_KIND( 12, 60 )
INTEGER,  PARAMETER  ::  SP = kind(1.0), DP = kind(1.0d0)
integer, parameter :: REAL_KIND = SP

integer, parameter :: nfin=10, nfout=11, nflog=12, nfres=13, nfrun=14
integer, parameter :: neumann(3,6) = reshape((/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /), (/3,6/))

integer, parameter :: OUTSIDE_TAG = -99999

type occupancy_type
	integer :: indx(2)
	real :: C
	type(boundary_type), pointer :: bdry
end type

type cell_type
	integer :: ID
	integer :: site(3)
	integer :: state
	logical :: exists
end type

type boundary_type
    integer :: site(3)
    type (boundary_type), pointer :: next
end type


type(occupancy_type), allocatable :: occupancy(:,:,:)
type(cell_type), allocatable :: cell_list(:)
!type(boundary_type), pointer :: bdrylist

integer :: NX, NY, NZ
integer, allocatable :: zdomain(:),zoffset(:)
integer :: blobrange(3,2)
real :: DELTA_X, DELTA_T, PI
real :: x0,y0,z0   ! centre in global coordinates (units = grids)
real :: blob_radius, Radius, Centre(3)
integer :: jumpvec(3,27)

integer :: max_nlist, nlist, Ncells, Ncells0, lastNcells, lastID
integer :: max_ngaps, ngaps, nadd_sites, Nsites
integer :: nbdry
integer :: istep, nsteps
integer :: Mnodes

character*(128) :: inputfile
character*(128) :: outputfile
character*(2048) :: logmsg
!TYPE(winsockport) :: awp_0, awp_1
logical :: use_TCP = .true.         ! turned off in para_main()
logical :: use_CPORT1 = .false.
!logical :: stopped, clear_to_send
logical :: simulation_start, par_zig_init, initialized
logical :: dbug = .false.

integer :: seed(2)

!DEC$ ATTRIBUTES DLLEXPORT :: nsteps

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine logger(msg)
character*(*) :: msg
integer :: error
logical :: isopen
character*(1) :: LF = char(94)

error = 0
if (use_TCP) then
!    if (awp_0%is_open) then
!        call winsock_send(awp_0,trim(msg)//LF,len_trim(msg)+1,error)
!    endif
else
	write(*,*) trim(msg)
endif
inquire(unit=nflog,OPENED=isopen)
if (isopen) then
	write(nflog,*) 'msg: ',trim(msg)
	if (error /= 0) then
	    write(nflog,'(a,i4)') 'winsock_send error: ',error
	    close(nflog)
	endif
endif
if (error /= 0) stop
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine make_jumpvec
integer :: k,ix,iy,iz

k = 0
do ix = -1,1
	do iy = -1,1
		do iz = -1,1
			k = k+1
			jumpvec(:,k) = (/ix,iy,iz/)
		enddo
	enddo
enddo
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine SetRadius(N)
integer :: N
Radius = (N*3.0/(4*PI))**(1./3)
end subroutine

!---------------------------------------------------------------------
! Uniform randomly generates an integer I: n1 <= I <= n2
!---------------------------------------------------------------------
integer function random_int(n1,n2,kpar)
integer :: n1,n2,kpar
integer :: k,R

if (n1 == n2) then
    random_int = n1
elseif (n1 > n2) then
    write(logmsg,*) 'ERROR: random_int: n1 > n2: ',n1,n2
    call logger(logmsg)
    stop
endif
R = par_shr3(kpar)
k = abs(R)
random_int = n1 + mod(k,(n2-n1+1))

end function

!--------------------------------------------------------------------------------
! Returns a permutation of the elements of a()
!--------------------------------------------------------------------------------
subroutine permute(a,n,kpar)
integer :: a(*),n,kpar
integer :: i,k,tmp

do i = 1,n
    k = random_int(1,n,kpar)
	tmp = a(i)
	a(i) = a(k)
	a(k) = tmp
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine waste_time(n,dummy)
integer :: k, n
real :: dummy
real(DP) :: rsum,R
integer :: kpar=0

rsum = 0
do k = 1,n
    R = par_uni(kpar)
    rsum = rsum + R
enddo
dummy = rsum
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real function norm(r)
real :: r(3)

norm = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real function norm2(r)
real :: r(3)

norm2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine normalize(r)
real :: r(3)

r = r/norm(r)
end subroutine

!-----------------------------------------------------------------------------------------
! Distance from the blob centre (units = grids)
!-----------------------------------------------------------------------------------------
real function cdistance(site)
integer :: site(3)
real :: r(3)

r = site - Centre
cdistance = norm(r)
end function

!--------------------------------------------------------------------------------
! Crude check for a site inside, just looking at the allowable ranges of x, y and z.
!--------------------------------------------------------------------------------
logical function inside_xyz(site)
integer :: site(3)

if (site(1) < 1 .or. site(1) > NX .or. site(2) < 1 .or. site(2) > NY .or. site(3) < 1 .or. site(3) > NZ) then
    inside_xyz = .false.
else
    inside_xyz = .true.
endif
end function


!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine get_vnorm(v,vnorm)
real :: v(3), vnorm(3)
real :: d

d = dot_product(v,v)
vnorm = v/sqrt(d)
end subroutine

end module
