! Global definitions

module global

use omp_lib
use par_zig_mod
use winsock
use, intrinsic :: ISO_C_BINDING

implicit none

!INTEGER,  PARAMETER  ::  DP=SELECTED_REAL_KIND( 12, 60 )
INTEGER,  PARAMETER  ::  SP = kind(1.0), DP = kind(1.0d0)
integer, parameter :: REAL_KIND = DP
integer, parameter :: TCP_PORT_0 = 5000		! main communication port (logging)
integer, parameter :: TCP_PORT_1 = 5001		! data transfer port (plotting)
integer, parameter :: NORMAL_DIST      = 1
integer, parameter :: LOGNORMAL_DIST   = 2
integer, parameter :: EXPONENTIAL_DIST = 3
integer, parameter :: CONSTANT_DIST    = 4

integer, parameter :: DIVIDING  = 1
integer, parameter :: QUIESCENT = 2
integer, parameter :: DEAD      = 3

integer, parameter :: OUTSIDE_TAG  = -1

integer, parameter :: DIVIDE_ALWAYS_PUSH  = 1
integer, parameter :: DIVIDE_USE_CLEAR_SITE  = 2

integer, parameter :: nfin=10, nfout=11, nflog=12, nfres=13, nfrun=14, nfcell=15, nftreatment=16
integer, parameter :: neumann(3,6) = reshape((/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /), (/3,6/))

!real(REAL_KIND), parameter :: DELTA_X = 0.002	! cm = 20 um
!real(REAL_KIND), parameter :: MM_THRESHOLD = 0.0001
integer, parameter :: MAX_CHEMO = 6
integer, parameter :: OXYGEN = 1
integer, parameter :: GLUCOSE = 2
integer, parameter :: DRUG_A = 3
integer, parameter :: DRUG_A_METAB = 4
integer, parameter :: DRUG_B = 5
integer, parameter :: DRUG_B_METAB = 6
integer, parameter :: EXTRA = 1
integer, parameter :: INTRA = 2
integer, parameter :: SN30000 = DRUG_A
integer, parameter :: SN30000_METAB = DRUG_A_METAB
logical, parameter :: use_ODE_diffusion = .true.
logical, parameter :: compute_concentrations = .true.
logical, parameter :: use_division = .true.
logical, parameter :: use_death = .true.
logical, parameter :: use_react = .true.
logical, parameter :: use_migration = .true.
logical, parameter :: use_medium_flux = .true.	! flux of constituents between spheroid and medium is accounted for.
logical, parameter :: use_metabolites = .true.

!integer, parameter :: MAX_RECEPTOR = 1
real(REAL_KIND), parameter :: PI = 4.0*atan(1.0)

type occupancy_type
	integer :: indx(2)
	real(REAL_KIND) :: C(MAX_CHEMO)
	type(boundary_type), pointer :: bdry
end type

type cell_type
	integer :: ID
	integer :: site(3)
	integer :: state
	real(REAL_KIND) :: conc(MAX_CHEMO)
!	real(REAL_KIND) :: oxygen
!	real(REAL_KIND) :: drug_A
!	real(REAL_KIND) :: drug_B
	real(REAL_KIND) :: dVdt
	real(REAL_KIND) :: volume
	real(REAL_KIND) :: divide_volume
	real(REAL_KIND) :: t_divide_last
	real(REAL_KIND) :: t_divide_next
	real(REAL_KIND) :: t_hypoxic
	real(REAL_KIND) :: t_anoxia_die
	real(REAL_KIND) :: M
	logical :: radiation_tag, drug_tag, anoxia_tag
	logical :: exists
end type

type SN30K_type
	real(REAL_KIND) :: diff_coef
	real(REAL_KIND) :: C1
	real(REAL_KIND) :: C2
	real(REAL_KIND) :: Kmet0
	real(REAL_KIND) :: KO2
	real(REAL_KIND) :: gamma
	real(REAL_KIND) :: Klesion
	real(REAL_KIND) :: halflife
	real(REAL_KIND) :: metabolite_halflife
	real(REAL_KIND) :: kill_O2
	real(REAL_KIND) :: kill_drug
	real(REAL_KIND) :: kill_duration
	real(REAL_KIND) :: kill_fraction
	real(REAL_KIND) :: Kd
end type

type boundary_type
    integer :: site(3)
    type (boundary_type), pointer :: next
end type

type dist_type
	integer :: class
	real(REAL_KIND) :: p1, p2, p3
end type

type, bind(C) :: field_data
	integer(c_int) :: site(3)
	integer(c_int) :: state
	real(c_double) :: dVdt
	real(c_double) :: volume
	real(c_double) :: conc(MAX_CHEMO+1)	! This needs to agree with the definition in field.h
end type

type treatment_type
	integer :: n
	character*(16) :: name
	real(REAL_KIND), allocatable :: tstart(:)
	real(REAL_KIND), allocatable :: tend(:)
	real(REAL_KIND), allocatable :: conc(:)
	real(REAL_KIND), allocatable :: dose(:)
	logical, allocatable :: started(:)
	logical, allocatable :: ended(:)
end type

type LQ_type
	real(REAL_KIND) :: OER_am, OER_bm
	real(REAL_KIND) :: alpha_H, beta_H
	real(REAL_KIND) :: K_ms
end type

type(dist_type) :: divide_dist
type(occupancy_type), allocatable :: occupancy(:,:,:)
type(cell_type), allocatable :: cell_list(:)
type(treatment_type), allocatable :: protocol(:)
!type(boundary_type), pointer :: bdrylist

integer :: NX, NY, NZ
integer :: initial_count
integer, allocatable :: zdomain(:),zoffset(:)
integer :: blobrange(3,2)
real(REAL_KIND) :: x0,y0,z0   ! centre in global coordinates (units = grids)
real(REAL_KIND) :: blob_radius, Radius, Centre(3)
integer :: jumpvec(3,27)

integer :: max_nlist, nlist, Ncells, Ncells0, lastNcells, lastID
integer :: max_ngaps, ngaps, nadd_sites, Nsites, Nreuse
integer :: Ndrug_tag, Nradiation_tag, Nanoxia_tag, Ndrug_dead, Nradiation_dead, Nanoxia_dead
integer :: nbdry
integer :: istep, nsteps, NT_CONC, NT_GUI_OUT
integer :: Mnodes
real(REAL_KIND) :: DELTA_T, DELTA_X, fluid_fraction, Vsite, Vextra, medium_volume
real(REAL_KIND) :: MM_THRESHOLD, ANOXIA_FACTOR, t_anoxic_limit, anoxia_death_delay, Vdivide0, dVdivide
real(REAL_KIND) :: divide_time_median, divide_time_shape, divide_time_mean
real(REAL_KIND) :: t_simulation
type(SN30K_type) :: SN30K
type(LQ_type) :: LQ
character*(128) :: inputfile
character*(128) :: treatmentfile
character*(128) :: outputfile
character*(2048) :: logmsg
logical :: test_case(4)

TYPE(winsockport) :: awp_0, awp_1
logical :: use_TCP = .true.         ! turned off in para_main()
logical :: use_CPORT1 = .false.
logical :: stopped, clear_to_send
logical :: simulation_start, par_zig_init, initialized
logical :: use_radiation, use_treatment
logical :: use_V_dependence
logical :: dbug = .false.

integer :: divide_option = DIVIDE_USE_CLEAR_SITE
!integer :: divide_option = DIVIDE_ALWAYS_PUSH
integer :: idbug = 0
integer :: seed(2)

!DEC$ ATTRIBUTES DLLEXPORT :: nsteps

contains

!-----------------------------------------------------------------------------------------
! WTIME returns a reading of the wall clock time.
!-----------------------------------------------------------------------------------------
real(DP) function wtime()
!DEC$ ATTRIBUTES DLLEXPORT :: wtime
  integer :: clock_max, clock_rate, clock_reading

  call system_clock ( clock_reading, clock_rate, clock_max )
  wtime = real(clock_reading,kind=DP)/clock_rate
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine logger(msg)
character*(*) :: msg
integer :: error
logical :: isopen
character*(1) :: LF = char(94)

error = 0
if (use_TCP) then
    if (awp_0%is_open) then
        call winsock_send(awp_0,trim(msg)//LF,len_trim(msg)+1,error)
    else
        write(99,*) trim(msg)
    endif
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
Radius = (3.0*N/(4.0*PI))**(1./3.)
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
real(REAL_KIND) :: dummy
real(REAL_KIND) :: rsum,R
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
real(REAL_KIND) function norm(r)
real(REAL_KIND) :: r(3)

norm = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function norm2(r)
real(REAL_KIND) :: r(3)

norm2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine normalize(r)
real(REAL_KIND) :: r(3)

r = r/norm(r)
end subroutine

!-----------------------------------------------------------------------------------------
! Distance from the blob centre (units = grids)
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function cdistance(site)
integer :: site(3)
real(REAL_KIND) :: r(3)

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
real(REAL_KIND) :: v(3), vnorm(3)
real(REAL_KIND) :: d

d = dot_product(v,v)
vnorm = v/sqrt(d)
end subroutine

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(REAL_KIND) function DivideTime()
real(REAL_KIND) :: p1, p2
integer :: kpar = 0

dividetime = 0
p1 = divide_dist%p1
p2 = divide_dist%p2
select case (divide_dist%class)
case (NORMAL_DIST)
	DivideTime = rv_normal(p1,p2,kpar)
case (LOGNORMAL_DIST)
	DivideTime = rv_lognormal(p1,p2,kpar)
case (CONSTANT_DIST)
	DivideTime = p1
end select
!write(*,'(a,f8.2)') 'DivideTime: ',DivideTime/3600
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real function DivisionTime()
integer :: kpar = 0
real(REAL_KIND), parameter :: rndfraction = 0.2

DivisionTime = rv_lognormal(divide_dist%p1,divide_dist%p2,kpar)
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_normal(p1,p2,kpar)
integer :: kpar
real(REAL_KIND) :: p1,p2
real(REAL_KIND) :: R

R = par_rnor(kpar)
rv_normal = p1+R*p2
end function

!--------------------------------------------------------------------------------------
! When Y is normal N(p1,p2) then X = exp(Y) is lognormal with
!   median = m = exp(p1)
!   shape  = s = exp(p2)
! Also median = m = mean/(s^2/2)
! kpar = parallel process number
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_lognormal(p1,p2,kpar)
integer :: kpar
real(REAL_KIND) :: p1,p2
real(REAL_KIND) :: R,z

R = par_rnor(kpar)
z = p1 + R*p2
rv_lognormal = exp(z)
end function

!--------------------------------------------------------------------------------------
! For testing.
!--------------------------------------------------------------------------------------
real(REAL_KIND) function my_rnor()
real(REAL_KIND) :: sum, R
integer :: k
integer :: kpar=0

sum = 0
do k = 1,12
    R = par_uni(kpar)
    sum = sum + R
enddo
my_rnor = sum - 6.0
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_exponential(p1)
real(REAL_KIND) :: p1
real(REAL_KIND) :: r
integer :: kpar = 0

r = par_rexp(kpar)
rv_exponential = p1*r
end function

!--------------------------------------------------------------------------------------
! Cumulative probability distribution for a lognormal variate with median m, shape s
! Computes Pr(X < a) where X = exp(Y) and Y is N(p1,p2), p1 = log(m), p2 = log(s)
! Pr(X < a) = Pr(Y < log(a)) = Pr(p1 + p2*R < log(a)) = Pr(R < (log(a)-p1)/p2)
! where R is N(0,1)
!--------------------------------------------------------------------------------------
real(REAL_KIND) function cum_prob_lognormal(a,p1,p2)
real(REAL_KIND) :: a, p1, p2
real(REAL_KIND) :: b, prob

b = (log(a) - p1)/p2
prob = 0.5 + 0.5*erf(b/sqrt(2.0))
cum_prob_lognormal = prob
end function

end module
