! Global definitions

module global

use omp_lib
use real_kind_mod
use par_zig_mod
use winsock
use, intrinsic :: ISO_C_BINDING

implicit none

#include "itsol_interface.f90"															+++add to spheroid

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
integer, parameter :: UNREACHABLE_TAG  = -2

integer, parameter :: DIVIDE_ALWAYS_PUSH  = 1
integer, parameter :: DIVIDE_USE_CLEAR_SITE  = 2
integer, parameter :: DIVIDE_USE_CLEAR_SITE_RANDOM  = 3

integer, parameter :: nfin=10, nfout=11, nflog=12, nfres=13, nfrun=14, nfcell=15, nftreatment=16, nfprofile=17, nfslice=18
integer, parameter :: neumann(3,6) = reshape((/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /), (/3,6/))

integer, parameter :: CFSE = 0
integer, parameter :: OXYGEN = 1
integer, parameter :: GLUCOSE = 2
integer, parameter :: TRACER = 3
integer, parameter :: DRUG_A = 4
integer, parameter :: TPZ_DRUG = DRUG_A
integer, parameter :: TPZ_DRUG_METAB_1 = TPZ_DRUG + 1
integer, parameter :: TPZ_DRUG_METAB_2 = TPZ_DRUG + 2
integer, parameter :: DRUG_B = DRUG_A + 3
integer, parameter :: DNB_DRUG = DRUG_B
integer, parameter :: DNB_DRUG_METAB_1 = DNB_DRUG + 1
integer, parameter :: DNB_DRUG_METAB_2 = DNB_DRUG + 2
integer, parameter :: MAX_CHEMO = DRUG_B + 2
integer, parameter :: GROWTH_RATE = MAX_CHEMO + 1	! (not used here, used in the GUI)
integer, parameter :: CELL_VOLUME = MAX_CHEMO + 2
integer, parameter :: O2_BY_VOL = MAX_CHEMO + 3

integer, parameter :: N_EXTRA = O2_BY_VOL - MAX_CHEMO + 1	! = 4 = total # of variables - MAX_CHEMO

integer, parameter :: TPZ_CLASS = 1
integer, parameter :: DNB_CLASS = 2
integer, parameter :: DRUG_EVENT = 1
integer, parameter :: RADIATION_EVENT = 2
integer, parameter :: MEDIUM_EVENT = 3

integer, parameter :: DIST_NV = 20

integer, parameter :: EXTRA = 1
integer, parameter :: INTRA = 2
integer, parameter :: MAX_CELLTYPES = 4
integer, parameter :: MAX_DRUGTYPES = 2
integer, parameter :: max_nlist = 200000
integer, parameter :: NRF = 4
integer, parameter :: LIMIT_THRESHOLD = 1500

logical, parameter :: use_ODE_diffusion = .true.
logical, parameter :: compute_concentrations = .true.
logical, parameter :: use_division = .true.
logical, parameter :: use_death = .true.
logical, parameter :: use_react = .true.
logical, parameter :: use_migration = .false.	! causing an error with vacant site becoming bdry 
logical, parameter :: use_medium_flux = .true.	! flux of constituents between spheroid and medium is accounted for.
logical, parameter :: use_metabolites = .true.
logical, parameter :: use_celltype_colour = .true.

logical, parameter :: use_Cex_Cin = .true.		! assume equilibrium to derive Cin from Cex
logical, parameter :: suppress_growth = .false.

logical, parameter :: OFF_LATTICE = .false.

real(REAL_KIND), parameter :: PI = 4.0*atan(1.0)
real(REAL_KIND), parameter :: CFSE_std = 0.05
real(REAL_KIND), parameter :: small_d = 0.1e-4          ! 0.1 um -> cm

type occupancy_type
	integer :: indx(2)
	real(REAL_KIND) :: C(MAX_CHEMO)
	type(boundary_type), pointer :: bdry
end type

type cell_type
	integer :: ID
	integer :: celltype
	integer :: site(3)
	integer :: ivin
	logical :: active
	integer :: state
	integer :: generation
	real(REAL_KIND) :: conc(MAX_CHEMO)
	real(REAL_KIND) :: Cex(MAX_CHEMO)
	real(REAL_KIND) :: dCdt(MAX_CHEMO)
	real(REAL_KIND) :: dMdt(MAX_CHEMO)      ! mumol/s
	real(REAL_KIND) :: CFSE
	real(REAL_KIND) :: dVdt
	real(REAL_KIND) :: volume			! fractional volume (fraction of nominal cell volume Vcell_cm3)
	real(REAL_KIND) :: divide_volume	! fractional divide volume (normalised)
	real(REAL_KIND) :: divide_time
	real(REAL_KIND) :: t_divide_last	! these two values are used for colony simulation
	real(REAL_KIND) :: t_divide_next
	real(REAL_KIND) :: t_hypoxic
	real(REAL_KIND) :: t_anoxia_die
	real(REAL_KIND) :: M
	real(REAL_KIND) :: p_rad_death
	real(REAL_KIND) :: p_drug_death(MAX_DRUGTYPES)
	logical :: growth_delay
	real(REAL_KIND) :: dt_delay
	real(REAL_KIND) :: t_growth_delay_end			! this is for suppression of growth before first division
	integer :: N_delayed_cycles_left		! decremented by 1 at each cell division
	logical :: radiation_tag, anoxia_tag	!, drugA_tag, drugB_tag
	logical :: drug_tag(MAX_DRUGTYPES)
	logical :: G2_M
	logical :: exists
	integer :: cnr(3,8)
	real(REAL_KIND) :: wt(8)
end type

type drug_type
	character*(3)   :: classname
	integer         :: drugclass
	character*(16)  :: name
	integer         :: nmetabolites
	logical         :: use_metabolites
	real(REAL_KIND) :: diff_coef(0:2)
	real(REAL_KIND) :: medium_diff_coef(0:2)
	real(REAL_KIND) :: membrane_diff_in(0:2)
	real(REAL_KIND) :: membrane_diff_out(0:2)
	real(REAL_KIND) :: halflife(0:2)
	logical         :: kills(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Kmet0(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: C2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: KO2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: n_O2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Vmax(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Km(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Klesion(MAX_CELLTYPES,0:2)
	integer         :: kill_model(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: death_prob(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Kd(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_O2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_drug(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_duration(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_fraction(MAX_CELLTYPES,0:2)
	logical         :: sensitises(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: SER_max(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: SER_Km(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: SER_KO2(MAX_CELLTYPES,0:2)
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
!	real(c_double) :: CFSE
!	real(c_double) :: dVdt
	real(c_double) :: volume
	real(c_double) :: conc(0:MAX_CHEMO+N_EXTRA)	! Must agree with the definition in field.h: 0 = CFSE, MAX_CHEMO+1 = dVdt, MAX_CHEMO+2 = cellvolume
end type

type, bind(C) :: dist_data
	logical(c_bool) :: used
	real(c_double) :: dv
	real(c_double) :: v0
	real(c_double) :: prob(DIST_NV)
end type

type treatment_type
	integer :: ichemo
	integer :: n
!	character*(16) :: name
	real(REAL_KIND), allocatable :: tstart(:)
	real(REAL_KIND), allocatable :: tend(:)
	real(REAL_KIND), allocatable :: conc(:)
	real(REAL_KIND), allocatable :: dose(:)
	logical, allocatable :: started(:)
	logical, allocatable :: ended(:)
end type

type event_type
	integer :: etype
	real(REAL_KIND) :: time
	integer :: idrug			! DRUG
	integer :: ichemo			! DRUG CHEMO INDEX
	real(REAL_KIND) :: volume	! DRUG MEDIUM
	real(REAL_KIND) :: conc		! DRUG
	real(REAL_KIND) :: O2conc		! DRUG
	real(REAL_KIND) :: O2flush		! DRUG
	real(REAL_KIND) :: dose		! RADIATION
	real(REAL_KIND) :: O2medium	! MEDIUM
	logical :: done
end type	

type LQ_type
	real(REAL_KIND) :: OER_am, OER_bm
	real(REAL_KIND) :: alpha_H, beta_H
	real(REAL_KIND) :: K_ms
	real(REAL_KIND) :: death_prob
	real(REAL_KIND) :: growth_delay_factor
	real(REAL_KIND) :: growth_delay_N
end type

type savedata_type
    logical :: active
    character*(128) :: filebase
    real(REAL_KIND) :: dt
    integer :: nt, it
end type

type(dist_type) :: divide_dist(MAX_CELLTYPES)
type(occupancy_type), allocatable :: occupancy(:,:,:)
type(cell_type), allocatable, target :: cell_list(:)
type(treatment_type), allocatable :: protocol(:)
type(event_type), allocatable :: event(:)
real(REAL_KIND), allocatable, target :: Cslice(:,:,:,:)

character*(12) :: dll_version, dll_run_version
character*(12) :: gui_version, gui_run_version
integer :: NX, NY, NZ, NXB, NYB, NZB
integer :: ixb0, iyb0, izb0
integer :: initial_count
integer, allocatable :: zdomain(:),zoffset(:)
integer :: blobrange(3,2)
!real(REAL_KIND) :: Radius, Centre(3)		! sphere radius and centre
real(REAL_KIND) :: x0,y0,z0					! initial sphere centre in lattice coordinates (units = sites)
real(REAL_KIND) :: Centre_b(3)              ! sphere centre in coarse grid axes
real(REAL_KIND) :: xb0,yb0,zb0              ! sphere centre in coarse grid axes
real(REAL_KIND) :: blob_volume, blob_area       ! blob volume, max z-slice area (units = sites)
real(REAL_KIND) :: blob_centre(3), blob_radius  ! blob centre, radius in lattice coordinates (units = sites)
real(REAL_KIND) :: Vex_min, Vex_max

logical :: use_dropper
integer :: Ndrop
real(REAL_KIND) :: alpha_shape, beta_shape	! squashed sphere shape parameters
real(REAL_KIND) :: adrop, bdrop, cdrop		! drop shape transformation parameters
integer :: zmin     						! drop lower bound at drop time = lower limit of blob thereafter
logical :: is_dropped

integer :: jumpvec(3,27)

integer :: nlist, Ncells, Ncells0, lastNcells, lastID, Ncelltypes, Ncells_type(MAX_CELLTYPES)
integer :: diam_count_limit
logical :: limit_stop
integer :: nadd_sites, Nsites, Nreuse
integer :: Ndrugs_used
integer :: Nradiation_tag(MAX_CELLTYPES), Nanoxia_tag(MAX_CELLTYPES)
integer :: Ndrug_tag(MAX_DRUGTYPES,MAX_CELLTYPES)
integer :: Nradiation_dead(MAX_CELLTYPES), Nanoxia_dead(MAX_CELLTYPES)
integer :: Ndrug_dead(MAX_DRUGTYPES,MAX_CELLTYPES)
logical :: use_radiation_growth_delay_all = .true.
!logical :: radiation_dosed

logical :: drug_gt_cthreshold(MAX_DRUGTYPES)
real(REAL_KIND), parameter :: Cthreshold = 1.0e-5

type(savedata_type) :: saveprofile, saveslice

integer :: istep, nsteps, it_solve, NT_CONC, NT_GUI_OUT, show_progeny, ichemo_curr
integer :: Mnodes, ncpu_input
integer :: Nevents
real(REAL_KIND) :: DELTA_T, DELTA_X, fluid_fraction, Vsite_cm3, Vextra_cm3, Vcell_cm3, Vcell_pL
real(REAL_KIND) :: dxb, dxb3, dxf, dx3
real(REAL_KIND) :: medium_volume0, total_volume, cell_radius, d_layer, t_lastmediumchange
real(REAL_KIND) :: celltype_fraction(MAX_CELLTYPES)
logical :: celltype_display(MAX_CELLTYPES)
real(REAL_KIND) :: MM_THRESHOLD, ANOXIA_THRESHOLD, t_anoxic_limit, anoxia_death_delay, Vdivide0, dVdivide
real(REAL_KIND) :: divide_time_median(MAX_CELLTYPES), divide_time_shape(MAX_CELLTYPES), divide_time_mean(MAX_CELLTYPES)
real(REAL_KIND) :: t_simulation, execute_t1
real(REAL_KIND) :: O2cutoff(3), hypoxia_threshold
real(REAL_KIND) :: growthcutoff(3)
real(REAL_KIND) :: spcrad_value
real(REAL_KIND) :: total_dMdt
real(REAL_KIND) :: total_flux_prev, medium_Cbnd_prev
real(REAL_KIND) :: start_wtime

type(drug_type), allocatable, target :: drug(:)

integer, allocatable :: gaplist(:)
integer :: ngaps
integer, parameter :: max_ngaps = 10000

logical :: bdry_changed
type(LQ_type) :: LQ(MAX_CELLTYPES)
character*(128) :: inputfile
character*(128) :: treatmentfile
character*(128) :: outputfile
character*(2048) :: logmsg
character*(1024) :: header
logical :: test_case(4)

TYPE(winsockport) :: awp_0, awp_1
logical :: use_TCP = .true.         ! turned off in para_main()
logical :: use_CPORT1 = .false.
logical :: stopped, clear_to_send
logical :: simulation_start, par_zig_init, initialized
logical :: use_radiation, use_treatment
!logical :: use_growth_suppression = .true.	! see usage in subroutine CellGrowth
logical :: use_extracellular_O2
logical :: use_V_dependence
logical :: use_divide_time_distribution = .true.
logical :: use_constant_divide_volume = .true.
logical :: randomise_initial_volume
logical :: use_FD = .true.
logical :: use_gaplist = .true.
logical :: relax
logical :: use_parallel
logical :: dbug = .false.
logical :: bdry_debug

logical :: use_events = .true.

real(REAL_KIND) :: ysave(100000),dCreactsave(100000)

integer :: divide_option = DIVIDE_USE_CLEAR_SITE
!integer :: divide_option = DIVIDE_ALWAYS_PUSH
integer :: idbug = 0
integer :: Nbnd
integer :: seed(2)

! Off-lattice parameters, in the input file but unused here
real(REAL_KIND) :: a_separation
real(REAL_KIND) :: a_force, b_force, c_force, x0_force, x1_force, kdrag, frandom

real(REAL_KIND), allocatable :: omp_x(:), omp_y(:), omp_z(:)

!DEC$ ATTRIBUTES DLLEXPORT :: nsteps, DELTA_T

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
logical :: logfile_isopen
character*(1) :: LF = char(94)

error = 0
inquire(unit=nflog,OPENED=logfile_isopen)
if (use_TCP) then
    if (awp_0%is_open) then
        call winsock_send(awp_0,trim(msg)//LF,len_trim(msg)+1,error)
    elseif (logfile_isopen) then
        write(nflog,*) trim(msg)
    else
        write(99,*) trim(msg)
    endif
else
	write(*,*) trim(msg)
endif
if (logfile_isopen) then
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
! This calculates the radius of the equivalent sphere.
! Note: this assumes that N is an accurate measure of the number of sites
! occupied by the blob.
!---------------------------------------------------------------------
subroutine SetRadius(N)
integer :: N
blob_radius = (3.0*N/(4.0*PI))**(1./3.)
if (is_dropped) then
	z0 = zmin + (bdrop-cdrop)*blob_radius
	blob_centre(3) = z0
endif
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
! Make a random choice of an integer from 1 - N on the basis of probabilities in
! the array p(:) (assumed to be normalized).
!--------------------------------------------------------------------------------
integer function random_choice(p,N,kpar)
integer :: N,kpar
real(REAL_KIND) :: p(:)
integer :: k
real(REAL_KIND) :: R, psum

R = par_uni(kpar)
psum = 0
do k = 1,N
    psum = psum + p(k)
    if (R <= psum) then
        random_choice = k
        return
    endif
enddo
write(logmsg,*) 'ERROR: random_choice: ',N,p
call logger(logmsg)
stop
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

r = site - blob_centre
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
real(REAL_KIND) function DivideTime(ityp)
integer :: ityp
real(REAL_KIND) :: p1, p2
integer :: kpar = 0

dividetime = 0
p1 = divide_dist(ityp)%p1
p2 = divide_dist(ityp)%p2
select case (divide_dist(ityp)%class)
case (NORMAL_DIST)
	DivideTime = rv_normal(p1,p2,kpar)
case (LOGNORMAL_DIST)
	DivideTime = rv_lognormal(p1,p2,kpar)
case (CONSTANT_DIST)
	DivideTime = p1
end select
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real function DivisionTime(ityp)
integer :: ityp
integer :: kpar = 0
real(REAL_KIND), parameter :: rndfraction = 0.2

DivisionTime = rv_lognormal(divide_dist(ityp)%p1,divide_dist(ityp)%p2,kpar)
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

!-----------------------------------------------------------------------------------------
! Generate a random value for CFSE from a distribution with mean = average
! In the simplest case we can allow a uniform distribution about the average.
! Multiplying factor in the range (1-a, 1+a)
! Better to make it a Gaussian distribution: 
!  = average*(1+s*R)
! where R = N(0,1), s = std deviation
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function generate_CFSE(average)
real(REAL_KIND) :: average, std
integer :: kpar = 0
real(REAL_KIND) :: R

! Uniform distribution
!R = par_uni(kpar)
!generate_CFSE = (1 - a + 2*a*R)*average
! Gaussian distribution
R = par_rnor(kpar)	! N(0,1)
generate_CFSE = (1 + CFSE_std*R)*average
end function

!-----------------------------------------------------------------------------------------
! ityp = cell type
! V0 = cell starting volume (after division) = %volume
! Two approaches:
! 1. Use Vdivide0 and dVdivide to generate a volume
! 2. Use the divide time log-normal distribution
!    (a) use_V_dependence = true
!    (b) use_V_dependence = false
! NOTE: %volume and %divide_volume are normalised.
!-----------------------------------------------------------------------------------------
function get_divide_volume(ityp,V0,Tdiv) result(Vdiv)
integer :: ityp
real(REAL_KIND) :: V0, Tdiv
real(REAL_KIND) :: Vdiv
real(REAL_KIND) :: Tmean, b, R
integer :: kpar=0

Tmean = divide_time_mean(ityp)
if (use_divide_time_distribution) then
	Tdiv = DivideTime(ityp)
	if (use_constant_divide_volume) then
		Vdiv = Vdivide0
	else
		if (use_V_dependence) then
			b = log(2.0)*(Tdiv/Tmean)
			Vdiv = V0*exp(b)
		else
			Vdiv = V0 + (Vdivide0/2)*(Tdiv/Tmean)
		endif
	endif
else
	if (use_constant_divide_volume) then
		Vdiv = Vdivide0
	else
		R = par_uni(kpar)
		Vdiv = Vdivide0 + dVdivide*(2*R-1)
	endif
	Tdiv = Tmean
endif
end function	

!--------------------------------------------------------------------------------------
! Determine real roots r(:) of the cubic equation:
! x^3 + a.x^2 + b.x + c = 0
! If there is one real root, n=1 and the root is r(1)
! If there are three distinct real roots, n=3 and the roots are r(1), r(2), r(3)
! If there is a repeated root, n=2 and the single root is r(1), the repeated root is r(2)
!--------------------------------------------------------------------------------------
subroutine cubic_roots(a, b, c, r, n)
real(REAL_KIND) :: a, b, c, r(3)
integer :: n
real(REAL_KIND) :: QQ, RR, theta, R2, Q3, AA, BB

QQ = (a*a - 3*b)/9
RR = (2*a*a*a - 9*a*b + 27*c)/54
Q3 = QQ*QQ*QQ
R2 = RR*RR
if (R2 < Q3) then
	n = 3
	theta = acos(RR/sqrt(Q3))
	r(1) = -2*sqrt(QQ)*cos(theta/3) - a/3
	r(2) = -2*sqrt(QQ)*cos((theta+2*PI)/3) - a/3
	r(3) = -2*sqrt(QQ)*cos((theta-2*PI)/3) - a/3
else
	n = 1
	AA = -sign(1.d0,RR)*(abs(RR) + sqrt(R2 - Q3))**(1.d0/3.d0)
	if (AA == 0) then
		BB = 0
	else
		BB = QQ/AA
	endif
	r(1) = AA + BB - a/3
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Squeeze gaps out of cellist array, adjusting occupancy array.
!-----------------------------------------------------------------------------------------
subroutine squeezer()
integer :: last, kcell, site(3), indx(2), i, j, idc, n, region

!write(*,*) 'squeezer'
!call logger('squeezer')
if (ngaps == 0) return
last = nlist
kcell = 0
n = 0
do
    kcell = kcell+1
    if (cell_list(kcell)%state == DEAD) then    ! a gap
        if (kcell == last) exit
        do
            if (last == 0) then
                write(nflog,*) 'last = 0: kcell: ',kcell
                stop
            endif
            if (cell_list(last)%state == DEAD) then
                last = last-1
                n = n+1
                if (n == ngaps) exit
            else
                exit
            endif
        enddo
        if (n == ngaps) exit
!        call copycell2cell(cell_list(last),cell_list(kcell),kcell) 
        cell_list(kcell) = cell_list(last)
        site = cell_list(last)%site
!        indx = occupancy(site(1),site(2),site(3))%indx
!        do i = 1,2
!            if (indx(i) == last) indx(i) = kcell
!        enddo
		if (occupancy(site(1),site(2),site(3))%indx(1) /= last) then
			write(nflog,*) 'squeezer: bad occupancy: ',occupancy(site(1),site(2),site(3))%indx(1),last
			stop
		endif
        occupancy(site(1),site(2),site(3))%indx(1) = kcell
        last = last-1
        n = n+1
    endif
    if (n == ngaps) exit
enddo
nlist = nlist - ngaps
ngaps = 0
if (dbug) write(nflog,*) 'squeezed: ',n,nlist

end subroutine

end module
