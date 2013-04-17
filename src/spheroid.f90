module spheroid_mod
use global
use boundary
!use behaviour 
!use ode_diffuse
!use FDC
use fields
use cellstate
use winsock  

IMPLICIT NONE

contains 

!-----------------------------------------------------------------------------------------
! This subroutine is called to initialize a simulation run.
! ncpu = the number of processors to use
! infile = file with the input data
! outfile = file to hold the output
! runfile = file to pass info to the master program (e.g. Python) as the program executes.
!-----------------------------------------------------------------------------------------
subroutine setup(ncpu,infile,outfile,ok)
integer :: ncpu
character*(*) :: infile, outfile
logical :: ok
character*(64) :: msg
integer :: error

ok = .true.
initialized = .false.
par_zig_init = .false.
Mnodes = ncpu
inputfile = infile
outputfile = outfile
write(logmsg,*) 'ncpu: ',ncpu 
call logger(logmsg)

#if defined(OPENMP) || defined(_OPENMP)
    call logger("OPENMP defined")
    call omp_initialisation(ok)
    if (.not.ok) return
#else
    call logger("OPENMP NOT defined")
    if (Mnodes > 1) then
        write(logmsg,'(a)') 'No OpenMP, using one thread only'
        call logger(logmsg)
        Mnodes = 1
    endif
#endif

call logger("read_cell_params")
call read_cell_params(ok)
if (.not.ok) return
call logger("did read_cell_params")

call array_initialisation(ok)
if (.not.ok) return
call logger('did array_initialisation')

call SetupChemo

call PlaceCells(ok)
call SetRadius(Nsites)
write(logmsg,*) 'did placeCells: Ncells: ',Ncells,Radius
call logger(logmsg)
if (.not.ok) return

!call test_get_path

call CreateBdryList

!chemo_N = 8
!call ChemoSetup

call make_split(.true.)

!call init_counters
!if (save_input) then
!    call save_inputfile(inputfile)
!    call save_parameters
!	call save_inputfile(treatmentfile)
!endif
!
!call AllocateConcArrays
!
!call ChemoSteadystate
!
!firstSummary = .true.
!initialized = .true.
!
!call checkcellcount(ok)

if (use_ODE_diffusion) then
!	call SetupODEDiff
	call InitConcs
	call AdjustMM
!	call TestODEDiffusion
!	call TestSolver
endif
call SetupODEdiff
Nradiation_tag = 0
Ndrug_tag = 0
Nanoxia_tag = 0
Nradiation_dead = 0
Ndrug_dead = 0
Nanoxia_dead = 0
t_simulation = 0
istep = 0
write(logmsg,'(a,i6)') 'Startup procedures have been executed: initial T cell count: ',Ncells0
call logger(logmsg)

!call check_Coxygen

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine omp_initialisation(ok)
logical :: ok
integer :: npr, nth

ok = .true.
if (Mnodes == 1) return
#if defined(OPENMP) || defined(_OPENMP)
write(logmsg,'(a,i2)') 'Requested Mnodes: ',Mnodes
call logger(logmsg)
npr = omp_get_num_procs()
write(logmsg,'(a,i2)') 'Machine processors: ',npr
call logger(logmsg)

nth = omp_get_max_threads()
write(logmsg,'(a,i2)') 'Max threads available: ',nth
call logger(logmsg)
if (nth < Mnodes) then
    Mnodes = nth
    write(logmsg,'(a,i2)') 'Setting Mnodes = max thread count: ',nth
	call logger(logmsg)
endif

call omp_set_num_threads(Mnodes)
!$omp parallel
nth = omp_get_num_threads()
write(logmsg,*) 'Threads, max: ',nth,omp_get_max_threads()
call logger(logmsg)
!$omp end parallel
#endif

call logger('did omp_initialisation')
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine array_initialisation(ok)
logical :: ok
integer :: x,y,z,k, ichemo
integer :: MAXX, z1, z2, nc0, inflow
integer :: cog_size
real(REAL_KIND) :: d, rr(3)

ok = .false.
call logger("call rng_initialisation")
call rng_initialisation
call logger("did rng_initialisation")

! These are deallocated here instead of in subroutine wrapup so that when a simulation run ends 
! it will still be possible to view the cell distributions and chemokine concentration fields.
if (allocated(occupancy)) deallocate(occupancy)
if (allocated(cell_list)) deallocate(cell_list)
if (allocated(allstate)) deallocate(allstate)
if (allocated(ODEdiff%ivar)) deallocate(ODEdiff%ivar)
!do ichemo = 1,MAX_CHEMO
!	if (allocated(chemo(ichemo)%coef)) deallocate(chemo(ichemo)%coef)
!	if (allocated(chemo(ichemo)%conc)) deallocate(chemo(ichemo)%conc)
!	if (allocated(chemo(ichemo)%grad)) deallocate(chemo(ichemo)%grad)
!enddo
call logger('did deallocation')

!nsteps_per_min = 1.0/DELTA_T
NY = NX
NZ = NX
ngaps = 0
max_ngaps = NY*NZ
nlist = 0

allocate(zoffset(0:2*Mnodes))
allocate(zdomain(NZ))
x0 = (NX + 1.0)/2.        ! global value
y0 = (NY + 1.0)/2.
z0 = (NZ + 1.0)/2.
Centre = (/x0,y0,z0/)   ! now, actually the global centre (units = grids)
Radius = blob_radius    ! starting value

nc0 = (4./3.)*PI*Radius**3
!max_nlist = 200*nc0
max_nlist = 200000
write(logmsg,*) 'Initial radius, nc0, max_nlist: ',Radius, nc0, max_nlist
call logger(logmsg)

allocate(cell_list(max_nlist))
allocate(occupancy(NX,NY,NZ))
!allocate(gaplist(max_ngaps))

call make_jumpvec

lastNcells = 0
nadd_sites = 0

ok = .true.

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine rng_initialisation
integer, allocatable :: zig_seed(:)
integer :: i
integer :: npar, grainsize = 32

npar = Mnodes
write(logmsg,*) 'npar = ',npar,seed
call logger(logmsg)
allocate(zig_seed(0:npar-1))
do i = 0,npar-1
    zig_seed(i) = seed(1)*seed(2)*(i+1)
enddo
call par_zigset(npar,zig_seed,grainsize)
par_zig_init = .true.
end subroutine

!----------------------------------------------------------------------------------------1123
!----------------------------------------------------------------------------------------
subroutine read_cell_params(ok)
logical :: ok
integer :: i, idrug, imetab, itestcase, ncpu_dummy, Nmm3, ichemo, itreatment, iuse_extra
integer :: iuse_drug, iuse_metab, idrug_decay, imetab_decay, iV_depend, iV_random
real(REAL_KIND) :: days, bdry_conc
real(REAL_KIND) :: sigma, DXmm, anoxia_tag_hours, anoxia_death_hours
character*(12) :: drug_name

ok = .true.
chemo(:)%used = .false.

open(nfcell,file=inputfile,status='old')

read(nfcell,*) NX							! rule of thumb: about 4*BLOB_RADIUS
read(nfcell,*) initial_count				! initial number of tumour cells
read(nfcell,*) divide_time_median
read(nfcell,*) divide_time_shape
read(nfcell,*) iV_depend
read(nfcell,*) iV_random
read(nfcell,*) days							! number of days to simulate
read(nfcell,*) DELTA_T						! time step size (sec)
read(nfcell,*) NT_CONC						! number of subdivisions of DELTA_T for diffusion computation
read(nfcell,*) Nmm3							! number of cells/mm^3
read(nfcell,*) fluid_fraction				! fraction of the (non-necrotic) tumour that is fluid
read(nfcell,*) medium_volume				! volume of medium that the spheroid is growing in (cm^3)
read(nfcell,*) Vdivide0						! nominal cell volume multiple for division
read(nfcell,*) dVdivide						! variation about nominal divide volume
read(nfcell,*) MM_THRESHOLD					! O2 concentration threshold Michaelis-Menten "soft-landing" (uM)
read(nfcell,*) ANOXIA_FACTOR			    ! multiplying factor for MM threshold for anoxia
read(nfcell,*) anoxia_tag_hours				! hypoxic time leading to tagging to die by anoxia (h)
read(nfcell,*) anoxia_death_hours			! time after tagging to death by anoxia (h)
read(nfcell,*) itestcase                    ! test case to simulate
read(nfcell,*) seed(1)						! seed vector(1) for the RNGs
read(nfcell,*) seed(2)						! seed vector(2) for the RNGs
read(nfcell,*) ncpu_dummy					! just a placeholder for ncpu, not used currently
read(nfcell,*) NT_GUI_OUT					! interval between GUI outputs (timesteps)
read(nfcell,*) chemo(OXYGEN)%used
read(nfcell,*) chemo(OXYGEN)%diff_coef
read(nfcell,*) chemo(OXYGEN)%cell_diff
read(nfcell,*) chemo(OXYGEN)%bdry_conc
read(nfcell,*) chemo(OXYGEN)%max_cell_rate
read(nfcell,*) chemo(GLUCOSE)%used
read(nfcell,*) chemo(GLUCOSE)%diff_coef
read(nfcell,*) chemo(GLUCOSE)%cell_diff
read(nfcell,*) chemo(GLUCOSE)%bdry_conc
read(nfcell,*) chemo(GLUCOSE)%max_cell_rate

do i = 1,2			! currently allowing for just two different drugs
	read(nfcell,*) iuse_drug
	read(nfcell,'(a12)') drug_name
	read(nfcell,*) bdry_conc
	read(nfcell,*) idrug_decay
	read(nfcell,*) iuse_metab
	read(nfcell,*) imetab_decay
!	if (iuse_drug == 0) cycle
	call get_indices(drug_name, idrug, imetab)
	if (idrug < 0 .and. iuse_drug /= 0) then
		write(logmsg,*) 'Unrecognized drug name: ',drug_name
		call logger(logmsg)
		ok = .false.
		return
	endif
	if (idrug == 0) cycle
	chemo(idrug)%used = (iuse_drug == 1)
	chemo(idrug)%bdry_conc = bdry_conc
	chemo(idrug)%decay = (idrug_decay == 1)
	if (chemo(idrug)%used) then
		chemo(imetab)%used = (iuse_metab == 1)
	else
		chemo(imetab)%used = .false.
	endif
	chemo(imetab)%decay = (imetab_decay == 1)
	if (idrug == SN30000) then
		read(nfcell,*) SN30K%diff_coef
		read(nfcell,*) SN30K%cell_diff
		read(nfcell,*) SN30K%Kmet0
		read(nfcell,*) SN30K%C1
		read(nfcell,*) SN30K%C2
		read(nfcell,*) SN30K%KO2
		read(nfcell,*) SN30K%gamma
		read(nfcell,*) SN30K%Klesion
		read(nfcell,*) SN30K%halflife
		read(nfcell,*) SN30K%metabolite_halflife
		read(nfcell,*) SN30K%kill_O2
		read(nfcell,*) SN30K%kill_drug
		read(nfcell,*) SN30K%kill_duration
		read(nfcell,*) SN30K%kill_fraction
		SN30K%KO2 = 1.0e-3*SN30K%KO2                    ! um -> mM
		SN30K%kill_duration = 60*SN30K%kill_duration    ! minutes -> seconds
		chemo(idrug)%halflife = SN30K%halflife
		chemo(imetab)%halflife = SN30K%metabolite_halflife
		chemo(idrug)%diff_coef = SN30K%diff_coef		! Note that metabolite is given the same diff_coef
		chemo(imetab)%diff_coef = SN30K%diff_coef		! and cell_diff as the drug.
		chemo(idrug)%cell_diff = SN30K%cell_diff
		chemo(imetab)%cell_diff = SN30K%cell_diff
	endif
	if (chemo(idrug)%used .and. chemo(idrug)%decay) then
		chemo(idrug)%decay_rate = DecayRate(chemo(idrug)%halflife)
	else
		chemo(idrug)%decay_rate = 0
	endif
	if (chemo(imetab)%used .and. chemo(imetab)%decay) then
		chemo(imetab)%decay_rate = DecayRate(chemo(imetab)%halflife)
	else
		chemo(imetab)%decay_rate = 0
	endif
enddo
read(nfcell,*) O2cutoff(1)
read(nfcell,*) O2cutoff(2)
read(nfcell,*) O2cutoff(3)
read(nfcell,*) spcrad_value
read(nfcell,*) iuse_extra
read(nfcell,*) itreatment
read(nfcell,*) treatmentfile						! file with treatment programme
close(nfcell)

MM_THRESHOLD = MM_THRESHOLD/1000					! uM -> mM
O2cutoff = O2cutoff/1000							! uM -> mM
blob_radius = (initial_count*3./(4.*PI))**(1./3)	! units = grids
divide_dist%class = LOGNORMAL_DIST
divide_time_median = 60*60*divide_time_median		! hours -> seconds
sigma = log(divide_time_shape)
!divide_dist%p1 = log(divide_time_mean/exp(sigma*sigma/2))	
divide_dist%p1 = log(divide_time_median)	
divide_dist%p2 = sigma
divide_time_mean = exp(divide_dist%p1 + 0.5*divide_dist%p2**2)	! mean = median.exp(sigma^2/2)
write(logmsg,'(a,2e12.4)') 'shape, sigma: ',divide_time_shape,sigma
call logger(logmsg)
write(logmsg,'(a,2e12.4)') 'Median, mean divide time: ',divide_time_median/3600,divide_time_mean/3600
call logger(logmsg)
use_V_dependence = (iV_depend == 1)
randomise_initial_volume = (iV_random == 1)
use_extracellular_O2 = (iuse_extra == 1)
t_anoxic_limit = 60*60*anoxia_tag_hours				! hours -> seconds
anoxia_death_delay = 60*60*anoxia_death_hours		! hours -> seconds
DXmm = 1.0/(Nmm3**(1./3))
DELTA_X = DXmm/10									! mm -> cm
Vsite = DELTA_X*DELTA_X*DELTA_X						! total site volume (cm^3)
Vextra = fluid_fraction*Vsite						! extracellular volume in a site

write(logmsg,'(a,e12.4)') 'DELTA_X: ',DELTA_X
call logger(logmsg)
write(logmsg,'(a,3e12.4)') 'Volumes: site, extra, cell (average): ',Vsite, Vextra, Vsite-Vextra
call logger(logmsg)

if (itreatment == 1) then
	use_treatment = .true.
	call read_treatment(ok)
	if (.not.ok) then
		write(logmsg,'(a,a)') 'Error reading treatment programme file: ',treatmentfile
		call logger(logmsg)
		return
	endif
else	
	use_treatment = .false.
!	do ichemo = DRUG_A,MAX_CHEMO
!		chemo(ichemo)%used = (iuse(ichemo) == 1)
!	enddo
!	do ichemo = DRUG_A,MAX_CHEMO
!		if (idecay(ichemo) == 1) then
!			chemo(ichemo)%decay = .true.
!!			chemo(ichemo)%bdry_decay_rate = DecayRate(chemo(ichemo)%bdry_halflife)
!			chemo(ichemo)%decay_rate = DecayRate(chemo(ichemo)%halflife)
!		else
!			chemo(ichemo)%decay = .false.
!			chemo(ichemo)%decay_rate = 0
!		endif
!	enddo
endif

!if (use_metabolites) then
!	do ichemo = DRUG_A,MAX_CHEMO
!		if (chemo(ichemo)%used .and. imetabolite(ichemo) == 1) then
!			chemo(ichemo+1)%used = .true.	! the metabolite is also simulated
!		endif
!	enddo
!endif

! Setup test_case
test_case = .false.
if (itestcase /= 0) then
    test_case(itestcase) = .true.
endif

if (mod(NX,2) /= 0) NX = NX+1					! ensure that NX is even
open(nfout,file=outputfile,status='replace')
write(logmsg,*) 'Opened nfout: ',outputfile
call logger(logmsg)

Nsteps = days*24*60*60/DELTA_T		! DELTA_T in seconds
write(logmsg,'(a,2i6,f6.0)') 'nsteps, NT_CONC, DELTA_T: ',nsteps,NT_CONC,DELTA_T
call logger(logmsg)

call determine_Kd
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_indices(drug_name, idrug, imetab)
character*(12) :: drug_name
integer :: idrug, imetab

if (drug_name == ' ') then
	idrug = 0
	imetab = 0
elseif (drug_name == 'SN30000') then
	idrug = SN30000
	imetab = SN30000_METAB
else
	idrug = -1
	imetab = 0
endif
end subroutine

!-----------------------------------------------------------------------------------------
! This overrides the drug usage specified by USE_DRUG_A etc
! RADIATION -> 0
! DRUG_A -> 1
! DRUG_B -> 2
! By default this assumes two drugs.
! Times are hours
! Drug concs are mM
! Radiation dose is Gy
!-----------------------------------------------------------------------------------------
subroutine read_treatment(ok)
logical :: ok
character*(64) :: line
integer :: ichemo, idrug, nmax, i
real(REAL_KIND) :: tstart,tend,conc,dose

allocate(protocol(0:2))

use_treatment = .true.
chemo(GLUCOSE+1:)%used = .false.
use_radiation = .false.
open(nftreatment,file=treatmentfile,status='old')
nmax = 0
do
	read(nftreatment,'(a)',end=99) line
	if (line(1:3) == 'END' .or. line(1:3) == 'end') exit
	if (line(1:6) == 'DRUG_A') then
		ichemo = DRUG_A
		idrug = 1
	elseif (line(1:6) == 'DRUG_B') then
		ichemo = DRUG_B
		idrug = 2
	elseif (line(1:9) == 'RADIATION') then
		ichemo = 0
		idrug = 0
	endif
	if (ichemo > 0) then
		read(nftreatment,'(a)') chemo(ichemo)%name
		read(nftreatment,*) protocol(idrug)%n
		nmax = max(nmax,protocol(idrug)%n)
		if (protocol(idrug)%n > 0) then
			do i = 1,protocol(idrug)%n
				read(nftreatment,*) tstart
				read(nftreatment,*) tend
				read(nftreatment,*) conc
			enddo
		endif
	else
		use_radiation = .true.
		read(nftreatment,*) protocol(idrug)%n
		nmax = max(nmax,protocol(idrug)%n)
		do i = 1,protocol(idrug)%n
			read(nftreatment,*) tstart
			read(nftreatment,*) dose
		enddo
	endif
	cycle
99	exit
enddo
rewind(nftreatment)
write(*,*) 'nmax: ',nmax
do idrug = 0,2
	allocate(protocol(idrug)%tstart(nmax))
	allocate(protocol(idrug)%tend(nmax))
	allocate(protocol(idrug)%conc(nmax))
	allocate(protocol(idrug)%dose(nmax))
	allocate(protocol(idrug)%started(nmax))
	allocate(protocol(idrug)%ended(nmax))
enddo
do
	read(nftreatment,'(a)',end=199) line
	if (line(1:3) == 'END' .or. line(1:3) == 'end') exit
	if (line(1:6) == 'DRUG_A') then
		ichemo = DRUG_A
		idrug = 1
	elseif (line(1:6) == 'DRUG_B') then
		ichemo = DRUG_B
		idrug = 2
	elseif (line(1:9) == 'RADIATION') then
		ichemo = 0
		idrug = 0
	endif
	if (ichemo > 0) then
		read(nftreatment,'(a)') chemo(ichemo)%name
		read(nftreatment,*) protocol(idrug)%n
		if (protocol(idrug)%n > 0) then
			chemo(ichemo)%used = .true.
			do i = 1,protocol(idrug)%n
				read(nftreatment,*) tstart
				protocol(idrug)%tstart(i) = 3600*tstart		! hours -> seconds
				read(nftreatment,*) tend
				protocol(idrug)%tend(i) = 3600*tend
				read(nftreatment,*) protocol(idrug)%conc(i)
			enddo
		endif
	else
		use_radiation = .true.
		read(nftreatment,*) protocol(idrug)%n
		do i = 1,protocol(idrug)%n
			read(nftreatment,*) tstart
			protocol(idrug)%tstart(i) = 3600*tstart
			read(nftreatment,*) protocol(idrug)%dose(i)
		enddo
	endif
	cycle
199	exit
enddo
close(nftreatment)
do i = 0,2
	protocol(i)%started = .false.
	protocol(i)%ended = .false.
enddo	
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
! The rate of cell killing, which becomes a cell death probability rate, is inferred from
! the cell kill experiment.
! The basic assumption is that the rate of killing is proportional to the drug metabolism rate:
! killing rate = c = Kd.dM/dt
! where dM/dt = F(O2).kmet0.Ci
! In the kill experiment both O2 and Ci are held constant:
! O2 = CkillO2, Ci = Ckill
! In this case c is constant and the cell population N(t) is given by:
! N(t) = N(0).exp(-ct), i.e.
! c = -log(N(T)/N(0))/T where T is the duration of the experiment
! N(T)/N(0) = 1 - f, where f = kill fraction
! c = Kd.F(CkillO2).kmet0.Ckill
! therefore
! Kd = -log(1-f)/(T.F(CkillO2).kmet0.Ckill)
!-----------------------------------------------------------------------------------------
subroutine determine_Kd
real(REAL_KIND) :: kmet

kmet = (SN30K%C1 + SN30K%C2*SN30K%KO2/(SN30K%KO2 + SN30K%kill_O2))*SN30K%Kmet0
SN30K%Kd = -log(1-SN30K%kill_fraction)/(SN30K%kill_duration*kmet*SN30K%kill_drug)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine placeCells(ok)
logical :: ok
integer :: x, y, z, k, site(3), ichemo, kpar=0
real(REAL_KIND) :: r2lim,r2,rad(3)
real(REAL_KIND) :: R, tpast, tdiv

occupancy(:,:,:)%indx(1) = 0
occupancy(:,:,:)%indx(2) = 0
r2lim = Radius*Radius
lastID = 0
k = 0
do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
			rad = (/x-x0,y-y0,z-z0/)
			r2 = dot_product(rad,rad)
			if (r2 < r2lim) then
				k = k+1
				lastID = lastID + 1
				site = (/x,y,z/)
				cell_list(k)%ID = lastID
				cell_list(k)%site = site
				cell_list(k)%state = 1
                cell_list(k)%drug_tag = .false.
                cell_list(k)%radiation_tag = .false.
                cell_list(k)%anoxia_tag = .false.
				cell_list(k)%exists = .true.
!				do
!					R = par_uni(kpar)
!					tpast = -R*divide_time_median
!					tdiv = DivideTime()
!					if (tdiv + tpast > 0) exit
!				enddo
!				cell_list(k)%divide_volume = Vdivide0
				R = par_uni(kpar)
				cell_list(k)%divide_volume = Vdivide0 + dVdivide*(2*R-1)
				R = par_uni(kpar)
				if (randomise_initial_volume) then
					cell_list(k)%volume = Vdivide0*0.5*(1 + R)
				else
					cell_list(k)%volume = 1.0
				endif
!				write(nflog,'(i6,2f8.4)') k,R,cell_list(k)%volume
				cell_list(k)%t_divide_last = 0		! not used
!				cell_list(k)%t_divide_next = tdiv + tpast
				cell_list(k)%t_hypoxic = 0
				cell_list(k)%conc = 0
				cell_list(k)%conc(OXYGEN) = chemo(OXYGEN)%bdry_conc
				cell_list(k)%conc(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
				cell_list(k)%M = 0
				occupancy(x,y,z)%indx(1) = k
			else
				occupancy(x,y,z)%indx(1) = OUTSIDE_TAG
			endif
		enddo
	enddo
enddo
do ichemo = 1,MAX_CHEMO
    occupancy(:,:,:)%C(ichemo) = 0
enddo
occupancy(:,:,:)%C(OXYGEN) = chemo(OXYGEN)%bdry_conc
occupancy(:,:,:)%C(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
nlist = k
Nsites = k
Ncells = k
Ncells0 = Ncells
Nreuse = 0	
ok = .true.
write(logmsg,*) 'idbug: ',idbug
call logger(logmsg)
end subroutine

!--------------------------------------------------------------------------------
! Generates the arrays zoffset() and zdomain().
! The domains (slices) are numbered 0,...,2*Mnodes-1
! wz(k) = width of the slice for kth domain
! zoffset(k) = offset of kth domain occupancy array in the occupancy array.
! zdomain(x) = domain that global z lies in.
! The kth domain (slice) extends from z = zoffset(k)+1 to z = zoffset(k+1)
! The idea is to set the domain boundaries such that each domain has roughly the
! same number of available sites.
! This is the initial split, which will continue to be OK if:
! not using a blob, or Mnodes <= 2
! blobrange(:,:) holds the info about the ranges of x, y and z that the blob occupies.
! blobrange(1,1) <= x <= blobrange(1,2)
! blobrange(2,1) <= y <= blobrange(2,2)
! blobrange(3,1) <= z <= blobrange(3,2)
!--------------------------------------------------------------------------------
subroutine make_split(force)
logical :: force
integer :: k, wsum, kdomain, nsum, Ntot, N, last, x, y, z
integer, allocatable :: scount(:)
integer, allocatable :: wz(:), ztotal(:)
integer :: Mslices
real(REAL_KIND) :: dNT, diff1, diff2
logical :: show = .false.

!write(*,*) 'make_split: istep,Mnodes: ',istep,Mnodes
if (Mnodes == 1) then
    Mslices = 1
    zdomain = 0
else
	Mslices = 2*Mnodes
endif
dNT = abs(Ncells - lastNcells)/real(lastNcells+1)
if (.not.force .and. dNT < 0.03) then
    return
endif
lastNcells = Ncells
if (Mslices > 1) then
	allocate(wz(0:Mslices))
	allocate(ztotal(0:Mslices))
	allocate(scount(NX))
endif
blobrange(:,1) = 99999
blobrange(:,2) = 0
nsum = 0
do z = 1,NZ
    k = 0
    do y = 1,NY
        do x = 1,NX
            if (occupancy(x,y,z)%indx(1) /= OUTSIDE_TAG) then
                k = k + 1
                blobrange(1,1) = min(blobrange(1,1),x)
                blobrange(1,2) = max(blobrange(1,2),x)
                blobrange(2,1) = min(blobrange(2,1),y)
                blobrange(2,2) = max(blobrange(2,2),y)
                blobrange(3,1) = min(blobrange(3,1),z)
                blobrange(3,2) = max(blobrange(3,2),z)
            endif
        enddo
    enddo
    if (Mslices > 1) then
	    scount(z) = k
	    nsum = nsum + scount(z)
	endif
enddo
if (Mslices == 1) return

Ntot = nsum
N = Ntot/Mslices
nsum = 0
last = 0
k = 0
do z = 1,NZ
    nsum = nsum + scount(z)
    if (nsum >= (k+1)*N) then
        diff1 = nsum - (k+1)*N
        diff2 = diff1 - scount(z)
        if (abs(diff1) < abs(diff2)) then
            wz(k) = z - last
            last = z
        else
            wz(k) = z - last - 1
            last = z - 1
        endif
        k = k+1
        if (k == Mslices-1) exit
    endif
enddo
wz(Mslices-1) = NZ - last
if (show) then
    write(*,*) 'Ntot, N: ',Ntot,N
    write(*,'(10i6)') scount
endif
zoffset(0) = 0
do k = 1,Mslices-1
    zoffset(k) = zoffset(k-1) + wz(k-1)
enddo
zoffset(Mslices) = NZ
z = 0
do kdomain = 0,Mslices-1
    do k = 1,wz(kdomain)
        z = z+1
        zdomain(z) = kdomain      ! = kpar with two sweeps
    enddo
enddo
if (show) then
    write(*,*) 'zoffset: ',zoffset
    write(*,*) 'wz:      ',wz
    write(*,*) 'zdomain: '
    write(*,'(10i4)') zdomain
endif
ztotal = 0
do k = 0,2*Mnodes-1
    do z = zoffset(k)+1,zoffset(k+1)
        ztotal(k) = ztotal(k) + scount(z)
    enddo
    if (show) write(*,*) k,ztotal(k)
enddo
deallocate(wz)
deallocate(ztotal)
deallocate(scount)
end subroutine

!--------------------------------------------------------------------------------
! Makes an approximate count of the number of sites of the spherical blob that
! are in the xth slice.  Uses the area of the slice.
! The blob centre is at (x0,y0,z0), and the blob radius is R = Radius%x
! NOT USED
!--------------------------------------------------------------------------------
integer function slice_count(x)
integer :: x
real(REAL_KIND) :: r2

r2 = Radius**2 - (x-x0)**2
if (r2 < 0) then
    slice_count = 0
else
    slice_count = PI*r2
endif
end function

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_dimensions(NX_dim, NY_dim, NZ_dim, nsteps_dim, deltat, maxchemo, cused) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dimensions
use, intrinsic :: iso_c_binding
integer(c_int) :: NX_dim,NY_dim,NZ_dim,nsteps_dim, maxchemo
real(c_double) :: deltat
logical(c_bool) :: cused(*)
integer :: ichemo

NX_dim = NX
NY_dim = NY
NZ_dim = NZ
nsteps_dim = nsteps
deltat = DELTA_T
maxchemo = MAX_CHEMO
do ichemo = 1,MAX_CHEMO
	cused(ichemo) = chemo(ichemo)%used
enddo
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine InitConcs
integer :: nextra, ic, ichemo, kcell, site(3)

write(logmsg,*) 'InitConcs: ',nchemo
call logger(logmsg)
allocate(allstate(MAX_VARS,MAX_CHEMO))
allocate(work_rkc(8+4*MAX_VARS,MAX_CHEMO))

do kcell = 1,Ncells
    site = cell_list(kcell)%site
    do ic = 1,nchemo
	    ichemo = chemomap(ic)
        occupancy(site(1),site(2),site(3))%C(ic) = chemo(ichemo)%bdry_conc
    enddo
enddo
!do ic = 1,nchemo
!	ichemo = chemomap(ic)
!	call InitState(ichemo,allstate(1:nextra,ichemo))
!enddo
end subroutine

!----------------------------------------------------------------------------------
! Radiation treatment is stored in protocol(0)
! 
!----------------------------------------------------------------------------------
subroutine treatment(radiation_dose)
real(REAL_KIND) :: radiation_dose
integer :: i, idrug, ichemo, ichemo_metab

radiation_dose = 0
do i = 1,protocol(0)%n
	if (t_simulation >= protocol(0)%tstart(i) .and. .not.protocol(0)%started(i)) then
		radiation_dose = protocol(0)%dose(i)
		protocol(0)%started(i) = .true.
		protocol(0)%ended(i) = .true.
		exit
	endif
enddo
do idrug = 1,2
	ichemo = idrug + GLUCOSE
	if (idrug == 1) then
		ichemo_metab = DRUG_A_METAB
	elseif (idrug == 2) then
		ichemo_metab = DRUG_B_METAB
	endif
	do i = 1,protocol(idrug)%n
		if (i == 1 .and. t_simulation < protocol(idrug)%tstart(i)) then
			chemo(ichemo)%bdry_conc = 0
			chemo(ichemo_metab)%bdry_conc = 0
			exit
		endif
		if (t_simulation >= protocol(idrug)%tstart(i) .and. .not.protocol(idrug)%started(i)) then
			chemo(ichemo)%bdry_conc = protocol(idrug)%conc(i)
			protocol(idrug)%started(i) = .true.
			protocol(idrug)%ended(i) = .false.
!			write(*,*) 'Started DRUG_A: ',chemo(ichemo)%bdry_conc
			exit
		endif
	enddo
	do i = 1,protocol(idrug)%n
		if (t_simulation >= protocol(idrug)%tend(i) .and. .not.protocol(idrug)%ended(i)) then
			chemo(ichemo)%bdry_conc = 0
			chemo(ichemo_metab)%bdry_conc = 0
			protocol(idrug)%ended(i) = .true.
!			write(*,*) 'Ended DRUG_A'
			exit
		endif
	enddo
enddo	
end subroutine

!-----------------------------------------------------------------------------------------
! Advance simulation through one big time step (DELTA_T)
! The concentration fields are first solved through NT_CONC subdivisions of DELTA_T,
! then the cell states are updated, which includes cell death and division.
! On either death or division cell positions are adjusted, and site concentrations
! (if necessary), and the ODE solver mappings in ODEdiff.
!-----------------------------------------------------------------------------------------
subroutine simulate_step(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: simulate_step  
use, intrinsic :: iso_c_binding
integer(c_int) :: res
integer :: kcell, site(3), hour, it, nthour, kpar=0
real(REAL_KIND) :: r(3), rmax, tstart, dt, radiation_dose
!integer, parameter :: NT_CONC = 6
integer :: nchemo

!call logger('simulate_step')
if (Ncells == 0) then
    res = 2
    return
endif
nthour = 3600/DELTA_T
dt = DELTA_T/NT_CONC
if (istep == 0) then
    call SetupODEdiff
endif
istep = istep + 1
if (mod(istep,nthour) == 0) then
	write(logmsg,*) 'istep, hour: ',istep,istep/nthour,nlist,ncells,nsites-ncells
	call logger(logmsg)
endif
t_simulation = (istep-1)*DELTA_T	! seconds
if (use_treatment) then
	call treatment(radiation_dose)
	if (radiation_dose > 0) then
		write(logmsg,'(a,f6.1)') 'Radiation dose: ',radiation_dose
		call logger(logmsg)
	endif
endif
call grow_cells(radiation_dose,DELTA_T)
call SetupODEdiff
call SiteCellToState
!call logger('solving')
do it = 1,NT_CONC
!	write(logmsg,*) 'it: ',it
!	call logger(logmsg)
	tstart = (it-1)*dt
	t_simulation = (istep-1)*DELTA_T + tstart
	call Solver(it,tstart,dt,Ncells)
enddo
!call logger('solved')
call StateToSiteCell
res = 0
if (mod(istep,60) == -1) then
	rmax = 0
	do kcell = 1,nlist
		r = cell_list(kcell)%site - Centre
		rmax = max(rmax,norm(r))
	enddo
	hour = istep/60
	write(logmsg,'(3i6,2f6.1)') istep, hour, Ncells, Radius, rmax
	call logger(logmsg)
!	call CheckBdryList
	call ShowConcs
!	call check_bdry
endif
!call test_cell_division
end subroutine

!--------------------------------------------------------------------------------  
! The GUI calls this subroutine to fetch the cell info needed to identify and render 
! the cells:
!   id			the cell's sequence number
!   position	(x,y,z)
!   state       this is translated into a colour
!
! The info is stored in integer arrays, one for B cells, one for FDCs,
! and one for cell-cell bonds (not used).
! As a quick-and-dirty measure, the first 7 B cells in the list are actually 
! markers to provide a visual indication of the extent of the follicular blob.
! Improving this:
! blobrange(:,:) holds the info about the ranges of x, y and z that the blob occupies.
! blobrange(1,1) <= x <= blobrange(1,2)
! blobrange(2,1) <= y <= blobrange(2,2)
! blobrange(3,1) <= z <= blobrange(3,2)
!--------------------------------------------------------------------------------
subroutine get_bcell_scene(nBC_list,BC_list,nFDCMRC_list,FDCMRC_list,nbond_list,bond_list) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_bcell_scene
use, intrinsic :: iso_c_binding
integer(c_int) :: nFDCMRC_list, nBC_list, nbond_list, FDCMRC_list(*), BC_list(*), bond_list(*)
integer :: k, kc, kcell, site(3), j, jb
integer :: col(3)
integer :: x, y, z
integer :: ifdcstate, ibcstate, ctype, stage, region
integer :: last_id1, last_id2
logical :: ok
integer, parameter :: axis_centre = -2	! identifies the ellipsoid centre
integer, parameter :: axis_end    = -3	! identifies the ellipsoid extent in 5 directions
integer, parameter :: axis_bottom = -4	! identifies the ellipsoid extent in the -Y direction, i.e. bottom surface
integer, parameter :: ninfo = 5			! the size of the info package for a cell (number of integers)
integer, parameter :: nax = 6			! number of points used to delineate the follicle

nBC_list = 0
nFDCMRC_list = 0
nbond_list = 0

! Need some markers to delineate the follicle extent.  These nax "cells" are used to convey (the follicle centre
! and) the approximate ellipsoidal blob limits in the 3 axis directions.
do k = 1,nax
	select case (k)
!	case (1)
!		x = Centre(1) + 0.5
!		y = Centre(2) + 0.5
!		z = Centre(3) + 0.5
!		site = (/x, y, z/)
!		ibcstate = axis_centre
	case (1)
!		x = Centre(1) - Radius%x - 2
		x = blobrange(1,1) - 1
		y = Centre(2) + 0.5
		z = Centre(3) + 0.5
		site = (/x, y, z/)
		ibcstate = axis_end
	case (2)
!		x = Centre(1) + Radius%x + 2
		x = blobrange(1,2) + 1
		y = Centre(2) + 0.5
		z = Centre(3) + 0.5
		site = (/x, y, z/)
		ibcstate = axis_end
	case (3)
		x = Centre(1) + 0.5
!		y = Centre(2) - Radius%y - 2
		y = blobrange(2,1) - 1
		z = Centre(3) + 0.5
		site = (/x, y, z/)
		ibcstate = axis_bottom
	case (4)
		x = Centre(1) + 0.5
!		y = Centre(2) + Radius%y + 2
		y = blobrange(2,2) + 1
		z = Centre(3) + 0.5
		site = (/x, y, z/)
		ibcstate = axis_end
	case (5)
		x = Centre(1) + 0.5
		y = Centre(2) + 0.5
!		z = Centre(3) - Radius%z - 2
		z = blobrange(3,1) - 1
		site = (/x, y, z/)
		ibcstate = axis_end
	case (6)
		x = Centre(1) + 0.5
		y = Centre(2) + 0.5
!		z = Centre(3) + Radius%z + 2
		z = blobrange(3,2) + 1
		site = (/x, y, z/)
		ibcstate = axis_end
	end select

	j = ninfo*(k-1)
	BC_list(j+1) = k-1
	BC_list(j+2:j+4) = site
	BC_list(j+5) = ibcstate
	last_id1 = k-1
enddo
k = last_id1 + 1

! Cells
do kcell = 1,nlist
	if (cell_list(kcell)%exists) then
		k = k+1
		j = ninfo*(k-1)
		site = cell_list(kcell)%site
		call cellColour(kcell,col)
		BC_list(j+1) = kcell + last_id1
		BC_list(j+2:j+4) = site
		BC_list(j+5) = rgb(col)
		last_id2 = kcell + last_id1
	endif
enddo
nBC_list = last_id2
end subroutine

!--------------------------------------------------------------------------------
! TC = tumour cell
!--------------------------------------------------------------------------------
subroutine get_scene(nTC_list,TC_list) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scene
use, intrinsic :: iso_c_binding
integer(c_int) :: nTC_list, TC_list(*)
integer :: k, kc, kcell, site(3), j, jb
integer :: col(3)
integer :: x, y, z
integer :: itcstate, ctype, stage, region
integer :: last_id1, last_id2
logical :: ok
integer, parameter :: axis_centre = -2	! identifies the spheroid centre
integer, parameter :: axis_end    = -3	! identifies the spheroid extent in 5 directions
integer, parameter :: axis_bottom = -4	! identifies the spheroid extent in the -Y direction, i.e. bottom surface
integer, parameter :: ninfo = 5			! the size of the info package for a cell (number of integers)
integer, parameter :: nax = 6			! number of points used to delineate the spheroid

nTC_list = 0

k = 0
last_id1 = 0
if (1 == 0) then
	! Need some markers to delineate the follicle extent.  These nax "cells" are used to convey (the follicle centre
	! and) the approximate ellipsoidal blob limits in the 3 axis directions.
	do k = 1,nax
		select case (k)
	!	case (1)
	!		x = Centre(1) + 0.5
	!		y = Centre(2) + 0.5
	!		z = Centre(3) + 0.5
	!		site = (/x, y, z/)
	!		itcstate = axis_centre
		case (1)
	!		x = Centre(1) - Radius%x - 2
			x = blobrange(1,1) - 1
			y = Centre(2) + 0.5
			z = Centre(3) + 0.5
			site = (/x, y, z/)
			itcstate = axis_end
		case (2)
	!		x = Centre(1) + Radius%x + 2
			x = blobrange(1,2) + 1
			y = Centre(2) + 0.5
			z = Centre(3) + 0.5
			site = (/x, y, z/)
			itcstate = axis_end
		case (3)
			x = Centre(1) + 0.5
	!		y = Centre(2) - Radius%y - 2
			y = blobrange(2,1) - 1
			z = Centre(3) + 0.5
			site = (/x, y, z/)
			itcstate = axis_bottom
		case (4)
			x = Centre(1) + 0.5
	!		y = Centre(2) + Radius%y + 2
			y = blobrange(2,2) + 1
			z = Centre(3) + 0.5
			site = (/x, y, z/)
			itcstate = axis_end
		case (5)
			x = Centre(1) + 0.5
			y = Centre(2) + 0.5
	!		z = Centre(3) - Radius%z - 2
			z = blobrange(3,1) - 1
			site = (/x, y, z/)
			itcstate = axis_end
		case (6)
			x = Centre(1) + 0.5
			y = Centre(2) + 0.5
	!		z = Centre(3) + Radius%z + 2
			z = blobrange(3,2) + 1
			site = (/x, y, z/)
			itcstate = axis_end
		end select

		j = ninfo*(k-1)
		TC_list(j+1) = k-1
		TC_list(j+2:j+4) = site
		TC_list(j+5) = itcstate
		last_id1 = k-1
	enddo
	k = last_id1 + 1
endif

! Cells
do kcell = 1,nlist
	if (idbug /= 0 .and. cell_list(kcell)%ID /= idbug) cycle
	if (cell_list(kcell)%exists) then
		k = k+1
		j = ninfo*(k-1)
		site = cell_list(kcell)%site
		call cellColour(kcell,col)
		TC_list(j+1) = kcell + last_id1
		TC_list(j+2:j+4) = site
		TC_list(j+5) = rgb(col)
		last_id2 = kcell + last_id1
	endif
enddo
!nTC_list = last_id2
nTC_list = k
end subroutine

!-----------------------------------------------------------------------------------------
! Rendered cognate B cell colour depends on stage, state, receptor expression level.
! col(:) = (r,g,b)
!-----------------------------------------------------------------------------------------
subroutine cellColour(kcell,col)
integer :: kcell, col(3)
integer :: stage, status
integer, parameter :: WHITE(3) = (/255,255,255/)
integer, parameter :: RED(3) = (/255,0,0/)
integer, parameter :: GREEN(3) = (/0,255,0/)
integer, parameter :: BLUE(3) = (/0,0,255/)
integer, parameter :: DEEPRED(3) = (/200,0,0/)
integer, parameter :: DEEPBLUE(3) = (/30,20,255/)
integer, parameter :: DEEPGREEN(3) = (/0,150,0/)
integer, parameter :: LIGHTRED(3) = (/255,70,90/)
integer, parameter :: LIGHTBLUE(3) = (/0,200,255/)
integer, parameter :: LIGHTGREEN(3) = (/50,255,150/)
integer, parameter :: DEEPORANGE(3) = (/240,70,0/)
integer, parameter :: LIGHTORANGE(3) = (/255,130,0/)
integer, parameter :: YELLOW(3) = (/255,255,0/)
integer, parameter :: DEEPPURPLE(3) = (/180,180,30/)
integer, parameter :: LIGHTPURPLE(3) = (/230,230,100/)
integer, parameter :: DEEPBROWN(3) = (/130,70,0/)
integer, parameter :: LIGHTBROWN(3) = (/200,100,0/)
integer, parameter :: GRAY(3) = (/128,128,128/)

integer, parameter :: Qt_white = 3
integer, parameter :: Qt_black = 2
integer, parameter :: Qt_red = 7
integer, parameter :: Qt_darkRed = 13
integer, parameter :: Qt_green = 8
integer, parameter :: Qt_darkGreen = 14
integer, parameter :: Qt_blue = 9
integer, parameter :: Qt_darkBlue = 15
integer, parameter :: Qt_cyan = 10
integer, parameter :: Qt_darkCyan = 16
integer, parameter :: Qt_magenta = 11
integer, parameter :: Qt_darkMagenta = 17
integer, parameter :: Qt_yellow = 12
integer, parameter :: Qt_darkYellow = 18
integer, parameter :: Qt_gray = 5
integer, parameter :: Qt_darkGray = 4
integer, parameter :: Qt_lightGray = 6

col = LIGHTORANGE
end subroutine

!-----------------------------------------------------------------------------------------
! Pack the colours (r,g,b) into an integer.
!-----------------------------------------------------------------------------------------
integer function rgb(col)
integer :: col(3)

rgb = ishft(col(1),16) + ishft(col(2),8) + col(3)
end function


!-----------------------------------------------------------------------------------------
! Live cells = Ncells
! Total deaths = Ndead
! Drug deaths = Ndrugdead
! Hypoxia deaths = Ndead - Ndrugdead
! Total tagged for drug death on division = Ndrug_tag
! Current tagged = Ntodie - Ntagdead 
!-----------------------------------------------------------------------------------------
subroutine get_summary(summaryData,icutoff) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*), icutoff
integer :: Ndead, Ntagged, Ntodie, Ntagdead, Ntagged_anoxia, diam_um, vol_mm3_1000, &
	nhypoxic(3), hypoxic_percent_10, necrotic_percent_10
real(REAL_KIND) :: vol_cm3, vol_mm3, hour

hour = istep*DELTA_T/3600.
vol_cm3 = Vsite*Nsites				! total volume in cm^3
vol_mm3 = vol_cm3*1000				! volume in mm^3
vol_mm3_1000 = vol_mm3*1000			! 1000 * volume in mm^3
diam_um = 2*DELTA_X*Radius*10000
Ntodie = Nradiation_tag + Ndrug_tag			! total that have been tagged by drug or radiation
Ntagdead = Nradiation_dead + Ndrug_dead		! total that died from drug or radiation
Ndead = Nsites + Nreuse - Ncells			! total that died from any cause
Ntagged = Ntodie - Ntagdead					! number currently tagged by drug or radiation
Ntagged_anoxia = Nanoxia_tag - Nanoxia_dead	! number currently tagged by anoxia
call get_hypoxic_count(nhypoxic)
hypoxic_percent_10 = (1000*nhypoxic(icutoff))/Ncells
necrotic_percent_10 = (1000*(Nsites-Ncells))/Nsites
summaryData(1:11) = (/ istep, Ncells, Nradiation_dead, Ndrug_dead, Ntagged, &
	diam_um, vol_mm3_1000, Nanoxia_dead, Ntagged_anoxia, hypoxic_percent_10, necrotic_percent_10 /)
write(nfres,'(i8,f8.2,f8.4,7i6,4f7.3)') istep, hour, vol_mm3, Ncells, Nradiation_dead, Ndrug_dead, Ntagged, &
	diam_um, Nanoxia_dead, Ntagged_anoxia, nhypoxic(:)/real(Ncells), (Nsites-Ncells)/real(Nsites)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_hypoxic_count(nhypoxic)
integer :: nhypoxic(3)
integer :: kcell, i

nhypoxic = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	do i = 1,3
		if (cell_list(kcell)%conc(OXYGEN) < O2cutoff(i)) nhypoxic(i) = nhypoxic(i) + 1
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
! Note: axis = 0,1,2
!--------------------------------------------------------------------------------
subroutine get_fieldinfo(nxx, axis, fraction, ns, nc, cused) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_fieldinfo
use, intrinsic :: iso_c_binding
integer(c_int) :: nxx, axis, ns, nc, cused(*)
real(c_double) :: fraction
integer rng(3,2), ichemo, kcell, x, y, z

nxx = NX
nc = MAX_CHEMO
do ichemo = 1,MAX_CHEMO
    if (chemo(ichemo)%used) then
        cused(ichemo) = 1
    else
        cused(ichemo) = 0
    endif
enddo
cused(MAX_CHEMO+1) = 1		! Growth rate
rng(:,1) = Centre(:) - (Radius + 2)
rng(:,2) = Centre(:) + (Radius + 2)
rng(axis,:) = Centre(axis) + fraction*Radius
write(logmsg,*) 'Centre, Radius, axis, fraction: ',Centre, Radius, axis, fraction
call logger(logmsg)
write(logmsg,*) 'rng: ',rng
call logger(logmsg)
ns = 0
do z = rng(3,1),rng(3,2)
    do y = rng(2,1),rng(2,2)
        do x = rng(1,1),rng(1,2)
            kcell = occupancy(x,y,z)%indx(1)
            if (kcell == OUTSIDE_TAG) cycle
            ns = ns + 1
        enddo
    enddo
enddo
write(logmsg,*) 'get_fieldinfo: ns: ',ns
call logger(logmsg)
end subroutine

!--------------------------------------------------------------------------------
! Need to transmit medium concentration data.  This could be a separate subroutine.
!--------------------------------------------------------------------------------
subroutine get_fielddata(axis, fraction, nfdata, nc, fdata) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_fielddata
use, intrinsic :: iso_c_binding
real(c_double) :: fraction
integer(c_int) :: axis, nc, nfdata
type(FIELD_DATA) :: fdata(*)
integer rng(3,2), kcell, x, y, z, i, ns
real(REAL_KIND) :: growthrate

write(logmsg,*) 'get_fielddata: nfdata, nc: ',nfdata, nc
call logger(logmsg)
if (nc > MAX_CHEMO) then
	write(logmsg,*) 'Error: get_fielddata: dimension of conc(MAX_CHEMO) not big enough!'
	call logger(logmsg)
	stop
endif
rng(:,1) = Centre(:) - (Radius + 2)
rng(:,2) = Centre(:) + (Radius + 2)
rng(axis,:) = Centre(axis) + fraction*Radius
ns = 0
do z = rng(3,1),rng(3,2)
    do y = rng(2,1),rng(2,2)
        do x = rng(1,1),rng(1,2)
            kcell = occupancy(x,y,z)%indx(1)
            if (kcell == OUTSIDE_TAG) cycle
            ns = ns + 1
	        i = ODEdiff%ivar(x,y,z)
            fdata(ns)%site = (/x,y,z/)
            fdata(ns)%state = 1
            if (kcell > 0) then
                fdata(ns)%volume = cell_list(kcell)%volume
                call get_growthrate(kcell,growthrate)
            else
                fdata(ns)%volume = 0
                growthrate = 0
            endif
            fdata(ns)%dVdt = growthrate
            fdata(ns)%conc(1:MAX_CHEMO) = allstate(i,1:MAX_CHEMO)
!            fdata(ns)%conc(MAX_CHEMO+1) = growthrate
        enddo
    enddo
enddo
write(logmsg,*) 'get_fielddata: ns: ',ns
call logger(logmsg)
if (ns /= nfdata) then
    write(logmsg,*) 'Error: inconsistent nsites: ',nfdata, ns
    call logger(logmsg)
    stop
endif

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_growthrate(kcell,growthrate)
integer :: kcell
real(REAL_KIND) :: growthrate

growthrate = cell_list(kcell)%dVdt
if (growthrate == 0) then
	write(nflog,'(a,2i6,e12.3)') 'growthrate: ',istep,kcell,growthrate
	stop
endif
end subroutine

!--------------------------------------------------------------------------------
! Returns all the extracellular concentrations along a line through the blob centre.
! Together with the growth rate dVdt
!--------------------------------------------------------------------------------
subroutine get_concdata(ns, dx, conc) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_concdata
use, intrinsic :: iso_c_binding
integer(c_int) :: ns
real(c_double) :: dx, conc(*)
real(REAL_KIND) :: cbnd, cmin = 1.0e-6
integer rng(3,2), i, k, ichemo, kcell, x, y, z

!call logger('get_concdata')
dx = DELTA_X
rng(:,1) = Centre(:) - (Radius + 2)
rng(:,2) = Centre(:) + (Radius + 2)
!rng(axis,:) = Centre(axis) + fraction*Radius
y = Centre(2) + 0.5
z = Centre(3) + 0.5
ns = 1
do x = rng(1,1),rng(1,2)
    kcell = occupancy(x,y,z)%indx(1)
    if (kcell == OUTSIDE_TAG) cycle
    ns = ns + 1
    do ichemo = 1,MAX_CHEMO+1
        i = ODEdiff%ivar(x,y,z)
        k = (ns-1)*(MAX_CHEMO+1) + ichemo
!        if (istep > 240) then
!	        write(logmsg,'(7i6)') x,kcell,ichemo,i,k
!		    call logger(logmsg)
!		endif
        if (ichemo <= MAX_CHEMO) then
			if (chemo(ichemo)%used) then
				if (i > 0) then
					conc(k) = allstate(i,ichemo)
				else
					conc(k) = 0
				endif
			else
				conc(k) = 0
			endif
        elseif (ichemo == MAX_CHEMO+1) then	! growth rate
			if (kcell > 0) then
				conc(k) = cell_list(kcell)%dVdt
			else
				conc(k) = 0
			endif
		endif
    enddo
enddo
! Add concentrations at the two boundaries
! At ns=1, at at ns=ns+1
ns = ns+1
do ichemo = 1,MAX_CHEMO
    if (chemo(ichemo)%used) then
        cbnd = BdryConc(ichemo,t_simulation)
        k = ichemo
        conc(k) = cbnd
        k = (ns-1)*(MAX_CHEMO+1) + ichemo
        conc(k) = cbnd
    else
        k = ichemo
        conc(k) = 0
        k = (ns-1)*(MAX_CHEMO+1) + ichemo
        conc(k) = 0
    endif
enddo
!do k = 1,(MAX_CHEMO+1)*ns
!    conc(k) = max(cmin,conc(k))
!enddo
!call logger('did get_concdata')
end subroutine

!--------------------------------------------------------------------------------
! Returns the distribution of cell volume.
! nv is passed from the GUI
! Min divide volume = Vdivide0 - dVdivide
! therefore Vmin = (Vdivide0 - dVdivide)/2
! Max divide volume = Vmax = Vdivide0 + dVdivide
! dv = (Vmax - Vmin)/nv
! v0 = Vmin + dv/2
!--------------------------------------------------------------------------------
subroutine get_volprob(nv, v0, dv, prob) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_volprob
use, intrinsic :: iso_c_binding
integer(c_int) :: nv
real(c_double) :: v0, dv, prob(*)
integer :: n, kcell, k
real(REAL_KIND) :: v, Vmin, Vmax

!call logger('get_volprob')
Vmin = (Vdivide0 - dVdivide)/2
Vmax = Vdivide0 + dVdivide
dv = (Vmax - Vmin)/nv
v0 = Vmin + dv/2
prob(1:nv) = 0
n = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	v = cell_list(kcell)%volume
	k = (v - Vmin)/dv + 1
	k = min(k,nv)
	prob(k) = prob(k) + 1
	n = n+1
enddo
prob(1:nv) = prob(1:nv)/n
end subroutine

!--------------------------------------------------------------------------------
! Returns the distribution of intracellular O2 level
! nv is passed from the GUI

!--------------------------------------------------------------------------------
subroutine get_oxyprob(nv, dv, prob) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_oxyprob
use, intrinsic :: iso_c_binding
integer(c_int) :: nv
real(c_double) :: v0, dv, prob(*)
integer :: n, kcell, k
real(REAL_KIND) :: v, Vmin, Vmax, O2max

!call logger('get_oxyprob')
Vmin = 0
Vmax = chemo(OXYGEN)%bdry_conc
dv = (Vmax - Vmin)/nv
prob(1:nv) = 0
n = 0
O2max = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	v = cell_list(kcell)%conc(OXYGEN)
	O2max = max(O2max,v)
	k = (v - Vmin)/dv + 1
	k = min(k,nv)
	prob(k) = prob(k) + 1
	n = n+1
enddo
prob(1:nv) = prob(1:nv)/n
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine execute(ncpu,infile_array,inbuflen,outfile_array,outbuflen) BIND(C)
!!DEC$ ATTRIBUTES DLLEXPORT :: EXECUTE
!!DEC$ ATTRIBUTES C, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"EXECUTE" :: execute
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char) :: infile_array(128), outfile_array(128)
integer(c_int) :: ncpu, inbuflen, outbuflen
character*(128) :: infile, outfile
logical :: ok, success
integer :: i, res

!use_CPORT1 = .false.	! DIRECT CALLING FROM C++
infile = ''
do i = 1,inbuflen
	infile(i:i) = infile_array(i)
enddo
outfile = ''
do i = 1,outbuflen
	outfile(i:i) = outfile_array(i)
enddo

open(nflog,file='spheroid.log',status='replace')
open(nfres,file='spheroid_ts.out',status='replace')
write(nfres,'(a)') 'istep hour vol_mm3 Ncells Nradiation_dead Ndrug_dead Ntagged diam_um Nanoxia_dead f_hypox f_necrot'
!awp_0%is_open = .false.
!awp_1%is_open = .false.

#ifdef GFORTRAN
    write(logmsg,'(a)') 'Built with GFORTRAN'
	call logger(logmsg)
#endif

logmsg = 'OS??'
#ifdef LINUX
    write(logmsg,'(a)') 'OS is Linux'
#endif
#ifdef OSX
    write(logmsg,'(a)') 'OS is OS-X'
#endif
#ifdef _WIN32
    write(logmsg,'(a)') 'OS is Windows'
#endif
#ifdef WINDOWS
    write(logmsg,'(a)') 'OS is Windows'
#endif
call logger(logmsg)

!#ifdef OPENMP
#if defined(OPENMP) || defined(_OPENMP)
    write(logmsg,'(a)') 'Executing with OpenMP'
	call logger(logmsg)
#endif

write(logmsg,*) 'inputfile:  ', infile
call logger(logmsg)
write(logmsg,*) 'outputfile: ', outfile 
call logger(logmsg)
if (use_TCP) then
!	call logger('call connector')
	call connecter(ok)
!	call logger('did connector')
	if (.not.ok) then
		call logger('Failed to make TCP connections')
		return
	endif
endif

istep = 0
NX=100
NY=100
NZ=100
DELTA_T = 600
nsteps = 100
res=0

call setup(ncpu,infile,outfile,ok)
if (ok) then
!	clear_to_send = .true.
!	simulation_start = .true.
else
	call logger('=== Setup failed ===')
endif
if (ok) then
	res = 0
else
	res = 1
endif

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine disableTCP
!DEC$ ATTRIBUTES DLLEXPORT :: disableTCP
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"DISABLETCP" :: disableTCP

use_TCP = .false.   ! because this is called from spheroid_main()	
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine connection(awp,port,error)
TYPE(winsockport) :: awp
integer :: port, error
integer :: address = 0
!!!character*(64) :: ip_address = "127.0.0.1"C      ! need a portable way to make a null-terminated C string
character*(64) :: host_name = "localhost"

if (.not.winsock_init(1)) then
    call logger("winsock_init failed")
    stop
endif

awp%handle = 0
awp%host_name = host_name
awp%ip_port = port
awp%protocol = IPPROTO_TCP
call Set_Winsock_Port (awp,error)

if (.not.awp%is_open) then
    write(nflog,*) 'Error: connection: awp not open: ',port
else
    write(nflog,*) 'connection: awp open: ',port, error
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine connecter(ok)
logical :: ok
integer :: error

! Main connection
ok = .true.
error = 0
call connection(awp_0,TCP_PORT_0,error)
if (awp_0%handle < 0 .or. error /= 0) then
    write(logmsg,'(a)') 'TCP connection to TCP_PORT_0 failed'
    call logger(logmsg)
    ok = .false.
    return
endif
if (.not.awp_0%is_open) then
	write(logmsg,'(a)') 'No connection to TCP_PORT_0'
    call logger(logmsg)
    ok = .false.
    return
endif
write(logmsg,'(a)') 'Connected to TCP_PORT_0  '
call logger(logmsg)

if (use_CPORT1) then
	call connection(awp_1,TCP_PORT_1,error)
	if (awp_1%handle < 0 .or. error /= 0) then
		write(logmsg,'(a)') 'TCP connection to TCP_PORT_1 failed'
		call logger(logmsg)
		ok = .false.
		return
	endif
	if (.not.awp_1%is_open) then
		write(logmsg,'(a)') 'No connection to TCP_PORT_1'
		call logger(logmsg)
		ok = .false.
		return
	endif
	write(logmsg,'(a)') 'Connected to TCP_PORT_1  '
	call logger(logmsg)
endif
! Allow time for completion of the connection
call sleeper(2)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine terminate_run(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: terminate_run 
use, intrinsic :: iso_c_binding
integer(c_int) :: res
character*(8), parameter :: quit = '__EXIT__'
integer :: error, i

call logger('terminate_run')
call wrapup

if (res == 0) then
	call logger(' Execution successful!')
elseif (res == 2) then
	call logger(' Execution successful, no live cells')
else
	call logger('  === Execution failed ===')
!	call sleeper(1)
endif
close(nflog)

if (use_TCP) then
	if (stopped) then
	    call winsock_close(awp_0)
	    if (use_CPORT1) call winsock_close(awp_1)
	else
	    call winsock_send(awp_0,quit,8,error)
	    call winsock_close(awp_0)
		if (use_CPORT1) then
			call winsock_send(awp_1,quit,8,error)
			call winsock_close(awp_1)
		endif
	endif
endif
write(logmsg,*) 'awp_0%is_open: ',awp_0%is_open
call logger(logmsg)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine wrapup
integer :: ierr, ichemo, idrug
logical :: isopen

call logger('doing wrapup ...')
ierr = 0
if (allocated(zoffset)) deallocate(zoffset)
if (allocated(zdomain)) deallocate(zdomain)
!if (allocated(occupancy)) deallocate(occupancy)
!if (allocated(cell_list)) deallocate(cell_list)
!if (allocated(allstate)) deallocate(allstate)
if (allocated(allstatep)) deallocate(allstatep)
if (allocated(work_rkc)) deallocate(work_rkc)
do ichemo = 1,MAX_CHEMO
	if (allocated(chemo(ichemo)%coef)) deallocate(chemo(ichemo)%coef)
	if (allocated(chemo(ichemo)%conc)) deallocate(chemo(ichemo)%conc)
	if (allocated(chemo(ichemo)%grad)) deallocate(chemo(ichemo)%grad)
enddo
!if (allocated(ODEdiff%ivar)) deallocate(ODEdiff%ivar)
if (allocated(ODEdiff%varsite)) deallocate(ODEdiff%varsite)
if (allocated(ODEdiff%icoef)) deallocate(ODEdiff%icoef)
if (allocated(protocol)) then
	do idrug = 0,2	!<------  change this to a variable
		if (allocated(protocol(idrug)%tstart)) deallocate(protocol(idrug)%tstart)
		if (allocated(protocol(idrug)%tend)) deallocate(protocol(idrug)%tend)
		if (allocated(protocol(idrug)%conc)) deallocate(protocol(idrug)%conc)
		if (allocated(protocol(idrug)%dose)) deallocate(protocol(idrug)%dose)
		if (allocated(protocol(idrug)%started)) deallocate(protocol(idrug)%started)
		if (allocated(protocol(idrug)%ended)) deallocate(protocol(idrug)%ended)
	enddo
	deallocate(protocol)
endif
call logger('deallocated all arrays')

! Close all open files
inquire(unit=nfout,OPENED=isopen)
if (isopen) then
	close(nfout)
	call logger('closed nfout')
endif
inquire(nfres,OPENED=isopen)
if (isopen) close(nfres)
call logger('closed files')

if (par_zig_init) then
	call par_zigfree
endif
call logger('freed par_zig')
end subroutine

end module
