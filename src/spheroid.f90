module spheroid_mod
use global
use boundary
use chemokine
use cellstate
use winsock  
use deform
use drop

#include "../src/version.h"

IMPLICIT NONE

contains 

!-----------------------------------------------------------------------------------------
! This subroutine is called to initialize a simulation run. 
! ncpu = the number of processors to use 
! infile = file with the input data
! outfile = file to hold the output
!-----------------------------------------------------------------------------------------
subroutine Setup(ncpu,infile,outfile,ok)
integer :: ncpu
character*(*) :: infile, outfile
logical :: ok
character*(64) :: msg
integer :: ichemo, error

ok = .true.
initialized = .false.
par_zig_init = .false.

inputfile = infile
outputfile = outfile
call logger("ReadCellParams")
call ReadCellParams(ok)
if (.not.ok) return
call logger("did ReadCellParams")

if (ncpu == 0) then
	ncpu = ncpu_input
endif
Mnodes = ncpu
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

call ArrayInitialisation(ok)
if (.not.ok) return
call logger('did ArrayInitialisation')

call SetupChemo

is_dropped = .false.
adrop = 1
bdrop = 1
cdrop = 0
zmin = 1
call PlaceCells(ok)
call SetRadius(Nsites)
write(logmsg,*) 'did PlaceCells: Ncells: ',Ncells,Radius
call logger(logmsg)
if (.not.ok) return

!call test_get_path

call CreateBdryList

!chemo_N = 8
!call ChemoSetup
!call MakeSplit(.true.)
!call init_counters
!if (save_input) then
!    call save_inputfile(inputfile)
!    call save_parameters
!	call save_inputfile(treatmentfile)
!endif
!call AllocateConcArrays
!call ChemoSteadystate
!firstSummary = .true.
!initialized = .true.
!call checkcellcount(ok)

!if (use_ODE_diffusion) then
!	call SetupODEDiff
!	call TestODEDiffusion
!	call TestSolver
!endif
istep = 0
call SetupODEdiff
allocate(allstate(MAX_VARS,MAX_CHEMO))
allocate(work_rkc(8+5*MAX_VARS))
do ichemo = 1,TRACER
	if (chemo(ichemo)%used) then
		call InitConcs(ichemo)
		call SetupMedium(ichemo)
	endif
enddo
call AdjustMM
call SetInitialGrowthRate
Nradiation_tag = 0
NdrugA_tag = 0
NdrugB_tag = 0
Nanoxia_tag = 0
Nradiation_dead = 0
NdrugA_dead = 0
NdrugB_dead = 0
Nanoxia_dead = 0
t_simulation = 0
it_saveprofiledata = 1
write(logmsg,'(a,i6)') 'Startup procedures have been executed: initial T cell count: ',Ncells0
call logger(logmsg)

end subroutine

!----------------------------------------------------------------------------------------- 
! Initialise medium concentrations, etc.
!-----------------------------------------------------------------------------------------
subroutine SetupMedium(ichemo)
integer :: ichemo
real(REAL_KIND) :: V, V0, R1

if (chemo(ichemo)%present) then
	call SetRadius(Nsites)
	R1 = Radius*DELTA_X			! cm
	V0 = total_volume			! cm3
	V = V0 - (4./3.)*PI*R1**3	! cm3
	chemo(ichemo)%medium_Cext = chemo(ichemo)%bdry_conc
	chemo(ichemo)%medium_Cbnd = chemo(ichemo)%bdry_conc
	chemo(ichemo)%medium_M = V*chemo(ichemo)%bdry_conc
	chemo(ichemo)%medium_U = 0
else
	chemo(ichemo)%medium_Cext = 0
	chemo(ichemo)%medium_Cbnd = 0
	chemo(ichemo)%medium_M = 0
	chemo(ichemo)%medium_U = 0
endif
end subroutine

!----------------------------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------
subroutine omp_initialisation(ok)
logical :: ok
integer :: npr, nth

ok = .true.
!if (Mnodes == 1) return
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
!call test_omp1

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_omp
integer, parameter :: n = 10
integer :: i

integer :: sum1, sum2
integer, allocatable :: y1(:)
integer :: y2(n)

allocate(y1(n))
y1 = 1
y2 = 1

sum1 = 0
sum2 = 0
!$omp parallel do
do i = 1,n
	sum1 = sum1 + y1(i)
	sum2 = sum2 + y2(i)
enddo
!$omp end parallel do
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_omp1
integer :: n = 1000000
!integer :: ncpu = 2
integer :: i, k, ith

real(REAL_KIND) :: x, y, z

real(REAL_KIND) :: sum1, sum2
real(REAL_KIND), allocatable :: y1(:), ysum(:)
!integer :: y2(100)

allocate(omp_x(n))
allocate(omp_y(n))
allocate(omp_z(n))
omp_x = 1.23456
omp_y = 9.87654

!call omp_set_num_threads(ncpu)

allocate(y1(n))
allocate(ysum(0:n-1))
y1 = 1
!y2 = 1
ysum = 0

sum1 = 0
sum2 = 0
!$omp parallel do private(x, y, z, k)
do i = 1,n
	z = 0
	omp_x(i) = omp_x(i)/i
	omp_y(i) = omp_y(i)/i + 1.0
	if (i < 4) then
		x = omp_x(i)**i
	elseif (i < 8) then
		x = 1.0/omp_x(i)**i
	else 
		x = 1./(1 + omp_x(i))
	endif
	omp_x(i) = x
	y = sqrt(i*omp_y(i))
	y = omp_y(i)
	do k = 1,1000
		y = 0.5*y
		call omp_sub(i,x,y)
		z = z + 1/omp_z(i)
	enddo
!	sum1 = sum1 + z
	ith = omp_get_thread_num()
	ysum(ith) = ysum(ith) + z
enddo
!$omp end parallel do
stop
end subroutine

subroutine omp_sub(i,x,y)
integer :: i
real(REAL_KIND) :: x, y
omp_z(i) = x/y + omp_x(i)/omp_y(i)
omp_y(i) = y
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ArrayInitialisation(ok)
logical :: ok
integer :: x,y,z,k, ichemo
integer :: MAXX, z1, z2, nc0, inflow
integer :: cog_size
real(REAL_KIND) :: d, rr(3)

ok = .false.
call RngInitialisation

! These are deallocated here instead of in subroutine wrapup so that when a simulation run ends 
! it will still be possible to view the cell distributions and chemokine concentration fields.
if (allocated(occupancy)) deallocate(occupancy)
if (allocated(cell_list)) deallocate(cell_list)
if (allocated(allstate)) deallocate(allstate)
if (allocated(ODEdiff%ivar)) deallocate(ODEdiff%ivar)
call logger('did deallocation')

!nsteps_per_min = 1.0/DELTA_T
NY = NX
NZ = NX
ngaps = 0
!max_ngaps = NY*NZ
nlist = 0

allocate(zoffset(0:2*Mnodes))
allocate(zdomain(NZ))
x0 = (NX + 1.0)/2.        ! global value
y0 = (NY + 1.0)/2.
z0 = (NZ + 1.0)/2.
Centre = (/x0,y0,z0/)   ! now, actually the global centre (units = grids)
call SetRadius(initial_count)
write(logmsg,*) 'Initial radius, nc0, max_nlist: ',Radius, initial_count, max_nlist
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
subroutine RngInitialisation
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
subroutine ReadCellParams(ok)
logical :: ok
integer :: i, idrug, imetab, nmetab, im, itestcase, Nmm3, ichemo, itreatment, iuse_extra, iuse_relax, iuse_par_relax
integer :: iuse_oxygen, iuse_glucose, iuse_tracer, iuse_drug, iuse_metab, iV_depend, iV_random
!integer ::  idrug_decay, imetab_decay
integer :: ictype, idisplay, isconstant
integer :: iuse_drop, iconstant, isaveprofiledata
logical :: use_metabolites
real(REAL_KIND) :: days, bdry_conc, percent
real(REAL_KIND) :: sigma, DXmm, anoxia_tag_hours, anoxia_death_hours
character*(12) :: drug_name
character*(1) :: numstr

ok = .true.
chemo(:)%used = .false.

open(nfcell,file=inputfile,status='old')

read(nfcell,*) gui_run_version				! program run version number
read(nfcell,*) dll_run_version				! DLL run version number
read(nfcell,*) NX							! size of grid
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
read(nfcell,*) medium_volume0				! initial total volume (medium + spheroid) (cm^3)
read(nfcell,*) d_layer						! thickness of the unstirred layer around the spheroid (cm)
read(nfcell,*) Vdivide0						! nominal cell volume multiple for division
read(nfcell,*) dVdivide						! variation about nominal divide volume
read(nfcell,*) MM_THRESHOLD					! O2 concentration threshold Michaelis-Menten "soft-landing" (uM)
read(nfcell,*) ANOXIA_THRESHOLD			    ! O2 threshold for anoxia (uM)
read(nfcell,*) anoxia_tag_hours				! hypoxic time leading to tagging to die by anoxia (h)
read(nfcell,*) anoxia_death_hours			! time after tagging to death by anoxia (h)
read(nfcell,*) itestcase                    ! test case to simulate
read(nfcell,*) seed(1)						! seed vector(1) for the RNGs
read(nfcell,*) seed(2)						! seed vector(2) for the RNGs
read(nfcell,*) ncpu_input					! for GUI just a placeholder for ncpu, used only when execute parameter ncpu = 0
read(nfcell,*) Ncelltypes					! maximum number of cell types in the spheroid
do ictype = 1,Ncelltypes
	read(nfcell,*) percent
	celltype_fraction(ictype) = percent/100
!	read(nfcell,*) idisplay
!	celltype_display(ictype) = (idisplay == 1)
enddo
read(nfcell,*) NT_GUI_OUT					! interval between GUI outputs (timesteps)
read(nfcell,*) show_progeny                 ! if != 0, the number of the cell to show descendents of
read(nfcell,*) iuse_oxygen		! chemo(OXYGEN)%used
read(nfcell,*) chemo(OXYGEN)%diff_coef
read(nfcell,*) chemo(OXYGEN)%medium_diff_coef
read(nfcell,*) chemo(OXYGEN)%membrane_diff_in
chemo(OXYGEN)%membrane_diff_in = chemo(OXYGEN)%membrane_diff_in/60		! /min -> /sec
chemo(OXYGEN)%membrane_diff_out = chemo(OXYGEN)%membrane_diff_in
read(nfcell,*) chemo(OXYGEN)%bdry_conc
read(nfcell,*) iconstant
chemo(OXYGEN)%constant = (iconstant == 1)
read(nfcell,*) chemo(OXYGEN)%max_cell_rate
read(nfcell,*) chemo(OXYGEN)%MM_C0
read(nfcell,*) chemo(OXYGEN)%Hill_N
read(nfcell,*) iuse_glucose		!chemo(GLUCOSE)%used
read(nfcell,*) chemo(GLUCOSE)%diff_coef
read(nfcell,*) chemo(GLUCOSE)%medium_diff_coef
read(nfcell,*) chemo(GLUCOSE)%membrane_diff_in
chemo(GLUCOSE)%membrane_diff_in = chemo(GLUCOSE)%membrane_diff_in/60	! /min -> /sec
chemo(GLUCOSE)%membrane_diff_out = chemo(GLUCOSE)%membrane_diff_in
read(nfcell,*) chemo(GLUCOSE)%bdry_conc
read(nfcell,*) iconstant
chemo(GLUCOSE)%constant = (iconstant == 1)
read(nfcell,*) chemo(GLUCOSE)%max_cell_rate
read(nfcell,*) chemo(GLUCOSE)%MM_C0
read(nfcell,*) chemo(GLUCOSE)%Hill_N
read(nfcell,*) iuse_tracer		!chemo(TRACER)%used
read(nfcell,*) chemo(TRACER)%diff_coef
read(nfcell,*) chemo(TRACER)%medium_diff_coef
read(nfcell,*) chemo(TRACER)%membrane_diff_in
chemo(TRACER)%membrane_diff_in = chemo(TRACER)%membrane_diff_in/60		! /min -> /sec
chemo(TRACER)%membrane_diff_out = chemo(TRACER)%membrane_diff_in
read(nfcell,*) chemo(TRACER)%bdry_conc
read(nfcell,*) iconstant
chemo(TRACER)%constant = (iconstant == 1)
read(nfcell,*) chemo(TRACER)%max_cell_rate
read(nfcell,*) chemo(TRACER)%MM_C0
read(nfcell,*) chemo(TRACER)%Hill_N
do i = 1,2			! currently allowing for just two different drugs: 1 = TPZ-type, 2 = DNB-type
	read(nfcell,*) iuse_drug
	read(nfcell,'(a12)') drug_name
	read(nfcell,*) bdry_conc
	read(nfcell,*) iconstant
	read(nfcell,*) iuse_metab
	call getIndices(i, idrug, nmetab)
	if (idrug < 0 .and. iuse_drug /= 0) then
		write(logmsg,*) 'Unrecognized drug name: ',drug_name
		call logger(logmsg)
		ok = .false.
		return
	endif
	if (idrug == 0) cycle
	chemo(idrug)%name = drug_name
	chemo(idrug)%used = (iuse_drug == 1)
	chemo(idrug)%bdry_conc = bdry_conc
	chemo(idrug)%constant = (iconstant == 1)
!	chemo(idrug)%decay = (idrug_decay == 1)
	if (chemo(idrug)%used) then
		use_metabolites = (iuse_metab == 1)
		write(nflog,*) 'drug: ',idrug,'  name: ',chemo(idrug)%name,' use_metabolites: ',use_metabolites
	else
		use_metabolites = .false.
	endif
!	chemo(imetab)%decay = (imetab_decay == 1)
	if (idrug == TPZ_DRUG) then
		TPZ%name = drug_name
		TPZ%nmetabolites = nmetab
		TPZ%use_metabolites = use_metabolites
		do im = 0,nmetab
			read(nfcell,*) TPZ%diff_coef(im)
			read(nfcell,*) TPZ%medium_diff_coef(im)
			read(nfcell,*) TPZ%membrane_diff_in(im)
			TPZ%membrane_diff_in(im) = TPZ%membrane_diff_in(im)/60		! /min -> /sec
			read(nfcell,*) TPZ%membrane_diff_out(im)
			TPZ%membrane_diff_out(im) = TPZ%membrane_diff_out(im)/60	! /min -> /sec
			read(nfcell,*) TPZ%halflife(im)
			ichemo = idrug + im
			if (im > 0) then
				write(numstr,'(i1)') im
				chemo(ichemo)%name = 'TPZ_metab_'//numstr
			endif
			chemo(ichemo)%diff_coef = TPZ%diff_coef(im)
			chemo(ichemo)%medium_diff_coef = TPZ%medium_diff_coef(im)
			chemo(ichemo)%membrane_diff_in = TPZ%membrane_diff_in(im)
			chemo(ichemo)%membrane_diff_out = TPZ%membrane_diff_out(im)
			chemo(ichemo)%halflife = TPZ%halflife(im)
		enddo
		do ictype = 1,Ncelltypes
			do im = 0,nmetab
				read(nfcell,*) TPZ%Kmet0(ictype,im)
				read(nfcell,*) TPZ%C2(ictype,im)
				read(nfcell,*) TPZ%KO2(ictype,im)
				read(nfcell,*) TPZ%Vmax(ictype,im)
				read(nfcell,*) TPZ%Km(ictype,im)
				read(nfcell,*) TPZ%Klesion(ictype,im)
				if (im == 0) then
					read(nfcell,*) TPZ%kill_model(ictype)
					read(nfcell,*) TPZ%kill_O2(ictype)
					read(nfcell,*) TPZ%kill_drug(ictype)
					read(nfcell,*) TPZ%kill_duration(ictype)
					read(nfcell,*) TPZ%kill_fraction(ictype)
				endif
			enddo
			TPZ%Kmet0(ictype,:) = TPZ%Kmet0(ictype,:)/60				! /min -> /sec
			TPZ%KO2(ictype,:) = 1.0e-3*TPZ%KO2(ictype,:)                ! um -> mM
			TPZ%kill_duration(ictype) = 60*TPZ%kill_duration(ictype)    ! minutes -> seconds
		enddo

	elseif (idrug == DNB_DRUG) then
		DNB%name = drug_name
		DNB%nmetabolites = nmetab
		DNB%use_metabolites = use_metabolites
		do im = 0,nmetab
			read(nfcell,*) DNB%diff_coef(im)
			read(nfcell,*) DNB%medium_diff_coef(im)
			read(nfcell,*) DNB%membrane_diff_in(im)
			DNB%membrane_diff_in(im) = DNB%membrane_diff_in(im)/60		! /min -> /sec
			read(nfcell,*) DNB%membrane_diff_out(im)
			DNB%membrane_diff_out(im) = DNB%membrane_diff_out(im)/60	! /min -> /sec
			read(nfcell,*) DNB%halflife(im)
			ichemo = idrug + im
			if (im > 0) then
				write(numstr,'(i1)') im
				chemo(ichemo)%name = 'DNB_metab_'//numstr
			endif
			chemo(ichemo)%diff_coef = DNB%diff_coef(im)
			chemo(ichemo)%medium_diff_coef = DNB%medium_diff_coef(im)
			chemo(ichemo)%membrane_diff_in = DNB%membrane_diff_in(im)
			chemo(ichemo)%membrane_diff_out = DNB%membrane_diff_out(im)
			chemo(ichemo)%halflife = DNB%halflife(im)	
		enddo
		do ictype = 1,Ncelltypes
			do im = 0,nmetab
				read(nfcell,*) DNB%Kmet0(ictype,im)
				read(nfcell,*) DNB%C2(ictype,im)
				read(nfcell,*) DNB%KO2(ictype,im)
				read(nfcell,*) DNB%Vmax(ictype,im)
				read(nfcell,*) DNB%Km(ictype,im)
				read(nfcell,*) DNB%Klesion(ictype,im)
				read(nfcell,*) DNB%kill_model(ictype,im)
				read(nfcell,*) DNB%kill_O2(ictype,im)
				read(nfcell,*) DNB%kill_drug(ictype,im)
				read(nfcell,*) DNB%kill_duration(ictype,im)
				read(nfcell,*) DNB%kill_fraction(ictype,im)
			enddo
			DNB%Kmet0(ictype,:) = DNB%Kmet0(ictype,:)/60					! /min -> /sec
			DNB%KO2(ictype,:) = 1.0e-3*DNB%KO2(ictype,:)					! um -> mM
			DNB%kill_duration(ictype,:) = 60*DNB%kill_duration(ictype,:)    ! minutes -> seconds
		enddo
	endif
	do im = 0,nmetab
		ichemo = idrug + im
		chemo(ichemo)%decay = (chemo(ichemo)%halflife > 0)
		if (chemo(ichemo)%used .and. chemo(ichemo)%decay) then
			chemo(ichemo)%decay_rate = DecayRate(chemo(ichemo)%halflife)
		else
			chemo(ichemo)%decay_rate = 0
		endif
!		if (chemo(imetab)%used .and. chemo(imetab)%decay) then
!			chemo(imetab)%decay_rate = DecayRate(chemo(imetab)%halflife)
!		else
!			chemo(imetab)%decay_rate = 0
!		endif
	enddo
enddo

read(nfcell,*) LQ%alpha_H
read(nfcell,*) LQ%beta_H
read(nfcell,*) LQ%OER_am
read(nfcell,*) LQ%OER_bm
read(nfcell,*) LQ%K_ms
read(nfcell,*) O2cutoff(1)
read(nfcell,*) O2cutoff(2)
read(nfcell,*) O2cutoff(3)
read(nfcell,*) growthcutoff(1)
read(nfcell,*) growthcutoff(2)
read(nfcell,*) growthcutoff(3)
read(nfcell,*) spcrad_value
read(nfcell,*) iuse_extra
read(nfcell,*) iuse_relax
read(nfcell,*) iuse_par_relax
!read(nfcell,*) iuse_radiation
!read(nfcell,*) itreatment
!read(nfcell,*) treatmentfile						! file with treatment programme
read(nfcell,*) iuse_drop
read(nfcell,*) Ndrop
read(nfcell,*) alpha_shape
read(nfcell,*) beta_shape
read(nfcell,*) isaveprofiledata
read(nfcell,*) profiledatafilebase
read(nfcell,*) dt_saveprofiledata
read(nfcell,*) nt_saveprofiledata

if (use_events) then
	call ReadProtocol(nfcell)
	use_treatment = .false.
endif

close(nfcell)

if (chemo(OXYGEN)%Hill_N /= 1 .and. chemo(OXYGEN)%Hill_N /= 2) then
	call logger('Error: OXYGEN_HILL_N must be 1 or 2')
	ok = .false.
	return
endif
if (chemo(GLUCOSE)%Hill_N /= 1 .and. chemo(GLUCOSE)%Hill_N /= 2) then
	call logger('Error: GLUCOSE_HILL_N must be 1 or 2')
	ok = .false.
	return
endif
MM_THRESHOLD = MM_THRESHOLD/1000					! uM -> mM
ANOXIA_THRESHOLD = ANOXIA_THRESHOLD/1000			! uM -> mM
O2cutoff = O2cutoff/1000							! uM -> mM
relax = (iuse_relax == 1)
use_parallel = (iuse_par_relax == 1)
chemo(OXYGEN)%used = (iuse_oxygen == 1)
chemo(GLUCOSE)%used = (iuse_glucose == 1)
chemo(TRACER)%used = (iuse_tracer == 1)
chemo(OXYGEN)%MM_C0 = chemo(OXYGEN)%MM_C0/1000		! uM -> mM
chemo(GLUCOSE)%MM_C0 = chemo(GLUCOSE)%MM_C0/1000	! uM -> mM
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
Vsite_cm3 = DELTA_X*DELTA_X*DELTA_X					! total site volume (cm^3)
Vextra_cm3 = fluid_fraction*Vsite_cm3				! extracellular volume in a site (cm^3)
cell_radius = (3*(1-fluid_fraction)*Vsite_cm3/(4*PI))**(1./3.)
! In a well-oxygenated tumour the average cell fractional volume is intermediate between vdivide0/2 and vdivide0.
! We assume that 0.75*vdivide0*Vcell_cm3 = (1 - fluid_fraction)*Vsite_cm3
Vcell_cm3 = (1 - fluid_fraction)*Vsite_cm3/(0.75*vdivide0)	! nominal cell volume in cm^3 (corresponds to %volume = 1)
															! fractional volume: cell%volume = (cell volume in cm^3)/Vcell_cm3, 
															! actual volume (cm^3) = Vcell_cm3*cell%volume
Vcell_pL = 1.0e9*Vcell_Cm3									! nominal cell volume in pL
															! actual volume (pL) = Vcell_pL*cell%volume
total_volume = medium_volume0

write(logmsg,'(a,3e12.4)') 'DELTA_X, cell_radius: ',DELTA_X,cell_radius
call logger(logmsg)
write(logmsg,'(a,4e12.4)') 'Volumes: site, extra, cell (average, base): ',Vsite_cm3, Vextra_cm3, Vsite_cm3-Vextra_cm3, Vcell_cm3
call logger(logmsg)

saveprofiledata = (isaveprofiledata == 1)
dt_saveprofiledata = 60*dt_saveprofiledata			! mins -> seconds
!if (.not.use_events .and. itreatment == 1) then
!	use_treatment = .true.
!	call ReadTreatment(ok)
!	if (.not.ok) then
!		write(logmsg,'(a,a)') 'Error reading treatment programme file: ',treatmentfile
!		call logger(logmsg)
!		return
!	endif
!else	
!	use_treatment = .false.
!endif
use_dropper = (iuse_drop == 1)

! Setup test_case
test_case = .false.
if (itestcase /= 0) then
    test_case(itestcase) = .true.
endif

if (mod(NX,2) /= 0) NX = NX+1					! ensure that NX is even

open(nfout,file=outputfile,status='replace')
write(nfout,'(a,a)') 'GUI version: ',gui_run_version
write(nfout,'(a,a)') 'DLL version: ',dll_run_version
write(nfout,*)

write(nflog,*)
write(nflog,'(a,a)') 'GUI version: ',gui_run_version
write(nflog,'(a,a)') 'DLL version: ',dll_run_version
write(nflog,*)

open(nfres,file='spheroid_ts.out',status='replace')
write(nfres,'(a,a)') 'GUI version: ',gui_run_version
write(nfres,'(a,a)') 'DLL version: ',dll_run_version
write(nfres,*)
write(nfres,'(a)') 'istep hour vol_mm3 diam_um Ncells &
Nanoxia_dead NdrugA_dead NdrugB_dead Nradiation_dead Ntagged_anoxia Ntagged_drugA Ntagged_drugB Ntagged_radiation &
f_hypox_1 f_hypox_2 f_hypox_3 f_growth_1 f_growth_2 f_growth_3 f_necrot plating_efficiency &
medium_oxygen medium_glucose medium_TPZ medium_DNB'

write(logmsg,*) 'Opened nfout: ',outputfile
call logger(logmsg)

Nsteps = days*24*60*60/DELTA_T		! DELTA_T in seconds
write(logmsg,'(a,2i6,f6.0)') 'nsteps, NT_CONC, DELTA_T: ',nsteps,NT_CONC,DELTA_T
call logger(logmsg)

call DetermineKd	! Kd is now set or computed in the GUI
ok = .true.

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getIndices(indx, idrug, nmetab)
integer :: indx, idrug, nmetab

if (indx == 1) then
	idrug = TPZ_DRUG
	nmetab = 2
elseif (indx == 2) then
	idrug = DNB_DRUG
	nmetab = 2
else
	idrug = 0
	nmetab = 0
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Skip lines until the 'PROTOCOL' line
!-----------------------------------------------------------------------------------------
subroutine ReadProtocol(nf)
integer :: nf
integer :: itime, ntimes, kevent, ichemo, im
character*(64) :: line
real(REAL_KIND) :: t, dt, vol, conc, dose
type(event_type) :: E

chemo(TRACER+1:)%used = .false.
do
	read(nf,'(a)') line
	if (trim(line) == 'PROTOCOL') exit
enddo
read(nf,*) ntimes
if (allocated(event)) deallocate(event)
allocate(event(2*ntimes))
kevent = 0
do itime = 1,ntimes
	read(nf,'(a)') line
	if (trim(line) == 'DRUG') then
		kevent = kevent + 1
		event(kevent)%etype = DRUG_EVENT
		read(nf,'(a)') line
		if (trim(line) == 'SN30000' .or. trim(line) == 'TPZ') then
			ichemo = TPZ_DRUG
		elseif (trim(line) == 'PR104A') then
			ichemo = DNB_DRUG
		endif
		read(nf,*) t
		read(nf,*) dt
		read(nf,*) vol
		read(nf,*) conc
		event(kevent)%time = t
		event(kevent)%ichemo = ichemo
		event(kevent)%volume = vol
		event(kevent)%conc = conc
		event(kevent)%dose = 0
		chemo(ichemo)%used = .true.
		if (ichemo == TPZ_DRUG .and. TPZ%use_metabolites) then
			do im = 1,TPZ%nmetabolites
				chemo(ichemo+im)%used = .true.
			enddo
		endif
		if (ichemo == DNB_DRUG .and. DNB%use_metabolites) then
			do im = 1,DNB%nmetabolites
				chemo(ichemo+im)%used = .true.
			enddo
		endif
		kevent = kevent + 1
		event(kevent)%etype = MEDIUM_EVENT
		event(kevent)%time = t + dt
		event(kevent)%ichemo = 0
		event(kevent)%volume = medium_volume0*1.0e3		! -> uL for consistency
		event(kevent)%conc = 0
		event(kevent)%dose = 0
	elseif (trim(line) == 'MEDIUM') then
		kevent = kevent + 1
		event(kevent)%etype = MEDIUM_EVENT
		read(nf,*) t
		read(nf,*) vol
		event(kevent)%time = t
		event(kevent)%volume = vol	
		event(kevent)%ichemo = 0
		event(kevent)%conc = 0
		event(kevent)%dose = 0
	elseif (trim(line) == 'RADIATION') then
		kevent = kevent + 1
		event(kevent)%etype = RADIATION_EVENT
		read(nf,*) t
		read(nf,*) dose
		event(kevent)%time = t
		event(kevent)%dose = dose	
		event(kevent)%ichemo = 0
		event(kevent)%volume = 0
		event(kevent)%conc = 0
	endif
enddo
Nevents = kevent
! Set events not done
! convert time from hours to seconds
! convert volume from uL to cm^3
do kevent = 1,Nevents
	event(kevent)%done = .false.
	event(kevent)%time = event(kevent)%time*60*60
	event(kevent)%volume = event(kevent)%volume*1.0e-3
	E = event(kevent)
!	write(*,'(a,i3,f8.0,2i3,3f8.4)') 'event: ',kevent,E%time,E%etype,E%ichemo,E%volume,E%conc,E%dose
enddo
! Check that events are sequential
do kevent = 1,Nevents-1
	if (event(kevent)%time >= event(kevent+1)%time) then
		write(logmsg,*) 'Error: non-sequential event: ',kevent,event(kevent)%time
		call logger(logmsg)
		stop
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! This overrides the drug usage specified by USE_DRUG_A etc
! RADIATION -> 0
! SN30000 -> 1
! PR104A  -> 2
! By default this assumes two drugs.
! Times are hours
! Drug concs are mM
! Radiation dose is Gy
!-----------------------------------------------------------------------------------------
subroutine ReadTreatment(ok)
logical :: ok
character*(64) :: line
character*(12) :: drug_name
integer :: ichemo, idrug, nmax, i, im
real(REAL_KIND) :: tstart,tend,conc,dose
logical :: use_it(0:2)
allocate(protocol(0:2))

!use_treatment = .true.
chemo(TRACER+1:)%used = .false.
!use_radiation = .false.
open(nftreatment,file=treatmentfile,status='old')
use_it = .false.
nmax = 0
do
	read(nftreatment,'(a)',end=99) line
	if (line(1:3) == 'END' .or. line(1:3) == 'end') exit
	if (line(1:6) == 'DRUG_A') then
		ichemo = -1	! ignore
	elseif (line(1:6) == 'DRUG_B') then
		ichemo = -1	! ignore
	elseif (trim(line) == 'SN30000') then
		ichemo = TPZ_DRUG
		chemo(ichemo)%name = 'SN30000'
		idrug = 1
	elseif (trim(line) == 'TPZ') then
		ichemo = TPZ_DRUG
		chemo(ichemo)%name = 'TPZ'
		idrug = 1
	elseif (trim(line) == 'PR104A') then
		ichemo = DNB_DRUG
		chemo(ichemo)%name = 'PR104A'
		idrug = 2
	elseif (line(1:9) == 'RADIATION') then
		ichemo = 0
		idrug = 0
	else
		write(logmsg,'(a,a)') 'Error: ReadTreatment: unknown drug: ',trim(line)
		ok = .false.
		return
	endif
	if (ichemo > 0) then
!		read(nftreatment,'(a)') chemo(ichemo)%name
		read(nftreatment,*) protocol(idrug)%n
		nmax = max(nmax,protocol(idrug)%n)
		if (protocol(idrug)%n > 0) then
			use_it(idrug) = .true.
			do i = 1,protocol(idrug)%n
				read(nftreatment,*) tstart
				read(nftreatment,*) tend
				read(nftreatment,*) conc
			enddo
		endif
	elseif (ichemo == 0) then
!		use_radiation = .true.
		read(nftreatment,*) protocol(idrug)%n
		nmax = max(nmax,protocol(idrug)%n)
		if (protocol(idrug)%n > 0) then
			use_it(idrug) = .true.
			do i = 1,protocol(idrug)%n
				read(nftreatment,*) tstart
				read(nftreatment,*) dose
			enddo
		endif
	endif
	cycle
99	exit
enddo
rewind(nftreatment)
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
		ichemo = -1
!		idrug = 1
	elseif (line(1:6) == 'DRUG_B') then
		ichemo = -1
!		idrug = 2
	elseif (trim(line) == 'SN30000') then
		ichemo = TPZ_DRUG
		chemo(ichemo)%name = 'SN30000'
		idrug = 1
	elseif (trim(line) == 'TPZ') then
		ichemo = TPZ_DRUG
		chemo(ichemo)%name = 'TPZ'
		idrug = 1
	elseif (trim(line) == 'PR104A') then
		ichemo = DNB_DRUG
		chemo(ichemo)%name = 'PR104A'
		idrug = 2
	elseif (line(1:9) == 'RADIATION') then
		ichemo = 0
		idrug = 0
	endif
	if (ichemo > 0) then
!		read(nftreatment,'(a)') chemo(ichemo)%name
		protocol(idrug)%ichemo = ichemo
		read(nftreatment,*) protocol(idrug)%n
		if (protocol(idrug)%n > 0) then
			chemo(ichemo)%used = .true.
			if (ichemo == TPZ_DRUG .and. TPZ%use_metabolites) then
				do im = 1,TPZ%nmetabolites
					chemo(ichemo+im)%used = .true.
				enddo
			endif
			if (ichemo == DNB_DRUG .and. DNB%use_metabolites) then
				do im = 1,DNB%nmetabolites
					chemo(ichemo+im)%used = .true.
				enddo
			endif
			chemo(ichemo)%bdry_conc = 0
			do i = 1,protocol(idrug)%n
				read(nftreatment,*) tstart
				protocol(idrug)%tstart(i) = 3600*tstart		! hours -> seconds
				read(nftreatment,*) tend
				protocol(idrug)%tend(i) = 3600*tend
				read(nftreatment,*) protocol(idrug)%conc(i)
				write(nflog,*) 'treatment: ',chemo(ichemo)%name,protocol(idrug)%n,i,tstart,tend,protocol(idrug)%conc(i)
			enddo
		endif
	elseif (ichemo == 0) then
		protocol(idrug)%ichemo = 0
		read(nftreatment,*) protocol(idrug)%n
		if (protocol(idrug)%n > 0) then
			do i = 1,protocol(idrug)%n
				read(nftreatment,*) tstart
				protocol(idrug)%tstart(i) = 3600*tstart
				read(nftreatment,*) protocol(idrug)%dose(i)
			enddo
		endif
	endif
	cycle
199	exit
enddo
close(nftreatment)
do i = 0,2
	protocol(i)%started = .false.
	protocol(i)%ended = .false.
enddo	
if (use_it(0)) then
	use_radiation = .true.
endif
if (use_it(0) .or. use_it(1) .or. use_it(2)) then
	use_treatment = .true.
endif
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
! SN30000 CT model
! ----------------
! The rate of cell killing, which becomes a cell death probability rate, is inferred from
! the cell kill experiment.
! The basic assumption is that the rate of killing depends on the drug metabolism rate.
! There are five models:
! kill_model = 1:
!   killing rate = c = Kd.dM/dt
! kill_model = 2:
!   killing rate = c = Kd.Ci.dM/dt
! kill_model = 3:
!   killing rate = c = Kd.(dM/dt)^2
! kill_model = 4:
!   killing rate = c = Kd.Ci
! kill_model = 5:
!   killing rate = c = Kd.Ci^2
! where dM/dt = F(O2).kmet0.Ci
! In the kill experiment both O2 and Ci are held constant:
! O2 = CkillO2, Ci = Ckill
! In this case c is constant and the cell population N(t) is given by:
! N(t) = N(0).exp(-ct), i.e.
! c = -log(N(T)/N(0))/T where T is the duration of the experiment
! N(T)/N(0) = 1 - f, where f = kill fraction
! kill_model = 1:
!   c = Kd.F(CkillO2).kmet0.Ckill => Kd = -log(1-f)/(T.F(CkillO2).kmet0.Ckill)
! kill_model = 2:
!   c = Kd.F(CkillO2).kmet0.Ckill^2 => Kd = -log(1-f)/(T.F(CkillO2).kmet0.Ckill^2)
! kill_model = 3:
!   c = Kd.(F(CkillO2).kmet0.Ckill)^2 => Kd = -log(1-f)/(T.(F(CkillO2).kmet0.Ckill)^2)
! kill_model = 4:
!   c = Kd.Ckill => Kd = -log(1-f)/(T.Ckill)
! kill_model = 5:
!   c = Kd.Ckill^2 => Kd = -log(1-f)/(T.Ckill^2)
!-----------------------------------------------------------------------------------------
subroutine DetermineKd
real(REAL_KIND) :: C2, KO2, Kmet0, kmet 
real(REAL_KIND) :: f, T, Ckill, Ckill_O2, Kd
integer :: idrug, i, im, kill_model

do idrug = 1,2
	if (idrug == 1 .and. .not.chemo(TPZ_DRUG)%used) cycle
	if (idrug == 2 .and. .not.chemo(DNB_DRUG)%used) cycle
	do i = 1,Ncelltypes
		do im = 0,2
			if (idrug == 1) then		! TPZ
				if (im /= 0) cycle
				C2 = TPZ%C2(i,im)
				KO2 = TPZ%KO2(i,im)
				Kmet0 = TPZ%Kmet0(i,im)
				kill_model = TPZ%kill_model(i)
				Ckill_O2 = TPZ%kill_O2(i)
				f = TPZ%kill_fraction(i)
				T = TPZ%kill_duration(i)
				Ckill = TPZ%kill_drug(i)
			elseif (idrug == 2) then	! DNB
				C2 = DNB%C2(i,im)
				KO2 = DNB%KO2(i,im)
				Kmet0 = TPZ%Kmet0(i,im)
				kill_model = DNB%kill_model(i,im)
				Ckill_O2 = DNB%kill_O2(i,im)
				f = DNB%kill_fraction(i,im)
				T = DNB%kill_duration(i,im)
				Ckill = DNB%kill_drug(i,im)
			endif
			kmet = (1 - C2 + C2*KO2/(KO2 + Ckill_O2))*Kmet0
			if (kill_model == 1) then
				Kd = -log(1-f)/(T*kmet*Ckill)
			elseif (kill_model == 2) then
				Kd = -log(1-f)/(T*kmet*Ckill**2)
			elseif (kill_model == 3) then
				Kd = -log(1-f)/(T*(kmet*Ckill)**2)
			elseif (kill_model == 4) then
				Kd = -log(1-f)/(T*Ckill)
			elseif (kill_model == 5) then
				Kd = -log(1-f)/(T*Ckill**2)
			endif
			if (idrug == 1) then
				TPZ%Kd(i) = Kd
			elseif (idrug == 2) then
				DNB%Kd(i,im) = Kd
			endif
		enddo
	enddo
enddo

!if (chemo(TPZ_DRUG)%used) then
!	do i = 1,Ncelltypes
!		kmet = (1 - TPZ%C2(i,0) + TPZ%C2(i,0)*TPZ%KO2(i,0)/(TPZ%KO2(i,0) + TPZ%kill_O2(i)))*TPZ%Kmet0(i,0)
!		kill_model = TPZ%kill_model(i)
!		if (kill_model == 1) then
!			Kd = -log(1-TPZ%kill_fraction(i))/(TPZ%kill_duration(i)*kmet*TPZ%kill_drug(i))
!		elseif (kill_model == 2) then
!			Kd = -log(1-TPZ%kill_fraction(i))/(TPZ%kill_duration(i)*kmet*TPZ%kill_drug(i)**2)
!		elseif (kill_model == 3) then
!			Kd = -log(1-TPZ%kill_fraction(i))/(TPZ%kill_duration(i)*(kmet*TPZ%kill_drug(i))**2)
!		elseif (kill_model == 4) then
!			Kd = -log(1-TPZ%kill_fraction(i))/(TPZ%kill_duration(i)*TPZ%kill_drug(i))
!		elseif (kill_model == 5) then
!			Kd = -log(1-TPZ%kill_fraction(i))/(TPZ%kill_duration(i)*TPZ%kill_drug(i)**2)
!		endif
!		TPZ%Kd(i) = Kd
!	enddo
!endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine PlaceCells(ok)
logical :: ok
integer :: x, y, z, kcell, site(3), ichemo
real(REAL_KIND) :: r2lim,r2,rad(3)

occupancy(:,:,:)%indx(1) = 0
occupancy(:,:,:)%indx(2) = 0
r2lim = 0.95*Radius*Radius
lastID = 0
kcell = 0
do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
			rad = (/x-x0,y-y0,z-z0/)
			r2 = dot_product(rad,rad)
			if (r2 < r2lim) then
				kcell = kcell+1
				site = (/x,y,z/)
				call AddCell(kcell,site)
				if (x == NX/2 .and. y == NY/2 .and. z == NZ/2) then
					idbug = kcell
					write(nfout,*) 'Mid-blob cell: idbug: ',idbug, x,y,z
				endif
			else
				occupancy(x,y,z)%indx(1) = OUTSIDE_TAG
			endif
		enddo
	enddo
enddo

if (kcell > initial_count) then
	write(logmsg,*) 'Cell count already exceeds specified number: ',kcell,initial_count
	call logger(logmsg)
	ok = .false.
	return
endif
! Now add cells to make the count up to the specified initial_count
if (kcell < initial_count) then
	call AddBdryCells(kcell)
	kcell = initial_count
endif
	
do ichemo = 1,MAX_CHEMO
    occupancy(:,:,:)%C(ichemo) = 0
enddo
occupancy(:,:,:)%C(OXYGEN) = chemo(OXYGEN)%bdry_conc
occupancy(:,:,:)%C(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
occupancy(:,:,:)%C(TRACER) = chemo(TRACER)%bdry_conc
nlist = kcell
Nsites = kcell
Ncells = kcell
Ncells0 = Ncells
Nreuse = 0	
ok = .true.
!write(logmsg,*) 'idbug: ',idbug
!call logger(logmsg)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine AddCell(k,site)
integer :: k, site(3)
integer :: kpar = 0
real(REAL_KIND) :: R

lastID = lastID + 1
cell_list(k)%ID = lastID
cell_list(k)%celltype = random_choice(celltype_fraction,Ncelltypes,kpar)
cell_list(k)%site = site
cell_list(k)%state = 1
cell_list(k)%drugA_tag = .false.
cell_list(k)%drugB_tag = .false.
cell_list(k)%radiation_tag = .false.
cell_list(k)%anoxia_tag = .false.
cell_list(k)%exists = .true.
cell_list(k)%active = .true.
R = par_uni(kpar)
cell_list(k)%divide_volume = Vdivide0 + dVdivide*(2*R-1)
R = par_uni(kpar)
if (randomise_initial_volume) then
	cell_list(k)%volume = cell_list(k)%divide_volume*0.5*(1 + R)
else
	cell_list(k)%volume = 1.0
endif
!write(nflog,'(a,i6,f6.2)') 'volume: ',k,cell_list(k)%volume
cell_list(k)%t_divide_last = 0		! not used
cell_list(k)%t_hypoxic = 0
cell_list(k)%conc = 0
cell_list(k)%conc(OXYGEN) = chemo(OXYGEN)%bdry_conc
cell_list(k)%conc(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
cell_list(k)%conc(TRACER) = chemo(TRACER)%bdry_conc
cell_list(k)%CFSE = generate_CFSE(1.d0)
cell_list(k)%M = 0
occupancy(site(1),site(2),site(3))%indx(1) = k
end subroutine

!--------------------------------------------------------------------------------
! Add cells at the boundary to bring the total count from k up to initial_count
! (1) Make a list of all boundary sites (sites in contact with an OUTSIDE site)
! (2) Iteratively traverse the list to select the adjacent OUTSIDE site closest 
! to the centre.
!--------------------------------------------------------------------------------
subroutine AddBdryCells(klast)
integer :: klast
integer :: kcell, i, kb, site(3), nbsite(3), nbt, kbmin, imin
integer, allocatable :: sitelist(:,:)
real(REAL_KIND) :: r2, r2min

nbt = 0
do kcell = 1,klast
	site = cell_list(kcell)%site
	do i = 1,27
		if (i == 14) cycle
		nbsite = site + jumpvec(:,i)
		if (occupancy(nbsite(1),nbsite(2),nbsite(3))%indx(1) == OUTSIDE_TAG) then
			nbt = nbt+1
			exit
		endif
	enddo
enddo

allocate(sitelist(3,nbt))

nbt = 0
do kcell = 1,klast
	site = cell_list(kcell)%site
	do i = 1,27
		if (i == 14) cycle
		nbsite = site + jumpvec(:,i)
		if (occupancy(nbsite(1),nbsite(2),nbsite(3))%indx(1) == OUTSIDE_TAG) then
			nbt = nbt+1
			sitelist(:,nbt) = site
			exit
		endif
	enddo
enddo
	

do kcell = klast+1,initial_count
	r2min = 1.0e10
	do kb = 1,nbt
		site = sitelist(:,kb)
		do i = 1,27
			if (i == 14) cycle
			nbsite = site + jumpvec(:,i)
			if (occupancy(nbsite(1),nbsite(2),nbsite(3))%indx(1) == OUTSIDE_TAG) then
				r2 = (nbsite(1) - Centre(1))**2 + (nbsite(2) - Centre(2))**2 + (nbsite(3) - Centre(3))**2
				if (r2 < r2min) then
					kbmin = kb
					imin = i
					r2min = r2
				endif
			endif
		enddo
	enddo
	site = sitelist(:,kbmin) + jumpvec(:,imin)
	call AddCell(kcell,site)
enddo

deallocate(sitelist)
		
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_dimensions(NX_dim, NY_dim, NZ_dim, nsteps_dim, deltat, maxchemo, nextra, cused, dfraction, deltax) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dimensions
use, intrinsic :: iso_c_binding
integer(c_int) :: NX_dim,NY_dim,NZ_dim,nsteps_dim, maxchemo, nextra
real(c_double) :: deltat, dfraction, deltax
logical(c_bool) :: cused(*)
integer :: ichemo

NX_dim = NX
NY_dim = NY
NZ_dim = NZ
nsteps_dim = nsteps
deltat = DELTA_T
deltax = DELTA_X
maxchemo = MAX_CHEMO
nextra = N_EXTRA
do ichemo = 1,MAX_CHEMO
	cused(ichemo+1) = chemo(ichemo)%used
enddo
cused(1) = .true.			! CFSE
cused(MAX_CHEMO+2) = .true.	! Growth rate
dfraction = 2*cell_radius/DELTA_X
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine InitConcs1
integer :: ic, ichemo, kcell, site(3), ntvars

write(logmsg,*) 'InitConcs: ',nchemo
call logger(logmsg)

do kcell = 1,Ncells
    site = cell_list(kcell)%site
    do ic = 1,nchemo
	    ichemo = chemomap(ic)
        occupancy(site(1),site(2),site(3))%C(ichemo) = chemo(ichemo)%bdry_conc
    enddo
enddo
ntvars = ODEdiff%nextra + ODEdiff%nintra
do ic = 1,nchemo
	ichemo = chemomap(ic)
	allstate(1:ntvars,ichemo) = chemo(ichemo)%bdry_conc
enddo
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine InitConcs(ichemo)
integer :: ichemo
integer :: i, kcell, site(3), ntvars
real(REAL_KIND) :: c0

ntvars = ODEdiff%nextra + ODEdiff%nintra
write(logmsg,*) 'InitConcs: ',ichemo,ntvars
call logger(logmsg)

if (istep == 0) then
	if (ichemo == OXYGEN .or. ichemo == GLUCOSE .or. ichemo == TRACER) then
		c0 = chemo(ichemo)%bdry_conc
	else	! drug or metabolite
		c0 = 0
	endif
	do kcell = 1,Ncells
		site = cell_list(kcell)%site
		occupancy(site(1),site(2),site(3))%C(ichemo) = c0
	enddo
	!allstate(1:ntvars,ichemo) = c0
	allstate(:,ichemo) = c0
else	! drug or metabolite
	do i = 1,ODEdiff%nvars
		site = ODEdiff%varsite(i,:)
		if (ODEdiff%vartype(i) == EXTRA) then
			occupancy(site(1),site(2),site(3))%C(ichemo) = 0
		endif
	enddo
endif
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine ProcessEvent(radiation_dose)
real(REAL_KIND) :: radiation_dose
integer :: kevent, ichemo, im, nmetab
real(REAL_KIND) :: V, C(MAX_CHEMO)
type(event_type) :: E

do kevent = 1,Nevents
	E = event(kevent)
	if (t_simulation >= E%time .and. .not.E%done) then
!	write(*,'(i3,2f8.0,i3,2f10.4)') E%etype,t_simulation,E%time,E%ichemo,E%volume,E%conc
		if (E%etype == RADIATION_EVENT) then
			radiation_dose = E%dose
!			write(*,*) 'radiation_dose: ',radiation_dose
		elseif (E%etype == MEDIUM_EVENT) then
			C = 0
			C(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
			V = E%volume
			call MediumChange(V,C)
		elseif (E%etype == DRUG_EVENT) then
			C = 0
			C(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
			ichemo = E%ichemo
			C(ichemo) = E%conc
			V = E%volume
			! set %present
			chemo(ichemo)%present = .true.
			chemo(ichemo)%bdry_conc = 0
			if (ichemo == TPZ_DRUG) then
				nmetab = TPZ%nmetabolites
			elseif (ichemo == DNB_DRUG) then
				nmetab = DNB%nmetabolites
			endif			
			do im = 1,nmetab
				if (chemo(ichemo + im)%used) then
					chemo(ichemo + im)%present = .true.
					chemo(ichemo + im)%bdry_conc = 0
				endif
			enddo
			call MediumChange(V,C)
		endif
		event(kevent)%done = .true.
	endif
enddo
end subroutine

!----------------------------------------------------------------------------------
! Radiation treatment is stored in protocol(0)
!----------------------------------------------------------------------------------
subroutine Treatment(radiation_dose)
real(REAL_KIND) :: radiation_dose
integer :: i, idrug, ichemo, nmetab, im	!, ichemo_metab

radiation_dose = 0
do i = 1,protocol(0)%n
	if (t_simulation >= protocol(0)%tstart(i) .and. .not.protocol(0)%started(i)) then
		radiation_dose = protocol(0)%dose(i)
		protocol(0)%started(i) = .true.
		protocol(0)%ended(i) = .true.
		write(nflog,*) 'Radiation started: dose: ',radiation_dose
		exit
	endif
enddo
do idrug = 1,2
!	ichemo = idrug + TRACER		!!!!!!!!!!!!!!!!!! wrong
	ichemo = protocol(idrug)%ichemo
	if (idrug == 1) then
!		ichemo_metab = DRUG_A_METAB
		nmetab = 2
	elseif (idrug == 2) then
!		ichemo_metab = DRUG_B_METAB
		nmetab = 2
	endif
	do i = 1,protocol(idrug)%n
		if (i == 1 .and. t_simulation < protocol(idrug)%tstart(i)) then
			chemo(ichemo)%bdry_conc = 0
			do im = 1,nmetab
				chemo(ichemo + im)%bdry_conc = 0
			enddo
			exit
		endif
		if (t_simulation >= protocol(idrug)%tstart(i) .and. .not.protocol(idrug)%started(i)) then
			chemo(ichemo)%bdry_conc = protocol(idrug)%conc(i)
			protocol(idrug)%started(i) = .true.
			protocol(idrug)%ended(i) = .false.
			chemo(ichemo)%present = .true.
			call InitConcs(ichemo)
			call SetupMedium(ichemo)
			do im = 1,nmetab
				if (chemo(ichemo + im)%used) then
					chemo(ichemo + im)%present = .true.
					call InitConcs(ichemo + im)
					call SetupMedium(ichemo + im)
				endif
			enddo
			write(nflog,*) 'Started DRUG: ',chemo(ichemo)%name,chemo(ichemo)%bdry_conc, i
			exit
		endif
	enddo
	do i = 1,protocol(idrug)%n
		if (t_simulation >= protocol(idrug)%tend(i) .and. .not.protocol(idrug)%ended(i)) then
			chemo(ichemo)%bdry_conc = 0
			protocol(idrug)%ended(i) = .true.
!			chemo(ichemo)%present = .false.
			call InitConcs(ichemo)
			call SetupMedium(ichemo)
			do im = 1,nmetab
				chemo(ichemo + im)%bdry_conc = 0
				if (chemo(ichemo + im)%used) then
	!				chemo(ichemo_metab)%present = .false.
					call InitConcs(ichemo + im)
					call SetupMedium(ichemo + im)
				endif
			enddo
			write(nflog,*) 'Ended DRUG: ',chemo(ichemo)%name,i
			exit
		endif
	enddo
enddo	
end subroutine

!-----------------------------------------------------------------------------------------
! If the volume removed is Vr, the fraction of constituent mass that is retained
! in the medium is (Vm - Vr)/Vm.  The calculation does not apply to oxygen.
! Usually Vr = Ve.
!-----------------------------------------------------------------------------------------
subroutine MediumChange(Ve,Ce)
real(REAL_KIND) :: Ve, Ce(:)
real(REAL_KIND) :: R, Vm, Vr, Vblob

call SetRadius(Nsites)
R = Radius*DELTA_X		! cm
Vblob = (4./3.)*PI*R**3	! cm3
Vm = total_volume - Vblob
Vr = min(Vm,Ve)
chemo(OXYGEN+1:)%medium_M = ((Vm - Vr)/Vm)*chemo(OXYGEN+1:)%medium_M + Vr*Ce(OXYGEN+1:)
total_volume = Vm - Vr + Ve + Vblob
chemo(OXYGEN+1:)%medium_Cext = chemo(OXYGEN+1:)%medium_M/(total_volume - Vblob)
chemo(OXYGEN)%medium_Cext = chemo(OXYGEN)%bdry_conc
call UpdateCbnd
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
integer :: kcell, site(3), hour, nthour, kpar=0
real(REAL_KIND) :: r(3), rmax, tstart, dt, radiation_dose, diam_um
!integer, parameter :: NT_CONC = 6
integer :: nchemo, i
logical :: ok = .true.
logical :: dbug

!call logger('simulate_step')
dbug = .false.
if (Ncells == 0) then
	call logger('Ncells = 0')
    res = 3
    return
endif
nthour = 3600/DELTA_T
dt = DELTA_T/NT_CONC

bdry_debug = (istep >= 25000)
if (use_dropper .and. Ncells >= Ndrop .and. .not.is_dropped) then
    call shaper
    call dropper
endif

!if (bdry_debug) call CheckBdryList('simulate_step a')
if (bdry_debug) write(*,*) 'istep, bdry_changed: ',istep,bdry_changed
if (dbug .or. mod(istep,nthour) == 0) then
	diam_um = 2*DELTA_X*Radius*10000
	write(logmsg,'(a,4i8,f8.1)') 'istep, hour: ',istep,istep/nthour,nlist,ncells,diam_um
	call logger(logmsg)
endif
	if (bdry_changed) then
		if (dbug) write(nflog,*) 'UpdateBdryList'
		call UpdateBdrylist
!		write(nflog,*) 'did UpdateBdryList'
	endif
	if (mod(istep,6*nthour) == 0) then
		if (dbug) write(nflog,*) 'CheckBdryList'
		call CheckBdryList('simulate_step')
!		write(nflog,*) 'did CheckBdryList'
	endif
istep = istep + 1
t_simulation = (istep-1)*DELTA_T	! seconds
radiation_dose = 0
if (use_treatment) then
	call treatment(radiation_dose)
endif
if (use_events) then
	call ProcessEvent(radiation_dose)
endif
if (radiation_dose > 0) then
	write(logmsg,'(a,f6.1)') 'Radiation dose: ',radiation_dose
	call logger(logmsg)
endif
if (bdry_debug) call CheckBdryList('simulate_step c')
if (dbug) write(nflog,*) 'GrowCells'
call GrowCells(radiation_dose,DELTA_T,ok)
if (bdry_debug) call CheckBdryList('simulate_step d')
!write(nflog,*) 'did GrowCells'
if (.not.ok) then
	res = 3
	return
endif

! Update Cbnd using current M, R1 and previous U, Cext
if (dbug) write(nflog,*) 'UpdateCbnd'
call UpdateCbnd
!write(nflog,*) 'did UpdateCbnd'

if (dbug) write(nflog,*) 'SetupODEdiff'
call SetupODEdiff

!call check_allstate('before SiteCellToState')
if (dbug) write(nflog,*) 'SiteCellToState'
call SiteCellToState
!call check_allstate('after SiteCellToState')

if (dbug) write(nflog,*) 'Solver'
do it_solve = 1,NT_CONC
	tstart = (it_solve-1)*dt
	t_simulation = (istep-1)*DELTA_T + tstart
	call Solver(it_solve,tstart,dt,Ncells)
enddo
!write(nflog,*) 'did Solver'
!call check_allstate('after Solver')

if (dbug) write(nflog,*) 'StateToSiteCell'
call StateToSiteCell
!call check_allstate('after StateToSiteCell')

res = 0

! Compute U and update M, Cext, Cbnd
if (dbug) write(nflog,*) 'UpdateMedium'
call UpdateMedium(DELTA_T)
!write(nflog,*) 'did UpdateMedium'

if (mod(istep,60) == -1) then
	rmax = 0
	do kcell = 1,nlist
		r = cell_list(kcell)%site - Centre
		rmax = max(rmax,norm(r))
	enddo
	hour = istep/60
	write(logmsg,'(3i6,2f6.1)') istep, hour, Ncells, Radius, rmax
	call logger(logmsg)
!	call ShowConcs
	call check_bdry
endif
!call test_CellDivision
if (saveprofiledata) then
	if (istep*DELTA_T >= it_saveprofiledata*dt_saveprofiledata) then
		call WriteProfileData
		it_saveprofiledata = it_saveprofiledata + 1
		if (it_saveprofiledata > nt_saveprofiledata) then
			saveprofiledata = .false.
		endif
	endif
endif
end subroutine

!--------------------------------------------------------------------------------
! TC = tumour cell
!--------------------------------------------------------------------------------
subroutine get_scene(nTC_list,TC_list) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scene
use, intrinsic :: iso_c_binding
integer(c_int) :: nTC_list, TC_list(*)
integer :: k, kc, kcell, site(3), j, jb, colour, nlive
integer :: col(3)
integer :: x, y, z
real(REAL_KIND) :: vol, r
integer :: itcstate, ctype, stage, region, highlight
integer :: last_id1, last_id2
logical :: ok, highlighting
integer, parameter :: axis_centre = -2	! identifies the spheroid centre
integer, parameter :: axis_end    = -3	! identifies the spheroid extent in 5 directions
integer, parameter :: axis_bottom = -4	! identifies the spheroid extent in the -Y direction, i.e. bottom surface
integer, parameter :: ninfo = 7			! the size of the info package for a cell (number of integers)
integer, parameter :: nax = 6			! number of points used to delineate the spheroid

highlighting = (show_progeny /= 0)
nTC_list = 0

k = 0
last_id1 = 0
if (.false.) then
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
		TC_list(j+6) = 1
		last_id1 = k-1
	enddo
	k = last_id1 + 1
endif

! Cells
nlive = 0
do kcell = 1,nlist
!	if (idbug /= 0 .and. cell_list(kcell)%ID /= idbug) cycle
	if (cell_list(kcell)%exists) then
		k = k+1
		j = ninfo*(k-1)
		site = cell_list(kcell)%site
	    if (highlighting .and. cell_list(kcell)%ID == show_progeny) then
	        highlight = 1
	    else
	        highlight = 0
	    endif
	    if (use_celltype_colour) then
			colour = cell_list(kcell)%celltype
		else
			call CellColour(kcell,highlight,col)
			colour = rgb(col)
		endif
		TC_list(j+1) = kcell + last_id1
		TC_list(j+2:j+4) = site
		TC_list(j+5) = colour
		vol = cell_list(kcell)%volume*Vcell_cm3		! cell volume in cm3
		r = (3*vol/(4*PI))**(1./3.)					! cell radius in cm
		TC_list(j+6) = 200*r/DELTA_X				! 100*diameter as fraction of DELTA_X
!		TC_list(j+6) = (cell_list(kcell)%volume)**(1./3.)	! diameter: 0.928 - 1.17
		TC_list(j+7) = highlight
		last_id2 = kcell + last_id1
		nlive = nlive + 1
	endif
enddo
!nTC_list = last_id2
nTC_list = k
!write(logmsg,*) 'nlive: ',nlive
!call logger(logmsg)
end subroutine

!-----------------------------------------------------------------------------------------
! Rendered cell colour may depend on stage, state, receptor expression level.
! col(:) = (r,g,b)
!-----------------------------------------------------------------------------------------
subroutine CellColour(kcell,highlight,col)
integer :: kcell, highlight, col(3)
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

if (highlight == 0) then
    col = LIGHTORANGE
else
    col = LIGHTRED
endif
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
! Total tagged for drug death on division = NdrugA_tag, NdrugB_tag
! Current tagged = Ntodie - Ntagdead 
!-----------------------------------------------------------------------------------------
subroutine get_summary(summaryData,i_hypoxia_cutoff,i_growth_cutoff) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*), i_hypoxia_cutoff,i_growth_cutoff
integer :: Ndead, Ntagged, Ntodie, Ntagdead, Ntagged_anoxia, &
	Ntagged_drugA, Ntagged_drugB, Ntagged_radiation, nviable, &
    diam_um, vol_mm3_1000, nhypoxic(3), ngrowth(3), &
    hypoxic_percent_10, growth_percent_10, necrotic_percent_10, plate_eff_10, &
    medium_oxygen_100, medium_glucose_100, medium_TPZ_drug_1000, medium_DNB_drug_1000
    
real(REAL_KIND) :: vol_cm3, vol_mm3, hour, plate_eff

hour = istep*DELTA_T/3600.
vol_cm3 = Vsite_cm3*Nsites			! total volume in cm^3
vol_mm3 = vol_cm3*1000				! volume in mm^3
vol_mm3_1000 = vol_mm3*1000			! 1000 * volume in mm^3
diam_um = 2*DELTA_X*Radius*10000
!Ntodie = Nradiation_tag + Ndrug_tag			! total that have been tagged by drug or radiation
!Ntagdead = Nradiation_dead + Ndrug_dead		! total that died from drug or radiation
!Ndead = Nsites + Nreuse - Ncells				! total that died from any cause
!Ntagged = Ntodie - Ntagdead					! number currently tagged by drug or radiation
Ntagged_anoxia = Nanoxia_tag - Nanoxia_dead				! number currently tagged by anoxia
Ntagged_drugA = NdrugA_tag - NdrugA_dead				! number currently tagged by drugA
Ntagged_drugB = NdrugB_tag - NdrugB_dead				! number currently tagged by drugB
Ntagged_radiation = Nradiation_tag - Nradiation_dead	! number currently tagged by radiation
call getHypoxicCount(nhypoxic)
hypoxic_percent_10 = (1000*nhypoxic(i_hypoxia_cutoff))/Ncells
call getGrowthCount(ngrowth)
growth_percent_10 = (1000*ngrowth(i_growth_cutoff))/Ncells
necrotic_percent_10 = (1000*(Nsites-Ncells))/Nsites
call getNviable(nviable)
plate_eff = real(nviable)/Ncells
plate_eff_10 = 1000*plate_eff
medium_oxygen_100 = 100*chemo(OXYGEN)%medium_Cext
medium_glucose_100 = 100*chemo(GLUCOSE)%medium_Cext
medium_TPZ_drug_1000 = 1000*chemo(TPZ_DRUG)%medium_Cext
medium_DNB_drug_1000 = 1000*chemo(DNB_DRUG)%medium_Cext
summaryData(1:20) = [ istep, Ncells, Nanoxia_dead, NdrugA_dead, NdrugB_dead, Nradiation_dead, &
    Ntagged_anoxia, Ntagged_drugA, Ntagged_drugB, Ntagged_radiation, &
	diam_um, vol_mm3_1000, hypoxic_percent_10, growth_percent_10, necrotic_percent_10, plate_eff_10, &
	medium_oxygen_100, medium_glucose_100, medium_TPZ_drug_1000, medium_DNB_drug_1000 ]
write(nfres,'(2a12,i8,f8.2,f8.4,10i7,12f7.3)') gui_run_version, dll_run_version, istep, hour, vol_mm3, diam_um, Ncells, &
    Nanoxia_dead, NdrugA_dead, NdrugB_dead, Nradiation_dead, Ntagged_anoxia, Ntagged_drugA, Ntagged_drugB, Ntagged_radiation, &
	nhypoxic(:)/real(Ncells), ngrowth(:)/real(Ncells), (Nsites-Ncells)/real(Nsites), plate_eff, &
	chemo(OXYGEN)%medium_Cext, chemo(GLUCOSE)%medium_Cext, chemo(TPZ_DRUG)%medium_Cext, chemo(DNB_DRUG)%medium_Cext
	
	
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getNviable(nviable)
integer :: nviable
integer :: kcell

nviable = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	if (cell_list(kcell)%anoxia_tag .or. &
	    cell_list(kcell)%drugA_tag .or. &
	    cell_list(kcell)%drugB_tag .or. &
	    cell_list(kcell)%radiation_tag) cycle
	nviable = nviable + 1
enddo	
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getHypoxicCount(nhypoxic)
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
! Need to compare growth rate with a fraction of average growth rate
!--------------------------------------------------------------------------------
subroutine getGrowthCount(ngrowth)
integer :: ngrowth(3)
integer :: kcell, i
real(REAL_KIND) :: r_mean

r_mean = Vdivide0/(2*divide_time_mean)
ngrowth = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	do i = 1,3
		if (cell_list(kcell)%dVdt < growthcutoff(i)*r_mean) ngrowth(i) = ngrowth(i) + 1
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
! Note: axis = 0,1,2
! Now CFSE is first in the list (0)
!--------------------------------------------------------------------------------
subroutine get_fieldinfo(nxx, axis, fraction, ns, nc, cused, res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_fieldinfo
use, intrinsic :: iso_c_binding
integer(c_int) :: nxx, axis, ns, nc, cused(0:*), res
real(c_double) :: fraction
integer rng(3,2), ichemo, kcell, x, y, z

!call logger('get_fieldinfo')
res = 0
nxx = NX
nc = MAX_CHEMO
do ichemo = 1,MAX_CHEMO
    if (chemo(ichemo)%used) then
        cused(ichemo) = 1
    else
        cused(ichemo) = 0
    endif
enddo
cused(CFSE) = 1
cused(GROWTH_RATE) = 1
cused(CELL_VOLUME) = 1
cused(O2_BY_VOL) = 1
rng(:,1) = Centre(:) - (adrop*Radius + 2)
rng(:,2) = Centre(:) + (adrop*Radius + 2)
rng(axis,:) = Centre(axis) + fraction*Radius
!write(logmsg,*) 'Centre, Radius, axis, fraction: ',Centre, Radius, axis, fraction
!call logger(logmsg)
!write(logmsg,*) 'rng: ',rng
!call logger(logmsg)
ns = 0
do z = rng(3,1),rng(3,2)
    do y = rng(2,1),rng(2,2)
        do x = rng(1,1),rng(1,2)
            kcell = occupancy(x,y,z)%indx(1)
            if (kcell <= OUTSIDE_TAG) cycle
            ns = ns + 1
        enddo
    enddo
enddo
!write(logmsg,*) 'get_fieldinfo: ns: ',ns
!call logger(logmsg)
end subroutine

!--------------------------------------------------------------------------------
! Need to transmit medium concentration data.  This could be a separate subroutine.
!--------------------------------------------------------------------------------
subroutine get_fielddata(axis, fraction, nfdata, nc, fdata, res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_fielddata
use, intrinsic :: iso_c_binding
real(c_double) :: fraction
integer(c_int) :: axis, nc, nfdata, res
type(FIELD_DATA) :: fdata(*)
integer rng(3,2), kcell, x, y, z, i, ns
real(REAL_KIND) :: growthrate, cellvolume, cfse

!write(logmsg,*) 'get_fielddata: nfdata, nc: ',nfdata, nc, MAX_CHEMO
!call logger(logmsg)
res = 0
if (nc > MAX_CHEMO) then
	write(logmsg,*) 'Error: get_fielddata: dimension of conc(MAX_CHEMO) not big enough!'
	call logger(logmsg)
	res = 1
	return
endif
rng(:,1) = Centre(:) - (adrop*Radius + 2)
rng(:,2) = Centre(:) + (adrop*Radius + 2)
rng(axis,:) = Centre(axis) + fraction*Radius
ns = 0
do z = rng(3,1),rng(3,2)
    do y = rng(2,1),rng(2,2)
        do x = rng(1,1),rng(1,2)
            kcell = occupancy(x,y,z)%indx(1)
            if (kcell <= OUTSIDE_TAG) cycle
!            write(nflog,*) x,y,z,kcell
            ns = ns + 1
	        i = ODEdiff%ivar(x,y,z)
            fdata(ns)%site = (/x,y,z/)
            fdata(ns)%state = 1
            if (kcell > 0) then
                fdata(ns)%volume = cell_list(kcell)%volume
                call getExtraVars(kcell,growthrate,cellvolume,cfse)
            else
                fdata(ns)%volume = 0
                growthrate = 0
                cellvolume = 0
                cfse = 0
            endif
!            fdata(ns)%dVdt = growthrate
			fdata(ns)%conc(0) = cfse
            fdata(ns)%conc(1:MAX_CHEMO) = allstate(i,1:MAX_CHEMO)	! cell_list(kcell)%conc(1:MAX_CHEMO)
            fdata(ns)%conc(GROWTH_RATE) = growthrate
            fdata(ns)%conc(CELL_VOLUME) = Vcell_pL*cellvolume		! Note: = fdata(ns)%volume, therefore redundant
            fdata(ns)%conc(O2_BY_VOL) = allstate(i,OXYGEN)*Vcell_pL*cellvolume	
        enddo
    enddo
enddo
!write(logmsg,*) 'get_fielddata: ns: ',ns
!call logger(logmsg)
if (ns /= nfdata) then
    write(logmsg,*) 'Error: inconsistent nsites: ',nfdata, ns
    call logger(logmsg)
    res = 2
    return
endif

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getExtraVars(kcell,growthrate,cellvolume,cfse)
integer :: kcell
real(REAL_KIND) :: growthrate, cellvolume, cfse

cfse = cell_list(kcell)%CFSE
growthrate = cell_list(kcell)%dVdt
cellvolume = cell_list(kcell)%volume
end subroutine

!--------------------------------------------------------------------------------
! Returns all the extracellular concentrations along a line through the blob centre.
! Together with CFSE, growth rate (dVdt), cell volume,...
! Store the constituent profiles one after the other.
! Drop the extra end points
!--------------------------------------------------------------------------------
subroutine get_concdata(nvars, ns, dx, ex_conc) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_concdata
use, intrinsic :: iso_c_binding
integer(c_int) :: nvars, ns
real(c_double) :: dx, ex_conc(0:*)
real(REAL_KIND) :: cbnd, cmin = 1.0e-6
integer :: rng(3,2), i, ic, k, ichemo, kcell, x, y, z, x1, x2, offset

!call logger('get_concdata')
nvars = 1 + MAX_CHEMO + N_EXTRA
dx = DELTA_X
rng(:,1) = Centre(:) - (adrop*Radius + 2)
rng(:,2) = Centre(:) + (adrop*Radius + 2)
!rng(axis,:) = Centre(axis) + fraction*Radius
y = Centre(2) + 0.5
z = Centre(3) + 0.5

! First need to establish the range of x that is inside the blob: (x1,x2)
x1 = 0
do x = rng(1,1),rng(1,2)
	kcell = occupancy(x,y,z)%indx(1)
	if (kcell <= OUTSIDE_TAG) then
		if (x1 == 0) then
			cycle
		else
			exit
		endif
	elseif (x1 == 0) then
		x1 = x
	endif
	x2 = x
enddo

ns = x2 - x1 + 1 
do ichemo = 0,nvars-1
	offset = ichemo*ns
	k = offset - 1
	do x = x1, x2
		k = k + 1
		kcell = occupancy(x,y,z)%indx(1)
		if (kcell <= OUTSIDE_TAG) then
			ex_conc(k) = 0
			cycle
		endif
		i = ODEdiff%ivar(x,y,z)
        if (ichemo == 0) then	! CFSE
			if (kcell > 0) then
				ex_conc(k) = cell_list(kcell)%CFSE
			else
				ex_conc(k) = 0
			endif
       elseif (ichemo <= MAX_CHEMO) then
			if (chemo(ichemo)%used) then
				if (i > 0) then
					ex_conc(k) = allstate(i,ichemo)	
				else
					ex_conc(k) = 0
				endif
			else
				ex_conc(k) = 0
			endif
        elseif (ichemo == GROWTH_RATE) then	
			if (kcell > 0) then
				ex_conc(k) = cell_list(kcell)%dVdt
			else
				ex_conc(k) = 0
			endif
        elseif (ichemo == CELL_VOLUME) then	
			if (kcell > 0) then
				ex_conc(k) = Vcell_pL*cell_list(kcell)%volume
			else
				ex_conc(k) = 0
			endif
        elseif (ichemo == O2_BY_VOL) then	
			if (kcell > 0) then
				ex_conc(k) = allstate(i,OXYGEN)*Vcell_pL*cell_list(kcell)%volume
			else
				ex_conc(k) = 0
			endif
		endif
    enddo
enddo
! Add concentrations at the two boundaries 
! At ns=1, at at ns=ns+1
!ns = ns+1
!do ic = 1,nvars
!	ichemo = ic - 1
!	if (ichemo == 0) then	! CFSE
!		k = ic	
!		ex_conc(k) = 0
!        k = (ns-1)*nvars + ic
!        ex_conc(k) = 0	
!	elseif (ichemo <= MAX_CHEMO) then
!		k = ic
!		if (chemo(ichemo)%used) then
!			ex_conc(k) = BdryConc(ichemo,t_simulation)
!		else
!			ex_conc(k) = 0
!		endif      
!		k = (ns-1)*nvars + ic
!		if (chemo(ichemo)%used) then
!			ex_conc(k) = BdryConc(ichemo,t_simulation)
!		else
!			ex_conc(k) = 0
!		endif      
!    elseif (ichemo == MAX_CHEMO+1) then	! growth rate
!		k = ic
!		ex_conc(k) = 0
!! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,1,i,k,ex_conc(k)
!        k = (ns-1)*nvars + ic
!        ex_conc(k) = 0
!! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,ns,i,k,ex_conc(k)
!    elseif (ichemo == MAX_CHEMO+2) then	! cell volume
!		k = ic
!		ex_conc(k) = 0
!! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,1,i,k,ex_conc(k)
!        k = (ns-1)*nvars + ic
!        ex_conc(k) = 0
!! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,ns,i,k,ex_conc(k)
!    endif
!enddo
end subroutine

!--------------------------------------------------------------------------------
! Returns all the extracellular concentrations along a line through the blob centre.
! Together with CFSE, growth rate (dVdt), cell volume,...
!--------------------------------------------------------------------------------
subroutine old_get_concdata(ns, dx, ex_conc) BIND(C)
!!DEC$ ATTRIBUTES DLLEXPORT :: get_concdata
use, intrinsic :: iso_c_binding
integer(c_int) :: ns
real(c_double) :: dx, ex_conc(*)
real(REAL_KIND) :: cbnd, cmin = 1.0e-6
integer :: rng(3,2), i, ic, k, ichemo, kcell, x, y, z, nvars

nvars = 1 + MAX_CHEMO + N_EXTRA
!call logger('get_concdata')
dx = DELTA_X
rng(:,1) = Centre(:) - (adrop*Radius + 2)
rng(:,2) = Centre(:) + (adrop*Radius + 2)
!rng(axis,:) = Centre(axis) + fraction*Radius
y = Centre(2) + 0.5
z = Centre(3) + 0.5
ns = 1
do x = rng(1,1),rng(1,2)
    kcell = occupancy(x,y,z)%indx(1)
    if (kcell <= OUTSIDE_TAG) cycle
	i = ODEdiff%ivar(x,y,z)
    ns = ns + 1
    do ichemo = 1,nvars
        k = (ns-1)*nvars + ic
        if (ichemo == 0) then	! CFSE
			if (kcell > 0) then
				ex_conc(k) = cell_list(kcell)%CFSE
			else
				ex_conc(k) = 0
			endif
       elseif (ichemo <= MAX_CHEMO) then
			if (chemo(ichemo)%used) then
				if (i > 0) then
					ex_conc(k) = allstate(i,ichemo)	
				else
					ex_conc(k) = 0
				endif
			else
				ex_conc(k) = 0
			endif
        elseif (ichemo == MAX_CHEMO+1) then	! growth rate
			if (kcell > 0) then
				ex_conc(k) = cell_list(kcell)%dVdt
			else
				ex_conc(k) = 0
			endif
! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,ns,kcell,k,ex_conc(k)
        elseif (ichemo == MAX_CHEMO+2) then	! cell volume
			if (kcell > 0) then
				ex_conc(k) = cell_list(kcell)%volume
			else
				ex_conc(k) = 0
			endif
 			write(nflog,'(a,4i6,f8.3)') 'Cell volume: ',x,ns,i,k,ex_conc(k)
		endif
    enddo
enddo
! Add concentrations at the two boundaries 
! At ns=1, at at ns=ns+1
ns = ns+1
do ic = 1,nvars
	ichemo = ic - 1
	if (ichemo == 0) then	! CFSE
		k = ic	
		ex_conc(k) = 0
        k = (ns-1)*nvars + ic
        ex_conc(k) = 0	
	elseif (ichemo <= MAX_CHEMO) then
		k = ic
		if (chemo(ichemo)%used) then
			ex_conc(k) = BdryConc(ichemo,t_simulation)
		else
			ex_conc(k) = 0
		endif      
		k = (ns-1)*nvars + ic
		if (chemo(ichemo)%used) then
			ex_conc(k) = BdryConc(ichemo,t_simulation)
		else
			ex_conc(k) = 0
		endif      
    elseif (ichemo == MAX_CHEMO+1) then	! growth rate
		k = ic
		ex_conc(k) = 0
! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,1,i,k,ex_conc(k)
        k = (ns-1)*nvars + ic
        ex_conc(k) = 0
! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,ns,i,k,ex_conc(k)
    elseif (ichemo == MAX_CHEMO+2) then	! cell volume
		k = ic
		ex_conc(k) = 0
! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,1,i,k,ex_conc(k)
        k = (ns-1)*nvars + ic
        ex_conc(k) = 0
! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,ns,i,k,ex_conc(k)
    endif
enddo
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
subroutine get_oxyprob(nv, v0, dv, prob) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_oxyprob
use, intrinsic :: iso_c_binding
integer(c_int) :: nv
real(c_double) :: v0, dv, prob(*)
integer :: n, kcell, k
real(REAL_KIND) :: v, Vmin, Vmax, O2max

!call logger('get_oxyprob')
Vmin = 0
Vmax = chemo(OXYGEN)%bdry_conc
v0 = Vmin
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
	k = max(k,1)
	prob(k) = prob(k) + 1
	n = n+1
enddo
prob(1:nv) = prob(1:nv)/n
end subroutine

!--------------------------------------------------------------------------------
! Returns the distribution of intracellular constituent concentration
! nv is passed from the GUI
!--------------------------------------------------------------------------------
subroutine get_concprob(ichemo, nv, v0, dv, prob) !BIND(C)
!!DEC$ ATTRIBUTES DLLEXPORT :: get_concprob
use, intrinsic :: iso_c_binding
integer(c_int) :: ichemo, nv
real(c_double) :: v0, dv, prob(*)
integer :: n, kcell, k
real(REAL_KIND) :: v, Vmin, Vmax

!call logger('get_oxyprob')
Vmin = 0
!Vmax = chemo(ichemo)%bdry_conc
Vmax = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	v = cell_list(kcell)%conc(ichemo)
	Vmax = max(Vmax,v)
enddo
VMax = 1.05*Vmax;
v0 = Vmin
dv = (Vmax - Vmin)/nv
prob(1:nv) = 0
n = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	v = cell_list(kcell)%conc(ichemo)
	k = (v - Vmin)/dv + 1
	k = min(k,nv)
	k = max(k,1)
	prob(k) = prob(k) + 1
	n = n+1
enddo
prob(1:nv) = prob(1:nv)/n
end subroutine

!--------------------------------------------------------------------------------
!subroutine get_distdata(nv, ddata) BIND(C)
!!DEC$ ATTRIBUTES DLLEXPORT :: get_distdata
!integer(c_int) :: nv
!type(dist_data) :: ddata(0:*)
!integer :: ichemo
!
!do ichemo = 0,MAX_CHEMO+N_EXTRA
!	if (.not.ddata(ichemo).used) continue
!	if (ichemo  == 0) then	! CFSE
!	elseif (ichemo <= MAX_CHEMO) then
!	elseif (ichemo == GROWTH_RATE) then
!	elseif (ichemo == CELL_VOLUME) then
!	endif
!enddo
!
!end subroutine

!-----------------------------------------------------------------------------------------
! Get number of live cells
!-----------------------------------------------------------------------------------------
subroutine get_nFACS(n) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nfacs
use, intrinsic :: iso_c_binding
integer(c_int) :: n
integer :: k, kcell

n = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	n = n+1
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_FACS(facs_data) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_facs
use, intrinsic :: iso_c_binding
real(c_double) :: val, facs_data(*)
integer :: k, kcell, iextra, ichemo, ivar, nvars, var_index(32)
real(REAL_KIND) :: cfse_min

nvars = 1	! CFSE
var_index(nvars) = 0
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	nvars = nvars + 1
	var_index(nvars) = ichemo
enddo
do iextra = 1,N_EXTRA
	nvars = nvars + 1
	var_index(nvars) = MAX_CHEMO + iextra
enddo
cfse_min = 1.0e20
k = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
!	k = k+1
!	facs_data(k) = cell_list(kcell)%CFSE
!	k = k+1
!	facs_data(k) = cell_list(kcell)%dVdt
!	k = k+1
!	facs_data(k) = cell_list(kcell)%conc(OXYGEN)
!	if (cell_list(kcell)%conc(OXYGEN) <= 0.00001 .or. cell_list(kcell)%dVdt < 2.0e-6) then
!		write(nflog,'(2i6,2e12.3)') istep,kcell,cell_list(kcell)%dVdt,cell_list(kcell)%conc(OXYGEN)
!	endif
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cell_list(kcell)%CFSE
			cfse_min = min(val,cfse_min)
		elseif (ichemo <= MAX_CHEMO) then
			val = cell_list(kcell)%conc(ichemo)
		elseif (ichemo == GROWTH_RATE) then
			val = cell_list(kcell)%dVdt
		elseif (ichemo == CELL_VOLUME) then
			val = cell_list(kcell)%volume
		elseif (ichemo == O2_BY_VOL) then
			val = cell_list(kcell)%volume*cell_list(kcell)%conc(OXYGEN)
		endif
		k = k+1
		facs_data(k) = val
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! nhisto is the number of histogram boxes
! vmax(ivar) is the maximum value for variable ivar
! Probably need to adjust vmax() to a roundish value
!
! Compute 3 distributions: 1 = both cell types
!                          2 = type 1
!                          3 = type 2
! Stack three cases in vmax() and histo_data()
!-----------------------------------------------------------------------------------------
subroutine get_histo(nhisto, histo_data, vmin, vmax, histo_data_log, vmin_log, vmax_log) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_histo
use, intrinsic :: iso_c_binding
integer(c_int),value :: nhisto
real(c_double) :: vmin(*), vmax(*), histo_data(*)
real(c_double) :: vmin_log(*), vmax_log(*), histo_data_log(*)
real(REAL_KIND) :: val, val_log
integer :: n(3), i, ih, k, kcell, ict, ichemo, ivar, nvars, var_index(32)
integer,allocatable :: cnt(:,:,:)
real(REAL_KIND),allocatable :: dv(:,:), valmin(:,:), valmax(:,:)
integer,allocatable :: cnt_log(:,:,:)
real(REAL_KIND),allocatable :: dv_log(:,:), valmin_log(:,:), valmax_log(:,:)
!real(REAL_KIND) :: vmin_log(100), vmax_log(100)
!real(REAL_KIND),allocatable :: histo_data_log(:)

nvars = 1	! CFSE
var_index(nvars) = 0
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	nvars = nvars + 1
	var_index(nvars) = ichemo
enddo
nvars = nvars + 1
var_index(nvars) = GROWTH_RATE
nvars = nvars + 1
var_index(nvars) = CELL_VOLUME
nvars = nvars + 1
var_index(nvars) = O2_BY_VOL

allocate(cnt(3,nvars,nhisto))
allocate(dv(3,nvars))
allocate(valmin(3,nvars))
allocate(valmax(3,nvars))
allocate(cnt_log(3,nvars,nhisto))
allocate(dv_log(3,nvars))
allocate(valmin_log(3,nvars))
allocate(valmax_log(3,nvars))
!allocate(histo_data_log(10000))
cnt = 0
valmin = 1.0e10
valmax = -1.0e10
cnt_log = 0
valmin_log = 1.0e10
valmax_log = -1.0e10
n = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	ict = cell_list(kcell)%celltype
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cell_list(kcell)%CFSE
		elseif (ichemo <= MAX_CHEMO) then
			val = cell_list(kcell)%conc(ichemo)
		elseif (ichemo == GROWTH_RATE) then
			val = cell_list(kcell)%dVdt
		elseif (ichemo == CELL_VOLUME) then
			val = Vcell_pL*cell_list(kcell)%volume
		elseif (ichemo == O2_BY_VOL) then
			val = cell_list(kcell)%conc(OXYGEN)*Vcell_pL*cell_list(kcell)%volume
		endif
		valmax(ict+1,ivar) = max(valmax(ict+1,ivar),val)	! cell type 1 or 2
		valmax(1,ivar) = max(valmax(1,ivar),val)			! both
		if (val <= 1.0e-8) then
			val_log = -8
		else
			val_log = log10(val)
		endif
		valmin_log(ict+1,ivar) = min(valmin_log(ict+1,ivar),val_log)	! cell type 1 or 2
		valmin_log(1,ivar) = min(valmin_log(1,ivar),val_log)			! both
		valmax_log(ict+1,ivar) = max(valmax_log(ict+1,ivar),val_log)	! cell type 1 or 2
		valmax_log(1,ivar) = max(valmax_log(1,ivar),val_log)			! both
	enddo
	n(ict+1) = n(ict+1) + 1
	n(1) = n(1) + 1
enddo
do ivar = 1,nvars
	ichemo = var_index(ivar)
	if (ichemo == CELL_VOLUME) then
		valmin(:,ivar) = Vcell_pL*0.8
		valmin_log(:,ivar) = log10(Vcell_pL*0.8)
	else
		valmin(:,ivar) = 0
	endif
enddo

dv = (valmax - valmin)/nhisto
!write(nflog,*) 'dv'
!write(nflog,'(e12.3)') dv
dv_log = (valmax_log - valmin_log)/nhisto
!write(nflog,*) 'dv_log'
!write(nflog,'(e12.3)') dv_log
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	ict = cell_list(kcell)%celltype
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cell_list(kcell)%CFSE
		elseif (ichemo <= MAX_CHEMO) then
			val = cell_list(kcell)%conc(ichemo)
		elseif (ichemo == GROWTH_RATE) then
			val = cell_list(kcell)%dVdt
		elseif (ichemo == CELL_VOLUME) then
			val = Vcell_pL*cell_list(kcell)%volume
		elseif (ichemo == O2_BY_VOL) then
			val = cell_list(kcell)%conc(OXYGEN)*Vcell_pL*cell_list(kcell)%volume
		endif
		k = (val-valmin(1,ivar))/dv(1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt(1,ivar,k) = cnt(1,ivar,k) + 1
		k = (val-valmin(ict+1,ivar))/dv(ict+1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt(ict+1,ivar,k) = cnt(ict+1,ivar,k) + 1
		if (val <= 1.0e-8) then
			val_log = -8
		else
			val_log = log10(val)
		endif
		k = (val_log-valmin_log(1,ivar))/dv_log(1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt_log(1,ivar,k) = cnt_log(1,ivar,k) + 1
		k = (val_log-valmin_log(ict+1,ivar))/dv_log(ict+1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt_log(ict+1,ivar,k) = cnt_log(ict+1,ivar,k) + 1
	enddo
enddo

do i = 1,3
	if (n(i) == 0) then
		vmin((i-1)*nvars+1:i*nvars) = 0
		vmax((i-1)*nvars+1:i*nvars) = 0
		histo_data((i-1)*nvars*nhisto+1:i*nhisto*nvars) = 0
		vmin_log((i-1)*nvars+1:i*nvars) = 0
		vmax_log((i-1)*nvars+1:i*nvars) = 0
		histo_data_log((i-1)*nvars*nhisto+1:i*nhisto*nvars) = 0
	else
		do ivar = 1,nvars
			vmin((i-1)*nvars+ivar) = valmin(i,ivar)
			vmax((i-1)*nvars+ivar) = valmax(i,ivar)
			do ih = 1,nhisto
				k = (i-1)*nvars*nhisto + (ivar-1)*nhisto + ih
				histo_data(k) = (100.*cnt(i,ivar,ih))/n(i)
			enddo
			vmin_log((i-1)*nvars+ivar) = valmin_log(i,ivar)
			vmax_log((i-1)*nvars+ivar) = valmax_log(i,ivar)
			do ih = 1,nhisto
				k = (i-1)*nvars*nhisto + (ivar-1)*nhisto + ih
				histo_data_log(k) = (100.*cnt_log(i,ivar,ih))/n(i)
			enddo
		enddo
	endif
enddo
deallocate(cnt)
deallocate(dv)
deallocate(valmin)
deallocate(valmax)
deallocate(cnt_log)
deallocate(dv_log)
deallocate(valmin_log)
deallocate(valmax_log)
!deallocate(histo_data_log)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine WriteProfileData
integer :: ns
real(REAL_KIND) :: dx
real(REAL_KIND), allocatable :: ex_conc(:,:)
real(REAL_KIND) :: cbnd, cmin = 1.0e-6
integer :: rng(3,2), i, k, ichemo, kcell, x, y, z, ic, nc, kmax
character*(16) :: title(1:MAX_CHEMO+N_EXTRA)
character*(128) :: filename
character*(6) :: mintag

dx = DELTA_X
rng(:,1) = Centre(:) - (adrop*Radius + 2)
rng(:,2) = Centre(:) + (adrop*Radius + 2)
kmax = rng(1,2)-rng(1,1)+ 3
allocate(ex_conc(1:MAX_CHEMO+N_EXTRA,kmax))
y = Centre(2) + 0.5
z = Centre(3) + 0.5
ic = 0
do ichemo = 0,MAX_CHEMO+N_EXTRA
	if (ichemo == CFSE) then
		if (.not.chemo(ichemo)%used) cycle
		ic = ic + 1
		title(ic) = 'CFSE'
	    cbnd = 0
	elseif (ichemo <= MAX_CHEMO) then
		if (.not.chemo(ichemo)%used) cycle
		ic = ic + 1
		title(ic) = chemo(ichemo)%name
	    cbnd = BdryConc(ichemo,t_simulation)
	elseif (ichemo == GROWTH_RATE) then
		ic = ic + 1
		title(ic) = 'Growth_rate'
		cbnd = 0
	elseif (ichemo == CELL_VOLUME) then
		ic = ic + 1
		title(ic) = 'Cell_volume'
		cbnd = 0
	elseif (ichemo == O2_BY_VOL) then
		ic = ic + 1
		title(ic) = 'Cell_O2xVol'
		cbnd = 0
	endif
	k = 1
    ex_conc(ic,k) = cbnd
	do x = rng(1,1),rng(1,2)
		kcell = occupancy(x,y,z)%indx(1)
		if (kcell <= OUTSIDE_TAG) cycle
		i = ODEdiff%ivar(x,y,z)
		k = k+1
		if (ichemo == CFSE) then
			if (kcell > 0) then
				ex_conc(ic,k) = cell_list(kcell)%cfse
			else
				ex_conc(ic,k) = 0
			endif
		elseif (ichemo <= MAX_CHEMO) then
			if (i > 0) then
				ex_conc(ic,k) = allstate(i,ichemo)
			else
				ex_conc(ic,k) = 0
			endif
		elseif (ichemo == GROWTH_RATE) then
			if (kcell > 0) then
				ex_conc(ic,k) = cell_list(kcell)%dVdt
			else
				ex_conc(ic,k) = 0
			endif
		elseif (ichemo == CELL_VOLUME) then
			if (kcell > 0) then
				ex_conc(ic,k) = Vcell_pL*cell_list(kcell)%volume
			else
				ex_conc(ic,k) = 0
			endif
		elseif (ichemo == O2_BY_VOL) then
			if (kcell > 0) then
				ex_conc(ic,k) = allstate(i,OXYGEN)*Vcell_pL*cell_list(kcell)%volume
			else
				ex_conc(ic,k) = 0
			endif
		endif
	enddo
	k = k+1
    ex_conc(ic,k) = cbnd	! Add concentration at the boundary  
enddo
ns = k
nc = ic
! Remove time tag from the filename for download from NeSI
!write(mintag,'(i6)') int(istep*DELTA_T/60)
filename = profiledatafilebase
!filename = trim(filename)//'_'
!filename = trim(filename)//trim(adjustl(mintag))
!filename = trim(filename)//'min.dat'
open(nfprofile,file=filename,status='replace')
write(nfprofile,'(a,a)') 'GUI version: ',gui_run_version
write(nfprofile,'(a,a)') 'DLL version: ',dll_run_version
write(nfprofile,'(i6,a)') int(istep*DELTA_T/60),' minutes'
write(nfprofile,'(i6,a)') ns,' sites'
write(nfprofile,'(f6.2,a)') 1000*DELTA_X,' dx (um)'
write(nfprofile,'(32a16)') title(1:nc)
do k = 1,ns
	write(nfprofile,'(32(e12.3,4x))') ex_conc(1:nc,k)
enddo
close(nfprofile)
deallocate(ex_conc)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_constituents(nvars,cvar_index,nvarlen,name_array,narraylen) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_constituents
use, intrinsic :: iso_c_binding
character(c_char) :: name_array(0:*)
integer(c_int) :: nvars, cvar_index(0:*), nvarlen, narraylen
integer :: ivar, k, ichemo
character*(24) :: name
character(c_char) :: c

write(nflog,*) 'get_constituents'
nvarlen = 24
ivar = 0
k = ivar*nvarlen
cvar_index(ivar) = 0	! CFSE
name = 'CFSE'
call copyname(name,name_array(k),nvarlen)
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	ivar = ivar + 1
	k = ivar*nvarlen
	cvar_index(ivar) = ichemo
	name = chemo(ichemo)%name
	write(nflog,*) 'get_constituents: ',ichemo,name
	call copyname(name,name_array(k),nvarlen)
enddo
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = GROWTH_RATE
name = 'Growth rate'
call copyname(name,name_array(k),nvarlen)
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = CELL_VOLUME
name = 'Cell volume'
call copyname(name,name_array(k),nvarlen)
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = O2_BY_VOL
name = 'Cell O2xVol'
call copyname(name,name_array(k),nvarlen)
nvars = ivar + 1
write(nflog,*) 'did get_constituents'
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine copyname(name,name_array,n)
character*(*) :: name
character :: name_array(*)
integer :: n
integer :: k

do k = 1,n
	name_array(k) = name(k:k)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_DLL_build_version(version_array,array_len) BIND(C) 
!DEC$ ATTRIBUTES DLLEXPORT :: get_dll_build_version
use, intrinsic :: iso_c_binding
character(c_char) :: version_array(12)
integer(c_int) :: array_len
integer :: k

dll_version = DLL_BUILD_VERSION
gui_version = GUI_BUILD_VERSION
!write(nflog,*) 'get_DLL_build_version: ',dll_version
do k = 1,12
	version_array(k) = dll_version(k:k)
!	write(nflog,'(i2,a,a)') k,' ',version_array(k)
	if (version_array(k) == ' ') then
		version_array(k) = char(0)
		array_len = k
		exit
	endif
enddo
!write(nflog,*) 'array_len: ',array_len
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Execute(ncpu,infile_array,inbuflen,outfile_array,outbuflen) BIND(C) 
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char) :: infile_array(128), outfile_array(128)
integer(c_int) :: ncpu, inbuflen, outbuflen
character*(128) :: infile, outfile
logical :: ok, success
integer :: i, res

infile = ''
do i = 1,inbuflen
	infile(i:i) = infile_array(i)
enddo
outfile = ''
do i = 1,outbuflen
	outfile(i:i) = outfile_array(i)
enddo

open(nflog,file='spheroid.log',status='replace')

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
	call connecter(ok)
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

call Setup(ncpu,infile,outfile,ok)
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
execute_t1 = wtime()

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine DisableTCP
!DEC$ ATTRIBUTES DLLEXPORT :: disableTCP
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"DISABLETCP" :: disableTCP

use_TCP = .false.   ! because this is called from spheroid_main()	
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Connection(awp,port,error)
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
subroutine Connecter(ok)
logical :: ok
integer :: error

! Main connection
ok = .true.
error = 0
call Connection(awp_0,TCP_PORT_0,error)
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
call Wrapup

if (res == 0) then
	call logger(' Execution successful!')
elseif (res == -1) then
	call logger(' Execution stopped')
else
	call logger('  === Execution failed ===')
endif
write(logmsg,'(a,f10.2)') 'Execution time (min): ',(wtime() - execute_t1)/60
call logger(logmsg)

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
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Wrapup
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
