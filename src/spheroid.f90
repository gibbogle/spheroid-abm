module spheroid_mod
use global
use boundary
use chemokine
use cellstate
use react_diff
use winsock  
use deform
use drop
use colony
use envelope
use transfer

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
call logger("ReadCellParams new")
call ReadCellParams(ok)
if (.not.ok) return
call logger("did ReadCellParams")

start_wtime = wtime()

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

! Set up grid alignment
NY = NX
NZ = NX
x0 = (NX + 1.0)/2.
y0 = (NY + 1.0)/2.
z0 = (NZ + 1.0)/2.
blob_centre = [x0,y0,z0]   ! (units = grids)

DXB = 1.0e-4*DXB	! um -> cm
ixb0 = (1 + NXB)/2
iyb0 = (1 + NYB)/2
izb0 = 6
xb0 = (ixb0-1)*DXB
yb0 = (iyb0-1)*DXB 
zb0 = (izb0-1)*DXB
centre_b = [xb0, yb0, zb0]
dxb3 = dxb*dxb*dxb

write(nflog,'(a,3e12.3)') 'x0,y0,z0: ',x0,y0,z0
write(nflog,'(a,3e12.3)') 'xb0,yb0,zb0: ',xb0,yb0,zb0
grid_offset(1) = (ixb0-1)*DXB - ((NX+1)/2)*DELTA_X
grid_offset(2) = (iyb0-1)*DXB - ((NY+1)/2)*DELTA_X
grid_offset(3) = (izb0-1)*DXB - ((NZ+1)/2)*DELTA_X
write(nflog,'(a,3e12.3)') 'grid_offset: ',grid_offset
write(nflog,'(a,3e12.3)') 'blob_centre: ',blob_centre*DELTA_X + grid_offset

call ArrayInitialisation(ok)
if (.not.ok) return
call logger('did ArrayInitialisation')

call make_lattice_grid_weights

call SetupChemo

is_dropped = .false.
adrop = 1
bdrop = 1
cdrop = 0
zmin = 1
call PlaceCells(ok)
!call show_volume_data
!call SetRadius(Nsites)
call getVolume(blob_volume,blob_area)
blob_radius = sqrt(blob_area/PI)
blob_centre = getCentre()
write(logmsg,*) 'did PlaceCells: Ncells: ',Ncells,blob_radius
call logger(logmsg)
if (.not.ok) return

call CreateBdryList

if (use_FD) then
	call setup_react_diff(ok)
	if (.not.ok) return
endif
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
Ndrug_tag = 0
Nanoxia_tag = 0
Naglucosia_tag = 0
Nradiation_dead = 0
Ndrug_dead = 0
Nanoxia_dead = 0
Naglucosia_dead = 0
!radiation_dosed = .false.
t_simulation = 0
total_dMdt = 0
chemo(:)%total_flux_prev = 0
t_lastmediumchange = 0
limit_stop = .false.
Vex_min = 1.0
Vex_max = 0
kcell_dbug = 0
medium_change_step = .false.
write(logmsg,'(a,i6)') 'Startup procedures have been executed: initial T cell count: ',Ncells0
call logger(logmsg)

end subroutine

!----------------------------------------------------------------------------------------- 
!----------------------------------------------------------------------------------------- 
subroutine show_volume_data
integer :: kcell
real(REAL_KIND) :: Vsum, Vdivsum

write(nfout,*) 'Volume data:'
write(nfout,'(a,L)') 'use_divide_time_distribution: ',use_divide_time_distribution
write(nfout,'(a,L)') 'use_V_dependence: ',use_V_dependence
Vsum = 0
Vdivsum = 0
do kcell = 1,nlist
	write(nfout,'(i6,2f6.2)') kcell,cell_list(kcell)%volume, cell_list(kcell)%divide_volume
	Vsum = Vsum + cell_list(kcell)%volume
	Vdivsum = Vdivsum + cell_list(kcell)%divide_volume
enddo
write(nfout,*)
write(nfout,'(a,f6.2)') 'Average initial volume: ',Vsum/nlist
write(nfout,'(a,f6.2)') 'Average divide volume: ',Vdivsum/nlist
stop
end subroutine

!----------------------------------------------------------------------------------------- 
! Initialise medium concentrations, etc.
!-----------------------------------------------------------------------------------------
subroutine SetupMedium(ichemo)
integer :: ichemo
real(REAL_KIND) :: V, V0, R1

if (chemo(ichemo)%present) then
!	call SetRadius(Nsites)
	R1 = blob_radius*DELTA_X			! cm
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
chemo(OXYGEN)%medium_Cbnd_prev = chemo(OXYGEN)%bdry_conc
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
if (allocated(gaplist)) deallocate(gaplist)
if (allocated(Cslice)) deallocate(Cslice)
call logger('did deallocation')

!nsteps_per_min = 1.0/DELTA_T
ngaps = 0
nlist = 0

allocate(zoffset(0:2*Mnodes))
allocate(zdomain(NZ))
!x0 = (NX + 1.0)/2.        ! global value
!y0 = (NY + 1.0)/2.
!z0 = (NZ + 1.0)/2.
!blob_centre = [x0,y0,z0]   ! now, actually the global centre (units = grids)
call SetRadius(initial_count)
write(logmsg,*) 'Initial radius, count, max_nlist: ',blob_radius, initial_count, max_nlist
call logger(logmsg)

allocate(cell_list(max_nlist))
allocate(occupancy(NX,NY,NZ))
allocate(gaplist(max_ngaps))
!allocate(Cslice(NX/2,NY/2,NZ/2,MAX_CHEMO))
allocate(Cslice(NX,NY,NZ,MAX_CHEMO))

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

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine ReadCellParams(ok)
logical :: ok
integer :: i, idrug, imetab, nmetab, im, itestcase, Nmm3, ichemo, itreatment, iuse_extra, iuse_relax, iuse_par_relax, iuse_FD
integer :: iuse_oxygen, iuse_glucose, iuse_tracer, iuse_drug, iuse_metab, iV_depend, iV_random, iuse_gd_all
!integer ::  idrug_decay, imetab_decay
integer :: ictype, idisplay, isconstant, ioxygengrowth, iglucosegrowth, ioxygendeath, iglucosedeath
integer :: iuse_drop, iconstant, isaveprofiledata, isaveslicedata, isaveFACSdata
logical :: use_metabolites
real(REAL_KIND) :: days, bdry_conc, percent, d_n_limit
real(REAL_KIND) :: sigma(2), DXmm, anoxia_tag_hours, anoxia_death_hours, aglucosia_tag_hours, aglucosia_death_hours
character*(12) :: drug_name
character*(1) :: numstr

ok = .true.
chemo(:)%used = .false.

open(nfcell,file=inputfile,status='old')
read(nfcell,'(a)') header
if (header(1:3) == 'GUI') then
	gui_run_version = header
	header = 'DD/MM/YYYY header_string'
else
	read(nfcell,*) gui_run_version				! program run version number
endif
read(nfcell,*) dll_run_version				! DLL run version number 
read(nfcell,*) NX							! size of grid
read(nfcell,*) initial_count				! initial number of tumour cells
read(nfcell,*) divide_time_median(1)
read(nfcell,*) divide_time_shape(1)
read(nfcell,*) divide_time_median(2)
read(nfcell,*) divide_time_shape(2)
read(nfcell,*) iV_depend
read(nfcell,*) iV_random
read(nfcell,*) days							! number of days to simulate
read(nfcell,*) d_n_limit					! possible limit on diameter or number of cells
diam_count_limit = d_n_limit
read(nfcell,*) DELTA_T						! time step size (sec)
read(nfcell,*) NXB							! size of coarse grid = NXB = NYB
read(nfcell,*) NZB							! size of coarse grid = NZB
read(nfcell,*) DXF							! fine grid spacing, off-lattice model (um)
read(nfcell,*) a_separation
read(nfcell,*) a_force
read(nfcell,*) c_force
read(nfcell,*) x0_force
read(nfcell,*) x1_force
read(nfcell,*) kdrag
read(nfcell,*) frandom
read(nfcell,*) NT_CONC						! number of subdivisions of DELTA_T for diffusion computation
read(nfcell,*) Nmm3							! number of cells/mm^3
DXmm = 1.0/(Nmm3**(1./3))
DELTA_X = DXmm/10							! mm -> cm
Vsite_cm3 = DELTA_X*DELTA_X*DELTA_X			! total site volume (cm^3)
read(nfcell,*) fluid_fraction				! fraction of the (non-necrotic) tumour that is fluid
read(nfcell,*) medium_volume0				! initial total volume (medium + spheroid) (cm^3)
read(nfcell,*) d_layer						! thickness of the unstirred layer around the spheroid (cm)
read(nfcell,*) Vdivide0						! nominal cell volume multiple for division
read(nfcell,*) dVdivide						! variation about nominal divide volume
read(nfcell,*) MM_THRESHOLD					! O2 concentration threshold Michaelis-Menten "soft-landing" (uM)
read(nfcell,*) anoxia_threshold			    ! O2 threshold for anoxia (uM)
read(nfcell,*) anoxia_tag_hours				! hypoxic time leading to tagging to die by anoxia (h)
read(nfcell,*) anoxia_death_hours			! time after tagging to death by anoxia (h)
read(nfcell,*) aglucosia_threshold			! O2 threshold for aglucosia (uM)
read(nfcell,*) aglucosia_tag_hours			! hypoxic time leading to tagging to die by aglucosia (h)
read(nfcell,*) aglucosia_death_hours		! time after tagging to death by aglucosia (h)
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
read(nfcell,*) ioxygengrowth
chemo(OXYGEN)%controls_growth = (ioxygengrowth == 1)
read(nfcell,*) ioxygendeath
chemo(OXYGEN)%controls_death = (ioxygendeath == 1)
read(nfcell,*) chemo(OXYGEN)%diff_coef
read(nfcell,*) chemo(OXYGEN)%medium_diff_coef
read(nfcell,*) chemo(OXYGEN)%membrane_diff_in
chemo(OXYGEN)%membrane_diff_in = chemo(OXYGEN)%membrane_diff_in*Vsite_cm3/60		! /min -> /sec
chemo(OXYGEN)%membrane_diff_out = chemo(OXYGEN)%membrane_diff_in
read(nfcell,*) chemo(OXYGEN)%bdry_conc
read(nfcell,*) iconstant
chemo(OXYGEN)%constant = (iconstant == 1)
read(nfcell,*) chemo(OXYGEN)%max_cell_rate
chemo(OXYGEN)%max_cell_rate = chemo(OXYGEN)%max_cell_rate*1.0e6					! mol/cell/s -> mumol/cell/s
read(nfcell,*) chemo(OXYGEN)%MM_C0
read(nfcell,*) chemo(OXYGEN)%Hill_N
read(nfcell,*) iuse_glucose		!chemo(GLUCOSE)%used
read(nfcell,*) iglucosegrowth
chemo(GLUCOSE)%controls_growth = (iglucosegrowth == 1)
read(nfcell,*) iglucosedeath
chemo(GLUCOSE)%controls_death = (iglucosedeath == 1)
read(nfcell,*) chemo(GLUCOSE)%diff_coef
read(nfcell,*) chemo(GLUCOSE)%medium_diff_coef
read(nfcell,*) chemo(GLUCOSE)%membrane_diff_in
chemo(GLUCOSE)%membrane_diff_in = chemo(GLUCOSE)%membrane_diff_in*Vsite_cm3/60	! /min -> /sec
chemo(GLUCOSE)%membrane_diff_out = chemo(GLUCOSE)%membrane_diff_in
read(nfcell,*) chemo(GLUCOSE)%bdry_conc
read(nfcell,*) iconstant
chemo(GLUCOSE)%constant = (iconstant == 1)
read(nfcell,*) chemo(GLUCOSE)%max_cell_rate
chemo(GLUCOSE)%max_cell_rate = chemo(GLUCOSE)%max_cell_rate*1.0e6					! mol/cell/s -> mumol/cell/s
read(nfcell,*) chemo(GLUCOSE)%MM_C0
read(nfcell,*) chemo(GLUCOSE)%Hill_N
read(nfcell,*) iuse_tracer		!chemo(TRACER)%used
read(nfcell,*) chemo(TRACER)%diff_coef
read(nfcell,*) chemo(TRACER)%medium_diff_coef
read(nfcell,*) chemo(TRACER)%membrane_diff_in
chemo(TRACER)%membrane_diff_in = chemo(TRACER)%membrane_diff_in*Vsite_cm3/60		! /min -> /sec
chemo(TRACER)%membrane_diff_out = chemo(TRACER)%membrane_diff_in
read(nfcell,*) chemo(TRACER)%bdry_conc
read(nfcell,*) iconstant
chemo(TRACER)%constant = (iconstant == 1)
read(nfcell,*) chemo(TRACER)%max_cell_rate
read(nfcell,*) chemo(TRACER)%MM_C0
read(nfcell,*) chemo(TRACER)%Hill_N

! removed old read of TPZ and DNB drug data

read(nfcell,*) LQ(1)%alpha_H
read(nfcell,*) LQ(1)%beta_H
read(nfcell,*) LQ(1)%OER_am
read(nfcell,*) LQ(1)%OER_bm
read(nfcell,*) LQ(1)%K_ms
read(nfcell,*) LQ(1)%death_prob
read(nfcell,*) LQ(1)%growth_delay_factor
read(nfcell,*) LQ(1)%growth_delay_N
read(nfcell,*) LQ(2)%alpha_H
read(nfcell,*) LQ(2)%beta_H
read(nfcell,*) LQ(2)%OER_am
read(nfcell,*) LQ(2)%OER_bm
read(nfcell,*) LQ(2)%K_ms
read(nfcell,*) LQ(2)%death_prob
read(nfcell,*) LQ(2)%growth_delay_factor
read(nfcell,*) LQ(2)%growth_delay_N
read(nfcell,*) iuse_gd_all
use_radiation_growth_delay_all = (iuse_gd_all == 1)
read(nfcell,*) O2cutoff(1)
read(nfcell,*) O2cutoff(2)
read(nfcell,*) O2cutoff(3)
read(nfcell,*) hypoxia_threshold
read(nfcell,*) growthcutoff(1)
read(nfcell,*) growthcutoff(2)
read(nfcell,*) growthcutoff(3)
read(nfcell,*) Cthreshold
read(nfcell,*) spcrad_value
read(nfcell,*) iuse_extra
read(nfcell,*) iuse_relax
read(nfcell,*) iuse_par_relax
read(nfcell,*) iuse_FD
read(nfcell,*) iuse_drop
read(nfcell,*) Ndrop
read(nfcell,*) alpha_shape
read(nfcell,*) beta_shape
read(nfcell,*) isaveprofiledata
read(nfcell,*) saveprofile%filebase
read(nfcell,*) saveprofile%dt
read(nfcell,*) saveprofile%nt
read(nfcell,*) isaveslicedata
read(nfcell,*) saveslice%filebase
read(nfcell,*) saveslice%dt
read(nfcell,*) saveslice%nt
read(nfcell,*) isaveFACSdata
read(nfcell,*) saveFACS%filebase
read(nfcell,*) saveFACS%dt
read(nfcell,*) saveFACS%nt

read(nfcell,*) Ndrugs_used
call ReadDrugData(nfcell)

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
DXB = 4*DXF
MM_THRESHOLD = MM_THRESHOLD/1000					! uM -> mM
anoxia_threshold = anoxia_threshold/1000			! uM -> mM
aglucosia_threshold = aglucosia_threshold/1000		! uM -> mM
O2cutoff = O2cutoff/1000							! uM -> mM
hypoxia_threshold = hypoxia_threshold/1000			! uM -> mM
relax = (iuse_relax == 1)
use_parallel = (iuse_par_relax == 1)
use_FD = (iuse_FD == 1)
chemo(OXYGEN)%used = (iuse_oxygen == 1)
chemo(GLUCOSE)%used = (iuse_glucose == 1)
chemo(TRACER)%used = (iuse_tracer == 1)
chemo(OXYGEN)%MM_C0 = chemo(OXYGEN)%MM_C0/1000		! uM -> mM
chemo(GLUCOSE)%MM_C0 = chemo(GLUCOSE)%MM_C0/1000	! uM -> mM
if (.not.chemo(OXYGEN)%used) then
    chemo(OXYGEN)%controls_growth = .false.
    chemo(OXYGEN)%controls_death = .false.
endif
if (.not.chemo(GLUCOSE)%used) then
    chemo(GLUCOSE)%controls_growth = .false.
    chemo(GLUCOSE)%controls_death = .false.
endif
LQ(:)%growth_delay_factor = 60*60*LQ(:)%growth_delay_factor	! hours -> seconds
divide_dist(1:2)%class = LOGNORMAL_DIST
divide_time_median(1:2) = 60*60*divide_time_median(1:2)			! hours -> seconds
sigma(1:2) = log(divide_time_shape(1:2))
!divide_dist%p1 = log(divide_time_mean/exp(sigma*sigma/2))	
divide_dist(1:2)%p1 = log(divide_time_median(1:2))	
divide_dist(1:2)%p2 = sigma
divide_time_mean(1:2) = exp(divide_dist(1:2)%p1 + 0.5*divide_dist(1:2)%p2**2)	! mean = median.exp(sigma^2/2)
write(logmsg,'(a,24e12.4)') 'shape, sigma: ',divide_time_shape(1:2),sigma(1:2)
call logger(logmsg)
write(logmsg,'(a,4e12.4)') 'Median, mean divide time: ',divide_time_median(1:2)/3600,divide_time_mean(1:2)/3600
call logger(logmsg)

use_V_dependence = (iV_depend == 1)
randomise_initial_volume = (iV_random == 1)
use_extracellular_O2 = (iuse_extra == 1)
t_anoxia_limit = 60*60*anoxia_tag_hours				! hours -> seconds
anoxia_death_delay = 60*60*anoxia_death_hours		! hours -> seconds
t_aglucosia_limit = 60*60*aglucosia_tag_hours		! hours -> seconds
aglucosia_death_delay = 60*60*aglucosia_death_hours	! hours -> seconds
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

saveprofile%active = (isaveprofiledata == 1)
saveprofile%it = 1
saveprofile%dt = 60*saveprofile%dt		! mins -> seconds
saveslice%active = (isaveslicedata == 1)
saveslice%it = 1
saveslice%dt = 60*saveslice%dt			! mins -> seconds
saveFACS%active = (isaveFACSdata == 1)
saveFACS%it = 1
saveFACS%dt = 60*saveFACS%dt			! mins -> seconds

use_dropper = (iuse_drop == 1)

! Setup test_case
test_case = .false.
if (itestcase /= 0) then
    test_case(itestcase) = .true.
endif

if (mod(NX,2) /= 0) NX = NX+1					! ensure that NX is even
NYB = NXB

open(nfout,file=outputfile,status='replace')
write(nfout,'(a,a)') 'GUI version: ',gui_run_version
write(nfout,'(a,a)') 'DLL version: ',dll_run_version
write(nfout,*)

write(nflog,*)
write(nflog,'(a,a)') 'GUI version: ',gui_run_version
write(nflog,'(a,a)') 'DLL version: ',dll_run_version
write(nflog,*)

open(nfres,file='spheroid_ts.out',status='replace')
!write(nfres,'(a,a)') 'GUI version: ',gui_run_version
!write(nfres,'(a,a)') 'DLL version: ',dll_run_version
!write(nfres,*)
write(nfres,'(a)') 'date info GUI_version DLL_version &
istep hour vol_mm3 diam_um Ncells(1) Ncells(2) &
Nanoxia_dead(1) Nanoxia_dead(2) Naglucosia_dead(1) Naglucosia_dead(2) NdrugA_dead(1) NdrugA_dead(2) &
NdrugB_dead(1) NdrugB_dead(2) Nradiation_dead(1) Nradiation_dead(2) &
Ntagged_anoxia(1) Ntagged_anoxia(2) Ntagged_aglucosia(1) Ntagged_aglucosia(2) Ntagged_drugA(1) Ntagged_drugA(2) &
Ntagged_drugB(1) Ntagged_drugB(2) Ntagged_radiation(1) Ntagged_radiation(2) &
f_hypox(1) f_hypox(2) f_hypox(3) &
f_clonohypox(1) f_clonohypox(2) f_clonohypox(3) &
f_growth(1) f_growth(2) f_growth(3) &
f_necrot plating_efficiency(1) plating_efficiency(2) &
medium_oxygen medium_glucose medium_drugA medium_drugB &
bdry_oxygen bdry_glucose bdry_drugA bdry_drugB'

write(logmsg,*) 'Opened nfout: ',outputfile
call logger(logmsg)

Nsteps = days*24*60*60/DELTA_T		! DELTA_T in seconds
write(logmsg,'(a,2i6,f6.0)') 'nsteps, NT_CONC, DELTA_T: ',nsteps,NT_CONC,DELTA_T
call logger(logmsg)

call DetermineKd
ok = .true.

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ReadDrugData(nf)
integer :: nf
integer :: idrug, im, ictyp, ival
character*(16) :: drugname
type(drug_type), pointer :: dp

write(logmsg,*) 'ReadDrugData'
call logger(logmsg)
if (allocated(drug)) then
	deallocate(drug)
endif
allocate(drug(Ndrugs_used))
do idrug = 1,Ndrugs_used
	dp => drug(idrug)
	read(nf,'(a)') drug(idrug)%classname
	if (drug(idrug)%classname == 'TPZ') then
		drug(idrug)%drugclass = TPZ_CLASS
	elseif (drug(idrug)%classname == 'DNB') then
		drug(idrug)%drugclass = DNB_CLASS
	endif
	drug(idrug)%nmetabolites = 2			! currently all drugs have 2 metabolites
	drug(idrug)%use_metabolites = .true.	! currently simulate metabolites
    do im = 0,2			! 0 = parent, 1 = metab_1, 2 = metab_2
		read(nf,'(a)') drugname
		if (im == 0) then
			drug(idrug)%name = drugname
		endif
		read(nf,*) drug(idrug)%diff_coef(im)
		read(nf,*) drug(idrug)%medium_diff_coef(im)
		read(nf,*) drug(idrug)%membrane_diff_in(im)
		read(nf,*) drug(idrug)%membrane_diff_out(im)
		read(nf,*) drug(idrug)%halflife(im)
		drug(idrug)%membrane_diff_in(im) = drug(idrug)%membrane_diff_in(im)*Vsite_cm3/60	! /min -> /sec
		drug(idrug)%membrane_diff_out(im) = drug(idrug)%membrane_diff_out(im)*Vsite_cm3/60	! /min -> /sec
		do ictyp = 1,ncelltypes
            read(nf,*) drug(idrug)%Kmet0(ictyp,im)
            read(nf,*) drug(idrug)%C2(ictyp,im)
            read(nf,*) drug(idrug)%KO2(ictyp,im)
            read(nf,*) drug(idrug)%Vmax(ictyp,im)
            read(nf,*) drug(idrug)%Km(ictyp,im)
            read(nf,*) drug(idrug)%Klesion(ictyp,im)
            read(nf,*) drug(idrug)%kill_O2(ictyp,im)
            read(nf,*) drug(idrug)%kill_drug(ictyp,im)
            read(nf,*) drug(idrug)%kill_duration(ictyp,im)
            read(nf,*) drug(idrug)%kill_fraction(ictyp,im)
            read(nf,*) drug(idrug)%SER_max(ictyp,im)
            read(nf,*) drug(idrug)%SER_Km(ictyp,im)
            read(nf,*) drug(idrug)%SER_KO2(ictyp,im)
            read(nf,*) drug(idrug)%n_O2(ictyp,im)
            read(nf,*) drug(idrug)%death_prob(ictyp,im)
            read(nf,*) ival
            drug(idrug)%kills(ictyp,im) = (ival == 1)
            read(nf,*) ival
            drug(idrug)%kill_model(ictyp,im) = ival
            read(nf,*) ival
            drug(idrug)%sensitises(ictyp,im) = (ival == 1)
			write(nflog,'(a,i2,7e12.3)') 'ReadDrug: im,C2,KO2,Ckill_O2,Kmet0,f,T: ',im, &
				dp%C2(ictyp,im),dp%KO2(ictyp,im),dp%kill_drug(ictyp,im),dp%Kmet0(ictyp,im),dp%kill_fraction(ictyp,im),dp%kill_duration(ictyp,im)
            drug(idrug)%Kmet0(ictyp,im) = drug(idrug)%Kmet0(ictyp,im)/60					! /min -> /sec
            drug(idrug)%KO2(ictyp,im) = 1.0e-3*drug(idrug)%KO2(ictyp,im)					! um -> mM
            drug(idrug)%kill_duration(ictyp,im) = 60*drug(idrug)%kill_duration(ictyp,im)	! min -> sec
		enddo
    enddo
    write(nflog,*) 'drug: ',idrug,drug(idrug)%classname,'  ',drug(idrug)%name
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Skip lines until the 'PROTOCOL' line
!-----------------------------------------------------------------------------------------
subroutine ReadProtocol(nf)
integer :: nf
integer :: itime, ntimes, kevent, ichemo, idrug, im
character*(64) :: line
character*(16) :: drugname
character*(1)  :: numstr
real(REAL_KIND) :: t, dt, vol, conc, O2conc, O2flush, dose, O2medium
type(event_type) :: E

write(logmsg,*) 'ReadProtocol'
call logger(logmsg)
chemo(TRACER+1:)%used = .false.
do
	read(nf,'(a)') line
	if (trim(line) == 'PROTOCOL') exit
enddo
read(nf,*) ntimes
if (ntimes == 0) then
	Nevents = 0
	return
endif
Nevents = ntimes
if (allocated(event)) deallocate(event)
allocate(event(2*ntimes))
kevent = 0
do itime = 1,ntimes
	read(nf,'(a)') line
	if (trim(line) == 'DRUG') then
		kevent = kevent + 1
		event(kevent)%etype = DRUG_EVENT
		read(nf,'(a)') line
		drugname = trim(line)
		do idrug = 1,ndrugs_used
			if (drugname == drug(idrug)%name) then
				ichemo = 4 + 3*(idrug-1)
				exit
			endif
		enddo
		! Need to copy drug(idrug) parameters to chemo(ichemo) 
		call CopyDrugParameters(idrug,ichemo)
		read(nf,*) t
		read(nf,*) dt
		read(nf,*) vol
		read(nf,*) O2conc
		read(nf,*) O2flush
		read(nf,*) conc
		event(kevent)%time = t
		event(kevent)%ichemo = ichemo
		event(kevent)%idrug = idrug
		event(kevent)%volume = vol
		event(kevent)%conc = conc
		event(kevent)%O2conc = O2conc
		event(kevent)%dose = 0
		chemo(ichemo)%used = .true.
		write(nflog,'(a,i3,2f8.3)') 'define DRUG_EVENT: volume, O2conc: ',kevent,event(kevent)%volume,event(kevent)%O2conc
		if (drug(idrug)%use_metabolites) then
			do im = 1,drug(idrug)%nmetabolites
				chemo(ichemo+im)%used = .true.
			enddo
		endif

		kevent = kevent + 1
		event(kevent)%etype = MEDIUM_EVENT
		event(kevent)%time = t + dt
		event(kevent)%ichemo = 0
		event(kevent)%volume = medium_volume0
		event(kevent)%conc = 0
		event(kevent)%O2medium = O2flush
		event(kevent)%dose = 0
		write(nflog,'(a,i3,2f8.3)') 'define MEDIUM_EVENT: volume: ',kevent,event(kevent)%volume,event(kevent)%O2medium
	elseif (trim(line) == 'MEDIUM') then
		kevent = kevent + 1
		event(kevent)%etype = MEDIUM_EVENT
		read(nf,*) t
		read(nf,*) vol
		read(nf,*) O2medium
		event(kevent)%time = t
		event(kevent)%volume = vol	
		event(kevent)%ichemo = 0
		event(kevent)%O2medium = O2medium
		event(kevent)%dose = 0
		write(nflog,'(a,i3,2f8.3)') 'define MEDIUM_EVENT: volume: ',kevent,event(kevent)%volume,event(kevent)%O2medium
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
! convert volume from uL to cm^3  NO LONGER - now using cm^3 everywhere
write(logmsg,*) 'nevents: ',nevents
call logger(logmsg)
do kevent = 1,Nevents
	event(kevent)%done = .false.
	event(kevent)%time = event(kevent)%time*60*60
!	event(kevent)%volume = event(kevent)%volume*1.0e-3
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
!-----------------------------------------------------------------------------------------
subroutine CopyDrugParameters(idrug,ichemo)
integer :: idrug,ichemo
integer :: im, im1, im2
character*(1) :: numstr

im1 = 0
chemo(ichemo)%name = drug(idrug)%name
if (drug(idrug)%use_metabolites) then
	do im = 1,drug(idrug)%nmetabolites
		chemo(ichemo+im)%used = .true.
		chemo(ichemo+im)%name = trim(chemo(ichemo)%name) // '_metab'
		write(numstr,'(i1)') im
		chemo(ichemo+im)%name = trim(chemo(ichemo+im)%name) // numstr
	enddo
	im2 = 2
else
	im2 = 0
endif
do im = im1, im2
	chemo(ichemo+im)%diff_coef = drug(idrug)%diff_coef(im)
	chemo(ichemo+im)%medium_diff_coef = drug(idrug)%medium_diff_coef(im)
	chemo(ichemo+im)%membrane_diff_in = drug(idrug)%membrane_diff_in(im)
	chemo(ichemo+im)%membrane_diff_out = drug(idrug)%membrane_diff_out(im)
	chemo(ichemo+im)%halflife = drug(idrug)%halflife(im)
	chemo(ichemo+im)%medium_dlayer = d_layer
	chemo(ichemo+im)%decay = (chemo(ichemo+im)%halflife > 0)
	if (chemo(ichemo+im)%decay) then
		chemo(ichemo+im)%decay_rate = DecayRate(chemo(ichemo+im)%halflife)
	else
		chemo(ichemo+im)%decay_rate = 0
	endif
enddo
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
real(REAL_KIND) :: C2, KO2, n_O2, Kmet0, kmet 
real(REAL_KIND) :: f, T, Ckill, Ckill_O2, Kd
integer :: idrug, ictyp, im, kill_model

do idrug = 1,ndrugs_used
	do ictyp = 1,Ncelltypes
		do im = 0,2
			if (drug(idrug)%kills(ictyp,im)) then
				C2 = drug(idrug)%C2(ictyp,im)
				KO2 = drug(idrug)%KO2(ictyp,im)
				n_O2 = drug(idrug)%n_O2(ictyp,im)
				Kmet0 = drug(idrug)%Kmet0(ictyp,im)
				kill_model = drug(idrug)%kill_model(ictyp,im)
				Ckill_O2 = drug(idrug)%kill_O2(ictyp,im)
				f = drug(idrug)%kill_fraction(ictyp,im)
				T = drug(idrug)%kill_duration(ictyp,im)
				Ckill = drug(idrug)%kill_drug(ictyp,im)
				kmet = (1 - C2 + C2*(KO2**n_O2)/(KO2**n_O2 + Ckill_O2**n_O2))*Kmet0
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
				write(nflog,'(a,i2,8e12.3)') 'DetermineKd: im,C2,KO2,Ckill,Kmet0,f,T,kmet,Kd: ',im,C2,KO2,Ckill,Kmet0,f,T,kmet,Kd
				drug(idrug)%Kd(ictyp,im) = Kd
			endif
		enddo
	enddo
enddo

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine PlaceCells(ok)
logical :: ok
integer :: x, y, z, kcell, site(3), ichemo, ityp
real(REAL_KIND) :: r2lim,r2,rad(3)

occupancy(:,:,:)%indx(1) = 0
occupancy(:,:,:)%indx(2) = 0
r2lim = 0.95*blob_radius*blob_radius
lastID = 0
kcell = 0
Ncells_type = 0
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
integer :: ityp, kpar = 0
real(REAL_KIND) :: V0, Tdiv, R

lastID = lastID + 1
cell_list(k)%ID = lastID
cell_list(k)%celltype = random_choice(celltype_fraction,Ncelltypes,kpar)
ityp = cell_list(k)%celltype
Ncells_type(ityp) = Ncells_type(ityp) + 1
cell_list(k)%site = site
cell_list(k)%state = 1
cell_list(k)%generation = 1
!cell_list(k)%drugA_tag = .false.
!cell_list(k)%drugB_tag = .false.
cell_list(k)%drug_tag = .false.
cell_list(k)%radiation_tag = .false.
cell_list(k)%anoxia_tag = .false.
cell_list(k)%aglucosia_tag = .false.
cell_list(k)%exists = .true.
cell_list(k)%active = .true.
cell_list(k)%growth_delay = .false.
cell_list(k)%G2_M = .false.
cell_list(k)%p_rad_death = 0
!R = par_uni(kpar)
!cell_list(k)%divide_volume = Vdivide0 + dVdivide*(2*R-1)
V0 = Vdivide0/2
cell_list(k)%divide_volume = get_divide_volume(ityp,V0, Tdiv)
cell_list(k)%divide_time = Tdiv
R = par_uni(kpar)
if (randomise_initial_volume) then
	cell_list(k)%volume = cell_list(k)%divide_volume*0.5*(1 + R)
else
	cell_list(k)%volume = 1.0
endif
!write(nflog,'(a,i6,f6.2)') 'volume: ',k,cell_list(k)%volume
cell_list(k)%t_divide_last = 0		! used in colony growth
cell_list(k)%t_anoxia = 0
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
				r2 = (nbsite(1) - blob_centre(1))**2 + (nbsite(2) - blob_centre(2))**2 + (nbsite(3) - blob_centre(3))**2
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
integer :: kevent, ichemo, idrug, im, nmetab
real(REAL_KIND) :: V, C(MAX_CHEMO)
type(event_type) :: E

!write(logmsg,*) 'ProcessEvent'
!call logger(logmsg)
do kevent = 1,Nevents
	E = event(kevent)
	if (t_simulation >= E%time .and. .not.E%done) then
		write(nflog,'(a,i3,2f8.0,i3,2f10.4)') 'Event: ',E%etype,t_simulation,E%time,E%ichemo,E%volume,E%conc
		if (E%etype == RADIATION_EVENT) then
			radiation_dose = E%dose
			write(logmsg,'(a,f8.0,f8.3)') 'RADIATION_EVENT: time, dose: ',t_simulation,E%dose
			call logger(logmsg)
		elseif (E%etype == MEDIUM_EVENT) then
			write(logmsg,'(a,f8.0,f8.3,2f8.4)') 'MEDIUM_EVENT: time, volume, O2medium: ',t_simulation,E%volume,E%O2medium
			call logger(logmsg)
			C = 0
			C(OXYGEN) = E%O2medium
			C(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
			V = E%volume
			call MediumChange(V,C)
		elseif (E%etype == DRUG_EVENT) then
			C = 0
			C(OXYGEN) = E%O2conc
			C(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
			ichemo = E%ichemo
			idrug = E%idrug
			C(ichemo) = E%conc
			V = E%volume
			write(logmsg,'(a,2f8.3)') 'DRUG_EVENT: volume, conc: ',E%volume,E%conc
			call logger(logmsg)
			! set %present
			chemo(ichemo)%present = .true.
			chemo(ichemo)%bdry_conc = 0
			nmetab = drug(idrug)%nmetabolites
			do im = 1,nmetab
				if (chemo(ichemo + im)%used) then
					chemo(ichemo + im)%present = .true.
					chemo(ichemo + im)%bdry_conc = 0
				endif
			enddo
			call MediumChange(V,C)
			call UpdateCbnd(0.0d0)
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
! Revised treatment of concentrations, to allow for setting O2 in medium change
!-----------------------------------------------------------------------------------------
subroutine MediumChange(Ve,Ce)
real(REAL_KIND) :: Ve, Ce(:)
real(REAL_KIND) :: R, Vm, Vr, Vblob
real(REAL_KIND) :: O2_drop_factor
integer :: ichemo, kcell, ix, iy, iz
type(cell_type), pointer :: cp
logical :: reset_O2_concs = .false.
logical :: new_fix = .true.
logical :: O2_drop

write(nflog,*) 'MediumChange:'
write(nflog,'(a,f8.4)') 'Ve: ',Ve
write(nflog,'(a,13f8.4)') 'Ce: ',Ce
write(nflog,'(a,13e12.3)')'medium_M: ',chemo(OXYGEN+1:)%medium_M
!call SetRadius(Nsites)
R = blob_radius*DELTA_X		! cm
Vblob = (4./3.)*PI*R**3	! cm3
Vm = total_volume - Vblob
Vr = min(Vm,Ve)
write(nflog,'(a,4f8.4)') 'total_volume, Vblob, Vm, Vr: ',total_volume, Vblob, Vm, Vr
chemo(:)%medium_M = ((Vm - Vr)/Vm)*chemo(:)%medium_M + Ve*Ce(:)
total_volume = Vm - Vr + Ve + Vblob
chemo(:)%medium_Cext = chemo(:)%medium_M/(total_volume - Vblob)
chemo(:)%medium_Cbnd = chemo(:)%medium_Cext
!chemo(OXYGEN+1:)%medium_Cext = chemo(OXYGEN+1:)%medium_M/(total_volume - Vblob)
!chemo(OXYGEN)%medium_Cext = chemo(OXYGEN)%bdry_conc
write(nflog,'(a,13e12.3)')'medium_M: ',chemo(OXYGEN+1:)%medium_M
write(nflog,'(a,13f8.4)') 'medium_Cext ',chemo(OXYGEN+1:)%medium_Cext
if (use_FD) then	! need to set medium concentrations in Cave
!	do ichemo = OXYGEN+1,MAX_CHEMO
	do ichemo = 1,MAX_CHEMO
		if (chemo(ichemo)%used) then
			chemo(ichemo)%Cave_b = chemo(ichemo)%medium_Cext
			chemo(ichemo)%Cprev_b = chemo(ichemo)%medium_Cext
			if (new_fix) then
				chemo(ichemo)%Fprev_b = 0
				chemo(ichemo)%Fcurr_b = 0
			endif
			write(nflog,'(a,i2,e12.3)') 'set Cave_b: ',ichemo,chemo(ichemo)%medium_Cext
		endif
	enddo
endif
O2_drop = (Ce(OXYGEN) < chemo(OXYGEN)%bdry_conc)
if (O2_drop) then
	O2_drop_factor = Ce(OXYGEN)/chemo(OXYGEN)%bdry_conc
endif
chemo(OXYGEN)%bdry_conc = Ce(OXYGEN)
!t_lastmediumchange = istep*DELTA_T
t_lastmediumchange = t_simulation

! This is new code
! Note: bdry_conc = 0 is a special case, no framp needed
if (chemo(OXYGEN)%bdry_conc > 0) then
	if (O2_drop .and. reset_O2_concs) then		! reset IC and EC for O2 only
		do kcell = 1,nlist
			cp => cell_list(kcell)
			if (cp%state == DEAD) cycle
			cp%conc(OXYGEN) = O2_drop_factor*cp%conc(OXYGEN)
			cp%Cex(OXYGEN) = O2_drop_factor*cp%Cex(OXYGEN)
		enddo
	endif
!	medium_change_step = .true.
endif
medium_change_step = .true.

if (.not.new_fix) return

! Reset concs at sites outside the blob
do ix = 1,NX
	do iy = 1,NY
		do iz = 1,NZ
			if (occupancy(ix,iy,iz)%indx(1) == OUTSIDE_TAG) then
				occupancy(ix,iy,iz)%C(:) = chemo(:)%medium_Cbnd
			endif
		enddo
	enddo
enddo
		
end subroutine

!-----------------------------------------------------------------------------------------
! If the volume removed is Vr, the fraction of constituent mass that is retained
! in the medium is (Vm - Vr)/Vm.  The calculation does not apply to oxygen.
! Usually Vr = Ve.
!-----------------------------------------------------------------------------------------
subroutine MediumChange1(Ve,Ce)
real(REAL_KIND) :: Ve, Ce(:)
real(REAL_KIND) :: R, Vm, Vr, Vblob
integer :: ichemo

write(nflog,*) 'MediumChange:'
write(nflog,'(a,f8.4)') 'Ve: ',Ve
write(nflog,'(a,13f8.4)') 'Ce: ',Ce
write(nflog,'(a,13e12.3)')'medium_M: ',chemo(OXYGEN+1:)%medium_M
!call SetRadius(Nsites)
R = blob_radius*DELTA_X		! cm
Vblob = (4./3.)*PI*R**3	! cm3
Vm = total_volume - Vblob
Vr = min(Vm,Ve)
write(nflog,'(a,4f8.4)') 'total_volume, Vblob, Vm, Vr: ',total_volume, Vblob, Vm, Vr
chemo(OXYGEN+1:)%medium_M = ((Vm - Vr)/Vm)*chemo(OXYGEN+1:)%medium_M + Ve*Ce(OXYGEN+1:)
total_volume = Vm - Vr + Ve + Vblob
chemo(OXYGEN+1:)%medium_Cext = chemo(OXYGEN+1:)%medium_M/(total_volume - Vblob)
chemo(OXYGEN)%medium_Cext = chemo(OXYGEN)%bdry_conc
write(nflog,'(a,13e12.3)')'medium_M: ',chemo(OXYGEN+1:)%medium_M
write(nflog,'(a,13f8.4)') 'medium_Cext ',chemo(OXYGEN+1:)%medium_Cext
if (use_FD) then	! need to set medium concentrations in Cave
	do ichemo = OXYGEN+1,MAX_CHEMO
		if (chemo(ichemo)%used) then
			chemo(ichemo)%Cave_b = chemo(ichemo)%medium_Cext
			chemo(ichemo)%Cprev_b = chemo(ichemo)%medium_Cext
			write(nflog,'(a,i2,e12.3)') 'set Cave_b: ',ichemo,chemo(ichemo)%medium_Cext
		endif
	enddo
endif

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
real(REAL_KIND) :: r(3), rmax, tstart, dt, radiation_dose, diam_um, framp, tnow, area, diam
!integer, parameter :: NT_CONC = 6
integer :: i, ic, ichemo, ndt, iz
integer :: nvars, ns
real(REAL_KIND) :: dxc, ex_conc(120*O2_BY_VOL+1)		! just for testing
real(REAL_KIND) :: DELTA_T_save, t_sim_0, t_ramp = 3600
integer :: ndiv, idiv
logical :: ok = .true.
logical :: dbug

!call logger('simulate_step')
!write(*,'(a,f8.3)') 'simulate_step: time: ',wtime()-start_wtime
!write(nflog,'(a,f8.3)') 'simulate_step: time: ',wtime()-start_wtime
dbug = .false.
if (Ncells == 0) then
	call logger('Ncells = 0')
    res = 2
    return
endif
if (limit_stop) then
	call logger('Spheroid size limit reached')
	res = 6
	return
endif

call getVolume(blob_volume,blob_area)
blob_radius = sqrt(blob_area/PI)		! number of sites
blob_centre = getCentre()   ! units = sites

nthour = 3600/DELTA_T
dt = DELTA_T/NT_CONC

if (ngaps > 200) then
	call squeezer
endif

bdry_debug = (istep < -250000)

istep = istep + 1
tnow = istep*DELTA_T
t_simulation = (istep-1)*DELTA_T	! seconds

if (bdry_debug) write(*,*) 'istep, bdry_changed: ',istep,bdry_changed
if (bdry_changed) then
	if (dbug) write(nflog,*) 'UpdateBdryList'
	call UpdateBdrylist
!	write(nflog,*) 'did UpdateBdryList'
endif
if (mod(istep,6*nthour) == 0) then
	if (dbug) write(nflog,*) 'CheckBdryList'
	call CheckBdryList('simulate_step')
!	write(nflog,*) 'did CheckBdryList'
endif

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

if (dbug) write(nflog,*) 'GrowCells'
!write(*,*) 'GrowCells'
call GrowCells(radiation_dose,DELTA_T,ok)
if (bdry_debug) call CheckBdryList('simulate_step d')
if (dbug) write(nflog,*) 'did GrowCells'
!write(*,*) 'did GrowCells'
if (.not.ok) then
	res = 3
	return
endif
radiation_dose = 0
if (bdry_debug) call CheckBdryList('simulate_step c')

if (medium_change_step) then
	ndiv = 6
else
	ndiv = 1
endif
t_ramp = DELTA_T
DELTA_T_save = DELTA_T
DELTA_T = DELTA_T/ndiv
t_sim_0 = t_simulation
do idiv = 0,ndiv-1
t_simulation = t_sim_0 + idiv*DELTA_T

if (use_FD) then
	! The change in flux driving the FD solver is ramped over one time step, using a sinusoid
	framp = 1
	if (medium_change_step .and. t_simulation - t_lastmediumchange < t_ramp) then
		framp = (t_simulation - t_lastmediumchange)/t_ramp
		framp = sin(framp*PI/2)
		framp = max(framp,0.1)
	else
		framp = 1
	endif
	if (istep == 1) framp = 0.5
	call diff_solver(DELTA_T, framp,ok)
	if (.not.ok) then
		res = 4
		return
	endif
	call UpdateCbnd(DELTA_T)
endif
if (dbug) write(nflog,*) 'SetupODEdiff'
call SetupODEdiff

!call check_allstate('before SiteCellToState')
if (dbug) write(nflog,*) 'SiteCellToState'
call SiteCellToState
!call check_allstate('after SiteCellToState')

if (dbug) write(nflog,*) 'Solver'
do it_solve = 1,NT_CONC
	tstart = (it_solve-1)*dt
!	t_simulation = (istep-1)*DELTA_T + tstart
	t_simulation = t_simulation + (it_solve-1)*dt
	call Solver(it_solve,tstart,dt,Ncells,ok)
	if (.not.ok) then
		res = 5
		return
	endif
enddo
!write(nflog,*) 'did Solver'
!call check_allstate('after Solver')

if (dbug) write(nflog,*) 'StateToSiteCell'
call StateToSiteCell
!call check_allstate('after StateToSiteCell')
call CheckDrugConcs

if (.not.use_FD) then
	call UpdateCbnd(DELTA_T)		! need to check placements of UpdateCbnd
!	write(nflog,*) 'did UpdateCbnd'
endif
call CheckDrugPresence

enddo
DELTA_T = DELTA_T_save
medium_change_step = .false.

res = 0

if (mod(istep,60) == -1) then
	rmax = 0
	do kcell = 1,nlist
		r = cell_list(kcell)%site - blob_centre
		rmax = max(rmax,norm(r))
	enddo
	hour = istep/60
	write(logmsg,'(3i6,2f6.1)') istep, hour, Ncells, blob_radius, rmax
	call logger(logmsg)
!	call ShowConcs
	call check_bdry
endif
!call test_CellDivision
if (.not.use_TCP .and. (mod(istep,6) == 0)) then
	call get_concdata(nvars, ns, dxc, ex_conc)
endif

if (saveprofile%active) then
	if (istep*DELTA_T >= saveprofile%it*saveprofile%dt) then
		call WriteProfileData
		saveprofile%it = saveprofile%it + 1
		if (saveprofile%it > saveprofile%nt) then
			saveprofile%active = .false.
		endif
	endif
endif
if (saveslice%active) then
	if (istep*DELTA_T >= saveslice%it*saveslice%dt) then
		call WriteSliceData
		saveslice%it = saveslice%it + 1
		if (saveslice%it > saveslice%nt) then
			saveslice%active = .false.
		endif
	endif
endif
if (saveFACS%active) then
	if (istep*DELTA_T >= saveFACS%it*saveFACS%dt) then
		call WriteFACSData
		saveFACS%it = saveFACS%it + 1
		if (saveFACS%it > saveFACS%nt) then
			saveFACS%active = .false.
		endif
	endif
endif

if (dbug .or. mod(istep,nthour) == 0) then
	diam_um = 2*DELTA_X*blob_radius*10000
	write(logmsg,'(a,2i6,a,2i6,a,f6.1,a,i2,a,2f6.3)') &
		'istep, hour: ',istep,istep/nthour,' nlist, ncells: ',nlist,ncells,' diam: ',diam_um,' nchemo: ',nchemo
	call logger(logmsg)
	write(nflog,*) 'Ndrug_tag: ',Ndrug_tag(1,1)
	call showcells
endif
! write(nflog,'(a,f8.3)') 'did simulate_step: time: ',wtime()-start_wtime 

if (istep == istep_output_cell_data) then
	call output_cell_data
	stop
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine output_cell_data
integer :: kcell
type(cell_type), pointer :: cp

do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	write(nfres,'(4e12.3)') cp%volume,cp%conc(TRACER+1:TRACER+3)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine showcell(kcell)
integer :: kcell
type(cell_type), pointer :: cp

cp => cell_list(kcell)
write(nflog,'(a,i6,4e12.3)') 'kcell, volume, divide_volume, dVdt, divide_time: ', &
                kcell, cp%volume, cp%divide_volume, cp%dVdt, cp%divide_time
end subroutine

!-----------------------------------------------------------------------------------------
! Average volumes etc
!-----------------------------------------------------------------------------------------
subroutine showcells
integer :: kcell, n
real(REAL_KIND) :: Vsum,divVsum,dVdtsum,divtsum
!real(REAL_KIND) :: Vn   ! to normalise volumes
type(cell_type), pointer :: cp

!Vn = Vdivide0/1.6
Vsum=0
divVsum=0
dVdtsum=0
divtsum=0
n=0
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%state == DEAD) cycle
    n = n+1
    Vsum = Vsum + cp%volume
    divVsum = divVsum + cp%divide_volume
    dVdtsum = dVdtsum + cp%dVdt
    divtsum = divtsum + cp%divide_time
enddo
!write(nflog,'(a,4e12.3)') 'ave volume, divide_volume, dVdt, divide_time: ', Vsum/n,divVsum/n,dVdtsum/n,divtsum/n
!write(*,'(a,4e12.3)') 'ave volume, divide_volume, dVdt, divide_time: ', Vsum/n,divVsum/n,dVdtsum/n,divtsum/n
end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Execute(ncpu,infile_array,inbuflen,outfile_array,outbuflen,res) BIND(C) 
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char) :: infile_array(128), outfile_array(128)
integer(c_int) :: ncpu, inbuflen, outbuflen, res
character*(128) :: infile, outfile
logical :: ok, success
integer :: i

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
elseif (res == 2) then
	call logger(' No more live cells')
elseif (res == 6) then
	call logger(' Spheroid size limit reached')
elseif (res == 3) then
	call logger(' === Execution failed === ERROR in GrowCells')
elseif (res == 4) then
	call logger(' === Execution failed === ERROR in diff_solver')
elseif (res == 5) then
	call logger(' === Execution failed === ERROR in Solver')
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
if (allocated(gaplist)) deallocate(gaplist,stat=ierr)
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
