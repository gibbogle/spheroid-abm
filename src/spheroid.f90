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
!	call save_inputfile(fixedfile)
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
	call SetupODEDiff
	call InitConcs
!	call TestODEDiffusion
!	call TestSolver
endif
istep = 0
write(logmsg,'(a,i6)') 'Startup procedures have been executed: initial T cell count: ',Ncells0
call logger(logmsg)

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
! it will still be possible to view the chemokine concentration gradient fields.
if (allocated(occupancy)) deallocate(occupancy)
!do ichemo = 1,MAX_CHEMO
!	if (allocated(chemo(ichemo)%coef)) deallocate(chemo(ichemo)%coef)
!	if (allocated(chemo(ichemo)%conc)) deallocate(chemo(ichemo)%conc)
!	if (allocated(chemo(ichemo)%grad)) deallocate(chemo(ichemo)%grad)
!enddo

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
max_nlist = 10*nc0
write(logmsg,*) 'Initial radius: ',Radius, nc0
call logger(logmsg)

allocate(cell_list(max_nlist))
allocate(occupancy(NX,NY,NZ))
!allocate(gaplist(max_ngaps))

!do k = 1,max_nlist
!	nullify(cell_list(k)%cptr)
!enddo

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
integer :: itestcase, ncpu_dummy
real(REAL_KIND) :: days

ok = .true.
open(nfcell,file=inputfile,status='old')

read(nfcell,*) NX							! rule of thumb: about 4*BLOB_RADIUS
read(nfcell,*) initial_count				! initial number of tumour cells
read(nfcell,*) divide_time_median
read(nfcell,*) divide_time_shape
read(nfcell,*) days							! number of days to simulate
read(nfcell,*) DELTA_T						! time step size (sec)
read(nfcell,*) NT_CONC						! number of subdivisions of DELTA_T for diffusion computation
read(nfcell,*) itestcase                    ! test case to simulate
read(nfcell,*) seed(1)						! seed vector(1) for the RNGs
read(nfcell,*) seed(2)						! seed vector(2) for the RNGs
read(nfcell,*) ncpu_dummy					! just a placeholder for ncpu, not used currently
read(nfcell,*) NT_GUI_OUT					! interval between GUI outputs (timesteps)
read(nfcell,*) fixedfile					! file with "fixed" parameter values
close(nfcell)

blob_radius = (initial_count*3./(4.*PI))**(1./3)

! Setup test_case
test_case = .false.
if (itestcase /= 0) then
    test_case(itestcase) = .true.
endif

call read_fixed_params(ok)
if (.not.ok) then
	write(logmsg,'(a,a)') 'Error reading fixed input data file: ',fixedfile
	call logger(logmsg)
	return
endif

if (mod(NX,2) /= 0) NX = NX+1					! ensure that NX is even
open(nfout,file=outputfile,status='replace')
write(logmsg,*) 'Opened nfout: ',outputfile
call logger(logmsg)

Nsteps = days*24*60*60/DELTA_T		! DELTA_T in seconds
write(logmsg,'(a,2i6,f6.0)') 'nsteps, NT_CONC, DELTA_T: ',nsteps,NT_CONC,DELTA_T
call logger(logmsg)
!call setup_dists

ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine read_fixed_params(ok)
logical :: ok

ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine placeCells(ok)
logical :: ok
integer :: x, y, z, k, site(3)
real(REAL_KIND) :: r2lim,r2,r(3)

occupancy(:,:,:)%indx(1) = 0
occupancy(:,:,:)%indx(2) = 0
r2lim = Radius*Radius
lastID = 0
k = 0
do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
			r = (/x-x0,y-y0,z-z0/)
			r2 = dot_product(r,r)
			if (r2 < r2lim) then
				k = k+1
				lastID = lastID + 1
				site = (/x,y,z/)
				cell_list(k)%ID = lastID
				cell_list(k)%site = site
				cell_list(k)%state = 1
				cell_list(k)%exists = .true.
				occupancy(x,y,z)%indx(1) = k
			else
				occupancy(x,y,z)%indx(1) = OUTSIDE_TAG
			endif
		enddo
	enddo
enddo
nlist = k
Nsites = k
Ncells = k
Ncells0 = Ncells		
ok = .true.
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

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine cell_division
integer :: kcell0, kpar=0
integer :: j, k, kcell1, site0(3), site1(3), site2(3), site(3), npath, path(3,200)
type (boundary_type), pointer :: bdry

kcell0 = random_int(1,nlist,kpar)
site0 = cell_list(kcell0)%site
if (dbug) write(*,*) 'cell_division: ',kcell0,site0,occupancy(site0(1),site0(2),site0(3))%indx
if (bdrylist_present(site0,bdrylist)) then	! site0 is on the boundary
	npath = 1
	site1 = site0
	path(:,1) = site0
else
	call choose_bdrysite(site0,site1)
	call get_path(site0,site1,path,npath)
endif
if (dbug) write(*,*) 'path: ',npath
do k = 1,npath
	if (dbug) write(*,'(i3,2x,3i4)') path(:,k)
enddo

!call get_nearbdrypath(site0,path,npath)
! Need to choose an outside site near site1
call get_outsidesite(site1,site2)
npath = npath+1
path(:,npath) = site2
if (dbug) write(*,*) 'outside site: ',site2,occupancy(site2(1),site2(2),site2(3))%indx

call push_path(path,npath)

if (dbug) write(*,*) 'did push_path'
call add_cell(kcell0,kcell1,site0)
! Now need to fix the bdrylist.  
! site1 was on the boundary, but may no longer be.
! site2 is now on the boundary
! First add site2
if (dbug) write(*,*) 'add site2 to bdrylist: ',site2
if (isbdry(site2)) then   ! add it to the bdrylist
    nbdry = nbdry + 1
    allocate(bdry)
    bdry%site = site2
!    bdry%chemo_influx = .false.
    nullify(bdry%next)
    call bdrylist_insert(bdry,bdrylist)
!    call AssignBdryRole(site,bdry)
    occupancy(site2(1),site2(2),site2(3))%bdry => bdry
    call SetBdryConcs(site2)
else
    write(logmsg,*) 'Added site is not bdry: ',site2,occupancy(site2(1),site2(2),site2(3))%indx
	call logger(logmsg)
    stop
endif
if (dbug) write(*,*) 'Check for changed boundary status'
! Now check sites near site2 that may have changed their boundary status (including site1)
do j = 1,6
	site = site2 + neumann(:,j)
	if (dbug) write(*,*) j,site
	if (isbdry(site)) then
		if (dbug) write(*,*) 'isbdry'
		if (.not.bdrylist_present(site,bdrylist)) then	! add it
			if (dbug) write(*,*) 'not present, add it'
			nbdry = nbdry + 1
			allocate(bdry)
			bdry%site = site
		!    bdry%chemo_influx = .false.
			nullify(bdry%next)
			call bdrylist_insert(bdry,bdrylist)
		!    call AssignBdryRole(site,bdry)
			occupancy(site(1),site(2),site(3))%bdry => bdry
		    call SetBdryConcs(site)
		endif
	else
		if (dbug) write(*,*) 'not isbdry'
		if (bdrylist_present(site,bdrylist)) then	! remove it
			if (dbug) write(*,*) 'present, remove it'
			call bdrylist_delete(site,bdrylist)
			nullify(occupancy(site(1),site(2),site(3))%bdry)
			nbdry = nbdry - 1
			call ResetConcs(site)
		endif
	endif
enddo
if (dbug) write(*,*) 'done!'
end subroutine

!-----------------------------------------------------------------------------------------
! Need to choose a site on the boundary in some sense near site0, and also to preserve
! the spherical shape
!-----------------------------------------------------------------------------------------
subroutine choose_bdrysite(site0,site1)
integer :: site0(3), site1(3)
integer :: site(3), sitemin(3)
real(REAL_KIND) :: vc(3), v(3), r, rfrac, d, alpha_max, cosa_min, dmin, cosa
logical :: hit
type (boundary_type), pointer :: bdry

if (dbug) write(*,*) 'choose_bdrysite: ',site0
vc = site0 - Centre
r = norm(vc)
vc = vc/r
rfrac = r/Radius
alpha_max = get_alpha_max(rfrac)
cosa_min = cos(alpha_max)
dmin = 1.0e10
hit = .false.
bdry => bdrylist
do while ( associated ( bdry )) 
    site = bdry%site
    v = site - Centre
    d = norm(v)
    cosa = dot_product(v,vc)/d
    if (cosa > cosa_min) then
		hit = .true.
		if (d < dmin) then
			dmin = d
			sitemin = site
		endif
	endif
    bdry => bdry%next
enddo
if (.not.hit) then
	write(logmsg,*) 'Error: choose_bdrysite: no candidate bdry site'
	call logger(logmsg)
	stop
endif
site1 = sitemin
if (dbug) write(*,*) 'site1: ',site1,occupancy(site1(1),site1(2),site1(3))%indx
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function get_alpha_max(r)
real(REAL_KIND) :: r
real(REAL_KIND) :: alf
real(REAL_KIND) :: r1 = 0.2, r2 = 0.8, alpha1 = PI, alpha2 = PI/6

if (r < r1) then
	get_alpha_max = alpha1
elseif (r > r2) then
	get_alpha_max = alpha2
else
	alf = (r-r1)/(r2-r1)
	get_alpha_max = (1-alf)*alpha1 + alf*alpha2
endif
end function

!-----------------------------------------------------------------------------------------
! For a given cell site site0(:), find a suitable free site on the boundary of
! the blob, and the path of sites leading from site0 to site1
!-----------------------------------------------------------------------------------------
subroutine get_nearbdrypath(site0,path,npath)
integer :: site0(3), path(3,200),npath
integer :: site1(3), jump(3), k, j, jmax, kpar=0
real(REAL_KIND) :: v(3), r, d, dmax, v_aim(3), v_try(3)
real(REAL_KIND) :: cosa, sina, d_try, del

v = site0 - Centre
do j = 1,3
	v(j) = v(j) + (par_uni(kpar) - 0.5)
enddo
r = norm(v)
v_aim = v/r
!write(*,'(a,3i4,3f6.3)') 'site0: ',site0,v_aim
k = 1
site1 = site0
path(:,k) = site1
do 
	dmax = 0
	do j = 1,6
		if (j == 14) cycle
		jump = neumann(:,j)
!		jump = jumpvec(:,j)
		v_try = site1 + jump - Centre
!		call get_vnorm(v,v_try)
		d_try = norm(v_try)
		d = dot_product(v_try,v_aim)
		cosa = d/d_try
		sina = sqrt(1 - cosa*cosa)
		del = d_try*sina
!		write(*,'(i6,2x,3i3,2x,3f6.3,2x,f8.4)') j,jump,v_try,d
		if (d-del > dmax) then
			dmax = d-del
			jmax = j
		endif
	enddo
	site1 = site1 + neumann(:,jmax)
!	site1 = site1 + jumpvec(:,jmax)
	k = k+1
	path(:,k) = site1
!	write(*,'(i3,2x,3i4,2x,i2,2x,3i3)') k,site1,jmax,neumann(:,jmax)
	if (occupancy(site1(1),site1(2),site1(3))%indx(1) == OUTSIDE_TAG) exit
enddo
npath = k

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_get_path
integer :: site1(3), site2(3), path(3,200), npath, kpar=0
integer :: x, y, z, it, k

do it = 1,10
	do
		x = random_int(1,NX,kpar)
		y = random_int(1,NY,kpar)
		z = random_int(1,NZ,kpar)
		if (occupancy(x,y,z)%indx(1) > 0) exit
	enddo
	site1 = (/x,y,z/)
	do
		x = random_int(1,NX,kpar)
		y = random_int(1,NY,kpar)
		z = random_int(1,NZ,kpar)
		if (occupancy(x,y,z)%indx(1) > 0) exit
	enddo
	site2 = (/x,y,z/)
	call get_path(site1,site2,path,npath)
	write(*,*) 'path: ',npath
	do k = 1,npath
		write(*,'(i3,2x,3i4)') path(:,k)
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Find a path of sites leading from site1 to site2
!-----------------------------------------------------------------------------------------
subroutine get_path(site1,site2,path,npath)
integer :: site1(3), site2(3), path(3,200), npath
integer :: v(3), jump(3), site(3), k, j, jmin
real(REAL_KIND) :: r, d2, d2min

if (dbug) write(*,'(a,3i4,2x,3i4)') 'site1, site2: ',site1, site2
k = 1
site = site1
path(:,k) = site
do 
	d2min = 1.0e10
	do j = 1,6
		jump = neumann(:,j)
		v = site + jump
		if (occupancy(v(1),v(2),v(3))%indx(1) < 0) cycle
		v = site2 - v
		d2 = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
		if (d2 < d2min) then
			d2min = d2
			jmin = j
		endif
	enddo
	site = site + neumann(:,jmin)
	k = k+1
	if (dbug) write(*,*) k,site,jmin,d2min
	path(:,k) = site
	if (site(1) == site2(1) .and. site(2) == site2(2) .and. site(3) == site2(3)) exit
enddo
npath = k

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine push_path(path,npath)
integer :: path(3,200),npath
integer :: k, site1(3), site2(3), kcell

do k = npath-1,1,-1
	site1 = path(:,k)
	kcell = occupancy(site1(1),site1(2),site1(3))%indx(1)
	if (dbug) write(*,*) k,' site1: ',site1,kcell
	site2 = path(:,k+1)
	if (dbug) write(*,*) 'site2: ',site2
	cell_list(kcell)%site = site2
	occupancy(site2(1),site2(2),site2(3))%indx(1) = kcell
enddo
occupancy(site1(1),site1(2),site1(3))%indx = 0
Nsites = Nsites + 1
call SetRadius(Nsites)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine add_cell(kcell0,kcell1,site)
integer :: kcell0, kcell1, site(3)

lastID = lastID + 1
nlist = nlist + 1
Ncells = Ncells + 1
kcell1 = nlist
cell_list(kcell1)%state = cell_list(kcell0)%state
cell_list(kcell1)%site = site
cell_list(kcell1)%ID = lastID
cell_list(kcell1)%exists = .true.
occupancy(site(1),site(2),site(3))%indx(1) = kcell1
end subroutine

!-----------------------------------------------------------------------------------------
! Look at all bdry sites (those with a neumann neighbour outside)
!-----------------------------------------------------------------------------------------
subroutine check_bdry
integer :: kcell, i, site(3), site1(3), minv(3), maxv(3)
real(REAL_KIND) :: v(3), r, rmin, rmax
logical :: bdry

rmin = 1.0e10
rmax = 0
do kcell = 1,nlist
	site = cell_list(kcell)%site
	bdry = .false.
	do i = 1,6
		site1 = site + neumann(:,i)
		if (occupancy(site1(1),site1(2),site1(3))%indx(1) == OUTSIDE_TAG) then
			bdry = .true.
			exit
		endif
	enddo
	if (bdry) then
		v = site - Centre
		r = norm(v)
		if (r < rmin) then
			rmin = r
			minv = site - Centre
		endif
		if (r > rmax) then
			rmax = r
			maxv = site - Centre
		endif
	endif
enddo
write(*,'(a,2(f7.1,3i4,2x))') 'rmin, rmax: ',rmin,minv,rmax,maxv
end subroutine

!-----------------------------------------------------------------------------------------
! No use
!-----------------------------------------------------------------------------------------
subroutine smoother
integer :: x, y, z, k, site(3), kpar=0
integer :: Nsmooth = 1000

do k = 1,Nsmooth
	x = NX/2 + 1
	y = random_int(1,NY,kpar)
	z = random_int(1,NZ,kpar)
	if (occupancy(x,y,z)%indx(1) >= 0) then
		do x = NX/2+2,NX
			if (occupancy(x,y,z)%indx(1) == OUTSIDE_TAG) exit
		enddo
		site = (/x,y,z/)
		call adjust(site)
	endif
	x = NX/2
	y = random_int(1,NY,kpar)
	z = random_int(1,NZ,kpar)
	if (occupancy(x,y,z)%indx(1) >= 0) then
		do x = NX/2-1,1,-1
			if (occupancy(x,y,z)%indx(1) == OUTSIDE_TAG) exit
		enddo
		site = (/x,y,z/)
		call adjust(site)
	endif
	x = random_int(1,NX,kpar)
	y = NY/2 + 1
	z = random_int(1,NZ,kpar)
	if (occupancy(x,y,z)%indx(1) >= 0) then
		do y = NY/2+2,NY
			if (occupancy(x,y,z)%indx(1) == OUTSIDE_TAG) exit
		enddo
		site = (/x,y,z/)
		call adjust(site)
	endif
	x = random_int(1,NX,kpar)
	y = NY/2
	z = random_int(1,NZ,kpar)
	if (occupancy(x,y,z)%indx(1) >= 0) then
		do y = NY/2-1,1,-1
			if (occupancy(x,y,z)%indx(1) == OUTSIDE_TAG) exit
		enddo
		site = (/x,y,z/)
		call adjust(site)
	endif
	x = random_int(1,NX,kpar)
	y = random_int(1,NY,kpar)
	z = NZ/2 + 1
	if (occupancy(x,y,z)%indx(1) >= 0) then
		do z = NZ/2+2,NZ
			if (occupancy(x,y,z)%indx(1) == OUTSIDE_TAG) exit
		enddo
		site = (/x,y,z/)
		call adjust(site)
	endif
	x = random_int(1,NX,kpar)
	y = random_int(1,NY,kpar)
	z = NZ/2
	if (occupancy(x,y,z)%indx(1) >= 0) then
		do z = NZ/2-1,1,-1
			if (occupancy(x,y,z)%indx(1) == OUTSIDE_TAG) exit
		enddo
		site = (/x,y,z/)
		call adjust(site)
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! The location site0 is just outside the blob.
! The aim is to find a cell near site0 that is further from the centre, and move it here.
!-----------------------------------------------------------------------------------------
subroutine adjust(site0)
integer :: site0(3)
integer :: i, site(3), kcell, sitemax(3)
real(REAL_KIND) :: r0, r, rmax

if (occupancy(site0(1),site0(2),site0(3))%indx(1) /= OUTSIDE_TAG) then
	write(*,*) 'Error: adjust: site is not OUTSIDE: ',site0
	stop
endif
r0 = cdistance(site0)
!write(*,'(a,3i4,f6.2)') 'adjust: ',site0,r0
rmax = 0
do i = 1,6
	site = site0 + neumann(:,i)
	kcell = occupancy(site(1),site(2),site(3))%indx(1)
	if (kcell > 0) then
		r = cdistance(site)
		if (r > r0 .and. r > rmax) then	! move the cell here
			rmax = r
			sitemax = site
		endif
	endif
enddo
if (rmax > 0) then
	kcell = occupancy(sitemax(1),sitemax(2),sitemax(3))%indx(1)
!	write(*,'(i6,2x,3i4,f6.2)') kcell,sitemax,rmax
	cell_list(kcell)%site = site0
	occupancy(site0(1),site0(2),site0(3))%indx(1) = kcell
	occupancy(sitemax(1),sitemax(2),sitemax(3))%indx(1) = OUTSIDE_TAG
endif
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_dimensions(NX_dim,NY_dim,NZ_dim,nsteps_dim) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dimensions
use, intrinsic :: iso_c_binding
integer(c_int) :: NX_dim,NY_dim,NZ_dim,nsteps_dim

NX_dim = NX
NY_dim = NY
NZ_dim = NZ
nsteps_dim = nsteps
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
integer :: kcell, site(3), hour, it, kpar=0
real(REAL_KIND) :: r(3), rmax, tstart, dt
!integer, parameter :: NT_CONC = 6
integer :: nchemo

dt = DELTA_T/NT_CONC
istep = istep + 1
write(logmsg,*) 'istep: ',istep
call logger(logmsg)
!if (istep == 6216) dbug = .true.
!call make_split(.true.)
!call UpdateFields(DELTA_T)
do it = 1,NT_CONC
	tstart = (it-1)*dt
!	call Solver(nchemo,tstart,dt)
	call Solver(tstart,dt)
enddo
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
!call cell_division
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

! Cells
do kcell = 1,nlist
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
nTC_list = last_id2
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
!-----------------------------------------------------------------------------------------
subroutine get_summary(summaryData) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*)

summaryData(1:20) = 0

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
	call connecter(ok)
	if (.not.ok) then
		call logger('Failed to make TCP connections')
		return
	endif
endif
call setup(ncpu,infile,outfile,ok)
if (ok) then
!	clear_to_send = .true.
!	simulation_start = .true.
	istep = 0
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
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine connecter(ok)
logical :: ok
integer :: error

! Main connection
ok = .true.
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
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine wrapup
integer :: ierr, ichemo
logical :: isopen

call logger('doing wrapup ...')
ierr = 0
if (allocated(zoffset)) deallocate(zoffset)
if (allocated(zdomain)) deallocate(zdomain)
!if (allocated(occupancy)) deallocate(occupancy)
if (allocated(cell_list)) deallocate(cell_list,stat=ierr)
if (allocated(allstate)) deallocate(allstate)
if (allocated(allstatep)) deallocate(allstatep)
if (allocated(work_rkc)) deallocate(work_rkc)
do ichemo = 1,MAX_CHEMO
	if (allocated(chemo(ichemo)%coef)) deallocate(chemo(ichemo)%coef)
	if (allocated(chemo(ichemo)%conc)) deallocate(chemo(ichemo)%conc)
	if (allocated(chemo(ichemo)%grad)) deallocate(chemo(ichemo)%grad)
enddo
if (allocated(ODEdiff%ivar)) deallocate(ODEdiff%ivar)
if (allocated(ODEdiff%varsite)) deallocate(ODEdiff%varsite)
if (allocated(ODEdiff%icoef)) deallocate(ODEdiff%icoef)
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
