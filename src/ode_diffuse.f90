! Solving 3D diffusion-decay eqtn as an ODE system.
! 
! Each chemokine has a list of sites with a specified rate of secretion into each site.
! The list needs to be updated when the boundary changes, or when FDCs move and/or the number
! changes.

!----------------------------------------------------------------------------------
! Note: The value of spcrad was determined by writing out the value computed in rkc
!----------------------------------------------------------------------------------
double precision function spcrad(neqn,t,y)
!DEC$ ATTRIBUTES DLLEXPORT :: spcrad
use global
integer :: neqn
double precision :: t, y(neqn)
!spcrad = 4d0*((nx+1)**2 + (ny+1)**2 + (nz+1)**2)
!spcrad = max(70.,4*BLOB_RADIUS)
spcrad = 75     ! maybe need to be increased with INTRA added to EXTRA variables
end function

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
module ode_diffuse

use chemokine
use rkf45s
use rksuite_90
use rksuite_90_prec, only:wp
use rkc_90

implicit none

integer, parameter :: MAX_VARS = 200000

integer, parameter :: RKF45_SOLVER = 1
integer, parameter :: RKSUITE_SOLVER = 2
integer, parameter :: RKC_SOLVER = 3
logical, parameter :: EXPLICIT_INTRA = .false.
real(REAL_KIND), allocatable :: allstate(:,:)
real(REAL_KIND), allocatable :: allstatep(:,:)
real(REAL_KIND), allocatable :: work_rkc(:,:)
!integer, allocatable :: cell_index(:)
!integer, allocatable :: intra_index(:)
!integer, allocatable :: extra_index(:)

integer :: nchemo, chemomap(MAX_CHEMO)
integer :: ivdbug
contains

!----------------------------------------------------------------------------------
! t in seconds
!----------------------------------------------------------------------------------
real(REAL_KIND) function BdryConc(ichemo,t)
integer :: ichemo
real(REAL_KIND) :: t

if ((use_medium_flux .and. medium_volume > 0) .or. .not.chemo(ichemo)%bdry_decay) then
    BdryConc = chemo(ichemo)%bdry_conc
else
    BdryConc = chemo(ichemo)%bdry_conc*exp(-chemo(ichemo)%bdry_decay_rate*t)
endif
!write(*,*) ichemo, BdryConc, exp(-chemo(ichemo)%bdry_decay_rate*30*60)
end function

!----------------------------------------------------------------------------------
! Solve for a test case with solute flux into the gridcells with x = 1.
! For each active site (x,y,z) mapping to i = ivar(x,y,z) there are
! 7 coefficients of the associated variables, i.e. those at:
! 1 <-- (x,y,z)
! 2 <-- (x-1,y,z)
! 3 <-- (x+1,y,z)
! 4 <-- (x,y-1,z)
! 5 <-- (x,y+1,z)
! 6 <-- (x,y,z-1)
! 7 <-- (x,y,z+1)
!
! PLUS for reactions, MAX_CHEMO 
! In each case, provided ivar(:,:,:) > 0
!
! If diffusion in an axis direction is suppressed (there is a boundary on one or
! the other neighbour sites), the current approximation used is to set the
! first derivative in that direction to zero, and compute the second derivative
! in the usual way.
! E.g. if we are considering the point (x,y,z) and ivar(x-1,y,z)=0,
! then the Laplacian d2C/dx2 + d2C/dy2 + d2C/dz2 that determines diffusion
! is reduced to d2C/dy2 + d2C/dz2.
! This is a crude simplification - in fact we should use a one-sided expression
! for d2C/dx2, derived by polynomial fitting.  The approach is as follows:
! Consider that x=0 is a boundary (dC/dx = 0), and let C(x) = a + bx + cx^2
! Then using the function values at the boundary and two grid points in,
! C0 = a
! C1 = a + b.dx + c.dx^2
! C2 = a + 2b.dx + 4c.dx^2
! solving these equations yields
! a = C0
! b = (-3C0 + 4C1 - C2)/(2dx)
! c = (C1 - C0 - b.dx)/(dx^2) = (C0 - 2C1 + C2)/(2dx^2)
! => d2C/dx2 = 2c = (C0 - 2C1 + C2)/dx^2
! In other words, the second derivative can be approximated by the second derivative
! at the neighbouring point (x+1,y,z), but this requires knowledge of concentration
! at more points.
!
! Note that now %coef(:,:) is not used.  Only %icoef(:,:) is needed.
!----------------------------------------------------------------------------------
subroutine SetupODEDiff
integer :: x, y, z, i, site(3), ifdc, ichemo, k, nc, nin, nex, kcell, ki, ie
!real(REAL_KIND) :: DX, DX2, c(MAX_CHEMO,7)
real(REAL_KIND) :: secretion
integer :: ierr
logical :: left, right

call logger('Set diffusion parameters')
if (.not.allocated(ODEdiff%ivar)) then
	allocate(ODEdiff%ivar(NX,NY,NZ))
endif
if (.not.allocated(ODEdiff%varsite)) then
	allocate(ODEdiff%varsite(MAX_VARS,3))
endif
if (.not.allocated(ODEdiff%icoef)) then
	allocate(ODEdiff%icoef(MAX_VARS,7))
endif
if (.not.allocated(ODEdiff%vartype)) then
	allocate(ODEdiff%vartype(MAX_VARS))
endif
if (.not.allocated(ODEdiff%cell_index)) then
	allocate(ODEdiff%cell_index(MAX_VARS))
endif
!DX = 1.0
!DX2 = DX*DX
ODEdiff%ivar = 0
i = 0
nex = 0
nin = 0
do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
		    kcell = occupancy(x,y,z)%indx(1)
			if (kcell /= OUTSIDE_TAG) then
				i = i+1
				nex = nex + 1
				ODEdiff%ivar(x,y,z) = i
				ODEdiff%varsite(i,:) = (/x,y,z/)
				ODEdiff%vartype(i) = EXTRA
				ODEdiff%cell_index(i) = kcell
!				if (y == 51 .and. z == 51) then
!				    write(logmsg,*) x,i,kcell
!				    call logger(logmsg)
!				endif
				if (kcell > 0) then ! with this numbering system the EXTRA variable corresponding to
				    i = i+1         ! an INTRA variable i is i-1
				    nin = nin + 1
    				ODEdiff%varsite(i,:) = (/x,y,z/)
	    			ODEdiff%vartype(i) = INTRA
    				ODEdiff%cell_index(i) = kcell
	    		endif
			endif
		enddo
	enddo
enddo
ODEdiff%nextra = nex
ODEdiff%nintra = nin
ODEdiff%nvars = i

if (.not.use_ODE_diffusion) return

!if (allocated(ODEdiff%icoef)) then
!	deallocate(ODEdiff%icoef)
!endif
!allocate(ODEdiff%icoef(ODEdiff%nextra,7))
!do ichemo = 1,MAX_CHEMO
!	if (.not.chemo(ichemo)%used) cycle
!	if (allocated(chemo(ichemo)%coef)) deallocate(chemo(ichemo)%coef)
!	allocate(chemo(ichemo)%coef(ODEdiff%nextra,7))
!enddo

ierr = 0
do i = 1,ODEdiff%nvars
!	c = 0
!	nc = 0
    if (ODEdiff%vartype(i) == INTRA) cycle
	site = ODEdiff%varsite(i,:)
	x = site(1)
	y = site(2)
	z = site(3)
	ODEdiff%icoef(i,:) = 0
	ODEdiff%icoef(i,1) = i
	
	left = .true.
	right = .true.
	if (x==1) then
!		left = .false.
		ierr = 1
		exit
	elseif (ODEdiff%ivar(x-1,y,z) == 0) then
		left = .false.
	endif
	if (x==NX) then 
!		right = .false.
		ierr = 2
		exit
	elseif (ODEdiff%ivar(x+1,y,z) == 0) then
		right = .false.
	endif
	if (left) then
		ODEdiff%icoef(i,2) = ODEdiff%ivar(x-1,y,z)
	endif
	if (right) then
		ODEdiff%icoef(i,3) = ODEdiff%ivar(x+1,y,z)
	endif
	left = .true.
	right = .true.
	if (y==1) then
!		left = .false.
		ierr = 3
		exit
	elseif (ODEdiff%ivar(x,y-1,z) == 0) then
		left = .false.
	endif
	if (y==NY) then 
!		right = .false.
		ierr = 4
		exit
	elseif (ODEdiff%ivar(x,y+1,z) == 0) then
		right = .false.
	endif
	if (left) then
		ODEdiff%icoef(i,4) = ODEdiff%ivar(x,y-1,z)
	endif
	if (right) then
		ODEdiff%icoef(i,5) = ODEdiff%ivar(x,y+1,z)
	endif
	left = .true.
	right = .true.
	if (z==1) then
!		left = .false.
		ierr = 5
		exit
	elseif (ODEdiff%ivar(x,y,z-1) == 0) then
		left = .false.
	endif
	if (z==NZ) then 
!		right = .false.
		ierr = 6
		exit
	elseif (ODEdiff%ivar(x,y,z+1) == 0) then
		right = .false.
	endif
	if (left) then
		ODEdiff%icoef(i,6) = ODEdiff%ivar(x,y,z-1)
	endif
	if (right) then
		ODEdiff%icoef(i,7) = ODEdiff%ivar(x,y,z+1)
	endif
enddo
if (ierr /= 0) then
	write(logmsg,*) 'Error: SetupODEDiff: lattice boundary reached: ierr: ',ierr
	call logger(logmsg)
	stop
endif
!if (allocated(extra_index)) deallocate(extra_index)
!if (allocated(intra_index)) deallocate(intra_index)
!if (allocated(cell_index)) deallocate(cell_index)
!allocate(extra_index(ODEdiff%nintra))   ! extra variable number associated with intra variable number
!allocate(intra_index(ODEdiff%nextra))   ! intra variable number associated with extra variable number
!allocate(cell_index(ODEdiff%nintra))    ! cell number associated with intra variable number

nchemo = 0
do ichemo = 1,MAX_CHEMO
	if (chemo(ichemo)%used) then
		nchemo = nchemo + 1
		chemomap(nchemo) = ichemo
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_ivar(kcell,msg)
integer :: kcell
character*(*) :: msg
integer :: site(3), iv

site = cell_list(kcell)%site
iv = ODEdiff%ivar(site(1),site(2),site(3))
write(*,*) 'check_ivar: ',msg,' ',kcell,site,iv
if (iv < 1) then
    stop
endif
end subroutine
!----------------------------------------------------------------------------------
! Makes a slight modification to the Michaelis-Menten function to create a
! "soft landing" as C -> 0
! f1(C) = (C-d)/(C0 + C-d)  is a shifted version of the MM curve
! f0(C) = kC^2             is a function with derivative -> 0 as C -> 0
! At C=e, we want f0(e) = f1(e), and f0'(e) = f1'(e)
! =>
! ke^2 = (e-d)/(C0 + e-d)
! 2ke = C0/(Co + e-d)^2
! => e/2 = (e-d)(C0 + e-d)/C0
! Set x = e-d
! C0(x+d) = 2x(C0+x)
! => x = e-d = (sqrt(Co^2 + 8dC0) - C0)/4
! k = (e-d)/(e^2(C0 + e-d))
! We fix d (as small as possible, by trial and error) then deduce e, k.
!----------------------------------------------------------------------------------
subroutine AdjustMM
real(REAL_KIND) :: deltaC, C0, C1

C0 = chemo(OXYGEN)%MM_C0
!deltaC = 5*CO2_DEATH_THRESHOLD
deltaC = 0.0005
C1 = deltaC + (sqrt(C0*C0 + 8*C0*deltaC) - C0)/4
ODEdiff%k = (C1-deltaC)/(C1*C1*(C0+C1-deltaC))
ODEdiff%C1 = C1
ODEdiff%deltaC = deltaC
write(logmsg,*) 'AdjustMM: C0, deltaC, C1, k: ',C0, ODEdiff%deltaC, ODEdiff%C1, ODEdiff%k
call logger(logmsg)
end subroutine

!----------------------------------------------------------------------------------
! A cell has moved to site(:), which was previously exterior (boundary), and the 
! variable-site mappings need to be updated, together with %icoef(:,:)
! The relevant neighbours are at x +/- 1, y +/- 1, z +/- 1
! %nextra is incremented, and allstate(nextra,:) is initialized.
!----------------------------------------------------------------------------------
subroutine ExtendODEDiff(site)
integer :: site(3)
integer :: x, y, z, x1, x2, y1, y2, z1, z2, n, kv, nb, ichemo
real(REAL_KIND) :: csum(MAX_CHEMO)

x = site(1)
y = site(2)
z = site(3)
n = ODEdiff%nextra + 1
if (n > MAX_VARS) then
	write(logmsg,*) 'Error: ExtendODEdiff: Too many variables: n > MAX_VARS: ',MAX_VARS
	call logger(logmsg)
	write(logmsg,*) 'Increase MAX_VARS and rebuild spheroid.DLL'
	call logger(logmsg)
	stop
endif
ODEdiff%ivar(x,y,z) = n
ODEdiff%varsite(n,:) = site
ODEdiff%nextra = n
ODEdiff%icoef(n,:) = 0
ODEdiff%icoef(n,1) = n
! See which neighbours of site are interior, and require to have 
! %icoef = 0 replaced by %icoef = n
nb = 0
csum = 0
x2 = x+1
kv = ODEdiff%ivar(x2,y,z) 
if (kv > 0) then
	if (ODEdiff%icoef(kv,2) /= 0) then
		write(logmsg,*) 'Error: ExtendODEDiff: icoef(kv,2): ',ODEdiff%icoef(kv,2)
		call logger(logmsg)
		stop
	endif
	ODEdiff%icoef(kv,2) = n
	ODEdiff%icoef(n,3) = kv
	nb = nb+1
	csum = csum + allstate(kv,:)
endif
x1 = x-1
kv = ODEdiff%ivar(x1,y,z) 
if (kv > 0) then
	if (ODEdiff%icoef(kv,3) /= 0) then
		write(logmsg,*) 'Error: ExtendODEDiff: icoef(kv,3): ',ODEdiff%icoef(kv,3)
		call logger(logmsg)
		stop
	endif
	ODEdiff%icoef(kv,3) = n
	ODEdiff%icoef(n,2) = kv
	nb = nb+1
	csum = csum + allstate(kv,:)
endif
y2 = y+1
kv = ODEdiff%ivar(x,y2,z) 
if (kv > 0) then
	if (ODEdiff%icoef(kv,4) /= 0) then
		write(logmsg,*) 'Error: ExtendODEDiff: icoef(kv,4): ',ODEdiff%icoef(kv,4)
		call logger(logmsg)
		stop
	endif
	ODEdiff%icoef(kv,4) = n
	ODEdiff%icoef(n,5) = kv
	nb = nb+1
	csum = csum + allstate(kv,:)
endif
y1 = y-1
kv = ODEdiff%ivar(x,y1,z) 
if (kv > 0) then
	if (ODEdiff%icoef(kv,5) /= 0) then
		write(logmsg,*) 'Error: ExtendODEDiff: icoef(kv,5): ',ODEdiff%icoef(kv,5)
		call logger(logmsg)
		stop
	endif
	ODEdiff%icoef(kv,5) = n
	ODEdiff%icoef(n,4) = kv
	nb = nb+1
	csum = csum + allstate(kv,:)
endif
z2 = z+1
kv = ODEdiff%ivar(x,y,z2) 
if (kv > 0) then
	if (ODEdiff%icoef(kv,6) /= 0) then
		write(logmsg,*) 'Error: ExtendODEDiff: icoef(kv,6): ',ODEdiff%icoef(kv,6)
		call logger(logmsg)
		stop
	endif
	ODEdiff%icoef(kv,6) = n
	ODEdiff%icoef(n,7) = kv
	nb = nb+1
	csum = csum + allstate(kv,:)
endif
z1 = z-1
kv = ODEdiff%ivar(x,y,z1) 
if (kv > 0) then
	if (ODEdiff%icoef(kv,7) /= 0) then
		write(logmsg,*) 'Error: ExtendODEDiff: icoef(kv,7): ',ODEdiff%icoef(kv,7)
		call logger(logmsg)
		stop
	endif
	ODEdiff%icoef(kv,7) = n
	ODEdiff%icoef(n,6) = kv
	nb = nb+1
	csum = csum + allstate(kv,:)
endif

!csum = csum + (6-nb)*chemo(:)%bdry_conc
do ichemo = 1,MAX_CHEMO
    csum(ichemo) = csum(ichemo) + (6-nb)*BdryConc(ichemo,t_simulation)
enddo
allstate(n,:) = csum/6
!allstatep(n,:) = 0
end subroutine

!----------------------------------------------------------------------------------
! Interpolate site and cell concentrations on cell division
!----------------------------------------------------------------------------------
subroutine InterpolateConc(site)
integer :: site(3)
integer :: x, y, z, x1, x2, y1, y2, z1, z2, n, kv, nb
real(REAL_KIND) :: csum(MAX_CHEMO)

end subroutine

!----------------------------------------------------------------------------------
! Copy site and cell concentrations to allstate(:,:)
!----------------------------------------------------------------------------------
subroutine SiteCellToState
integer :: i, site(3), kcell

do i = 1,ODEdiff%nvars
    site = ODEdiff%varsite(i,:)
    if (ODEdiff%vartype(i) == EXTRA) then
        allstate(i,:) = occupancy(site(1),site(2),site(3))%C(:)
!        write(*,'(4i6,4f8.3)') i,site,allstate(i,:)
    else
!        kcell = occupancy(site(1),site(2),site(3))%indx(1)
        kcell = ODEdiff%cell_index(i)
        allstate(i,:) = cell_list(kcell)%conc(:)
    endif
enddo
end subroutine

!----------------------------------------------------------------------------------
! Copy allstate(:,:) to site and cell concentrations 
!----------------------------------------------------------------------------------
subroutine StateToSiteCell
integer :: i, site(3), kcell

do i = 1,ODEdiff%nvars
    site = ODEdiff%varsite(i,:)
    if (ODEdiff%vartype(i) == EXTRA) then
        occupancy(site(1),site(2),site(3))%C(:) = allstate(i,:)
    else
!        kcell = occupancy(site(1),site(2),site(3))%indx(1)
        kcell = ODEdiff%cell_index(i)
        cell_list(kcell)%conc(:) = allstate(i,:)
    endif
enddo

end subroutine

!----------------------------------------------------------------------------------
! The rate of change of mass of constituent ichemo as a result of reactions is computed,
! with the assumption that the concentrations of the other constituents are fixed
! at the values from the previous time step, stored in allstate(:,:)
! isite = current site index
! C = current vector of concentrations
! dMdt = rate of change of mass of C(ichemo) at a site as a result of reactions (mol/sec)
! Note: the rate of consumption (of O2 or G) is in general a function of both the
! cell state and the concentration.  The consumption rate must in any case be less
! than some fraction of the total mass in the site (to avoid negative concentrations).
! We can use a Michaelis-Menten function to take the rate to zero as C -> 0.
! Note that currently Vextra is fixed - no accounting for cell death, gaps etc.
!----------------------------------------------------------------------------------
subroutine extra_react(ichemo,iv,C,vol,dCreact)
integer :: ichemo, iv
integer :: site(3), kcell
real(REAL_KIND) :: C(:), vol, dCreact
real(REAL_KIND) :: metab, dMdt
real(REAL_KIND) :: metab1, metab2, alpha

! Check for necrotic site - for now, no reactions
if (use_death) then
	site = ODEdiff%varsite(iv,:)
	if (occupancy(site(1),site(2),site(3))%indx(1) < 0) then
		dCreact = 0
		return
	endif
endif

dCreact = 0
if (ichemo == OXYGEN .or. ichemo == GLUCOSE) then
!	metab = max(0.0,C(OXYGEN))/(chemo(OXYGEN)%MM_C0 + C(OXYGEN))
!	metab = max(CO2_DEATH_THRESHOLD,C(OXYGEN))/(chemo(OXYGEN)%MM_C0 + C(OXYGEN))
	if (C(OXYGEN) > ODEdiff%C1) then
		metab = (C(OXYGEN)-ODEdiff%deltaC)/(chemo(OXYGEN)%MM_C0 + C(OXYGEN) - ODEdiff%deltaC)
	elseif (C(OXYGEN) > 0) then
		metab = ODEdiff%k*C(OXYGEN)*C(OXYGEN)
	else
		metab = 0
	endif
!	dMdt = -metab*chemo(ichemo)%max_cell_rate	! mol/s
!	dCreact = dMdt*1.0e6/Vextra	! convert mass rate (mol/s) to concentration rate (mM/s)
	dCreact = -metab*chemo(ichemo)%max_cell_rate*1.0e6/vol	! convert mass rate (mol/s) to concentration rate (mM/s)
endif
end subroutine


!----------------------------------------------------------------------------------
! In the original version, the neqn variables (for any constituent) are associated
! with the lattice sites encompassed by the blob.  Not all of these sites contain
! live cells.
! This can be generalized by making the ODEdiff variables apply to the extracellular
! concentration fields.
! One approach is to extend the vector of extracellular variables, for a given 
! constituent, to include the corresponding vector of intracellular variables for
! that constituent.  The approach would then be to treat all other constituents as
! constant for the purposes of solving for one constituent.
! Need a mapping between extra_index <-> intra_index
!                        1,...,nextra    1,...,nintra
! Say:
! intra_index(ie) = ki, 0 if no cell at site for extracellular variable ie
! extra_index(ki) = ie
! ki = 0
! do ie = 1, nextra
!   site = ODEdiff%varsite(ie,:)
!	kcell = occupancy(site(1),site(2),site(3))%indx(1)
!   if (kcell > 0) then
!       ki = ki + 1
!       intra_index(ie) = ki
!       extra_index(ki) = ie
!       cell_index(ki) = kcell
!       allstate(nextra+ki,:) = cell_list(kcell)%conc
!   else
!       intra_index(ie) = 0
!   endif
! enddo 
!
! Then kcell is recovered from ki: ki -> ie -> site -> indx(1)
! or from cell_index(ki)
!.....................................................................................
! This version has the intracellular variables interleaved with the extracellular.
! Volumes:
! If cell volume did not change, and there was no cell death, every site would contain
! fixed extra- and intracellular volumes, Vextra + Vcell = Vsite
! When there is no cell the extracellular volume is Vsite.
!----------------------------------------------------------------------------------
subroutine f_rkc(neqn,t,v,dvdt,icase)
integer :: neqn, icase
real(REAL_KIND) :: t, v(neqn), dvdt(neqn)
integer :: i, k, ie, ki, kv, nextra, nintra, ichemo, site(3), kcell
real(REAL_KIND) :: dCsum, dCdiff, dCreact,  DX2, DX3, vol, val, Cin(MAX_CHEMO), Cex
real(REAL_KIND) :: decay_rate, dc1, dc6, cbnd
logical :: bnd, dbug
real(REAL_KIND) :: metab, dMdt
logical :: use_compartments = .true.
logical :: intracellular, cell_exists

!write(*,*) 'neqn: ',neqn
ichemo = icase
DX2 = DELTA_X*DELTA_X
decay_rate = chemo(ichemo)%decay_rate
dc1 = chemo(ichemo)%diff_coef/DX2
dc6 = 6*dc1 + decay_rate
!cbnd = chemo(ichemo)%bdry_conc
cbnd = BdryConc(ichemo,t_simulation)
!if (t < 1.0) write(*,*) icase,t
do i = 1,neqn
    if (ODEdiff%vartype(i) == EXTRA) then
        intracellular = .false.
		vol = Vextra
        Cex = v(i)
        cell_exists = .false.
        if (i < neqn) then
            if (ODEdiff%vartype(i+1) == INTRA) then
                cell_exists = .true.
                vol = Vsite
	            Cin = allstate(i+1,:)
	            Cin(ichemo) = v(i+1)
	        endif
	    endif
	else
        intracellular = .true.
		vol = Vsite - Vextra	! for now, ignoring cell volume change!!!!!
        Cex = v(i-1)
	    Cin = allstate(i,:)
	    Cin(ichemo) = v(i)
	endif
	if (.not.intracellular) then
	    dCsum = 0
	    do k = 1,7
		    kv = ODEdiff%icoef(i,k)
		    if (k == 1) then
			    dCdiff = -dc6
		    else
			    dCdiff = dc1
		    endif
		    if (kv == 0) then
			    val = cbnd
		    else
			    val = v(kv)
		    endif
		    dCsum = dCsum + dCdiff*val
	    enddo
	    if (cell_exists) then
!		    call extra_react(ichemo,i,Cin,vol,dCreact)
		    dCreact = -chemo(ichemo)%cell_diff*(Cex - Cin(ichemo))
		else
            dCreact=0
		endif
    	dvdt(i) = dCsum + dCreact
	else
	    call intra_react(ichemo,Cin,Cex,vol,dCreact)
	    dvdt(i) = dCreact
	endif
enddo
end subroutine

!----------------------------------------------------------------------------------
! Reactions here + cross-membrane diffusion
! Should metab depend on cell volume, stage in cell cycle?
!----------------------------------------------------------------------------------
subroutine intra_react(ichemo,Cin,Cex,vol,dCreact)
integer :: ichemo
real(REAL_KIND) :: Cin(:), Cex, vol, dCreact
real(REAL_KIND) :: metab

dCreact = 0
if (ichemo == OXYGEN .or. ichemo == GLUCOSE) then
    if (Cin(OXYGEN) > ODEdiff%C1) then
	    metab = (Cin(OXYGEN)-ODEdiff%deltaC)/(chemo(OXYGEN)%MM_C0 + Cin(OXYGEN) - ODEdiff%deltaC)
    elseif (Cin(OXYGEN) > 0) then
	    metab = ODEdiff%k*Cin(OXYGEN)*Cin(OXYGEN)
    else
	    metab = 0
	    write(*,*) 'metab = 0'
    endif
    dCreact = -metab*chemo(ichemo)%max_cell_rate*1.0e6/vol	! convert mass rate (mol/s) to concentration rate (mM/s)
elseif (ichemo == SN30000) then
    dCreact = -(SN30K%C1 + SN30K%C2*SN30K%KO2/(SN30K%KO2 + Cin(OXYGEN)))*SN30K%Kmet0*Cin(ichemo)
elseif (ichemo == SN30000_METAB) then
    dCreact = (SN30K%C1 + SN30K%C2*SN30K%KO2/(SN30K%KO2 + Cin(OXYGEN)))*SN30K%Kmet0*Cin(SN30000)
endif
dCreact = dCreact + chemo(ichemo)%cell_diff*(Cex - Cin(ichemo)) 
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine showresults(v)
real(REAL_KIND) :: v(:), Cin(100), drugM(100)
integer :: x,y,z,z1,z2,i,ki,kcell

Cin = 0
drugM = 0
x = NX/2
y = NY/2
z1 = NZ/2
do z = z1,NZ/2+9
	i = ODEdiff%ivar(x,y,z)
	if (i == 0) exit
	z2 = z
	if (ODEdiff%cell_index(1) > 0) then
	    Cin(z) = v(i+1)
!	    kcell = occupancy(x,y,z)%indx(1)
!	    drugM(z) = cell_list(kcell)%M
	else
	    Cin(z) = 0
	endif
!	ki = intra_index(i)
!	if (ki > 0) then
!    	Ce(z) = v(ki+ODEdiff%nextra)
!    else
!        Ce(z) = 0
!    endif
enddo
write(logmsg,'(a,i6,10f7.4)') 'E:',ODEdiff%nextra,(v(ODEdiff%ivar(x,y,z)),z=z1,z2)	
call logger(logmsg)
write(logmsg,'(a,i6,10f7.4)') 'I:',ODEdiff%nintra,(Cin(z),z=z1,z2)	
call logger(logmsg)
!write(logmsg,'(a,i6,10f7.4)') 'M:',ODEdiff%nintra,(drugM(z),z=z1,z2)	
!call logger(logmsg)
!write(logmsg,'(a,10f7.4)') 'Diff:   ',((v(ODEdiff%ivar(x,y,z))-Cin(z))/Cin(z),z=z1,z2)	
!call logger(logmsg)
end subroutine

!----------------------------------------------------------------------------------
! In this version the diffusion/decay of each constituent is solved by a separate
! OMP thread.  Obviously this requires at least as many CPUs as there are constituents.
! Note that this required modifications to the way the ODE solver handles SAVEd variables,
! to avoid collisions between different threads.
! The ODE solver is RKC
! The subroutine assumes that:
!   * ODEdiff has been set up with the correct mappings
!   * allstate(:,:) holds the most recent solution, including estimates when the 
!     blob changes (i.e. grows or shrinks).  This could entail variable renumbering.
!   * work(:,:) is correctly sized (ODEdiff%nextra)
!----------------------------------------------------------------------------------
subroutine Solver(it,tstart,dt,ncells)
integer :: it, ncells
real(REAL_KIND) :: tstart, dt
integer :: ichemo, nvars, ntvars, ic, kcell, site(3), iv
integer :: ie, ki, i
real(REAL_KIND) :: t, tend
real(REAL_KIND), allocatable :: state(:,:)
real(REAL_KIND) :: C(MAX_CHEMO), Ce(MAX_CHEMO), dCreact(MAX_CHEMO)
real(REAL_KIND) :: dCsum, dC
real(REAL_KIND) :: timer1, timer2
logical :: ok
! Variables for RKC
integer :: info(4), idid
real(REAL_KIND) :: rtol, atol(1), sprad_ratio
type(rkc_comm) :: comm_rkc(MAX_CHEMO)

ODEdiff%nintra = ncells
ntvars = ODEdiff%nextra + ODEdiff%nintra
if (EXPLICIT_INTRA) then
    nvars = ODEdiff%nextra
else
    nvars = ntvars
endif
allocate(state(ntvars,MAX_CHEMO))

state(:,:) = allstate(1:ntvars,1:MAX_CHEMO)
if (it == 1) then
	call showresults(state(:,OXYGEN))
endif

info(1) = 1
info(2) = 1		! 1 = use spcrad() to estimate spectral radius, 2 = let rkc do it 
info(3) = 1
info(4) = 0
rtol = 1d-2
atol = rtol

!$omp parallel do private(t, tend, idid, ichemo)
do ic = 1,nchemo
	ichemo = chemomap(ic)
!	if (.not.chemo(ichemo)%used) cycle
	idid = 0
	t = tstart
	tend = t + dt
	call rkc(comm_rkc(ichemo),nvars,f_rkc,state(:,ichemo),t,tend,rtol,atol,info,work_rkc(:,ichemo),idid,ichemo)
	if (idid /= 1) then
		write(logmsg,*) ' Failed at t = ',t,' with idid = ',idid
		call logger(logmsg)
		stop
	endif
	if (info(2) == 2 .and. ichemo == OXYGEN) then
		sprad_ratio = rkc_sprad/blob_radius
		write(logmsg,'(a,2f8.4)') 'sprad_ratio: ',blob_radius,sprad_ratio
		call logger(logmsg)
	endif
enddo
!$omp end parallel do

! This doesn't work.  Because the O2 level goes very low in hypoxic cells, the explicit step easily
! generates negative O2 concentrations.  It seems that the intracellular computation needs to be treated
! implicitly together with the extracellular.
!if (EXPLICIT_INTRA) then
!    do ki = 1,ODEdiff%nintra
!        i = ki + ODEdiff%nextra
!        ie = extra_index(ki)
!	    C = state(i,:)     ! intracellular
!	    Ce = state(ie,:)   ! extracellular
!!        call intra_react(nchemo,C,Ce,dCreact)
!	    do ichemo = 1,nchemo
!	        state(i,ichemo) = C(ichemo) + dCreact(ichemo)*dt
!!	        if (ichemo == 1) write(*,*) ki,ichemo,dC
!	    enddo
!	    if (state(i,OXYGEN) > 0.2) then
!	        write(*,*) 'C(OXYGEN) too big: ',ki,ie,Ce(OXYGEN),C(OXYGEN),state(i,OXYGEN)
!	        stop
!	    endif
!    enddo
!endif

if (use_medium_flux .and. medium_volume > 0) then
	call update_medium(ntvars,state,dt)
endif

allstate(1:nvars,1:MAX_CHEMO) = state(:,:)
! Note: some time we need to copy the state values to the cell_list() array.

deallocate(state)
!deallocate(extra_index)
!deallocate(intra_index)
!deallocate(cell_index)

end subroutine

!----------------------------------------------------------------------------------
! The medium concentrations are updated explicitly, assuming a sphere with boundary
! concentrations equal to the mean extracellular concentrations of boundary sites.
! Note that concentrations of O2 and glucose are not varied.
!----------------------------------------------------------------------------------
subroutine update_medium(ntvars,state,dt)
integer :: ntvars
real(REAL_KIND) :: dt, state(:,:)
integer :: nb, i, k
real(REAL_KIND) :: Cbnd(MAX_CHEMO), F(MAX_CHEMO)
logical :: bnd

! First need the spheroid radius
call SetRadius(Nsites)
! Now compute the mean boundary site concentrations Cbnd(:)
nb = 0
Cbnd = 0
do i = 1,ntvars
	if (ODEdiff%vartype(i) /= EXTRA) cycle
	bnd = .false.
	do k = 1,7
		if (ODEdiff%icoef(i,k) == 0) then
			bnd = .true.
			exit
		endif
	enddo
	if (bnd) then
		nb = nb + 1
		Cbnd = Cbnd + state(i,:)
	endif
enddo
Cbnd = Cbnd/nb
F(:) = 4*PI*Radius*Radius*DELTA_X*chemo(:)%diff_coef*(Cbnd(:) - chemo(:)%bdry_conc*dt)
!write(*,*) 'Radius, C, F: ',Radius, chemo(DRUG_A)%bdry_conc, F(DRUG_A)
chemo(DRUG_A:MAX_CHEMO)%bdry_conc = (chemo(DRUG_A:MAX_CHEMO)%bdry_conc*medium_volume + F(DRUG_A:MAX_CHEMO)*dt)/medium_volume
end subroutine

!----------------------------------------------------------------------------------
! Update intracellular concentrations for cell kcell by time dt.
! The extracellular concentrations are Cextra(:)
! The mass rate of transport is as for the extracellular calculation (they must 
! balance) but the conversion to concentration change is different (volume is different).
!----------------------------------------------------------------------------------
subroutine intracellular(kcell,Cextra,dt)
integer :: kcell
real(REAL_KIND) :: Cextra(:), dt

end subroutine

!----------------------------------------------------------------------------------
! Initialise the state vector to the current concentrations
!----------------------------------------------------------------------------------
subroutine InitState(ichemo,state)
integer :: ichemo
real(REAL_KIND) :: state(:)
integer :: x, y, z, i, site(3), nz
real(REAL_KIND) :: smin, smax

write(logmsg,*) 'InitState: ',chemo(ichemo)%name
call logger(logmsg)
smin = 1.0e10
smax = -smin
do i = 1,ODEdiff%nextra
!	site = ODEdiff%varsite(i,:)
	state(i) = chemo(ichemo)%bdry_conc
!	state(i) = chemo(ichemo)%conc(site(1),site(2),site(3))
!	if (state(i) < smin) smin = state(i)
!	if (state(i) > smax) smax = state(i)
enddo
end subroutine

!----------------------------------------------------------------------------------
! Initialise the state vector to the current concentrations
!----------------------------------------------------------------------------------
subroutine InitStates(ichemo,n,allstate)
integer :: ichemo,n
real(REAL_KIND) :: allstate(n,*)
integer :: x, y, z, i, site(3), nz
real(REAL_KIND) :: smin, smax

write(logmsg,*) 'InitState: ',chemo(ichemo)%name
call logger(logmsg)
smin = 1.0e10
smax = -smin
do i = 1,ODEdiff%nextra
	site = ODEdiff%varsite(i,:)
	allstate(i,ichemo) = chemo(ichemo)%conc(site(1),site(2),site(3))
enddo
end subroutine

!---------------------------------------------------------------------------------- 
!----------------------------------------------------------------------------------
subroutine TestSolver
integer :: nvars, nt, ichemo, it, ncells
real(REAL_KIND) :: tstart, dt
real(REAL_KIND) :: timer1, timer2
real(REAL_KIND) :: trun = 400

nt = 20
ncells = 1000
dt = trun/nt
timer1 = wtime()
tstart = 0
do it = 1,nt
!	write(*,'(a,i4,f8.1)') 'it, tstart: ',it,tstart
	tstart = (it-1)*dt
	call Solver(it,tstart,dt,ncells)
enddo
timer2 = wtime()
write(*,'(a,f10.1)') 'Time: ',timer2-timer1
stop
end subroutine

!----------------------------------------------------------------------------------
! In this version the diffusion/decay of each constituent is solved by a separate
! OMP thread.  Obviously this requires at least as many CPUs as there are constituents.
! Note that this required modifications to the way r4_rkf45 handles SAVEd variables,
! to avoid collisions between different threads.
!----------------------------------------------------------------------------------
subroutine TestODEDiffusion
integer flag, i, j, k, ichemo
real(REAL_KIND) :: tstart, tend, t, relerr, abserr
!real(REAL_KIND), allocatable :: allstate(:,:), allstatep(:,:)
real(REAL_KIND), allocatable :: state(:), statep(:)
real(REAL_KIND), parameter :: dt = 10.0
real(REAL_KIND) :: res(10)
integer :: x, y, z, dx, dy, dz, xmid, ymid, zmid, ibnd, k0, site(3), nchemo, imax, nvars
real(REAL_KIND) :: grad(3)
real(REAL_KIND) :: amp, r, ctemp(10), dvmax
real(REAL_KIND) :: t1, t2
integer :: nt = 60
logical :: ok
integer :: flag_par(MAX_CHEMO)
real(REAL_KIND) :: abserr_par(MAX_CHEMO), relerr_par(MAX_CHEMO)
integer, parameter :: SOLVER = RKC_SOLVER

! Variables for RKSUITE
type(rk_comm_real_1d) :: comm
real(kind=wp) :: t_start, t_end, tolerance, t_want, t_inc, t_got
real(kind=wp), allocatable :: thresholds(:)

! Variables for RKC
real(REAL_KIND), allocatable :: work(:,:)
integer :: info(4), idid
real(REAL_KIND) :: rtol, atol(1)
type(rkc_comm) :: comm_rkc(MAX_CHEMO)

write(logmsg,*) 'TestODEDiffusion: nvars: ',ODEdiff%nextra
call logger(logmsg)
nchemo = 2	!MAX_CHEMO
nvars = ODEdiff%nextra
!allocate(prev_state(ODEdiff%nextra))
!allocate(allstate(MAX_VARS,MAX_CHEMO))
!allocate(allstatep(MAX_VARS,MAX_CHEMO))
allocate(allstate(nvars,MAX_CHEMO))
allocate(allstatep(nvars,MAX_CHEMO))
allocate(work(8+4*nvars,MAX_CHEMO))
xmid = NX/2
ymid = NY/2
zmid = NZ/2
k0 = ODEdiff%ivar(xmid,ymid,zmid)
tstart = 0
do ichemo = 1,nchemo
	if (.not.chemo(ichemo)%used) cycle
	call InitState(ichemo,allstate(1:nvars,ichemo))
	allstatep(:,ichemo) = 0
	if (SOLVER == RKF45_SOLVER) then
		call deriv(tstart,allstate(1:nvars,ichemo),allstatep(1:nvars,ichemo),ichemo)
	elseif (SOLVER == RKSUITE_SOLVER) then
		t_end = 10000
		tolerance = 1.0e-4
		allocate(thresholds(nvars))
		thresholds = 0.1
		call rk_setup(comm, tstart, allstate(1:nvars,ichemo), t_end,  tolerance, thresholds)
	elseif (SOLVER == RKC_SOLVER) then
		info(1) = 1
		info(2) = 1
		info(3) = 1
		info(4) = 0
		rtol = 1d-2
		atol = rtol
		idid = 0
	endif
!	write(*,'(10f7.3)') allstate(k0:k0+9,ichemo)
enddo

t1 = wtime()
flag_par = 1
do k = 1,nt

	!$omp parallel do private(tstart, tend, state, statep, flag, relerr, abserr, site, idid)
	do ichemo = 1,nchemo
		if (.not.chemo(ichemo)%used) cycle
		allocate(state(nvars))
		allocate(statep(nvars))
		state = allstate(1:nvars,ichemo)
		statep = allstatep(1:nvars,ichemo)
		if (ichemo == 1) then
!			write(*,'(i6,10f7.3)') k,state(k0:k0+9)
			call showresults(state)
		endif
		if (k == 1) then
			abserr = 10*sqrt ( epsilon ( abserr ) )
			abserr_par(ichemo) = abserr
			relerr = 10*sqrt ( epsilon ( relerr ) )
			relerr_par(ichemo) = relerr
		else
			abserr = abserr_par(ichemo)
			relerr = relerr_par(ichemo)
		endif
		tstart = (k-1)*dt
		tend = tstart + dt
		flag = flag_par(ichemo)
		if (SOLVER == RKF45_SOLVER) then
			if (REAL_KIND == SP) then	
!				call r4_rkf45 ( deriv, nvars, state, statep, tstart, tend, relerr, abserr, flag, ichemo )
			else
				call r8_rkf45 ( deriv, nvars, state, statep, tstart, tend, relerr, abserr, flag, ichemo )
			endif
			if (flag /= 2) then
				write(logmsg,*) 'Bad flag: ',flag
				call logger(logmsg)
				stop
				call r8_rkf45 ( deriv, nvars, state, statep, tstart, tend, relerr, abserr, flag, ichemo )
				flag = 2
			endif
			flag_par(ichemo) = flag
			abserr_par(ichemo) = abserr
			relerr_par(ichemo) = relerr
		elseif (SOLVER == RKSUITE_SOLVER) then
			call range_integrate(comm, f_deriv, tend, t_got, state, statep, flag)
			if (flag /= 1 .and. flag /= 3 .and. flag /= 4) then
				write(logmsg,*) 'Bad flag: ',flag
				call logger(logmsg)
				stop
			endif
			write(*,'(a,3f10.3)') 'tstart: ',tstart,tend,t_got
		elseif (SOLVER == RKC_SOLVER) then
			idid = 0
			call rkc(comm_rkc(ichemo),nvars,f_rkc,state,tstart,tend,rtol,atol,info,work(:,ichemo),idid,ichemo)
			if (idid /= 1) then
				write(*,*) ' Failed at t = ',tstart,' with idid = ',idid
				stop
			endif
		endif
!		do i = 1,ODEdiff%nextra
!			site = ODEdiff%varsite(i,:)
!			chemo(ichemo)%conc(site(1),site(2),site(3)) = state(i)
!		enddo
		allstate(1:nvars,ichemo) = state
		allstatep(1:nvars,ichemo) = statep
! To test derivs
!		call deriv(tend,state,statep,ichemo)
!		dvmax = 0
!		do i = 1,ODEdiff%nextra
!			if (abs(statep(i)) > dvmax) then
!				dvmax = abs(statep(i))
!				imax = i
!			endif
!		enddo
!		write(*,*) 'Max dv: ',imax,dvmax
!		write(*,'(7i6)') ODEdiff%icoef(imax,:)
		deallocate(state)
		deallocate(statep)
	enddo
!	call TestAddSite
!   nvars = ODEdiff%nextra
!		prev_state = state
enddo

t2 = wtime()
write(*,'(a,f10.1)') 'Time: ',t2-t1

!	do x = 1,NX
!		do y = 1,NY
!			do z = 1,NZ
!				if (occupancy(x,y,z)%indx(1) < 0) cycle
!				i = ODEdiff%ivar(x,y,z)
!				chemo(ichemo)%conc(x,y,z) = state(i)
!				call compute_gradient(state,x,y,z,grad)
!				chemo(ichemo)%grad(:,x,y,z) = grad
!			enddo
!		enddo
!	enddo
!enddo
!deallocate(state)
!deallocate(prev_state)
!deallocate(statep)
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine TestAddSite
integer :: i, j, x, y, z, site(3), kpar=0

do 
	i = random_int(1,ODEdiff%nextra,kpar)
	do j = 1,7
		if (ODEdiff%icoef(i,j) == 0) then
			site = ODEdiff%varsite(i,:)
			x = site(1)
			y = site(2)
			z = site(3)
			if (j == 2) then
				x = x-1
			elseif (j == 3) then
				x = x+1
			elseif (j == 4) then
				y = y-1
			elseif (j == 5) then
				y = y+1
			elseif (j == 6) then
				z = z-1
			elseif (j == 7) then
				z = z+1
			endif
			if (ODEdiff%ivar(x,y,z) /= 0) then
				write(*,*) 'Error: TestAddSite: ',ODEdiff%ivar(x,y,z)
				stop
			endif
			site = (/x,y,z/)
			write(*,*) 'Add site at bdry: ',site
			call ExtendODEDiff(site)
			return
		endif
	enddo
enddo
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine SolveSteadystate_B
integer flag, i, j, k, ichemo
real(REAL_KIND) :: tstart, tend, relerr, abserr
real(REAL_KIND), allocatable :: state(:), prev_state(:), statep(:)
real(REAL_KIND), parameter :: dt = 100.0
real(REAL_KIND) :: res(10)
integer :: x, y, z, dx, dy, dz, xmid, ymid, zmid, ibnd
real(REAL_KIND) :: grad(3)
real(REAL_KIND) :: amp, r, ctemp(10)
integer :: nt(2) = (/20,100/)	! Number of rkf45 iterations with bdry concentration (1) and bdry secretion (2)
logical :: ok

write(logmsg,*) 'SolveSteadystate_B: nextra: ',ODEdiff%nextra
call logger(logmsg)
allocate(state(ODEdiff%nextra))
allocate(prev_state(ODEdiff%nextra))
allocate(statep(ODEdiff%nextra))
xmid = NX/2
ymid = NY/2
zmid = NZ/2
do ichemo = 1,MAX_CHEMO
	ODEdiff%ichemo = ichemo
	if (.not.chemo(ichemo)%used) cycle
	if (chemo(ichemo)%use_secretion) then
		ibnd = 2
	else
		ibnd = 1
	endif
	call InitState(ichemo,state)
	statep = 0
	tstart = 0
	call deriv(tstart,state,statep,ichemo)

	abserr = sqrt ( epsilon ( abserr ) )
	relerr = sqrt ( epsilon ( relerr ) )

	prev_state = 0
	flag = 1
	do k = 1,nt(ibnd)
		tstart = (k-1)*dt
		tend = tstart + dt
		if (REAL_KIND == SP) then	
!			call r4_rkf45 ( deriv, ODEdiff%nextra, state, statep, tstart, tend, relerr, abserr, flag )
		else
!			call r8_rkf45 ( deriv, ODEdiff%nextra, state, statep, tstart, tend, relerr, abserr, flag )
		endif
		if (flag /= 2) then
			write(logmsg,*) 'Bad flag: ',flag
			call logger(logmsg)
		endif
		flag = 2
		call CheckConvergence(state,prev_state,ok)
		if (ok) exit
		prev_state = state
	enddo
	do x = 1,NX
		do y = 1,NY
			do z = 1,NZ
				if (occupancy(x,y,z)%indx(1) /= OUTSIDE_TAG) cycle
				i = ODEdiff%ivar(x,y,z)
				chemo(ichemo)%conc(x,y,z) = state(i)
				call compute_gradient(state,x,y,z,grad)
				chemo(ichemo)%grad(:,x,y,z) = grad
			enddo
		enddo
	enddo
enddo
deallocate(state)
deallocate(prev_state)
deallocate(statep)
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine CheckConvergence(s1,s2,ok)
real(REAL_KIND) :: s1(:), s2(:)
logical :: ok
integer :: i, imax
real(REAL_KIND) :: ds, df, dfmax, dsmax
real(REAL_KIND), parameter :: tol = 0.001

dfmax = 0
dsmax = 0
do i = 1,ODEdiff%nextra
	if (s1(i) > 0) then
		ds = (s1(i) - s2(i))
		df = ds/s1(i)
		if (abs(df) > abs(dfmax)) then
			imax = i
			dfmax = df
			dsmax = ds
		endif
	endif
enddo
if (abs(dfmax) < tol) then
	ok = .true.
else
	ok = .false.
endif
end subroutine

!----------------------------------------------------------------------------------
! Compute the gradient vector for chemokine concentration at (x,y,z).
! The gradient determination is 2-sided if possible, 1-sided if only one adjacent
! value is available, and the gradient is set to 0 if neither is possible.
!----------------------------------------------------------------------------------
subroutine compute_gradient(conc,x,y,z,grad)
real(REAL_KIND) :: conc(:)
integer :: x, y, z
real(REAL_KIND) :: grad(3)
integer :: i0, i1, i2
real(REAL_KIND) :: c0, c1, c2, del

i0 = ODEdiff%ivar(x,y,z)
c0 = conc(i0)
del = 2
if (x > 1) then
	i1 = ODEdiff%ivar(x-1,y,z)
else
	i1 = 0
endif
if (x < NX) then
	i2 = ODEdiff%ivar(x+1,y,z)
else
	i2 = 0
endif
if (i1 > 0) then
	c1 = conc(i1)
else
	c1 = c0
	del = del - 1
endif			
if (i2 > 0) then
	c2 = conc(i2)
else
	c2 = c0
	del = del - 1
endif
if (del > 0) then
	grad(1) = (c2 - c1)/del
else
	grad(1) = 0
endif			

del = 2
if (y > 1) then
	i1 = ODEdiff%ivar(x,y-1,z)
else
	i1 = 0
endif
if (y < NY) then
	i2 = ODEdiff%ivar(x,y+1,z)
else
	i2 = 0
endif
if (i1 > 0) then
	c1 = conc(i1)
else
	c1 = c0
	del = del - 1
endif			
if (i2 > 0) then
	c2 = conc(i2)
else
	c2 = c0
	del = del - 1
endif
if (del > 0) then
	grad(2) = (c2 - c1)/del
else
	grad(2) = 0
endif

del = 2
if (z > 1) then
	i1 = ODEdiff%ivar(x,y,z-1)
else
	i1 = 0
endif
if (z < NZ) then
	i2 = ODEdiff%ivar(x,y,z+1)
else
	i2 = 0
endif
if (i1 > 0) then
	c1 = conc(i1)
else
	c1 = c0
	del = del - 1
endif			
if (i2 > 0) then
	c2 = conc(i2)
else
	c2 = c0
	del = del - 1
endif
if (del > 0) then
	grad(3) = (c2 - c1)/del
else
	grad(3) = 0
endif
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
function f_deriv(t,v)
use rksuite_90_prec, only:wp
real(kind=wp), intent(in) :: t
real(kind=wp), dimension(:), intent(in) :: v
real(kind=wp), dimension(size(v,1)) :: f_deriv   

integer :: icase
integer :: i, k, kv, n, ichemo
real(REAL_KIND) :: dCsum, dCdiff, dCreact,  DX2, DX3, vol, val, C(MAX_CHEMO)
real(REAL_KIND) :: decay_rate, dc1, dc6, cbnd, dmax
logical :: bnd, dbug

!ichemo = icase
ichemo = 1
DX2 = DELTA_X*DELTA_X
decay_rate = chemo(ichemo)%decay_rate
dc1 = chemo(ichemo)%diff_coef/DX2
dc6 = 6*dc1 + decay_rate
vol = Vsite
!cbnd = chemo(ichemo)%bdry_conc
cbnd = BdryConc(ichemo,t_simulation)
n = ODEdiff%nextra
!if (t < 1.0) write(*,*) icase,t
dmax = 0
do i = 1,n
	C = allstate(i,:)
	C(ichemo) = v(i)
	dCsum = 0
	do k = 1,7
		kv = ODEdiff%icoef(i,k)
!		if (ODEdiff%icoef(i,k) /= 0) then	! interior
!			bnd = .false.
!		else								! boundary
!			bnd = .true.
!		endif
		if (k == 1) then
			dCdiff = -dc6
		else
			dCdiff = dc1
		endif
		if (kv == 0) then
			val = cbnd
		else
			val = v(kv)
		endif
		dCsum = dCsum + dCdiff*val
	enddo
	dCreact = 0
	call extra_react(ichemo,i,C,vol,dCreact)
!	if (i == 6822) then
!		write(*,'(5f8.4)') C(ichemo), dCsum, dCreact
!	endif
	f_deriv(i) = dCsum + dCreact
!	dmax = max(dmax,f_deriv(i))
enddo
end function

!----------------------------------------------------------------------------------
! Simple diffusion-decay for a single constituent (ODEdiff%ichemo)
! Now all chemo constituents share the same FD template, and there are always
! 6 neighbours (some may be boundary, with specified concentrations).
! For now assume uniform diffusion coefficient for a constituent
!----------------------------------------------------------------------------------
subroutine deriv(t,v,dv,icase)
real(REAL_KIND) :: t, v(*), dv(*)
integer :: icase
integer :: i, k, kv, n, ichemo
real(REAL_KIND) :: dCsum, dCdiff, dCreact,  DX2, DX3, vol, val, C(MAX_CHEMO)
real(REAL_KIND) :: decay_rate, dc1, dc6, cbnd
logical :: bnd, dbug

ichemo = icase
DX2 = DELTA_X*DELTA_X
decay_rate = chemo(ichemo)%decay_rate
dc1 = chemo(ichemo)%diff_coef/DX2
dc6 = 6*dc1 + decay_rate
vol = Vsite
!cbnd = chemo(ichemo)%bdry_conc
cbnd = BdryConc(ichemo,t_simulation)
n = ODEdiff%nextra
!if (t < 1.0) write(*,*) icase,t
do i = 1,n
	C = allstate(i,:)
	C(ichemo) = v(i)
	dCsum = 0
	do k = 1,7
		kv = ODEdiff%icoef(i,k)
!		if (ODEdiff%icoef(i,k) /= 0) then	! interior
!			bnd = .false.
!		else								! boundary
!			bnd = .true.
!		endif
		if (k == 1) then
			dCdiff = -dc6
		else
			dCdiff = dc1
		endif
		if (kv == 0) then
			val = cbnd
		else
			val = v(kv)
		endif
		dCsum = dCsum + dCdiff*val
	enddo
	call extra_react(ichemo,i,C,vol,dCreact)
!	if (i == 6822) then
!		write(*,'(5f8.4)') C(ichemo), dCsum, dCreact
!	endif
	dv(i) = dCsum + dCreact
enddo
end subroutine

!----------------------------------------------------------------------------------
! Now all chemo constituents share the same FD template, and there are always
! 6 neighbours (some may be boundary, with specified concentrations).
!----------------------------------------------------------------------------------
subroutine deriv_all(t,v,dv)
real(REAL_KIND) :: t, v(*), dv(*)
integer :: i, k, kv, n, x, y, z, ichemo, idbug, ifdc, site(3), dx, dy, dz, ic, nf_FDC, nf_MRC
real(REAL_KIND) :: csum(MAX_CHEMO), s, vtemp, ctemp, dc, DX2, val
logical :: bnd, dbug

DX2 = DELTA_X*DELTA_X
n = ODEdiff%nextra
do i = 1,n
	csum = 0
	do k = 1,7
!		if (chemo(ODEdiff%ichemo)%coef(i,k) /= 0) then
		if (ODEdiff%icoef(i,k) /= 0) then	! interior
			bnd = .false.
		else								! boundary
			bnd = .true.
		endif
		do ichemo = 1,MAX_CHEMO
			if (k == 1) then
				dc = -chemo(ichemo)%decay_rate - 6*chemo(ichemo)%diff_coef/DX2	! ODEdiff%ncoef(i)
			else
				dc = chemo(ichemo)%diff_coef/DX2
			endif
			if (bnd) then
!				val = chemo(ichemo)%bdry_conc
                val = BdryConc(ichemo,t_simulation)
			else
				kv = (ODEdiff%icoef(i,k)-1)*MAX_CHEMO + ichemo
				val = v(kv)
			endif
			csum(ichemo) = csum(ichemo) + dc*val
		enddo
	enddo
	! Need to add in reactions here
	do ichemo = 1,MAX_CHEMO
		kv = (i-1)*MAX_CHEMO + ichemo
		dv(kv) = csum(ichemo)
	enddo
enddo

end subroutine

!----------------------------------------------------------------------------------
! In the original version, the neqn variables (for any constituent) are associated
! with the lattice sites encompassed by the blob.  Not all of these sites contain
! live cells.
! This can be generalized by making the ODEdiff variables apply to the extracellular
! concentration fields.
! One approach is to extend the vector of extracellular variables, for a given 
! constituent, to include the corresponding vector of intracellular variables for
! that constituent.  The approach would then be to treat all other constituents as
! constant for the purposes of solving for one constituent.
! Need a mapping between extra_index <-> intra_index
!                        1,...,nextra    1,...,nintra
! Say:
! intra_index(ie) = ki, 0 if no cell at site for extracellular variable ie
! extra_index(ki) = ie
! ki = 0
! do ie = 1, nextra
!   site = ODEdiff%varsite(ie,:)
!	kcell = occupancy(site(1),site(2),site(3))%indx(1)
!   if (kcell > 0) then
!       ki = ki + 1
!       intra_index(ie) = ki
!       extra_index(ki) = ie
!       cell_index(ki) = kcell
!       allstate(nextra+ki,:) = cell_list(kcell)%conc
!   else
!       intra_index(ie) = 0
!   endif
! enddo 
!
! Then kcell is recovered from ki: ki -> ie -> site -> indx(1)
! or from cell_index(ki)
!----------------------------------------------------------------------------------
subroutine f_rkc1(neqn,t,v,dvdt,icase)
integer :: neqn, icase
real(REAL_KIND) :: t, v(neqn), dvdt(neqn)
integer :: i, k, ie, ki, kv, nextra, nintra, ichemo, site(3), kcell
real(REAL_KIND) :: dCsum, dCdiff, dCreact,  DX2, DX3, vol, val, C(MAX_CHEMO), Ce, Ci
real(REAL_KIND) :: decay_rate, dc1, dc6, cbnd
logical :: bnd, dbug
real(REAL_KIND) :: metab, dMdt, vs
logical :: use_compartments = .true.

!write(*,*) 'neqn: ',neqn
ichemo = icase
DX2 = DELTA_X*DELTA_X
decay_rate = chemo(ichemo)%decay_rate
dc1 = chemo(ichemo)%diff_coef/DX2
dc6 = 6*dc1 + decay_rate
!cbnd = chemo(ichemo)%bdry_conc
cbnd = BdryConc(ichemo,t_simulation)
!if (t < 1.0) write(*,*) icase,t
vs = 1.0e6/Vextra
do i = 1,ODEdiff%nextra
	C = allstate(i,:)
	C(ichemo) = v(i)
	dCsum = 0
	do k = 1,7
		kv = ODEdiff%icoef(i,k)
		if (k == 1) then
			dCdiff = -dc6
		else
			dCdiff = dc1
		endif
		if (kv == 0) then
			val = cbnd
		else
			val = v(kv)
		endif
		dCsum = dCsum + dCdiff*val
	enddo
    dCreact=0
!	if (use_react) then
!		call extra_react(ichemo,i,C,dCreact)
!	if (use_compartments) then
!!    	site = ODEdiff%varsite(i,:)
!!	    kcell = occupancy(site(1),site(2),site(3))%indx(1)  ! this can be stored directly as a link to i
!        ki = intra_index(i)
!	    if (ki > 0) then     ! there is an intracellular compartment
!	        if (EXPLICIT_INTRA) then
!	            Ce = allstate(ki+ODEdiff%nextra,ichemo)
!	        else
!	            Ce = v(ki+ODEdiff%nextra)
!    	    endif
!   	        dCreact = chemo(ichemo)%cell_diff*(Ce - C(ichemo))    !*vs
!	        !dCreact = chemo(ichemo)%cell_diff*(cell_list(kcell)%conc(ichemo) - C(ichemo))*vs
!	        ! Note: %cell_diff is the constant coef that determines the rate of diffusion across the cell membrane
!	        !       the mass rate (mol/s) must be converted to concentration rate (mM/s)
!	        !       the conversion factor is different for the extracellular and intracellular compartments
!	        !       how to handle changing cell volume during growth?
!	    endif
!	else
!        !GS - made some changes to reduce call to react timing
!        if (ichemo==OXYGEN .or. ichemo==GLUCOSE) then
!!GS            metab=allstate(i,ichemo)
!!GS	        dMdt = -(metab/(chemo(ichemo)%MM_C0 + metab))*chemo(ichemo)%max_cell_rate	! mol/s
!			metab = max(0.0,C(OXYGEN))/(chemo(OXYGEN)%MM_C0 + C(OXYGEN))
!!			metab = C(OXYGEN)/(chemo(OXYGEN)%MM_C0 + C(OXYGEN))
!!			dMdt = -metab*chemo(ichemo)%max_cell_rate	! mol/s
!!           dCreact = dMdt*vs ! convert mass rate (mol/s) to concentration rate (mM/s)
!            dCreact =  -metab*chemo(ichemo)%max_cell_rate*vs ! convert mass rate (mol/s) to concentration rate (mM/s)
!        endif
!    endif
	dvdt(i) = dCsum + dCreact
!	write(*,*) i,dCsum,dCreact
enddo
!if (dCreact > 0) stop
if (EXPLICIT_INTRA) return
!do ki = 1,ODEdiff%nintra
!    dCsum = 0
!    i = ki + ODEdiff%nextra
!    ie = extra_index(ki)
!	C = allstate(i,:)   ! intracellular
!	C(ichemo) = v(i)
!	! reactions here + cross-membrane diffusion
!!	call react(ichemo,i,C,dCreact)
!	dCreact = chemo(ichemo)%cell_diff*(v(ki+ODEdiff%nextra) - C(ichemo))    !*vs
!	dCsum = -1.0*C(ichemo)
!	dvdt(i) = dCsum + dCreact
!enddo
end subroutine

!----------------------------------------------------------------------------------
! Reactions here + cross-membrane diffusion
!----------------------------------------------------------------------------------
subroutine intra_react1(nchemo,C,Ce,dCreact)
integer :: nchemo
real(REAL_KIND) :: C(:), Ce(:), dCreact(:)
integer :: ichemo
real(REAL_KIND) :: metab

if (C(OXYGEN) > ODEdiff%C1) then
	metab = (C(OXYGEN)-ODEdiff%deltaC)/(chemo(OXYGEN)%MM_C0 + C(OXYGEN) - ODEdiff%deltaC)
elseif (C(OXYGEN) > 0) then
	metab = ODEdiff%k*C(OXYGEN)*C(OXYGEN)
else
	metab = 0
endif

dCreact = 0
dCreact(1:2) = -metab*chemo(1:2)%max_cell_rate*1.0e6/Vextra	! convert mass rate (mol/s) to concentration rate (mM/s)

do ichemo = 1,nchemo
    dCreact(ichemo) = dCreact(ichemo) + chemo(ichemo)%cell_diff*(Ce(ichemo) - C(ichemo)) 
enddo
end subroutine

!----------------------------------------------------------------------------------
! For all living cells, update the total amount of SN30000 metabolized.
! The rate depends on the oxygen concentration C_ox and the SN30000 level, C_drug
! NOT USED
!----------------------------------------------------------------------------------
subroutine update_M(dt,state)
real(REAL_KIND) :: dt
real(REAL_KIND) :: state(:,:)
real(REAL_KIND) :: C_ox, C_drug, dCdt
integer :: i, kcell
real(REAL_KIND) :: Vcell = 1.0

do i = 1,ODEdiff%nvars
    if (ODEdiff%vartype(i) == INTRA) then
        kcell = ODEdiff%cell_index(i)
        C_ox = state(i,OXYGEN)
        C_drug = state(i,DRUG_A)
        dCdt = (SN30K%C1 + SN30K%C2*SN30K%KO2/(SN30K%KO2 + C_ox))*SN30K%Kmet0*C_drug
        cell_list(kcell)%M = cell_list(kcell)%M + Vcell*dCdt*dt
!        write(*,'(2i6,4e12.4)') i,kcell,C_ox,C_drug,Vcell*dCdt*dt,cell_list(kcell)%M
    endif
enddo
end subroutine



end module

