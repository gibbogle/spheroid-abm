! Solving 3D diffusion-decay eqtn as an ODE system.
!
! Each chemokine has a list of sites with a specified rate of secretion into each site.
! The list needs to be updated when the boundary changes, or when FDCs move and/or the number
! changes.

!----------------------------------------------------------------------------------
! Note: The value of spcrad was first determined by writing out the value computed in rkc.
! Later it was just determined by trial, then made into a run parameter.
!----------------------------------------------------------------------------------
double precision function spcrad(neqn,t,y)
!DEC$ ATTRIBUTES DLLEXPORT :: spcrad
use global
integer :: neqn
double precision :: t, y(neqn)
!spcrad = 4d0*((nx+1)**2 + (ny+1)**2 + (nz+1)**2)
!spcrad = max(70.,4*BLOB_RADIUS)
!spcrad = 600	!180	!75     ! maybe need to be increased with INTRA added to EXTRA variables
spcrad = spcrad_value
end function

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
module ode_diffuse

use chemokine
!use rkf45s
!use rksuite_90
!use rksuite_90_prec, only:wp
use rkc_90

implicit none

integer, parameter :: MAX_VARS = 2*max_nlist

integer, parameter :: RKF45_SOLVE = 1
integer, parameter :: RKSUITE_SOLVE = 2
integer, parameter :: RKC_SOLVE = 3
integer, parameter :: IRKC_SOLVE = 4
logical, parameter :: EXPLICIT_INTRA = .false.
real(REAL_KIND), allocatable :: allstate(:,:)
real(REAL_KIND), allocatable :: allstatep(:,:)
real(REAL_KIND), allocatable :: work_rkc(:)

integer :: nchemo, chemomap(MAX_CHEMO)
integer :: ivdbug

real(REAL_KIND) :: ichemo_sol
integer, parameter :: RK_SOLVER = RKC_SOLVE

contains

!----------------------------------------------------------------------------------
! t in seconds
! %bdry_decay = .true. for the case that a simple exponential decay of bdry conc
! is specified.
! Now uses %medium_Cbnd, no t dependence here.
!----------------------------------------------------------------------------------
real(REAL_KIND) function BdryConc(ichemo,t)
integer :: ichemo
real(REAL_KIND) :: t

BdryConc = chemo(ichemo)%medium_Cbnd
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
! With dropping, need to account for boundary at zmin-1, the surface that the blob
! is sitting on.  
!----------------------------------------------------------------------------------
subroutine SetupODEDiff
integer :: x, y, z, i, site(3), ifdc, ichemo, k, nc, nin, nex, kcell, ki, ie	!, isite
integer :: kextra, j
real(REAL_KIND) :: secretion
integer, allocatable :: exmap(:)
integer :: ierr
logical :: left, right

!call logger('Set diffusion parameters')
if (.not.allocated(ODEdiff%ivar)) then
	allocate(ODEdiff%ivar(NX,NY,NZ))
endif
if (.not.allocated(ODEdiff%varsite)) then
	allocate(ODEdiff%varsite(MAX_VARS,3))
endif
if (.not.allocated(ODEdiff%icoef)) then
	allocate(ODEdiff%icoef(MAX_VARS,7))
endif
if (.not.allocated(ODEdiff%iexcoef)) then
	allocate(ODEdiff%iexcoef(MAX_VARS,7))
endif
if (.not.allocated(ODEdiff%vartype)) then
	allocate(ODEdiff%vartype(MAX_VARS))
endif
if (.not.allocated(ODEdiff%cell_index)) then
	allocate(ODEdiff%cell_index(MAX_VARS))
endif
!if (.not.allocated(ODEdiff%isite_extra)) then
!	allocate(ODEdiff%isite_extra(MAX_VARS))
!endif
!if (.not.allocated(ODEdiff%isite_intra)) then
!	allocate(ODEdiff%isite_intra(MAX_VARS))
!endif
!if (.not.allocated(ODEdiff%extra_isite)) then
!	allocate(ODEdiff%extra_isite(MAX_VARS))
!endif


! Now %ivar = variable index, OUTSIDE_TAG, or UNREACHABLE_TAG

!ODEdiff%ivar = OUTSIDE_TAG
i = 0
!isite = 0
nex = 0
nin = 0
do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
		    kcell = occupancy(x,y,z)%indx(1)
			if (kcell > OUTSIDE_TAG) then
				i = i + 1
!				isite = isite + 1
				nex = nex + 1
				ODEdiff%ivar(x,y,z) = i			! extracellular variable index
				ODEdiff%varsite(i,:) = (/x,y,z/)
				ODEdiff%vartype(i) = EXTRA
				ODEdiff%cell_index(i) = kcell
!				ODEdiff%isite_extra(isite) = i	! extracellular index corresponding to site index isite
!				ODEdiff%extra_isite(i) = isite	! site index corresponding to extracellular index i
!				ODEdiff%isite_intra(isite) = 0	! extracellular index corresponding to site index isite
				if (kcell > 0) then 
					if (cell_list(kcell)%active) then
						i = i+1  ! with this numbering system the EXTRA variable corresponding to an INTRA variable i is i-1
						nin = nin + 1
    					ODEdiff%varsite(i,:) = (/x,y,z/)
	    				ODEdiff%vartype(i) = INTRA
    					ODEdiff%cell_index(i) = kcell
    					cell_list(kcell)%ivin = i		! (intracellular) variable index corresponding to kcell
!						ODEdiff%isite_intra(isite) = i	! intracellular index corresponding to site index isite
	    			endif
	    		endif
	    	else
	    		ODEdiff%ivar(x,y,z) = kcell
			endif
		enddo
	enddo
enddo
ODEdiff%nextra = nex
ODEdiff%nintra = nin
ODEdiff%nvars = i
if (.not.use_ODE_diffusion) return

! Set up mapping for EXTRA variables from variable index to EXTRA sequence number
allocate(exmap(ODEdiff%nvars))
kextra = 0
do i = 1,ODEdiff%nvars
	if (ODEdiff%vartype(i) == EXTRA) then
		kextra = kextra + 1
		exmap(i) = kextra
	endif
enddo

ierr = 0
do i = 1,ODEdiff%nvars
    if (ODEdiff%vartype(i) == INTRA) cycle
	site = ODEdiff%varsite(i,:)
	x = site(1)
	y = site(2)
	z = site(3)
	ODEdiff%icoef(i,1) = i

	if (x==1) then
		ODEdiff%icoef(i,2) = UNREACHABLE_TAG
	else
		ODEdiff%icoef(i,2) = ODEdiff%ivar(x-1,y,z)
	endif
	if (x==NX) then
		ODEdiff%icoef(i,3) = UNREACHABLE_TAG
	else
		ODEdiff%icoef(i,3) = ODEdiff%ivar(x+1,y,z)
	endif
	if (y==1) then
		ODEdiff%icoef(i,4) = UNREACHABLE_TAG
	else
		ODEdiff%icoef(i,4) = ODEdiff%ivar(x,y-1,z)
	endif
	if (y==NY) then
		ODEdiff%icoef(i,5) = UNREACHABLE_TAG
	else
		ODEdiff%icoef(i,5) = ODEdiff%ivar(x,y+1,z)
	endif
	if (z==1) then
		ODEdiff%icoef(i,6) = UNREACHABLE_TAG
	else
		ODEdiff%icoef(i,6) = ODEdiff%ivar(x,y,z-1)
	endif
	if (z==NZ) then
		ODEdiff%icoef(i,7) = UNREACHABLE_TAG
	else
		ODEdiff%icoef(i,7) = ODEdiff%ivar(x,y,z+1)
	endif
enddo

nchemo = 0
do ichemo = 1,MAX_CHEMO
!	if (chemo(ichemo)%used) then
	if (chemo(ichemo)%present) then
		nchemo = nchemo + 1
		chemomap(nchemo) = ichemo
	endif
enddo
deallocate(exmap)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_ivar(kcell,msg)
integer :: kcell
character*(*) :: msg
integer :: site(3), iv

site = cell_list(kcell)%site
iv = ODEdiff%ivar(site(1),site(2),site(3))
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
! Note: This has really been superceded by the option of a Hill function with N = 2.
!----------------------------------------------------------------------------------
subroutine AdjustMM
real(REAL_KIND) :: deltaC, C0, C1

C0 = chemo(OXYGEN)%MM_C0
if (MM_THRESHOLD > 0) then
	deltaC = MM_THRESHOLD
	C1 = deltaC + (sqrt(C0*C0 + 8*C0*deltaC) - C0)/4
	ODEdiff%k_soft = (C1-deltaC)/(C1*C1*(C0+C1-deltaC))
	ODEdiff%C1_soft = C1
	ODEdiff%deltaC_soft = deltaC
else
	ODEdiff%k_soft = 0
	ODEdiff%C1_soft = 0
	ODEdiff%deltaC_soft = 0
endif
!write(logmsg,'(a,4e12.4)') 'AdjustMM: C0, deltaC, C1, k: ',C0, ODEdiff%deltaC_soft, ODEdiff%C1_soft, ODEdiff%k_soft
!call logger(logmsg)
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
integer :: i, site(3), kcell, ichemo
logical :: ZERO_CONCS = .false.

do i = 1,ODEdiff%nvars
    site = ODEdiff%varsite(i,:)
    if (ODEdiff%vartype(i) == EXTRA) then
		if (ZERO_CONCS) then
			do ichemo = 1,MAX_CHEMO
				if (chemo(ichemo)%used) then
					if (chemo(ichemo)%present) then
						allstate(i,ichemo) = occupancy(site(1),site(2),site(3))%C(ichemo)
					else
						allstate(i,ichemo) = 0
					endif
				else
					allstate(i,ichemo) = 0
				endif
			enddo
		else
			allstate(i,:) = occupancy(site(1),site(2),site(3))%C(:)
		endif
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
!----------------------------------------------------------------------------------
subroutine InitialiseDrug(ichemo)
integer :: ichemo
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
! fixed extra- and intracellular volumes, Vextra_cm3 + Vcell_cm3 = Vsite_cm3
! When there is no cell the extracellular volume is Vsite_cm3.
!----------------------------------------------------------------------------------
subroutine f_rkc(neqn,t,y,dydt,icase)
integer :: neqn, icase
real(REAL_KIND) :: t, y(neqn), dydt(neqn)
integer :: i, k, ie, ki, kv, nextra, nintra, ichemo, idrug, im, site(3), kcell, ict, ith, Ng
real(REAL_KIND) :: dCsum, dCdiff, dCreact,  DX2, DX3, vol_cm3, val, Cin(MAX_CHEMO), Cex
real(REAL_KIND) :: decay_rate, dc1, dc6, cbnd, yy, C, membrane_kin, membrane_kout, membrane_flux, area_factor
logical :: bnd, dbug
!logical :: TPZ_metabolised(MAX_CELLTYPES,0:2), DNB_metabolised(MAX_CELLTYPES,0:2)
logical :: metabolised(MAX_CELLTYPES,0:2)
real(REAL_KIND) :: metab, dMdt, KmetC
logical :: intracellular, cell_exists
logical :: use_actual_cell_volume = .false.
type(drug_type), pointer :: dp

ichemo = icase
if (ichemo == GLUCOSE) then
	Ng = chemo(GLUCOSE)%Hill_N
endif
if (ichemo > TRACER) then
	idrug = (ichemo - TRACER - 1)/3 + 1
	im = ichemo - TRACER - 1 - 3*(idrug-1)		! 0 = drug, 1 = metab1, 2 = metab2
	dp => drug(idrug)
	metabolised(:,:) = (dp%Kmet0(:,:) > 0)
	if (im == 0) then
		write(*,*) 'f_rkc: ichemo,idrug,im: ',ichemo,idrug,im
		write(*,*) 'metabolised: ',metabolised(:,:)
	endif
else
	idrug = 0
	metabolised(:,:) = .false.
endif

DX2 = DELTA_X*DELTA_X
decay_rate = chemo(ichemo)%decay_rate
dc1 = chemo(ichemo)%diff_coef/DX2
dc6 = 6*dc1 + decay_rate
cbnd = BdryConc(ichemo,t_simulation)
membrane_kin = chemo(ichemo)%membrane_diff_in
membrane_kout = chemo(ichemo)%membrane_diff_out
!TPZ_metabolised(:,:) = (TPZ%Kmet0(:,:) > 0)	
!DNB_metabolised(:,:) = (DNB%Kmet0(:,:) > 0)	

!write(*,*) 'f_rkc: neqn: ',neqn
!$omp parallel do private(intracellular, vol_cm3, Cex, cell_exists, Cin, dCsum, k, kv, dCdiff, val, dCreact, yy, C, metab, kcell, ict, membrane_flux, KmetC, area_factor) default(shared) schedule(static)
do i = 1,neqn
	yy = y(i)
	if (isnan(yy)) then
		write(nflog,*) 'f_rkc: isnan: ',i,ichemo,yy
		write(*,*) 'f_rkc: isnan: ',i,ichemo,yy
		stop
	endif
    if (ODEdiff%vartype(i) == EXTRA) then
        intracellular = .false.
		vol_cm3 = Vsite_cm3
        Cex = yy
        cell_exists = .false.
        if (i < neqn) then
            if (ODEdiff%vartype(i+1) == INTRA) then
                cell_exists = .true.
				kcell = ODEdiff%cell_index(i+1)		! for access to cell-specific parameters (note that the intra variable follows the extra variable)
				if (use_actual_cell_volume) then
					vol_cm3 = Vsite_cm3 - Vcell_cm3*cell_list(kcell)%volume	! accounting for cell volume change
				else
	                vol_cm3 = Vextra_cm3									! for now, ignoring cell volume change!!!!!
	            endif
				area_factor = (cell_list(kcell)%volume)**(2./3.)
	            Cin = allstate(i+1,:)
	            Cin(ichemo) = y(i+1)
	        endif
	    endif
	else
        intracellular = .true.
	    kcell = ODEdiff%cell_index(i)					! for access to cell-specific parameters
	    if (use_actual_cell_volume) then
			vol_cm3 = Vcell_cm3*cell_list(kcell)%volume	! accounting for cell volume change
		else
			vol_cm3 = Vsite_cm3 - Vextra_cm3			! for now, ignoring cell volume change!!!!!
		endif
		area_factor = (cell_list(kcell)%volume)**(2./3.)
        Cex = y(i-1)
	    Cin = allstate(i,:)
	    Cin(ichemo) = yy
	    ict = cell_list(kcell)%celltype
	endif
	if (.not.intracellular) then
		! Need to check diffusion eqtn. when Vextra_cm3 < Vsite_cm3 = DX^3 !!!!!!!!!!!!!!!!!!!!!!
	    dCsum = 0
	    do k = 1,7
		    kv = ODEdiff%icoef(i,k)
		    if (k == 1) then
			    dCdiff = -dc6
		    else
			    dCdiff = dc1
		    endif
		    if (kv == OUTSIDE_TAG) then
			    val = cbnd
		    elseif (kv == UNREACHABLE_TAG) then
				val = y(ODEdiff%icoef(i,1))		! reflect concentration --> no flux
		    else
			    val = y(kv)
		    endif
		    dCsum = dCsum + dCdiff*val
	    enddo
	    if (cell_exists) then
!			membrane_flux = chemo(ichemo)%membrane_diff*(Cex - Cin(ichemo))	! We should account for change in surface area
			membrane_flux = area_factor*(membrane_kin*Cex - membrane_kout*Cin(ichemo))
!			if (ichemo == OXYGEN .and. kcell == 1) write(*,'(a,2e12.3)') 'O2, O2 flux: ',Cex,membrane_flux
		else
		    membrane_flux = 0
		endif
    	dydt(i) = dCsum - membrane_flux/vol_cm3
	else
		C = Cin(ichemo)
!		membrane_flux = chemo(ichemo)%membrane_diff*(Cex - C)		! We should account for change in surface area
		membrane_flux = area_factor*(membrane_kin*Cex - membrane_kout*C)
		dCreact = 0
		if (ichemo == OXYGEN) then
			metab = O2_metab(C)
			dCreact = (-metab*chemo(ichemo)%max_cell_rate + membrane_flux)/vol_cm3	! convert mass rate (mol/s) to concentration rate (mM/s)
		elseif (ichemo == GLUCOSE) then
			metab = C**Ng/(chemo(ichemo)%MM_C0**Ng + C**Ng)
			dCreact = (-metab*chemo(ichemo)%max_cell_rate + membrane_flux)/vol_cm3	! convert mass rate (mol/s) to concentration rate (mM/s)
		elseif (ichemo == TRACER) then
			dCreact = membrane_flux/vol_cm3
		elseif (im == 0) then
		    if (metabolised(ict,0) .and. C > 0) then
				KmetC = dp%Kmet0(ict,0)*C
				if (dp%Vmax(ict,0) > 0) then
					KmetC = KmetC + dp%Vmax(ict,0)*C/(dp%Km(ict,0) + C)
				endif
				dCreact = -(1 - dp%C2(ict,0) + dp%C2(ict,0)*dp%KO2(ict,0)/(dp%KO2(ict,0) + Cin(OXYGEN)))*KmetC
			endif
!			write(*,'(a,2i6,4e12.3)') 'dCreact: ',im,i,dCreact,membrane_flux/vol_cm3,yy,decay_rate
			dCreact = dCreact + membrane_flux/vol_cm3
		elseif (im == 1) then	! ichemo-1 is the PARENT drug
			if (metabolised(ict,0) .and. Cin(ichemo-1) > 0) then
				dCreact = (1 - dp%C2(ict,0) + dp%C2(ict,0)*dp%KO2(ict,0)/(dp%KO2(ict,0) + Cin(OXYGEN)))*dp%Kmet0(ict,0)*Cin(ichemo-1)
			endif
			if (metabolised(ict,1) .and. C > 0) then
				dCreact = dCreact - (1 - dp%C2(ict,1) + dp%C2(ict,1)*dp%KO2(ict,1)/(dp%KO2(ict,1) + Cin(OXYGEN)))*dp%Kmet0(ict,1)*C
			endif
!			write(*,'(a,2i6,2e12.3)') 'dCreact: ',im,i,dCreact,membrane_flux/vol_cm3
			dCreact = dCreact + membrane_flux/vol_cm3
		elseif (im == 2) then	! ichemo-1 is the METAB1
			if (metabolised(ict,1) .and. Cin(ichemo-1) > 0) then
				dCreact = (1 - dp%C2(ict,1) + dp%C2(ict,1)*dp%KO2(ict,1)/(dp%KO2(ict,1) + Cin(OXYGEN)))*dp%Kmet0(ict,1)*Cin(ichemo-1)
			endif
			if (metabolised(ict,2) .and. C > 0) then
				dCreact = dCreact - (1 - dp%C2(ict,2) + dp%C2(ict,2)*dp%KO2(ict,2)/(dp%KO2(ict,2) + Cin(OXYGEN)))*dp%Kmet0(ict,2)*C
			endif
!			write(*,'(a,2i6,2e12.3)') 'dCreact: ',im,i,dCreact,membrane_flux/vol_cm3
			dCreact = dCreact + membrane_flux/vol_cm3
		endif
!		select case (ichemo)
!		case (OXYGEN)
!			metab = O2_metab(C)
!			dCreact = (-metab*chemo(ichemo)%max_cell_rate + membrane_flux)/vol_cm3	! convert mass rate (mol/s) to concentration rate (mM/s)
!		case (GLUCOSE)
!			metab = C**Ng/(chemo(ichemo)%MM_C0**Ng + C**Ng)
!			dCreact = (-metab*chemo(ichemo)%max_cell_rate + membrane_flux)/vol_cm3	! convert mass rate (mol/s) to concentration rate (mM/s)
!		case (TRACER)
!			dCreact = membrane_flux/vol_cm3
!		case (TPZ_DRUG)
!		    if (TPZ_metabolised(ict,0) .and. C > 0) then
!				KmetC = TPZ%Kmet0(ict,0)*C
!				if (TPZ%Vmax(ict,0) > 0) then
!					KmetC = KmetC + TPZ%Vmax(ict,0)*C/(TPZ%Km(ict,0) + C)
!					KmetC = KmetC + DNB%Vmax(ict,0)*C/(DNB%Km(ict,0) + C)
!				endif
!				dCreact = -(1 - TPZ%C2(ict,0) + TPZ%C2(ict,0)*TPZ%KO2(ict,0)/(TPZ%KO2(ict,0) + Cin(OXYGEN)))*KmetC
!				dCreact = -(1 - DNB%C2(ict,0) + DNB%C2(ict,0)*DNB%KO2(ict,0)/(DNB%KO2(ict,0) + Cin(OXYGEN)))*KmetC
!			endif
!			dCreact = dCreact + membrane_flux/vol_cm3
!		case (TPZ_DRUG_METAB_1)
!			if (TPZ_metabolised(ict,0) .and. Cin(TPZ_DRUG) > 0) then
!				dCreact = (1 - TPZ%C2(ict,0) + TPZ%C2(ict,0)*TPZ%KO2(ict,0)/(TPZ%KO2(ict,0) + Cin(OXYGEN)))*TPZ%Kmet0(ict,0)*Cin(TPZ_DRUG)
!				dCreact = (1 - DNB%C2(ict,0) + DNB%C2(ict,0)*DNB%KO2(ict,0)/(DNB%KO2(ict,0) + Cin(OXYGEN)))*DNB%Kmet0(ict,0)*Cin(DNB_DRUG)
!			endif
!			if (TPZ_metabolised(ict,1) .and. C > 0) then
!				dCreact = dCreact - (1 - TPZ%C2(ict,1) + TPZ%C2(ict,1)*TPZ%KO2(ict,1)/(TPZ%KO2(ict,1) + Cin(OXYGEN)))*TPZ%Kmet0(ict,1)*C
!				dCreact = dCreact - (1 - DNB%C2(ict,1) + DNB%C2(ict,1)*DNB%KO2(ict,1)/(DNB%KO2(ict,1) + Cin(OXYGEN)))*DNB%Kmet0(ict,1)*C
!			endif
!			dCreact = dCreact + membrane_flux/vol_cm3
!		case (TPZ_DRUG_METAB_2)
!			if (TPZ_metabolised(ict,1) .and. Cin(TPZ_DRUG_METAB_1) > 0) then
!				dCreact = (1 - TPZ%C2(ict,1) + TPZ%C2(ict,1)*TPZ%KO2(ict,1)/(TPZ%KO2(ict,1) + Cin(OXYGEN)))*TPZ%Kmet0(ict,1)*Cin(TPZ_DRUG_METAB_1)
!				dCreact = (1 - DNB%C2(ict,1) + DNB%C2(ict,1)*DNB%KO2(ict,1)/(DNB%KO2(ict,1) + Cin(OXYGEN)))*DNB%Kmet0(ict,1)*Cin(DNB_DRUG_METAB_1)
!			endif
!			if (TPZ_metabolised(ict,2) .and. C > 0) then
!				dCreact = dCreact - (1 - TPZ%C2(ict,2) + TPZ%C2(ict,2)*TPZ%KO2(ict,2)/(TPZ%KO2(ict,2) + Cin(OXYGEN)))*TPZ%Kmet0(ict,2)*C
!				dCreact = dCreact - (1 - DNB%C2(ict,2) + DNB%C2(ict,2)*DNB%KO2(ict,2)/(DNB%KO2(ict,2) + Cin(OXYGEN)))*DNB%Kmet0(ict,2)*C
!			endif
!			dCreact = dCreact + membrane_flux/vol_cm3
!		case (DNB_DRUG)
!		    if (DNB_metabolised(ict,0) .and. C > 0) then
!				KmetC = DNB%Kmet0(ict,0)*C
!				if (DNB%Vmax(ict,0) > 0) then
!					KmetC = KmetC + DNB%Vmax(ict,0)*C/(DNB%Km(ict,0) + C)
!				endif
!				dCreact = -(1 - DNB%C2(ict,0) + DNB%C2(ict,0)*DNB%KO2(ict,0)/(DNB%KO2(ict,0) + Cin(OXYGEN)))*KmetC
!			endif
!			dCreact = dCreact + membrane_flux/vol_cm3
!		case (DNB_DRUG_METAB_1)
!			if (DNB_metabolised(ict,0) .and. Cin(DNB_DRUG) > 0) then
!				dCreact = (1 - DNB%C2(ict,0) + DNB%C2(ict,0)*DNB%KO2(ict,0)/(DNB%KO2(ict,0) + Cin(OXYGEN)))*DNB%Kmet0(ict,0)*Cin(DNB_DRUG)
!			endif
!			if (DNB_metabolised(ict,1) .and. C > 0) then
!				dCreact = dCreact - (1 - DNB%C2(ict,1) + DNB%C2(ict,1)*DNB%KO2(ict,1)/(DNB%KO2(ict,1) + Cin(OXYGEN)))*DNB%Kmet0(ict,1)*C
!			endif
!			dCreact = dCreact + membrane_flux/vol_cm3
!		case (DNB_DRUG_METAB_2)
!			if (DNB_metabolised(ict,1) .and. Cin(DNB_DRUG_METAB_1) > 0) then
!				dCreact = (1 - DNB%C2(ict,1) + DNB%C2(ict,1)*DNB%KO2(ict,1)/(DNB%KO2(ict,1) + Cin(OXYGEN)))*DNB%Kmet0(ict,1)*Cin(DNB_DRUG_METAB_1)
!			endif
!			if (DNB_metabolised(ict,2) .and. C > 0) then
!				dCreact = dCreact - (1 - DNB%C2(ict,2) + DNB%C2(ict,2)*DNB%KO2(ict,2)/(DNB%KO2(ict,2) + Cin(OXYGEN)))*DNB%Kmet0(ict,2)*C
!			endif
!			dCreact = dCreact + membrane_flux/vol_cm3
!		end select
	    dydt(i) = dCreact - yy*decay_rate
		if (isnan(dydt(i))) then
			write(nflog,*) 'f_rkc: dydt isnan: ',i,ichemo,dydt(i)
			write(*,*) 'f_rkc: dydt isnan: ',i,ichemo,dydt(i)
			stop
		endif
	endif
enddo
!$omp end parallel do
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine ShowResults(v)
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
!write(logmsg,'(a,2f7.4)') 'Medium: ',chemo(DRUG_A)%bdry_conc,chemo(DRUG_METAB_A)%bdry_conc
!call logger(logmsg)
write(logmsg,'(a,i6,10f7.4)') 'E:',ODEdiff%nextra,(v(ODEdiff%ivar(x,y,z)),z=z1,z2)
call logger(logmsg)
write(logmsg,'(a,i6,10f7.4)') 'I:',ODEdiff%nintra,(Cin(z),z=z1,z2)
call logger(logmsg)
!write(logmsg,'(a,i6,10f7.4)') 'M:',ODEdiff%nintra,(drugM(z),z=z1,z2)
!call logger(logmsg)
!write(logmsg,'(a,10f7.4)') 'Diff:   ',((v(ODEdiff%ivar(x,y,z))-Cin(z))/Cin(z),z=z1,z2)
!call logger(logmsg)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_allstate(msg)
character*(*) :: msg
integer :: i

if (chemo(DRUG_A)%present) return
do i = 1,MAX_VARS
	if (allstate(i,DRUG_A) /= 0) then
		stop
	endif
enddo
end subroutine


!----------------------------------------------------------------------------------
! In this version the diffusion/decay of each constituent is solved by a separate
! OMP thread.  Obviously this requires at least as many CPUs as there are constituents.
! THIS HAS BEEN CHANGED
! Note that this required modifications to the way the ODE solver handles SAVEd variables,
! to avoid collisions between different threads.
! The ODE solver is RKC
! The subroutine assumes that:
!   * ODEdiff has been set up with the correct mappings
!   * allstate(:,:) holds the most recent solution, including estimates when the
!     blob changes (i.e. grows or shrinks).  This could entail variable renumbering.
!   * work(:,:) is correctly sized (ODEdiff%nextra)
! Currently the relaxation method is used only for OXYGEN
!----------------------------------------------------------------------------------
subroutine Solver(it,tstart,dt,nc)
integer :: it, nc
real(REAL_KIND) :: tstart, dt
integer :: ichemo, nvars, ntvars, ic, kcell, site(3), iv, nth
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

!if (RK_SOLVER == IRKC_SOLVE) then	! The performance is very poor
!	call irkc_solver(it,tstart,dt,nc)
!	return
!endif

ODEdiff%nintra = nc
ntvars = ODEdiff%nextra + ODEdiff%nintra
if (EXPLICIT_INTRA) then
    nvars = ODEdiff%nextra
else
    nvars = ntvars
endif
allocate(state(ntvars,MAX_CHEMO))
state(:,:) = allstate(1:ntvars,1:MAX_CHEMO)

info(1) = 1
info(2) = 1		! = 1 => use spcrad() to estimate spectral radius, != 1 => let rkc do it
info(3) = 1
info(4) = 0
rtol = 1d-2
atol = rtol

do ic = 1,nchemo
	ichemo = chemomap(ic)
	if (relax .and. ichemo == OXYGEN) cycle
	if (chemo(ichemo)%constant) cycle
	
	idid = 0
	t = tstart
	tend = t + dt
	call rkc(comm_rkc(ichemo),nvars,f_rkc,state(:,ichemo),t,tend,rtol,atol,info,work_rkc,idid,ichemo)
	if (idid /= 1) then
		write(logmsg,*) ' Failed at t = ',t,' with idid = ',idid
		call logger(logmsg)
		stop
	endif
enddo

ichemo = OXYGEN
if (relax .and. .not.chemo(ichemo)%constant) then
	if (use_parallel) then
	    call ParRelaxSolver(ichemo,state(:,ichemo))
	else
	    call RelaxSolver(ichemo,state(:,ichemo))
	endif
endif

allstate(1:nvars,1:MAX_CHEMO) = state(:,:)
! Note: some time we need to copy the state values to the cell_list() array.
deallocate(state)

end subroutine

!----------------------------------------------------------------------------------
! Solve for ichemo using the over- under-relaxation method to solve the Poisson's
! equation that corresponds to the steady state equation with all other concentrations
! held constant.
! This was suggested by Tim Secomb.
! dC/dt = K.grad^2C - f(C,...) = 0
! where:
!   C = C(i,j,k) = concentration of the constituent at the grid-cell we are solving for,
!   K = diffusion coeff.
!   f(C,...) = rate of consumption of the constituent as a function of C and other 
!   constituents, all held fixed.
! Tim's suggested procedure:
! (1) In the first iteration round the concentration of C used in f(C,...) is held
!     fixed at the "old" value.  The "diff" values are updated by a number of iterations
!     (e.g. up to 20) of the over-relaxation procedure (w_over = 1.9)
! (2) The "old" C values are now replaced by the "under-relaxed" values:
!     Cnew = Cold + w_under*(Cdiff - Cold)
! The cycle is repeated until convergence.
! The situation is complicated by the fact that the constituent values are normally
! stored in an array that has extra- and intracellular concentrations interleaved.
! This will need to be addressed later to improve execution speed.
! For now we need a mapping derived from ODEdiff%icoef(:,:) to convey the
! extracellular array indices of the neighbours of extracellular variable ie.
! This currently applies only to ichemo = OXYGEN
! (no decay, and because its diffusivity is much greater than all other constituents,
! O2 equilibrates rapidly)
! NOTE: Solving only for EXTRA concentrations (INTRA are inferred from these with getCin())
! NOTE: Currently used only for OXYGEN
!----------------------------------------------------------------------------------
subroutine RelaxSolver(ichemo,y)
integer :: ichemo
real(REAL_KIND) :: y(:)
integer :: nvar, ie, ia, k, je, ja, k_under, k_over, n, site(3)
integer, allocatable :: all_index(:), extra_index(:), icoef(:,:)
real(REAL_KIND), allocatable :: y0(:), ydiff(:), uptake(:)
real(REAL_KIND) :: DX2, Kdiff, Csum, dCreact, val, cbnd, esum2, sum, Cnew, y0temp, dy, ymin
real(REAL_KIND), parameter :: w_over = 1.6, w_under = 0.05
real(REAL_KIND), parameter :: tol1_over = 1.0e-5, tol1_under = 1.0e-6, tol2 = 1.0e-10
integer, parameter :: n_over = 1000, n_under = 1000
!integer, static :: iemin = 500
integer :: iemin = 500

DX2 = DELTA_X*DELTA_X
nvar = ODEdiff%nextra
Kdiff = chemo(ichemo)%diff_coef
cbnd = BdryConc(ichemo,t_simulation)

! Allocate y0(:), ydiff(:) and copy y(:) to y0(:), ydiff(:)
allocate(y0(nvar))
allocate(ydiff(nvar))
allocate(uptake(nvar))
allocate(extra_index(ODEdiff%nvars))
allocate(all_index(nvar))
allocate(icoef(nvar,6))
! This recapitulates SiteCellToState()
ie = 0
do ia = 1,ODEdiff%nvars
    if (ODEdiff%vartype(ia) == EXTRA) then	! Copy EXTRA variable values into y0()
        ie = ie+1
        y0(ie) = y(ia)
        extra_index(ia) = ie
        all_index(ie) = ia
        site = ODEdiff%varsite(ie,:)
    endif
enddo
ydiff = y0
! Set up icoef(:,:)
do ie = 1,nvar
    ia = all_index(ie)
    do k = 1,6
        ja = ODEdiff%icoef(ia,k+1)
        if (ja < 0) then
            je = ja
        else
            je = extra_index(ja)
        endif
        icoef(ie,k) = je
    enddo
enddo

! Loop over under-relaxation until convergence

do k_under = 1,n_under

    do ie = 1,nvar
        uptake(ie) = UptakeRate(ichemo,y0(ie))     ! rate of decrease of constituent concentration by reactions
    enddo

    ! Loop over diffusion by over-relaxation until convergence
    do k_over = 1,n_over
		sum = 0
		n = 0
        do ie = 1,nvar
	        Csum = 0
	        do k = 1,6
		        je = icoef(ie,k)
		        if (je == OUTSIDE_TAG) then
			        val = cbnd
		        elseif (je == UNREACHABLE_TAG) then
					val = ydiff(ie)
		        else
			        val = ydiff(je)
		        endif
		        Csum = Csum + val
	        enddo
            dCreact = uptake(ie)

            Cnew = (Csum - DX2*dCreact/Kdiff)/6
            val = ydiff(ie)
            ydiff(ie) = (1-w_over)*ydiff(ie)+ w_over*Cnew
            ydiff(ie) = max(ydiff(ie),0.000001)
            if (ydiff(ie) < 0) then
				write(nflog,*) 'ydiff < 0'
				stop
			endif
			dy = abs(val-ydiff(ie))/ydiff(ie)
			sum = sum + dy
			if (dy > tol1_over) then
				n = n+1
			endif
        enddo
        if (n == 0) exit
    enddo
    esum2 = 0
    sum = 0
    n = 0
    do ie = 1,nvar
        val = y0(ie)
        y0(ie) = y0(ie) + w_under*(ydiff(ie) - y0(ie))
        if (y0(ie) < 0) then
			write(logmsg,*) 'y0 < 0'
			call logger(logmsg)
			stop 
		endif
		dy = abs(val-y0(ie))/y0(ie)
!		sum = sum + dy
		if (dy > tol1_under) n = n+1
    enddo
    if (n == 0) exit
!    if (sum/nvar < tol2) exit
enddo
!call write_solution('end',y0,nvar)
!ysave(1:nvar) = y0(1:nvar)

! Copy y0(:) back to y(:)
ymin = 1.0e10
ie = 0
do k = 1,ODEdiff%nvars
    if (ODEdiff%vartype(k) == EXTRA) then
        ie = ie+1
        y(k) = y0(ie)
        if (y0(ie) < ymin) then
			ymin = y0(ie)
			iemin = ie
		endif
!		if (ODEdiff%cell_index(k+1) == idbug) then
!			write(nfout,*) 'k, idbug: ',k,idbug
!		endif

    else
!        y(k) = y(k-1)	! put Cin here
		y(k) = getCin(ichemo,y(k-1))
!		if (ODEdiff%cell_index(k) == idbug) then
!			write(nfout,'(a,2f10.4)') 'RelaxSolver: ',y(k-1),y(k)
!	        site = ODEdiff%varsite(ie,:)
!	        write(nfout,*) 'site: ',site
!		endif
!		y(k) = getCinO2(y(k-1))
    endif
enddo

deallocate(y0)
deallocate(ydiff)
deallocate(uptake)
deallocate(extra_index)
deallocate(all_index)
deallocate(icoef)

!write(nfout,'(a,2i6,e12.4)') 'Min y: ',istep,iemin,ymin

end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine ParRelaxSolver(ichemo,y)
integer :: ichemo
real(REAL_KIND) :: y(:)
integer :: nexvar, ie, ia, k, je, ja, k_under, k_over
integer :: site(3), x, z, i, nt, ie1, ie2, n1, n2, nrange, ierange(100,2), n(100), nn
integer, allocatable :: all_index(:), extra_index(:), icoef(:,:)
real(REAL_KIND), allocatable :: y0(:), ydiff(:), uptake(:)
real(REAL_KIND) :: DX2, Kdiff, Csum, dCreact, val, cbnd, esum2, sum, Cnew, y0temp, dy, ymin
real(REAL_KIND), parameter :: w_over = 1.7, w_under = 0.05
real(REAL_KIND), parameter :: tol1_over = 1.0e-5, tol1_under = 1.0e-6, tol2 = 1.0e-10
integer, parameter :: n_over = 1000, n_under = 1000
!integer, static :: iemin = 500
integer :: iemin = 500
logical, parameter :: interleave = .false.

DX2 = DELTA_X*DELTA_X
nexvar = ODEdiff%nextra
Kdiff = chemo(ichemo)%diff_coef
cbnd = BdryConc(ichemo,t_simulation)

! Allocate y0(:), ydiff(:) and copy y(:) to y0(:), ydiff(:)
allocate(y0(nexvar))
allocate(ydiff(nexvar))
allocate(uptake(nexvar))
allocate(extra_index(ODEdiff%nvars))
allocate(all_index(nexvar))
allocate(icoef(nexvar,6))
! This recapitulates SiteCellToState()
x = 0
z = 0
nrange = 0
ie = 0
do ia = 1,ODEdiff%nvars
    if (ODEdiff%vartype(ia) == EXTRA) then
        ie = ie+1
        y0(ie) = y(ia)
        extra_index(ia) = ie
        all_index(ie) = ia
        site = ODEdiff%varsite(ie,:)
        if (site(1) /= x) then
			x = site(1)
!        if (site(3) /= z) then
!			z = site(3)
			nrange = nrange + 1
			if (nrange > 1) then
				ierange(nrange-1,2) = ie-1
			endif
			ierange(nrange,1) = ie
		endif
    endif
enddo
ierange(nrange,2) = ie
n2 = nrange/2
n2 = 2*n2
if (n2 == nrange) then
	n1 = nrange - 1
else
	n1 = nrange
endif

ydiff = y0

! Set up icoef(:,:)
do ie = 1,nexvar
    ia = all_index(ie)
    do k = 1,6
        ja = ODEdiff%icoef(ia,k+1)
        if (ja < 0) then
            je = ja
        else
            je = extra_index(ja)
        endif
        icoef(ie,k) = je
    enddo
enddo

! Loop over under-relaxation until convergence
do k_under = 1,n_under
    do ie = 1,nexvar
        uptake(ie) = UptakeRate(ichemo,y0(ie))     ! rate of decrease of intracellular constituent concentration by reactions
    enddo
    ! Loop over diffusion by over-relaxation until convergence
    do k_over = 1,n_over
		
		if (interleave) then
! Interleave the sweeps to reduce memory access delays
!$omp parallel do private(ie1,ie2) default(shared) schedule(static)
		do i = 1,n1,2
			ie1 = ierange(i,1)
			ie2 = ierange(i,2)
			call PlaneSolver(ie1,ie2,icoef,uptake,cbnd,Kdiff,ydiff,n(i))
		enddo
!$omp end parallel do

!$omp parallel do private(ie1,ie2) default(shared) schedule(static)
		do i = 2,n2,2
			ie1 = ierange(i,1)
			ie2 = ierange(i,2)
			call PlaneSolver(ie1,ie2,icoef,uptake,cbnd,Kdiff,ydiff,n(i))
		enddo
!$omp end parallel do

		else
!$omp parallel do private(ie1,ie2) !default(shared) schedule(dynamic) !schedule(static)
		do i = 1,nrange
			ie1 = ierange(i,1)
			ie2 = ierange(i,2)
			call PlaneSolver(ie1,ie2,icoef,uptake,cbnd,Kdiff,ydiff,n(i))
		enddo
!$omp end parallel do
		endif

		nt = sum(n(1:nrange))
        if (nt == 0) exit
    enddo
!    sum = 0
    nt = 0
    do ie = 1,nexvar
        val = y0(ie)
        y0(ie) = y0(ie) + w_under*(ydiff(ie) - y0(ie))	! under-relax the extracellular solution
        if (y0(ie) < 0) then
			write(logmsg,*) 'y0 < 0'
			call logger(logmsg)
			stop
	    endif
		dy = abs(val-y0(ie))/y0(ie)
!		sum = sum + dy
        if (dy > tol1_under) nt = nt+1
    enddo
    if (nt == 0) exit
!    if (sum/nexvar < tol2) exit
enddo

!call write_solution('end',y0,nexvar)
!ysave(1:nexvar) = y0(1:nexvar)

! Copy y0(:) back to y(:)
ymin = 1.0e10
ie = 0
do k = 1,ODEdiff%nvars
    if (ODEdiff%vartype(k) == EXTRA) then
        ie = ie+1
        y(k) = y0(ie)
        if (y0(ie) < ymin) then
			ymin = y0(ie)
			iemin = ie
		endif
    else
!        y(k) = y(k-1)
		y(k) = getCin(ichemo,y(k-1))
!		y(k) = getCinO2(y(k-1))
    endif
enddo

deallocate(y0)
deallocate(ydiff)
deallocate(uptake)
deallocate(extra_index)
deallocate(all_index)
deallocate(icoef)

!write(nfout,'(a,2i6,e12.4)') 'Min y: ',istep,iemin,ymin

end subroutine

!----------------------------------------------------------------------------------
! No decay!!!
!----------------------------------------------------------------------------------
subroutine PlaneSolver(ie1,ie2,icoef,uptake,cbnd,Kdiff,ydiff,n)
integer :: ie1, ie2, n
integer :: icoef(:,:)
real(REAL_KIND) :: cbnd, Kdiff, uptake(:), ydiff(:)
integer :: ie, k, je
real(REAL_KIND) :: DX2, Csum, dCreact, val, Cnew, dy
real(REAL_KIND), parameter :: w_over = 1.7
real(REAL_KIND), parameter :: tol1_over = 1.0e-5

DX2 = DELTA_X*DELTA_X
n = 0
do ie = ie1,ie2
    Csum = 0
    do k = 1,6
        je = icoef(ie,k)
        if (je == OUTSIDE_TAG) then
	        val = cbnd
        elseif (je == UNREACHABLE_TAG) then
			val = ydiff(ie)
        else
	        val = ydiff(je)
        endif
        Csum = Csum + val
    enddo
    dCreact = uptake(ie)
    Cnew = (Csum - DX2*dCreact/Kdiff)/6
    val = ydiff(ie)
    ydiff(ie) = (1-w_over)*ydiff(ie)+ w_over*Cnew
    ydiff(ie) = max(ydiff(ie),0.000001)
    if (ydiff(ie) < 0) then
		write(nflog,*) 'ydiff < 0'
		stop
	endif
	dy = abs(val-ydiff(ie))/ydiff(ie)
	if (dy > tol1_over) then
		n = n+1
	endif
enddo
end subroutine

!----------------------------------------------------------------------------------
! Note: This computes a rate of change of concentration! mM/s
! Currently only for O2!!!
! There are two options: use_Cex_Cin = true/false
!
! use_Cex_Cin = true
! ------------------
! The idea is that the speed of the intracellular reactions, compared with other
! processes, is so fast that effectively the intracellular concentration is always
! in equilibrium with the extracellular value.  This means that the rate of consumption
! in the cell matches the rate of transport across the cell membrane: both these rates 
! depend on Cin, therefore we can solve for Cin given Cex then deduce uptake rate
!
! use_Cex_Cin = false
! -------------------
! In this case we just use Cin = Cex to calculate the consumption rate - no
! dependence on chemo(OXYGEN)%membrane_diff
!----------------------------------------------------------------------------------
real(REAL_KIND) function UptakeRate(ichemo,Cex)
integer :: ichemo
real(REAL_KIND) :: Cex
real(REAL_KIND) :: vol, K1, Cin, flux
integer :: n, i

if (ichemo == OXYGEN) then
!    vol = Vsite_cm3
!    vol = Vsite_cm3 - Vextra_cm3	! this was used in the RKC solution
    vol = Vextra_cm3	! the current extracellular volume should be used I think !!!!!!!!!!!!!!!
	if (use_Cex_Cin) then
		Cin = getCin(ichemo,Cex)
!		flux = chemo(ichemo)%membrane_diff*(Cex - Cin)
		flux = (chemo(ichemo)%membrane_diff_in*Cex - chemo(ichemo)%membrane_diff_out*Cin)
	else	! 
		flux = O2_metab(Cex)*chemo(ichemo)%max_cell_rate
	endif
	if (dbug) write(nfout,'(a,2e12.4)') 'Cex, flux: ',Cex,flux
	UptakeRate = flux/vol	! concentration rate (mM/s)
else
	write(logmsg,*) 'ERROR: UptakeRate: currently only for OXYGEN'
	call logger(logmsg)
	stop
endif
end function

!----------------------------------------------------------------------------------
! Computes intracellular O2 concentration as a function of the extracellular level C,
! assuming equilibrium, i.e. rate of consumption = rate of membrane transport.
! Note that the cell's O2 uptake rate is taken to be independent of any other factors,
! e.g. independent of cell size.
! NOTE: Currently only for OXYGEN - OK because membrane_diff_in = membrane_diff_out
!----------------------------------------------------------------------------------
!real(REAL_KIND) function getCinO2(C)
real(REAL_KIND) function getCin(ichemo,C)
integer :: ichemo
real(REAL_KIND) :: C
real(REAL_KIND) :: K1, K2, K2K1, C0, a, b, cc, D, r(3), Cin
integer :: i, n

if (ichemo /= OXYGEN) then
	write(logmsg,*) 'ERROR: getCin: currently only for OXYGEN'
	call logger(logmsg)
	stop
endif
!ichemo = OXYGEN
!K1 = chemo(OXYGEN)%membrane_diff*(Vsite_cm3 - Vextra_cm3)
K1 = chemo(ichemo)%membrane_diff_in
K2 = chemo(ichemo)%max_cell_rate
K2K1 = K2/K1
C0 = chemo(ichemo)%MM_C0
if (chemo(ichemo)%Hill_N == 2) then
	a = K2K1 - C
	b = C0*C0
	cc = -b*C
	call cubic_roots(a,b,cc,r,n)
	if (n == 1) then
		Cin = r(1)
	else
		n = 0
		do i = 1,3
			if (r(i) > 0) then
				n = n+1
				Cin = r(i)
			endif
		enddo
		if (n > 1) then
			write(nflog,*) 'getCin: two roots > 0: ',r
			stop
		endif
	endif
elseif (chemo(ichemo)%Hill_N == 1) then
	b = K2K1 + C0 - C
	cc = -C0*C
	D = sqrt(b*b - 4*cc)
	Cin = (D - b)/2
endif
getCin = Cin
end function

!----------------------------------------------------------------------------------
! Computes metabolism rate as a fraction of the maximum cell rate
! Use the "soft landing" option for Hill_N = 1 if MM_threshold = 0
!----------------------------------------------------------------------------------
real(REAL_KIND) function O2_metab(C)
!real(REAL_KIND) function metabolic_rate(ichemo,C)
integer :: ichemo
real(REAL_KIND) :: C
real(REAL_KIND) :: metab

ichemo = OXYGEN
if (ichemo == OXYGEN) then
	if (chemo(ichemo)%Hill_N == 2) then
		if (C > 0) then
			metab = C*C/(chemo(ichemo)%MM_C0*chemo(ichemo)%MM_C0 + C*C)
		else
			metab = 0
		endif
	else
		if (MM_THRESHOLD > 0) then
			if (C > ODEdiff%C1_soft) then
				metab = (C-ODEdiff%deltaC_soft)/(chemo(ichemo)%MM_C0 + C - ODEdiff%deltaC_soft)
			elseif (C > 0) then
				metab = ODEdiff%k_soft*C*C
			else
				metab = 0
			endif
		else
			if (C > 0) then
				metab = C/(chemo(ichemo)%MM_C0 + C)
			else
				metab = 0
			endif
		endif
	endif
endif
O2_metab = metab
!metabolic_rate = metab
end function

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine WriteSolution(msg,y,n)
character*(*) :: msg
real(REAL_KIND) :: y(:)
integer :: n

write(nfout,'(a)') msg
write(nfout,'(10f7.4)') y(1:n)
end subroutine

!----------------------------------------------------------------------------------
! Update Cbnd using current M, R1 and previous U, Cext
! No longer needed, with revised UpdateMedium
!----------------------------------------------------------------------------------
subroutine UpdateCbnd
integer :: ichemo
real(REAL_KIND) :: R1, R2, a, Vblob, Vm

return

call SetRadius(Nsites)
R1 = Radius*DELTA_X		! cm
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%present) cycle
	if (chemo(ichemo)%constant) then
		chemo(ichemo)%medium_Cbnd = chemo(ichemo)%bdry_conc
	else
		R2 = R1 + chemo(ichemo)%medium_dlayer
		a = (1/R2 - 1/R1)/(4*PI*chemo(ichemo)%medium_diff_coef)
		chemo(ichemo)%medium_Cbnd = chemo(ichemo)%medium_Cext + a*chemo(ichemo)%medium_U
	endif
	if (ichemo == OXYGEN .and. chemo(ichemo)%medium_Cbnd < 0) then
		write(logmsg,'(a,2e12.3,a,e12.3)') 'UpdateCbnd: O2 < 0: Cext: ',chemo(ichemo)%medium_Cbnd,chemo(ichemo)%medium_Cext,' U: ',chemo(ichemo)%medium_U
		call logger(logmsg)
		stop
	endif
enddo
Vblob = (4./3.)*PI*R1**3	! cm3
Vm = total_volume - Vblob
!write(*,'(a,4e12.3)') 'UpdateCbnd: ext glucose conc, mass: ', &
!chemo(GLUCOSE)%medium_Cext,chemo(GLUCOSE)%medium_Cbnd,chemo(GLUCOSE)%medium_Cext*Vm,chemo(GLUCOSE)%medium_M
!write(nflog,'(a,i6,2f10.6)') 'UpdateCbnd: istep,R1,R2: ',istep,R1,R2
!write(nflog,'(a,4e12.3)') 'UpdateCbnd: medium_Cext,medium_U: ',chemo(OXYGEN)%medium_Cext,chemo(GLUCOSE)%medium_Cext, &
!																chemo(OXYGEN)%medium_U,chemo(GLUCOSE)%medium_U
!write(nflog,'(a,2e12.3)') 'UpdateCbnd: medium_Cbnd: ',chemo(OXYGEN)%medium_Cbnd,chemo(GLUCOSE)%medium_Cbnd
end subroutine

!----------------------------------------------------------------------------------
! The medium concentrations are updated explicitly, assuming a sphere with 
! known total uptake rate U.
! Compute U and update M, Cext, Cbnd.
! Note that for O2 Cext is fixed.
! Need to include decay
! PROBLEM: For glucose, U ia about 20x the total cell uptake flux
! !! Nbnd > Ncells, but r resulting Csum/Nbnd is  reasonable estimate of the conc.
!----------------------------------------------------------------------------------
subroutine UpdateMedium1(dt)
real(REAL_KIND) :: dt
integer :: i, k, ichemo, ntvars
integer :: ichemo_log
real(REAL_KIND) :: dA, R1, R2, V0, Csum(MAX_CHEMO), U(MAX_CHEMO), tracer_C, tracer_N
real(REAL_KIND) :: a(MAX_CHEMO), b(MAX_CHEMO), Rlayer(MAX_CHEMO)
logical :: bnd

ichemo_log = OXYGEN
! First need the spheroid radius
call SetRadius(Nsites)
dA = DELTA_X*DELTA_X	! cm2
R1 = Radius*DELTA_X		! cm
V0 = total_volume		! cm3
!V0 = medium_volume0		! cm3
Rlayer(:) = R1 + chemo(:)%medium_dlayer
b = dA*chemo(:)%medium_diff_coef/DELTA_X
a = (1/Rlayer(:) - 1/R1)/(4*PI*chemo(:)%medium_diff_coef)
U = 0
! This could/should be done using sites in bdrylist
ntvars = ODEdiff%nextra + ODEdiff%nintra
tracer_C = 0
tracer_N = 0
Nbnd = 0
Csum = 0
do i = 1,ntvars
	if (ODEdiff%vartype(i) /= EXTRA) cycle
	bnd = .false.
	do k = 1,7
		if (ODEdiff%icoef(i,k) < 0) then
			bnd = .true.
			exit
		endif
	enddo
	if (.not.bnd) cycle
	do k = 1,7
		if (ODEdiff%icoef(i,k) < 0) then	! boundary with medium ???????????????????????????????????????????????????????????????????
			Nbnd = Nbnd + 1
			do ichemo = 1,MAX_CHEMO
				if (.not.chemo(ichemo)%present) cycle
!				U(ichemo) = U(ichemo) + dA*chemo(ichemo)%diff_coef*(chemo(ichemo)%medium_Cbnd - allstate(i,ichemo))/DELTA_X
				Csum(ichemo) = Csum(ichemo) + allstate(i,ichemo)
				if (ichemo == TRACER) then
					tracer_C = tracer_C + allstate(i,ichemo)
					tracer_N = tracer_N + 1
				endif
			enddo
		endif
	enddo
enddo

!!!U = (dA*chemo(:)%diff_coef/DELTA_X)*(Nbnd*chemo(:)%medium_Cbnd - Csum(:))	! Note: this is an approximation

chemo(:)%medium_Cbnd = (chemo(:)%medium_Cext - Csum(:)*b*a(:))/(1 - Nbnd*b*a(:))

!write(*,'(a,2i6,2e12.3)') 'UpdateMedium: glucose Csum/Nbnd: ',Nbnd,Ncells,Csum(GLUCOSE)/Nbnd,chemo(GLUCOSE)%medium_Cbnd
U(:) = (chemo(:)%medium_Cbnd - chemo(:)%medium_Cext)/a(:)

!write(nflog,'(a,2i6,4e12.3)') 'updateMedium: ',istep,Nbnd,chemo(ichemo_log)%medium_Cext,Csum(ichemo_log),U(ichemo_log),chemo(ichemo_log)%medium_Cbnd

do ichemo = 1,MAX_CHEMO
!	if (.not.chemo(ichemo)%used) cycle
	if (.not.chemo(ichemo)%present) cycle
	if (chemo(ichemo)%constant) then
		chemo(ichemo)%medium_Cext = chemo(ichemo)%bdry_conc
	else
		R2 = Rlayer(ichemo)
		if (ichemo /= OXYGEN) then
			chemo(ichemo)%medium_M = chemo(ichemo)%medium_M*(1 - chemo(ichemo)%decay_rate*dt) - U(ichemo)*dt
			chemo(ichemo)%medium_Cext = (chemo(ichemo)%medium_M - (U(ichemo)/(6*chemo(ichemo)%medium_diff_coef)*R2) &
				*(R1*R1*(3*R2 - 2*R1) - R2*R2*R2))/(V0 - 4*PI*R1*R1*R1/3.)
		endif
	endif
enddo	
!!!U = b(:)*(Nbnd*chemo(:)%medium_Cext - Csum(:))/(1 - b(:)*Nbnd*a(:))
!!!!chemo(:)%medium_Cbnd = chemo(:)%medium_Cext + (U(:)/(4*PI*chemo(:)%medium_diff_coef))*(1/Rlayer(:) - 1/R1)
chemo(:)%medium_U = U(:)
!write(nflog,'(a,i6,2f10.6)') 'UpdateMedium: istep,R1,R2: ',istep,R1,R2
!write(nflog,'(a,4e12.3)') 'UpdateMedium: medium_Cext,medium_U: ',chemo(OXYGEN)%medium_Cext,chemo(GLUCOSE)%medium_Cext, &
!																chemo(OXYGEN)%medium_U,chemo(GLUCOSE)%medium_U
!write(nflog,'(a,2e12.3)') 'UpdateMedium: medium_Cbnd: ',chemo(OXYGEN)%medium_Cbnd,chemo(GLUCOSE)%medium_Cbnd
!write(*,'(a,2e12.3)') 'UpdateMedium: glucose U: ',U(GLUCOSE),U(GLUCOSE)*dt
!!!!!!!!! Problem: U is about 20x sum_dMdt
end subroutine

!----------------------------------------------------------------------------------
! The medium concentrations are updated explicitly, assuming a sphere with 
! known total uptake rate U.
! Compute U and update M, Cext, Cbnd.
! Note that for O2 Cext is fixed.
! Need to include decay
! Need to estimate U in a different way.
! Use the total cell uptake and decay within the blob.
! The change in total mass in the medium M is minus the sum of decay and U.dt
! The far-field concentration Cext is inferred from M.
! If the blob radius is R1 and the radius of the boundary layer is R2 = R1 + d_layer,
! then the concentration distribution in the layer at radius r is:
! C(r) = Cext + U/(4.pi.K).(1/R2 - 1/r)
! where K = diffusion coefficient in the unstirred layer
! By integrating C(r) and adding the portion at Cext we find the total mass is:
! M = Cext.Vmedium + (U/K)[(R2^3-R1^3)/(3R2) + (R2^2 - R1^2)/2)]
! which can be solved for Cext when U, M, R1 and R2 are known.
! Note that Vmedium = total_volume - Vblob = total_volume - (4/3)pi.R1^3
! From Cext, Cbnd = Cext + U/(4.pi.K).(1/R2 - 1/R1)
!----------------------------------------------------------------------------------
subroutine UpdateMedium(dt)
real(REAL_KIND) :: dt
integer :: i, k, ichemo, ntvars
real(REAL_KIND) :: R1, R2, U(MAX_CHEMO)
real(REAL_KIND) :: a, Rlayer(MAX_CHEMO)
integer :: kcell, Nh, Nc
real(REAL_KIND) :: C, metab, dMdt, asum
real(REAL_KIND) :: Kin, Kout, decay_rate, vol_cm3, Cin, Cex

!write(*,*) 'UpdateMedium'
! Start by looking at a conservative constituent (e.g. glucose)
! Contribution from cell uptake
U = 0
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%present) cycle
	if (chemo(ichemo)%constant) then
		chemo(ichemo)%medium_Cbnd = chemo(ichemo)%bdry_conc
		cycle
	endif
	if (ichemo == OXYGEN .or. ichemo == GLUCOSE) then
		Nh = chemo(ichemo)%Hill_N
		asum = 0
		Nc = 0
		do kcell = 1,nlist
			if (cell_list(kcell)%state == DEAD) cycle
			Nc = Nc + 1
			C = cell_list(kcell)%conc(ichemo)
			metab = C**Nh/(chemo(ichemo)%MM_C0**Nh + C**Nh)
			dMdt = metab*chemo(ichemo)%max_cell_rate 
			asum = asum + dMdt
		enddo
		U(ichemo) = asum
	else
		! need to sum cell uptake and decay, and extracellular decay
		decay_rate = chemo(ichemo)%decay_rate
		Kin = chemo(ichemo)%membrane_diff_in
		Kout = chemo(ichemo)%membrane_diff_out
		asum = 0	
		ntvars = ODEdiff%nextra + ODEdiff%nintra
		do i = 1,ntvars
			if (ODEdiff%vartype(i) == EXTRA) then
				if (decay_rate == 0) cycle	! currently there is no distinction between intra- and extracellular decay rate
				vol_cm3 = Vsite_cm3
				Cex = allstate(i,ichemo)
				if (i < ntvars) then
					if (ODEdiff%vartype(i+1) == INTRA) then
						kcell = ODEdiff%cell_index(i+1)		! for access to cell-specific parameters (note that the intra variable follows the extra variable)
						vol_cm3 = Vsite_cm3 - Vcell_cm3*cell_list(kcell)%volume	! accounting for cell volume change
					endif
				endif
				asum = asum + vol_cm3*Cex*decay_rate
			else
				Cin = allstate(i,ichemo)
				kcell = ODEdiff%cell_index(i)
				vol_cm3 = Vcell_cm3*cell_list(kcell)%volume
				asum = asum + vol_cm3*Cin*decay_rate
				! add cell uptake rate
				Cex = allstate(i-1,ichemo)
				asum = asum + Kin*Cex - Kout*Cin
			endif
		enddo
		U(ichemo) = asum
!		if (ichemo == DRUG_A) write(*,*) 'UpdateMedium: ',ichemo,asum
	endif
enddo
chemo(:)%medium_U = U(:)

! First need the spheroid radius
call SetRadius(Nsites)
R1 = Radius*DELTA_X		! cm
Rlayer(:) = R1 + chemo(:)%medium_dlayer
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%present) cycle
	if (chemo(ichemo)%constant) cycle
	R2 = Rlayer(ichemo)
	if (ichemo /= OXYGEN) then	! update %medium_M, then %medium_Cext
!		if (ichemo == DRUG_A) write(*,*) 'UpdateMedium: medium_M: ',ichemo,chemo(ichemo)%medium_M
		chemo(ichemo)%medium_M = chemo(ichemo)%medium_M*(1 - chemo(ichemo)%decay_rate*dt) - U(ichemo)*dt
!		if (ichemo == DRUG_A) write(*,*) 'UpdateMedium: medium_M: ',ichemo,chemo(ichemo)%medium_M
		chemo(ichemo)%medium_Cext = (chemo(ichemo)%medium_M - (U(ichemo)/(6*chemo(ichemo)%medium_diff_coef)*R2) &
			*(R1*R1*(3*R2 - 2*R1) - R2*R2*R2))/(total_volume - 4*PI*R1*R1*R1/3.)
!		if (ichemo == DRUG_A) write(*,*) 'UpdateMedium: medium_Cext: ',ichemo,chemo(ichemo)%medium_Cext,(total_volume - 4*PI*R1*R1*R1/3.)
	endif
	a = (1/R2 - 1/R1)/(4*PI*chemo(ichemo)%medium_diff_coef)
	chemo(ichemo)%medium_Cbnd = chemo(ichemo)%medium_Cext + a*chemo(ichemo)%medium_U
enddo

!write(*,'(a,10e12.3)') 'UpdateMedium: ',chemo(:)%medium_Cbnd
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
do i = 1,ODEdiff%nvars
	state(i) = chemo(ichemo)%bdry_conc
enddo
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
integer, parameter :: SOLVER = RKC_SOLVE

! Variables for RKSUITE
!type(rk_comm_real_1d) :: comm
real(REAL_KIND) :: t_start, t_end, tolerance, t_want, t_inc, t_got
real(REAL_KIND), allocatable :: thresholds(:)

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
allocate(work(8+5*nvars,MAX_CHEMO))
xmid = NX/2
ymid = NY/2
zmid = NZ/2
k0 = ODEdiff%ivar(xmid,ymid,zmid)
tstart = 0
do ichemo = 1,nchemo
	if (.not.chemo(ichemo)%used) cycle
	call InitState(ichemo,allstate(1:nvars,ichemo))
	allstatep(:,ichemo) = 0
	if (SOLVER == RKF45_SOLVE) then
!		call deriv(tstart,allstate(1:nvars,ichemo),allstatep(1:nvars,ichemo),ichemo)
	elseif (SOLVER == RKSUITE_SOLVE) then
		t_end = 10000
		tolerance = 1.0e-4
		allocate(thresholds(nvars))
		thresholds = 0.1
!		call rk_setup(comm, tstart, allstate(1:nvars,ichemo), t_end,  tolerance, thresholds)
	elseif (SOLVER == RKC_SOLVE) then
		info(1) = 1
		info(2) = 1
		info(3) = 1
		info(4) = 0
		rtol = 1d-2
		atol = rtol
		idid = 0
	endif
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
		if (SOLVER == RKF45_SOLVE) then
			if (REAL_KIND == SP) then
!				call r4_rkf45 ( deriv, nvars, state, statep, tstart, tend, relerr, abserr, flag, ichemo )
			else
!				call r8_rkf45 ( deriv, nvars, state, statep, tstart, tend, relerr, abserr, flag, ichemo )
			endif
			if (flag /= 2) then
				write(logmsg,*) 'Bad flag: ',flag
				call logger(logmsg)
				stop
!				call r8_rkf45 ( deriv, nvars, state, statep, tstart, tend, relerr, abserr, flag, ichemo )
				flag = 2
			endif
			flag_par(ichemo) = flag
			abserr_par(ichemo) = abserr
			relerr_par(ichemo) = relerr
		elseif (SOLVER == RKSUITE_SOLVE) then
!			call range_integrate(comm, f_deriv, tend, t_got, state, statep, flag)
			if (flag /= 1 .and. flag /= 3 .and. flag /= 4) then
				write(logmsg,*) 'Bad flag: ',flag
				call logger(logmsg)
				stop
			endif
		elseif (SOLVER == RKC_SOLVE) then
			idid = 0
			call rkc(comm_rkc(ichemo),nvars,f_rkc,state,tstart,tend,rtol,atol,info,work(:,ichemo),idid,ichemo)
			if (idid /= 1) then
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
		deallocate(state)
		deallocate(statep)
	enddo
!	call TestAddSite
!   nvars = ODEdiff%nextra
!		prev_state = state
enddo

t2 = wtime()

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
				write(nflog,*) 'Error: TestAddSite: ',ODEdiff%ivar(x,y,z)
				stop
			endif
			site = (/x,y,z/)
			call ExtendODEDiff(site)
			return
		endif
	enddo
enddo
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
! NOT USED
!----------------------------------------------------------------------------------
subroutine ComputeGradient(conc,x,y,z,grad)
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
subroutine CheckCoxygen
real(REAL_KIND) :: MM0, dC, Cs, Ks, rmax, r, Kd, M, V, Ci, Ce
real(REAL_KIND) :: a, b, c, d, C1, C2
logical :: from_Ce = .false.

Cs = ODEdiff%C1_soft
dC = ODEdiff%deltaC_soft
MM0 = chemo(OXYGEN)%MM_C0
rmax = chemo(OXYGEN)%max_cell_rate
V = Vsite_cm3 - Vextra_cm3
r = rmax/V
Ks = ODEdiff%k_soft
Kd = chemo(OXYGEN)%membrane_diff_in

if (from_Ce) then
	! Specifying Ce -> Ci
	! If Ci > Cs
	Ce = 0.01

	a = 1
	b = r/Kd + MM0 - dC - Ce
	c = Ce*(dC - MM0) - (r/Kd)*dC

	d = sqrt(b*b - 4*a*c)
	C1 = (-b + d)/(2*a)
	C2 = (-b - d)/(2*a)
else
	! Specifying Ci -> Ce

	Ci = 0.00005
	do
		if (Ci > Cs) then
			M = (Ci - dC)/(MM0 + Ci - dC)
		else
			M = Ks*Ci*Ci
		endif
		Ce = Ci + M*r/Kd
		Ci = 2*Ci
		if (Ci > 0.15) exit
	enddo
endif
stop
end subroutine

!----------------------------------------------------------------------------------
! The idea is to split the solver into two parts.
!
! Extra:
! -----
! This solves the diffusion equation over all sites in the blob, 
! for each constituent separately.  For each constituent, nvar = nextra.
! Function: f_rkc_extra
!
! Intra:
! -----
! This solves the membrane diffusion, decay and intracellular reactions for all 
! constituents simultaneously, at each occupied site in turn.  The variables are
! both extra- and intracellular.  For each cell, nvar = 2*nchemo.
! Function: f_rkc_intra
!
! Variable indexing: 
! For i = 1,...,ntvars (= nextra + nintra):
!     ODEdiff%vartype(i) = EXTRA for an extracellular concentration
!                        = INTRA for an intracellular concentration
! NOT USED
!----------------------------------------------------------------------------------
subroutine NogoodSolver(it,tstart,dt,nc)
integer :: it, nc
real(REAL_KIND) :: tstart, dt
integer :: ichemo, ncvars, ntvars, ic, kcell, site(3), iv, nth
integer :: ie, ki, i, kextra, kintra
real(REAL_KIND) :: t, tend
real(REAL_KIND), allocatable :: state_ex(:,:)
real(REAL_KIND), allocatable :: state_in(:,:)
integer, allocatable :: exindex(:)
real(REAL_KIND), target :: Creact(2*MAX_CHEMO)
real(REAL_KIND) :: C(MAX_CHEMO), Ce(MAX_CHEMO), dCreact(MAX_CHEMO)
real(REAL_KIND) :: dCsum, dC
real(REAL_KIND) :: timer1, timer2
logical :: ok
logical :: dbug = .false.
! Variables for RKC
integer :: info(4), idid
real(REAL_KIND) :: rtol, atol(1), sprad_ratio
type(rkc_comm), allocatable :: comm_rkc(:)
logical :: solve_intra = .true.

ODEdiff%nintra = nc
ntvars = ODEdiff%nextra + ODEdiff%nintra
allocate(state_ex(ODEdiff%nextra,MAX_CHEMO))
allocate(state_in(ODEdiff%nintra,MAX_CHEMO))
allocate (exindex(ODEdiff%nintra))
allocate(comm_rkc(MAX_VARS))
kextra = 0
kintra = 0
do i = 1,ntvars
	if (ODEdiff%vartype(i) == EXTRA) then
		kextra = kextra + 1
		state_ex(kextra,1:nchemo) = allstate(i,1:nchemo)
	else
		kintra = kintra + 1
		exindex(kintra) = kextra	! the index into state_ex(:) to the corresponding extracellular variables
	endif
enddo

!if (it == 1) then
!	call showresults(state(:,DRUG_A))
!	call showresults(state(:,DRUG_METAB_A))
!endif

info(1) = 1
info(2) = 1		! = 1 => use spcrad() to estimate spectral radius, != 1 => let rkc do it
info(3) = 1
info(4) = 0
rtol = 1d-2
atol = rtol

if (dbug) write(nflog,'(a,2e12.3)') 'start: ',state_ex(1,1:nchemo)
! Solve diffusion equation for each constituent
do ic = 1,nchemo
	ichemo = chemomap(ic)
	idid = 0
	t = tstart
	tend = t + dt
	call rkc(comm_rkc(ichemo),ODEdiff%nextra,f_rkc_extra,state_ex(:,ichemo),t,tend,rtol,atol,info,work_rkc,idid,ichemo)
	if (idid /= 1) then
		write(logmsg,*) ' Failed at t = ',t,' with idid = ',idid
		call logger(logmsg)
		stop
	endif
enddo

if (solve_intra) then
ncvars = 2*nchemo
kintra = 0
! Solve reactions, decay, and membrane diffusion
! Note: the order could (should?) be randomised using a permutation
do i = 1,ntvars
	if (ODEdiff%vartype(i) == EXTRA) cycle
	kintra = kintra + 1
	kextra = exindex(kintra)
	! initialize Creact
	do ic = 1,nchemo
		Creact(ic) = state_ex(kextra,ic)	! extracellular concentrations first
		Creact(ic+nchemo) = allstate(i,ic)	! intracellular concentrations next
	enddo
	if (dbug .and. kintra == 1) then
		write(nflog,'(a,i4,4e12.3)') 'before: ',nchemo,Creact(1:4)
	endif
	idid = 0
	t = tstart
	tend = t + dt
	! what do do with comm_rkc?
!	call rkc(comm_rkc(kintra),ncvars,f_rkc_intra,Creact,t,tend,rtol,atol,info,work_rkc,idid,kintra)
	if (idid /= 1) then
		write(logmsg,*) ' Failed at t = ',t,' with idid = ',idid
		call logger(logmsg)
		stop
	endif
	if (dbug .and. kintra == 1) then
		write(nflog,'(a,i4,4e12.3)') 'after: ',nchemo,Creact(1:4)
	endif
	! transfer results
	do ic = 1,nchemo
		state_ex(kextra,ic) = Creact(ic)		! extracellular concentrations first
		state_in(kintra,ic) = Creact(ic+nchemo) ! intracellular concentrations next
	enddo
enddo
endif

! transfer results to allstate(:,:)
kextra = 0
kintra = 0
do i = 1,ntvars
	if (ODEdiff%vartype(i) == EXTRA) then
		kextra = kextra + 1
		allstate(i,1:nchemo) = state_ex(kextra,1:nchemo) 
	elseif (solve_intra) then
		kintra = kintra + 1
		allstate(i,1:nchemo) = state_in(kintra,1:nchemo) 
	endif
enddo
if (dbug) write(nflog,'(a,4e12.3)') 'allstate (a): ',allstate(1,1:2),allstate(2,1:2)
!allstate(1:nvars,1:MAX_CHEMO) = state(:,:)
! Note: some time we need to copy the state values to the cell_list() array.
deallocate(state_ex)
deallocate(state_in)
deallocate(exindex)
deallocate(comm_rkc)
end subroutine

!----------------------------------------------------------------------------------
! Simple extracellular diffusion of a single constituent
! The site volumes should enter into the calculation of dydt
! NOT USED
!----------------------------------------------------------------------------------
subroutine f_rkc_extra(neqn,t,y,dydt,icase)
integer :: neqn, icase
real(REAL_KIND) :: t, y(neqn), dydt(neqn)
integer :: i, k, ie, ki, kv, ichemo
real(REAL_KIND) :: dCsum, dCdiff,  DX2, DX3, vol, val
real(REAL_KIND) :: dc1, dc6, cbnd

ichemo = icase
DX2 = DELTA_X*DELTA_X
dc1 = chemo(ichemo)%diff_coef/DX2
dc6 = 6*dc1 !+ decay_rate
cbnd = BdryConc(ichemo,t_simulation)
!$omp parallel do private(vol, dCsum, k, kv, dCdiff, val) default(shared) schedule(static)
do ie = 1,neqn
	vol = Vextra_cm3			!!!!!! need to worry about volume !!!!
    dCsum = 0
    do k = 1,7
		kv = ODEdiff%iexcoef(ie,k)
	    if (k == 1) then
		    dCdiff = -dc6
	    else
		    dCdiff = dc1
	    endif
	    if (kv == 0) then
		    val = cbnd
	    else
		    val = y(kv)
	    endif
	    dCsum = dCsum + dCdiff*val
    enddo
	dydt(ie) = dCsum		!+ dCreact
enddo
!$omp end parallel do
end subroutine

!----------------------------------------------------------------------------------
! A cell has moved to site(:), which was previously exterior (boundary), and the
! variable-site mappings need to be updated, together with %icoef(:,:)
! The relevant neighbours are at x +/- 1, y +/- 1, z +/- 1
! %nextra is incremented, and allstate(nextra,:) is initialized.
! NOT USED
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

do ichemo = 1,MAX_CHEMO
    csum(ichemo) = csum(ichemo) + (6-nb)*BdryConc(ichemo,t_simulation)
enddo
allstate(n,:) = csum/6
end subroutine

!----------------------------------------------------------------------------------
! Initialise the state vector to the current concentrations
! NOT USED
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

end module

