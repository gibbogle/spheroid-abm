!----------------------------------------------------------------------------------------
! Identify S1P, CCL21 and OXY influx sites.
! Identify sites where cells can enter and exit.
! Candidate B cell entry sites are at the lower surface of the ellipsoidal blob,
! in the patch delineated by the plane y = -ENTRY_ALPHA*Rb, where:
!    Rb = Ra/ELLIPSE_RATIO
! and Ra and Rb are the major and minor axis radii of the defining ellipse.
!----------------------------------------------------------------------------------------
subroutine AssignBdryRole(site,bdry)
integer :: site(3), dy
type (boundary_type), pointer :: bdry
real, parameter :: S1Pfraction = 0.25

bdry%chemo_influx = .false.
site = bdry%site
dy = site(2) - Centre(2)
if (dy > 0) then
	bdry%chemo_influx(S1P) = .true.
	bdry%chemo_influx(OXY) = .true.
else
    if (dy > -S1Pfraction*Radius%y) then
		bdry%chemo_influx(S1P) = .true.
    else
		bdry%chemo_influx(CCL21) = .true.
    endif
endif
if (dy > -ENTRY_ALPHA*Radius%y) then
    bdry%entry_ok = .false.
else
    bdry%entry_ok = .true.
endif
if (dy > -EXIT_ALPHA*Radius%y) then
    bdry%exit_ok = .false.
else
    bdry%exit_ok = .true.
endif 
end subroutine

!-----------------------------------------------------------------------------------------
! The added site needs approximate values for chemokine concentrations and gradients,
! as an interim measure until the next steady-state computation.  Actually only the
! gradient is used (for chemotaxis) but the concentration serves as the starting value
! when the steady-state solution is computed.
!-----------------------------------------------------------------------------------------
subroutine SetBdryConcs(site)
integer :: site(3)
integer :: ic, i, nsum, nbsite(3)
type (boundary_type), pointer :: bdry
real :: csum
logical :: set(MAX_CHEMO)

bdry => occupancy(site(1),site(2),site(3))%bdry
if (.not.associated(bdry)) then
	write(logmsg,*) 'Error: SetBdryConcs: not a bdry site'
	call logger(logmsg)
	stop
endif

set = .false.
do ic = 1,MAX_CHEMO
	if (chemo(ic)%used .and. .not.chemo(ic)%use_secretion) then
		if (bdry%chemo_influx(ic)) then
			chemo(ic)%conc(site(1),site(2),site(3)) = chemo(ic)%bdry_conc
			set(ic) = .true.
		endif
	endif
enddo
do ic = 1,MAX_CHEMO
	if (chemo(ic)%used .and. .not.set(ic)) then
		! set the conc to the average of neighbour sites
		csum = 0
		nsum = 0
		do i = 1,27
			if (i==14) cycle
			nbsite = site + jumpvec(:,i)
			if (ODEdiff%ivar(nbsite(1),nbsite(2),nbsite(3)) == 0) cycle
			nsum = nsum + 1
			csum = csum + chemo(ic)%conc(nbsite(1),nbsite(2),nbsite(3))
		enddo
		if (nsum > 0) then
			chemo(ic)%conc(site(1),site(2),site(3)) = csum/nsum
		endif
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! This subroutine is for the case of a site that was boundary and is now in the blob
! interior.  The chemokine gradients are adjusted to reduce the inaccuracy in chemotactic 
! effects until the new steady-state is computed, while the aim in adjusting the concs
! is to speed up convergence in the steady-state computation.  In fact it doesn't have
! much effect on the convergence.  There might be a better way to do the adjustment than
! the simple neighbourhood averaging that is used.
!-----------------------------------------------------------------------------------------
subroutine SetConcs(site)
integer :: site(3)
integer :: i

do i = 1,MAX_CHEMO
	if (chemo(i)%used) then
		call AverageConc(chemo(i)%conc,chemo(i)%grad,site)
	endif
enddo
end subroutine



!-----------------------------------------------------------------------------------------
! Using bdrylist, randomly select bdry sites, check that the site is an entrysite and
! suitable for cell ingress.  In fact we look at all neighbouring 'inside' sites.
!-----------------------------------------------------------------------------------------
subroutine GetEntrySite(site,ok)
integer :: site(3)
logical :: ok
integer :: ibdry, k, it, bsite(3), indx(2), nt = 1000, kpar=0
type (boundary_type), pointer :: bdry

it = 0
do
    it = it+1
    if (it == nt) then
        ok = .false.
        return
    endif
    ibdry = random_int(1,nbdry,kpar)
    k = 0
    bdry => bdrylist
    do while ( associated ( bdry )) 
        k = k+1
        if (k == ibdry) exit
        bdry => bdry%next
    enddo
    if (.not.bdry%entry_ok) cycle
    do k = 1,27
	    site = bdry%site + jumpvec(:,k)
	    indx = occupancy(site(1),site(2),site(3))%indx
        if (indx(1) < 0) cycle                      ! outside or DC
        if (indx(1) == 0 .or. indx(2) == 0) then    ! free slot
            ok = .true.
            return
        endif
    enddo     
enddo
end subroutine

!----------------------------------------------------------------------------------------
! If the site is a source boundary for the chemokine, cbnd = boundary concentration,
! otherwise cbnd = -1, which is a flag to compute the estimate of concentration by
! averaging the over neighbour sites.
! For the gradient, averaging is always used.
!----------------------------------------------------------------------------------------
subroutine AverageBdryConc(bdry,C,G,site,cbnd)
integer :: site(3)
real :: C(:,:,:), G(:,:,:,:)
real :: cbnd
integer :: x, y, z, k, nave
real :: cave, gave(3)
type (boundary_type), pointer :: bdry

cave = 0
gave = 0
nave = 0
do k = 1,27
	if (k == 14) cycle
    x = site(1) + jumpvec(1,k)
    y = site(2) + jumpvec(2,k)
    z = site(3) + jumpvec(3,k)
    if (outside_xyz(x,y,z)) cycle
!    if (occupancy(x,y,z)%indx(1) < 0) cycle	! outside or DC
    if (.not.ChemoRegion(occupancy(x,y,z)%indx)) cycle	! site not accessible by chemokines
    nave = nave + 1
    cave = cave + C(x,y,z)
    gave = gave + G(:,x,y,z)
enddo
if (cbnd >= 0) then
	C(site(1),site(2),site(3)) = cbnd
else
	if (nave > 0) then
		C(site(1),site(2),site(3)) = cave/nave
	else
		C(site(1),site(2),site(3)) = 0
	endif
endif
if (nave > 0) then
	G(:,site(1),site(2),site(3)) = gave/nave
else
	G(:,site(1),site(2),site(3)) = 0
endif
end subroutine

!----------------------------------------------------------------------------------------
! The concentration and gradient are estimated by averaging the over neighbour sites.
!----------------------------------------------------------------------------------------
subroutine AverageConc(C,G,site)
integer :: site(3)
real :: C(:,:,:), G(:,:,:,:)
integer :: x, y, z, k, nave
real :: cave, gave(3)

cave = 0
gave = 0
nave = 0
do k = 1,27
	if (k == 14) cycle
    x = site(1) + jumpvec(1,k)
    y = site(2) + jumpvec(2,k)
    z = site(3) + jumpvec(3,k)
    if (outside_xyz(x,y,z)) cycle
!    if (occupancy(x,y,z)%indx(1) < 0) cycle	! outside or DC
    if (.not.ChemoRegion(occupancy(x,y,z)%indx)) cycle	! site not accessible by chemokines
    nave = nave + 1
    cave = cave + C(x,y,z)
    gave = gave + G(:,x,y,z)
enddo
if (nave > 0) then
	C(site(1),site(2),site(3)) = cave/nave
	G(:,site(1),site(2),site(3)) = gave/nave
else
	C(site(1),site(2),site(3)) = 0
	G(:,site(1),site(2),site(3)) = 0
endif
end subroutine

!----------------------------------------------------------------------------------------
! We can solve for a steady-state concentration field, compute the gradient field from it,
! and use this for all chemotaxis calculations until the ellipsoid size changes.  At this
! point we need to resolve for the concentration field.
! There are two possible approaches.  Either to solve for the steady-state field again,
! using the previous solution as the starting condition, or to allow the field to evolve.
! In the latter case the concentration in some number of sites in the neighbourhood of the
! added/removed site is averaged in a way that ensures mass conservation.
! This all works, but there is some doubt as to whether there is such a thing as
! S1P chemotaxis.  According to Irina Grigorova, S1P1 on a T cell enables it to cross
! the endothelial boundary into a sinus (high-S1P), but the cell must reach the sinus
! by random motion.
!----------------------------------------------------------------------------------------
subroutine AllocateConcArrays
integer :: ic

write(logmsg,*) 'AllocateConcArrays'
call logger(logmsg)
do ic = 1,MAX_CHEMO
	if (chemo(ic)%used) then
		allocate(chemo(ic)%conc(NX,NY,NZ))
		allocate(chemo(ic)%grad(3,NX,NY,NZ))
		chemo(ic)%conc = 0	
	endif
enddo
! Note: the arrays ivar and varsite are generally useful, not just for ODE diffusion
allocate(ODEdiff%ivar(NX,NY,NZ))
allocate(ODEdiff%varsite(NX*NY*NZ,3))
if (use_ODE_diffusion) then
	allocate(ODEdiff%icoef(NX*NY*NZ,7))
endif

end subroutine

!----------------------------------------------------------------------------------------
! Set up boundary concentrations for the chemokines.
! These values do not change as long as the boundaries or FDCs/MRCs do not move.
! Note that the concentrations produced by FDCs and MRCs are currently treated as the same.
!----------------------------------------------------------------------------------------
subroutine BdryConcentrations
integer :: ichemo, i, site(3)
type (boundary_type), pointer :: bdry

write(logmsg,*) 'BdryConcentrations'
call logger(logmsg)
call CheckBdryList
do ichemo = 1,MAX_CHEMO
	if (chemo(ichemo)%used .and. .not.chemo(ichemo)%use_secretion) then
		if (ichemo /= CXCL13) then
			bdry => bdrylist
			do while ( associated ( bdry )) 
				if (bdry%chemo_influx(ichemo)) then
					site = bdry%site
					chemo(ichemo)%conc(site(1),site(2),site(3)) = chemo(ichemo)%bdry_conc
				endif
				bdry => bdry%next
			enddo
		else	! special treatment needed for CXCL13, which is secreted by FDCs and MRCs
			if (USE_CELL_SITES) then
				do i = 1,NFDC
					site = FDC_list(i)%site
					chemo(ichemo)%conc(site(1),site(2),site(3)) = chemo(ichemo)%bdry_conc
				enddo
				do i = 1,NMRC
					site = MRC_list(i)%site
					chemo(ichemo)%conc(site(1),site(2),site(3)) = chemo(ichemo)%bdry_conc
				enddo
			else
				! need to maintain a list of FDC neighbour sites - for now use occupancy()%FDC_nbdry, OK if FDCs do not move
				do i = 1,ODEdiff%nvars
					site = ODEdiff%varsite(i,:)
					if (occupancy(site(1),site(2),site(3))%FDC_nbdry > 0 &
					.or. occupancy(site(1),site(2),site(3))%MRC_nbdry > 0) then
						chemo(ichemo)%conc(site(1),site(2),site(3)) = chemo(ichemo)%bdry_conc
					endif
				enddo
			endif
		endif
	endif
enddo
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine ChemoSteadystate
integer :: ichemo, x, y, z
real :: g(3), gamp, gmax(MAX_CHEMO)
logical, save :: first = .true.

write(logmsg,*) 'ChemoSteadystate'
call logger(logmsg)
call SetupODEDiffusion
call BdryConcentrations
if (use_ODE_diffusion) then
	call SolveSteadystate_B
else
	call SolveSteadystate_A
endif
gmax = 0
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	do x = 1,NX
		do y = 1,NY
			do z = 1,NZ
				g = chemo(ichemo)%grad(:,x,y,z)
				gamp = sqrt(dot_product(g,g))
				gmax(ichemo) = max(gamp,gmax(ichemo))
			enddo
		enddo
	enddo
enddo
write(logmsg,'(a,4(i3,f8.3))') 'Max gradients: ',(ichemo,gmax(ichemo),ichemo=1,MAX_CHEMO)
call logger(logmsg)
if (first) then
	call ShowConcs
endif
first = .false.

end subroutine

!----------------------------------------------------------------------------------------
! Recompute steady-state concentration fields if there has been a significant change in
! the B cell population.
!----------------------------------------------------------------------------------------
subroutine UpdateSSFields
integer, save :: NBlast = 0
integer :: x, y, z, site(3)
real :: delNB
real :: cmin, cmax, dc, dcmax
real, allocatable :: S1P_old(:,:,:)
logical :: doit

if (NBlast == 0) then
	NBlast = NBcells0
endif
doit = .false.
if (inflammation_level == 0) then
	if (mod(istep,4*60*6) == 0) then	! every 6 hours in no inflammation case
		doit = .true.
	endif
else
	delNB = abs(NBcells - NBlast)
	if (delNB/NBcells > 0.05) then
		doit = .true.
	endif
endif
if (.not.doit) return

call ChemoSteadystate
NBlast = NBcells
end subroutine

!----------------------------------------------------------------------------------------
! Advance the dynamic solution for chemokine concentration fields through time interval dtstep.
! For now only treat the case of NOT use_ODE_diffusion, because use_ODE_diffusion is a bit
! more complicated (should store derivatives).
!----------------------------------------------------------------------------------------
subroutine UpdateFields(dtstep)
real :: dtstep
type(chemokine_type), pointer :: Cptr
integer :: ichemo, it, nt = 4
real :: Kdiffusion, Kdecay, dt
integer :: z1, z2, n, kpar
real, allocatable :: C_par(:,:,:)

!write(*,*) 'UpdateFields: zoffset: ',zoffset
dt = dtstep/nt
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
!	write(*,*) 'UpdateFields: ichemo: ',ichemo
	Cptr => chemo(ichemo)
	Kdiffusion = Cptr%diff_coef
	Kdecay = Cptr%decay_rate
	do it = 1,nt
		!$omp parallel do private(z1,z2,n,C_par,dt)
		do kpar = 0,Mnodes-1
			if (Mnodes == 1) then
				z1 = blobrange(3,1)
				z2 = blobrange(3,2)
			else
    	        z1 = max(zoffset(2*kpar) + 1,blobrange(3,1))
	            z2 = min(zoffset(2*kpar+2),blobrange(3,2))
			endif
			n = z2 - z1 + 1
			dt = dtstep/nt
			allocate(C_par(NX,NY,n))
			call par_evolve_A(ichemo,Cptr%conc,Kdiffusion,Kdecay,C_par,z1,z2,dt,kpar)
			Cptr%conc(:,:,z1:z2) = C_par(:,:,1:n)
			deallocate(C_par)
		enddo
	enddo
!	write(*,*) 'call gradient'
	call gradient(Cptr%conc,Cptr%grad)
enddo
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine ShowConcs
integer :: x, y, z, i
real :: gamp

x = NX/2
z = 45
do i = 1,MAX_CHEMO
	if (chemo(i)%used) then
		write(nfout,'(a,a)') chemo(i)%name,'  conc        gradient' 
		do y = 1,NY
!			if (occupancy(x,y,z)%indx(1) < 0) cycle	! outside or DC
			if (.not.ChemoRegion(occupancy(x,y,z)%indx)) cycle	! site not accessible by chemokines
			gamp = norm(chemo(i)%grad(:,x,y,z))
		    write(nfout,'(3i4,4e12.4,4x,f8.4)') x,y,z,chemo(i)%conc(x,y,z),chemo(i)%grad(:,x,y,z),gamp
		enddo
	endif
enddo
end subroutine


!----------------------------------------------------------------------------------------
! Solve for steady-state chemokine concentrations, given levels at source sites, using
! Method A.
! On entry C contains the starting concentration field, which may be 0.
!----------------------------------------------------------------------------------------
subroutine SolveSteadystate_A
type(chemokine_type), pointer :: Cptr
integer :: ichemo
real :: Kdiffusion, Kdecay
real :: dx2diff, total, maxchange, maxchange_par, total_par
real, parameter :: alpha = 0.5
real, parameter :: tol = 1.0e-5		! max change in C at any site as fraction of average C
integer :: nc, nc_par, k, it, n, kpar
integer :: nsweeps, sweep, slice
integer :: z1, z2, zfr, zto	
real, allocatable :: C_par(:,:,:)
integer :: nt = 10000
character*(12) :: msg

if (Mnodes == 1) then
    nsweeps = 1
else
    nsweeps = 2
endif

do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	Cptr => chemo(ichemo)
	Kdiffusion = Cptr%diff_coef
	Kdecay = Cptr%decay_rate
	do it = 1,nt
		maxchange = 0
		total = 0
		nc = 0
	do sweep = 0,nsweeps-1
		!$omp parallel do private(slice,z1,z2,n,C_par,maxchange_par,total_par,nc_par)
		do kpar = 0,Mnodes-1
			slice = sweep + 2*kpar
	        if (Mnodes == 1) then
	            z1 = blobrange(3,1)
	            z2 = blobrange(3,2)
	        else
    	        z1 = max(zoffset(slice) + 1,blobrange(3,1))
	            z2 = min(zoffset(slice+1),blobrange(3,2))
	        endif
			n = z2 - z1 + 1
			allocate(C_par(NX,NY,n))
			call par_steadystate_A(ichemo,Cptr%conc,Kdiffusion,Kdecay,C_par,z1,z2,maxchange_par,total_par,nc_par,kpar,sweep)
			Cptr%conc(:,:,z1:z2) = C_par(:,:,1:n)
			deallocate(C_par)
			nc = nc + nc_par
			total = total + total_par
			maxchange = max(maxchange,maxchange_par)
		enddo
	enddo
		if (maxchange < tol*total/nc) then
			write(logmsg,'(a,a,a,i4)') 'Convergence reached for: ',chemo(ichemo)%name,' # of iterations: ',it
			call logger(logmsg)
			exit
		endif
	enddo
	call gradient(Cptr%conc,Cptr%grad)
enddo
end subroutine	

!----------------------------------------------------------------------------------------
subroutine SolveSteadystate_A_original
type(chemokine_type), pointer :: Cptr
integer :: ichemo
real :: Kdiffusion, Kdecay
real :: dx2diff, total, maxchange, maxchange_par, total_par
real, parameter :: alpha = 0.5
real, parameter :: tol = 1.0e-5		! max change in C at any site as fraction of average C
integer :: nc, nc_par, k, it, n, kpar
integer :: z1, z2, zfr, zto	
real, allocatable :: C_par(:,:,:)
integer :: nt = 10000
integer :: sweep = 0

do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	Cptr => chemo(ichemo)
	Kdiffusion = Cptr%diff_coef
	Kdecay = Cptr%decay_rate
	do it = 1,nt
		maxchange = 0
		total = 0
		nc = 0
		!$omp parallel do private(z1,z2,n,C_par,maxchange_par,total_par,nc_par)
		do kpar = 0,Mnodes-1
	        if (Mnodes == 1) then
	            z1 = blobrange(3,1)
	            z2 = blobrange(3,2)
	        else
    	        z1 = max(zoffset(2*kpar) + 1,blobrange(3,1))
	            z2 = min(zoffset(2*kpar+2),blobrange(3,2))
	        endif
			n = z2 - z1 + 1
			allocate(C_par(NX,NY,n))
			call par_steadystate_A(ichemo,Cptr%conc,Kdiffusion,Kdecay,C_par,z1,z2,maxchange_par,total_par,nc_par,kpar,sweep)
			Cptr%conc(:,:,z1:z2) = C_par(:,:,1:n)
			deallocate(C_par)
			nc = nc + nc_par
			total = total + total_par
			maxchange = max(maxchange,maxchange_par)
		enddo
		if (maxchange < tol*total/nc) then
			write(logmsg,'(a,a,a,i4)') 'Convergence reached for: ',chemo(ichemo)%name,' # of iterations: ',it
			call logger(logmsg)
			exit
		endif
	enddo
	call gradient(Cptr%conc,Cptr%grad)
enddo
end subroutine	

!----------------------------------------------------------------------------------------
! A different approach from that used in bone-abm.
! The bdry gridcells that are adjacent to a chemokine-rich region (i.e. there are 
! chemokine-secreting gridcells near the boundary) are given a fixed concentration.
! ***All boundaries are no-flux***
! This is a simplified approach.  We need only keep track of the boundary.
! The concentrations in a specified chemokine influx bdry site are input parameters.
! Note:  It is only the ratio Kdecay/Kdiffusion that matters.
! Method A solves in an iterative fashion dC/dt = 0, with specified concentrations
! in the source gridcells.
!----------------------------------------------------------------------------------------
subroutine par_steadystate_A(ichemo,C,Kdiffusion,Kdecay,Ctemp,z1,z2,maxchange,total,nc,kpar,sweep)
integer :: ichemo, sweep
real :: C(:,:,:), Ctemp(:,:,:)
integer :: z1, z2, kpar
real :: Kdiffusion, Kdecay
real :: dx2diff, total, maxchange, dC, sum, dV
real, parameter :: alpha = 0.5		!0.7
integer :: x, y, z, xx, yy, zz, nb, nc, k, zpar, indx(2), i, maxsite(3)
logical :: source_site

dx2diff = DELTA_X**2/Kdiffusion
dV = DELTA_X**3
maxchange = 0
total = 0
nc = 0
do zpar = 1,z2-z1+1
	z = zpar + z1-1
    do y = blobrange(2,1),blobrange(2,2)
        do x = blobrange(3,1),blobrange(3,2)
		    indx = occupancy(x,y,z)%indx
			if (.not.ChemoRegion(indx)) cycle	! site not accessible by chemokines
            source_site = .false.
            if (associated(occupancy(x,y,z)%bdry)) then
                ! Check for chemo bdry site - no change to the concentration at such a site
                do i = 1,MAX_CHEMO
	                if (ichemo == i .and. occupancy(x,y,z)%bdry%chemo_influx(i)) then
		                C(x,y,z) = chemo(i)%bdry_conc
			            Ctemp(x,y,zpar) = C(x,y,z)
				        source_site = .true.
					endif
				enddo
			elseif (ichemo == CXCL13 .and. &
			  (occupancy(x,y,z)%FDC_nbdry > 0 .or. occupancy(x,y,z)%MRC_nbdry > 0)) then
                C(x,y,z) = chemo(ichemo)%bdry_conc
	            Ctemp(x,y,zpar) = C(x,y,z)
		        source_site = .true.
            endif
            if (.not.source_site) then
			    sum = 0
			    nb = 0
			    do k = 1,6
				    xx = x + neumann(1,k)
				    yy = y + neumann(2,k)
				    zz = z + neumann(3,k)
				    if (outside_xyz(xx,yy,zz)) cycle
!    				if (occupancy(xx,yy,zz)%indx(1) < 0) cycle	! outside or DC
					indx = occupancy(xx,yy,zz)%indx
					if (.not.ChemoRegion(indx)) cycle	! site not accessible by chemokines
				    nb = nb + 1
				    sum = sum + C(xx,yy,zz)
			    enddo
			    Ctemp(x,y,zpar) = alpha*(DELTA_X*Kdiffusion*sum)/(Kdecay*dV + nb*DELTA_X*Kdiffusion) + (1-alpha)*C(x,y,z)
			    dC = abs(Ctemp(x,y,zpar) - C(x,y,z))
!			    maxchange = max(dC,maxchange)
				if (dC > maxchange) then
					maxchange = DC
					maxsite = (/x,y,z/)
				endif
			endif
			nc = nc + 1
			total = total + Ctemp(x,y,zpar)
		enddo
	enddo
enddo
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine debug_grad(msg)
character*(12) :: msg
integer :: x, y, z, indx(2)
real :: cbnd, cprev
logical :: bnd

write(*,'(a,$)') msg
y = Centre(2)
z = Centre(3)
cprev = -1
do x = 65,72
	indx = occupancy(x,y,z)%indx
	if (indx(1) == OUTSIDE_TAG) then
		write(*,*)
		if (cprev == 0) stop
		exit
	endif
	write(*,'(f6.1,$)') chemo(1)%conc(x,y,z)
	cprev = chemo(1)%conc(x,y,z)
	if (associated(occupancy(x,y,z)%bdry)) then
		bnd = .true.
		if (occupancy(x,y,z)%bdry%chemo_influx(1)) then
			cbnd = chemo(1)%bdry_conc
		else
			cbnd = 0
		endif
		write(*,'(a,f6.1,$)') '  BND:',cbnd
	else
		bnd = .false.
	endif
enddo
end subroutine

!----------------------------------------------------------------------------------------
! Updates concentrations through a timestep dt, solving the diffusion-decay eqtn by a 
! simple explicit method.  This assumes concentration boundary conditions.
!----------------------------------------------------------------------------------------
subroutine par_evolve_A(ichemo,C,Kdiffusion,Kdecay,Ctemp,z1,z2,dt,kpar)
integer :: ichemo, z1, z2, kpar
real :: C(:,:,:), Ctemp(:,:,:)
real :: Kdiffusion, Kdecay, dt
real :: C0, sum, dV, dMdt
integer :: x, y, z, zpar, xx, yy, zz, nb, k, indx(2), i
logical :: source_site

!write(*,*) 'par_evolve_A: ',ichemo,kpar,z1,z2,dt
dV = DELTA_X**3
do zpar = 1,z2-z1+1
	z = zpar + z1-1
    do y = blobrange(2,1),blobrange(2,2)
        do x = blobrange(1,1),blobrange(1,2)
		    indx = occupancy(x,y,z)%indx
!            if (indx(1) < 0) cycle      ! outside or DC
			if (.not.ChemoRegion(indx)) cycle	! site not accessible by chemokines
            source_site = .false.
            if (associated(occupancy(x,y,z)%bdry)) then
                ! Check for chemo bdry site - no change to the concentration at such a site
                do i = 1,MAX_CHEMO
	                if (ichemo == i .and. occupancy(x,y,z)%bdry%chemo_influx(i)) then
		                C(x,y,z) = chemo(i)%bdry_conc
			            Ctemp(x,y,zpar) = C(x,y,z)
				        source_site = .true.
					endif
				enddo
			elseif (ichemo == CXCL13 .and. &
			(occupancy(x,y,z)%FDC_nbdry > 0 .or. occupancy(x,y,z)%MRC_nbdry > 0)) then
                C(x,y,z) = chemo(ichemo)%bdry_conc
	            Ctemp(x,y,zpar) = C(x,y,z)
		        source_site = .true.
            endif
            if (.not.source_site) then
				C0 = C(x,y,z)
			    sum = 0
			    nb = 0
			    do k = 1,6
				    xx = x + neumann(1,k)
				    yy = y + neumann(2,k)
				    zz = z + neumann(3,k)
				    if (outside_xyz(xx,yy,zz)) cycle
!				    if (occupancy(xx,yy,zz)%indx(1) < 0) cycle	! outside or DC
					indx = occupancy(xx,yy,zz)%indx
					if (.not.ChemoRegion(indx)) cycle	! site not accessible by chemokines
				    nb = nb + 1
				    sum = sum + C(xx,yy,zz)
			    enddo
			    dMdt = Kdiffusion*DELTA_X*(sum - nb*C0) - Kdecay*C0*dV ! + influx(x,y,z)
			    Ctemp(x,y,zpar) = (C0*dV + dMdt*dt)/dV
			endif
		enddo
	enddo
enddo
end subroutine

!----------------------------------------------------------------------------------------
! Updates concentrations through a timestep dt, solving the diffusion-decay eqtn by a 
! simple explicit method.  This assumes concentration boundary conditions.
! Could be parallelized?
!----------------------------------------------------------------------------------------
subroutine evolve(ichemo,Kdiffusion,Kdecay,C,dt)
integer :: ichemo
real :: C(:,:,:)
real :: Kdiffusion, Kdecay, dt
real :: dx2diff, total, maxchange, C0, dC, sum, dV, dMdt
real, parameter :: alpha = 0.99
integer :: x, y, z, xx, yy, zz, nb, nc, k, it, i, indx(2)
real, allocatable :: Ctemp(:,:,:)
logical :: bdry_conc

dV = DELTA_X**3
allocate(Ctemp(NX,NY,NZ))
do z = blobrange(3,1),blobrange(3,2)
    do y = blobrange(2,1),blobrange(2,2)
        do x = blobrange(1,1),blobrange(1,2)
!			if (occupancy(x,y,z)%indx(1) < 0) cycle	! outside or DC
			indx = occupancy(x,y,z)%indx
			if (.not.ChemoRegion(indx)) cycle	! site not accessible by chemokines
			C0 = C(x,y,z)
            bdry_conc = .false.
            if (associated(occupancy(x,y,z)%bdry)) then
                ! Check for S1P bdry site - no change to the concentration at such a site
                do i = 1,MAX_CHEMO
	                if (ichemo == i .and. occupancy(x,y,z)%bdry%chemo_influx(i)) then
		                C(x,y,z) = chemo(i)%bdry_conc
			            Ctemp(x,y,z) = C(x,y,z)
				        bdry_conc = .true.
					endif
				enddo
            endif
            if (.not.bdry_conc) then
			    sum = 0
			    nb = 0
			    do k = 1,6
				    xx = x + neumann(1,k)
				    yy = y + neumann(2,k)
				    zz = z + neumann(3,k)
				    if (outside_xyz(xx,yy,zz)) cycle
!				    if (occupancy(xx,yy,zz)%indx(1) < 0) cycle	! outside or DC
					indx = occupancy(xx,yy,zz)%indx
					if (.not.ChemoRegion(indx)) cycle	! site not accessible by chemokines
				    nb = nb + 1
				    sum = sum + C(xx,yy,zz)
			    enddo
			    dMdt = Kdiffusion*DELTA_X*(sum - nb*C0) - Kdecay*C0*dV ! + influx(x,y,z)
			    Ctemp(x,y,z) = (C0*dV + dMdt*dt)/dV
			endif
		enddo
	enddo
enddo
C = Ctemp
deallocate(Ctemp) 
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine gradient(C,grad)
real :: C(:,:,:), grad(:,:,:,:)
integer :: x, y, z, xx, yy, zz, x1, x2, y1, y2, z1, z2, i, k, indx(2)
real :: g(3)
logical :: missed
real, parameter :: MISSING_VAL = 1.0e10

grad = 0
do z = blobrange(3,1),blobrange(3,2)
    do y = blobrange(2,1),blobrange(2,2)
        do x = blobrange(1,1),blobrange(1,2)
!			if (occupancy(x,y,z)%indx(1) < 0) cycle	! outside or DC
			indx = occupancy(x,y,z)%indx
			if (.not.ChemoRegion(indx)) cycle	! site not accessible by chemokines
			x1 = x - 1
			x2 = x + 1
			if (x1 < 1 .or. x2 > NX) then
				g(1) = 0
!			elseif (occupancy(x1,y,z)%indx(1) >= 0 .and. occupancy(x2,y,z)%indx(1) >= 0) then
			elseif (ChemoRegion(occupancy(x1,y,z)%indx) .and. ChemoRegion(occupancy(x2,y,z)%indx)) then
				g(1) = (C(x2,y,z) - C(x1,y,z))/(2*DELTA_X)
			else
				g(1) = MISSING_VAL
			endif
			y1 = y - 1
			y2 = y + 1
			if (y1 < 1 .or. y2 > NY) then
				g(2) = 0
!			elseif (occupancy(x,y1,z)%indx(1) >= 0 .and. occupancy(x,y2,z)%indx(1) >= 0) then
			elseif (ChemoRegion(occupancy(x,y1,z)%indx) .and. ChemoRegion(occupancy(x,y2,z)%indx)) then
				g(2) = (C(x,y2,z) - C(x,y1,z))/(2*DELTA_X)
			else
				g(2) = MISSING_VAL
			endif
			z1 = z - 1
			z2 = z + 1
			if (z1 < 1 .or. z2 > NZ) then
				g(3) = 0
!			elseif (occupancy(x,y,z1)%indx(1) >= 0 .and. occupancy(x,y,z2)%indx(1) >= 0) then
			elseif (ChemoRegion(occupancy(x,y,z1)%indx) .and. ChemoRegion(occupancy(x,y,z2)%indx)) then
				g(3) = (C(x,y,z2) - C(x,y,z1))/(2*DELTA_X)
			else
				g(3) = MISSING_VAL
			endif
			grad(:,x,y,z) = g
		enddo
	enddo
enddo
do z = blobrange(3,1),blobrange(3,2)
    do y = blobrange(2,1),blobrange(2,2)
        do x = blobrange(1,1),blobrange(1,2)
!			if (occupancy(x,y,z)%indx(1) < 0) cycle	! outside or DC
			if (.not.ChemoRegion(occupancy(x,y,z)%indx)) cycle
			do i = 1,3
				if (grad(i,x,y,z) == MISSING_VAL) then
					missed = .true.
					grad(i,x,y,z) = 0
					do k = 1,6
						xx = x + neumann(1,k)
						yy = y + neumann(2,k)
						zz = z + neumann(3,k)
						if (outside_xyz(xx,yy,zz)) cycle
!						if (occupancy(xx,yy,zz)%indx(1) < 0) cycle	! outside or DC
						if (.not.ChemoRegion(occupancy(xx,yy,zz)%indx)) cycle
						if (grad(i,xx,yy,zz) /= MISSING_VAL) then
							grad(i,x,y,z) = grad(i,xx,yy,zz)
							missed = .false.
							exit
						endif
					enddo
					if (missed) then
!						write(*,*) 'Missing gradient at: ',x,y,z,grad(:,x,y,z)
					endif
				endif
			enddo
		enddo
	enddo
enddo
end subroutine

!----------------------------------------------------------------------------------------
! This subroutine is called from the GUI, and it passes back the chemokine gradient info
! needed to size the gradient_array and interpret the data.
!----------------------------------------------------------------------------------------
subroutine get_gradient_info(chem_used, ntsites) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_gradient_info
use, intrinsic :: iso_c_binding
integer(c_int) :: chem_used(*), ntsites
integer :: i, x, y, z, ns
logical :: halve = .true.

do i = 1,4
    if (chemo(i)%used) then
        chem_used(i) = 1
    else
        chem_used(i) = 0
    endif
enddo
!ntsites = nsites
ns = 0
do z = blobrange(3,1),blobrange(3,2)
    if (halve .and. mod(z,2) == 0) cycle
    do y = blobrange(2,1),blobrange(2,2)
        if (halve .and. mod(y,2) == 0) cycle
        do x = blobrange(1,1),blobrange(1,2)
            if (halve .and. mod(x,2) == 0) cycle
!			if (occupancy(x,y,z)%indx(1) < 0) cycle	! outside or DC
			if (.not.ChemoRegion(occupancy(x,y,z)%indx)) cycle	! outside or DC
			ns = ns+1
        enddo
    enddo
enddo
ntsites = ns
end subroutine

!----------------------------------------------------------------------------------------
! The gradients are stored in a 1-D array of size = ntsites*(3 + nchem_used*3).
! Here we can check the ntsites value.
!----------------------------------------------------------------------------------------
subroutine get_gradients(chem_used, ntsites, gradient_array, use_strength) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_gradients
use, intrinsic :: iso_c_binding
integer(c_int) :: chem_used(*), ntsites, use_strength
real(c_float) :: gradient_array(*)
integer :: x, y, z, i, j, k, ic, ns
real :: strength
logical :: halve = .true.

ns = 0
k = 0
do z = blobrange(3,1),blobrange(3,2)
    if (halve .and. mod(z,2) == 0) cycle
    do y = blobrange(2,1),blobrange(2,2)
        if (halve .and. mod(y,2) == 0) cycle
        do x = blobrange(1,1),blobrange(1,2)
            if (halve .and. mod(x,2) == 0) cycle
!			if (occupancy(x,y,z)%indx(1) < 0) cycle	! outside or DC
			if (.not.ChemoRegion(occupancy(x,y,z)%indx)) cycle	! outside or DC
			ns = ns+1
			k = k+1
			gradient_array(k) = x
			k = k+1
			gradient_array(k) = y
			k = k+1
			gradient_array(k) = z
			do ic = 1,4
		        if (use_strength == 1) then
		            strength = receptor(ic)%strength
		        else
		            strength = 1
		        endif
		        do j = 1,3
		            k = k+1
    			    if (chemo(ic)%used) then
                        gradient_array(k) = strength*chemo(ic)%grad(j,x,y,z)
                    else
                        gradient_array(k) = 0
                    endif
                enddo
            enddo
        enddo
    enddo
enddo
if (ns /= ntsites) then
    write(logmsg,*) 'Error: get_gradients: inconsistent site count: ',ns,ntsites
    call logger(logmsg)
    stop
endif
!write(nflog,'(3f6.1,12f8.4)') gradient_array(1:ntsites*15)
end subroutine

!----------------------------------------------------------------------------------------
! This subroutine is called from the GUI, and it passes back the chemokine gradient info
! needed to size the gradient_array and interpret the data.
! The argument axis takes values 1,2,3
! 1 = Y-Z plane (normal to X axis)
! 2 = X-Z plane (normal to Y axis)
! 3 = X-Y plane (normal to Z axis)
! The argument fraction is the fractional distance along the blob radius (-1,1)
!----------------------------------------------------------------------------------------
subroutine get_gradient2D_info(chem_used, ntsites, axis, fraction) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_gradient2d_info
use, intrinsic :: iso_c_binding
integer(c_int) :: chem_used(*), ntsites, axis
real(c_float) :: fraction
integer :: i, x, y, z, ns
integer rng(3,2), rad(3)
logical :: halve = .false.

rng = blobrange
rad(1) = Radius%x
rad(2) = Radius%y
rad(3) = Radius%z
rng(axis,:) = Centre(axis) + fraction*rad(axis)
do i = 1,4
    if (chemo(i)%used) then
        chem_used(i) = 1
    else
        chem_used(i) = 0
    endif
enddo
ns = 0
do z = rng(3,1),rng(3,2)
    if (halve .and. axis /= 3 .and. mod(z,2) == 0) cycle
    do y = rng(2,1),rng(2,2)
        if (halve .and. axis /= 2 .and. mod(y,2) == 0) cycle
        do x = rng(1,1),rng(1,2)
            if (halve .and. axis /= 1 .and. mod(x,2) == 0) cycle
!		    if (occupancy(x,y,z)%indx(1) < 0) cycle	! outside or DC
		    if (.not.ChemoRegion(occupancy(x,y,z)%indx)) cycle	! outside or DC
		    ns = ns+1
        enddo
    enddo
enddo
ntsites = ns
end subroutine

!----------------------------------------------------------------------------------------
! The gradients are stored in a 1-D array of size = ntsites*(3 + nchem_used*3).
! Here we can check the ntsites value.
!----------------------------------------------------------------------------------------
subroutine get_gradients2D(chem_used, ntsites, gradient_array, axis, fraction, use_strength) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_gradients2d
use, intrinsic :: iso_c_binding
integer(c_int) :: chem_used(*), ntsites, axis, use_strength
real(c_float) :: gradient_array(*), fraction
integer :: x, y, z, i, j, k, ic, ns
integer rng(3,2), rad(3)
real :: strength
logical :: halve = .false.

rng = blobrange
rad(1) = Radius%x
rad(2) = Radius%y
rad(3) = Radius%z
rng(axis,:) = Centre(axis) + fraction*rad(axis)
ns = 0
k = 0
do z = rng(3,1),rng(3,2)
    if (halve .and. axis /= 3 .and. mod(z,2) == 0) cycle
    do y = rng(2,1),rng(2,2)
        if (halve .and. axis /= 2 .and. mod(y,2) == 0) cycle
        do x = rng(1,1),rng(1,2)
            if (halve .and. axis /= 1 .and. mod(x,2) == 0) cycle
!		    if (occupancy(x,y,z)%indx(1) < 0) cycle	! outside or DC
		    if (.not.ChemoRegion(occupancy(x,y,z)%indx)) cycle	! outside or DC
		    ns = ns+1
		    k = k+1
		    gradient_array(k) = x
		    k = k+1
		    gradient_array(k) = y
		    k = k+1
		    gradient_array(k) = z
		    do ic = 1,4
		        if (use_strength == 1) then
		            strength = receptor(ic)%strength
		        else
		            strength = 1
		        endif
	            do j = 1,3
	                k = k+1
			        if (chemo(ic)%used) then
                        gradient_array(k) = strength*chemo(ic)%grad(j,x,y,z)
                    else
                        gradient_array(k) = 0
                    endif
                enddo
            enddo
        enddo
    enddo
enddo
if (ns /= ntsites) then
    write(logmsg,*) 'Error: get_gradients: inconsistent site count: ',ns,ntsites
    call logger(logmsg)
    stop
endif
!write(nflog,'(3f6.1,12f8.4)') gradient_array(1:ntsites*15) 
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
logical function ChemoRegion(indx)
integer :: indx(2)

if (indx(1) == OUTSIDE_TAG .or. (.not.USE_CELL_SITES .and. indx(1) < 0)) then
	ChemoRegion = .false.
else
	ChemoRegion = .true.
endif
end function

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine CheckGradient
integer :: x, y, z
real :: g(3)
logical :: err = .false.

y = Centre(2)
z = Centre(3)
do x = Centre(1)+1,NX
	if (occupancy(x,y,z)%indx(1) == OUTSIDE_TAG) cycle
	g = chemo(1)%grad(:,x,y,z)
	if (g(1) < 0) then
		write(*,*) 'CheckGradient: ',x,g(1)
		err = .true.
	endif
enddo
if (err) stop
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
