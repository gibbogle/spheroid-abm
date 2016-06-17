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

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
function f_deriv(t,v)
!use rksuite_90_prec, only:wp
real(REAL_KIND), intent(in) :: t
real(REAL_KIND), dimension(:), intent(in) :: v
real(REAL_KIND), dimension(size(v,1)) :: f_deriv

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
	f_deriv(i) = dCsum + dCreact
enddo
end function

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
! Reactions here + cross-membrane diffusion
!----------------------------------------------------------------------------------
subroutine intra_react1(nchemo,C,Ce,dCreact)
integer :: nchemo
real(REAL_KIND) :: C(:), Ce(:), dCreact(:)
integer :: ichemo
real(REAL_KIND) :: metab

if (C(OXYGEN) > ODEdiff%C1_soft) then
	metab = (C(OXYGEN)-ODEdiff%deltaC_soft)/(chemo(OXYGEN)%MM_C0 + C(OXYGEN) - ODEdiff%deltaC_soft)
elseif (C(OXYGEN) > 0) then
	metab = ODEdiff%k_soft*C(OXYGEN)*C(OXYGEN)
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
if (ichemo == OXYGEN) then
	if (C(OXYGEN) > ODEdiff%C1_soft) then
		metab = (C(OXYGEN)-ODEdiff%deltaC_soft)/(chemo(OXYGEN)%MM_C0 + C(OXYGEN) - ODEdiff%deltaC_soft)
	elseif (C(OXYGEN) > 0) then
		metab = ODEdiff%k_soft*C(OXYGEN)*C(OXYGEN)
	else
		metab = 0
	endif
	dCreact = -metab*chemo(ichemo)%max_cell_rate*1.0e6/vol	! convert mass rate (mol/s) to concentration rate (mM/s)
elseif (ichemo == GLUCOSE) then
	dCreact = -chemo(ichemo)%max_cell_rate*1.0e6/vol	! convert mass rate (mol/s) to concentration rate (mM/s)
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Note: updating of concentrations is now done in extendODEdiff()
!-----------------------------------------------------------------------------------------
subroutine cell_divider1(kcell0)
integer :: kcell0
integer :: kpar=0
integer :: j, k, kcell1, site0(3), site1(3), site2(3), site(3), npath, path(3,200)
real(REAL_KIND) :: tnow, R
type (boundary_type), pointer :: bdry

tnow = istep*DELTA_T
site0 = cell_list(kcell0)%site
if (dbug) write(*,*) 'cell_divider: ',kcell0,site0,occupancy(site0(1),site0(2),site0(3))%indx
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
Nsites = Nsites + 1
call SetRadius(Nsites)
call extendODEdiff(site2)

if (dbug) write(*,*) 'did push_path'
cell_list(kcell0)%t_divide_last = tnow
!cell_list(kcell0)%t_divide_next = tnow + DivideTime()
cell_list(kcell0)%volume = cell_list(kcell0)%volume/2
R = par_uni(kpar)
cell_list(kcell0)%divide_volume = Vdivide0 + dVdivide*(2*R-1)
call add_cell(kcell0,kcell1,site0)	! daughter cell is placed at site0

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
    write(logmsg,'(a,3i4,i6)') 'Added site is not bdry: ',site2,occupancy(site2(1),site2(2),site2(3))%indx(1)
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
!		    call SetBdryConcs(site)
		endif
	else
		if (dbug) write(*,*) 'not isbdry'
		if (bdrylist_present(site,bdrylist)) then	! remove it
			if (dbug) write(*,*) 'present, remove it'
			call bdrylist_delete(site,bdrylist)
			nullify(occupancy(site(1),site(2),site(3))%bdry)
			nbdry = nbdry - 1
!			call ResetConcs(site)
		endif
	endif
enddo
if (dbug) write(*,*) 'done!'
end subroutine

!-----------------------------------------------------------------------------------------
! For a given cell site site0(:), find a suitable free site on the boundary of
! the blob, and the path of sites leading from site0 to site1
! NOT USED
!-----------------------------------------------------------------------------------------
subroutine get_nearbdrypath(site0,path,npath)
integer :: site0(3), path(3,200),npath
integer :: site1(3), jump(3), site(3), k, j, jmax, kpar=0
real(REAL_KIND) :: v(3), r, d, dmax, v_aim(3), v_try(3)
real(REAL_KIND) :: cosa, sina, d_try, del

v = site0 - Centre
do j = 1,3
	v(j) = v(j) + (par_uni(kpar) - 0.5)
enddo
r = norm(v)
v_aim = v/r
k = 1
site1 = site0
path(:,k) = site1
do 
	dmax = 0
	do j = 1,27
		if (j == 14) cycle
		jump = jumpvec(:,j)
		site = site1 + jump
		if (occupancy(site(1),site(2),site(3))%indx(1) < -100) cycle
		v_try = site - Centre
!		call get_vnorm(v,v_try)
		d_try = norm(v_try)
		d = dot_product(v_try,v_aim)
		cosa = d/d_try
		sina = sqrt(1 - cosa*cosa)
		del = d_try*sina
		if (d-del > dmax) then
			dmax = d-del
			jmax = j
		endif
	enddo
	site1 = site1 + jumpvec(:,jmax)
	k = k+1
	path(:,k) = site1
	if (occupancy(site1(1),site1(2),site1(3))%indx(1) == OUTSIDE_TAG) exit
enddo
npath = k

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
subroutine MakeSplit(force)
logical :: force
integer :: k, wsum, kdomain, nsum, Ntot, N, last, x, y, z
integer, allocatable :: scount(:)
integer, allocatable :: wz(:), ztotal(:)
integer :: Mslices
real(REAL_KIND) :: dNT, diff1, diff2
logical :: show = .false.

!write(*,*) 'MakeSplit: istep,Mnodes: ',istep,Mnodes
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
integer :: ifdcstate, ibcstate, ctype, stage, region, highlight=0
integer :: last_id1, last_id2
logical :: ok
integer, parameter :: axis_centre = -2	! identifies the ellipsoid centre
integer, parameter :: axis_end    = -3	! identifies the ellipsoid extent in 5 directions
integer, parameter :: axis_bottom = -4	! identifies the ellipsoid extent in the -Y direction, i.e. bottom surface
integer, parameter :: ninfo = 6			! the size of the info package for a cell (number of integers)
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
		call cellColour(kcell,highlight,col)
		BC_list(j+1) = kcell + last_id1
		BC_list(j+2:j+4) = site
		BC_list(j+5) = rgb(col)
		last_id2 = kcell + last_id1
	endif
enddo
nBC_list = last_id2
end subroutine

subroutine ComputeCexCin
real(REAL_KIND) :: K1, K2, K2K1, C0, a, b, c, Cex, Cin, D2, D
integer :: i

K1 = chemo(OXYGEN)%membrane_diff*Vsite
K2 = chemo(OXYGEN)%max_cell_rate*1.0e6
K2K1 = K2/K1
C0 = chemo(OXYGEN)%MM_C0
do i = 1,37
	if (i < 10) then
		Cex = i*0.0001
	elseif (i < 19) then
		Cex = i*0.001
	else
		Cex = (i-9)*0.01
	endif
	b = K2K1 + C0 - Cex
	c = -C0*Cex
	D2 = b*b - 4*c
	D = sqrt(D2)
	Cin = (D - b)/2
	write(nfout,'(3e12.4)') Cex,Cin,Cex-Cin
enddo	
end subroutine

!----------------------------------------------------------------------------------
! To implement IRKC solver.
! The basic idea is that at every grid point (site) there are Nc extracellular
! variables and Nc intracellular, i.e. 2*Nc in all.  At a vacant site the Nc
! intracellular concentrations will all be zero, and dy/dt = 0.
! The PDE system is split into two parts.  For the diffusion part, which is solved
! by explicit RK, the derivatives are given by F_E, while the reaction part, which
! involves no coupling between variables at different grid points, is given by F_I.
! We just need to take care of the mapping of variable indices.
! For now, stick with the approximate method of holding all other constituents fixed
! while solving for one constituent.
! The numbering of variables for IRKC (for a single constituent) needs to be:
! site  extra  intra
!   1     1      2
!   2     3      4
!   3     5      6
! ....
!   k   2k-1    2k
! ....
!   Ns  2Ns-1   2Ns
!
! where Ns is the number of sites.  Effectively, NPDES = 2.
! From the site index (1,..,Ns) we need to be able to determine:
!   if the site is vacant
!   if vacant, extracellular concentrations
!   if occupied, extra- and intracellular concentrations.
!----------------------------------------------------------------------------------
subroutine irkc_solver(it,tstart,dt,nc)
use irkc_m
integer :: it, nc
real(REAL_KIND) :: tstart, dt
integer :: ic, ichemo, isite, ivar, nextra
real(REAL_KIND) :: t, tend
real(REAL_KIND), allocatable :: y(:)
real(REAL_KIND), parameter :: rtol = 1.d-3, atol = rtol/10
TYPE(IRKC_SOL) :: SOL

nextra = ODEdiff%nextra
allocate(y(2*nextra))
do ic = 1,nchemo
	ichemo = chemomap(ic)
	! initialize y(:) from allstate(:,:)
	! this is inefficient because we are using the variable storage mechanism that was developed for RKC
	! it could be worth simplifying this (and changing ODEdiff) to store intracellular values even for
	! vacant sites
!	write(*,*) 'irk_solver: ',ichemo,ODEdiff%nvars,2*nextra
!	write(*,*) 'allstate'
!	write(*,'(10f7.3)') allstate(1:ODEdiff%nvars,ichemo)
	y = 0
	do ivar = 1,ODEdiff%nvars
		if (ODEdiff%vartype(ivar) == EXTRA) then
			isite = ODEdiff%extra_isite(ivar)
			y(2*isite-1) = allstate(ivar,ichemo)
		else
			isite = ODEdiff%extra_isite(ivar-1)
			y(2*isite) = allstate(ivar,ichemo)
		endif
	enddo
	ichemo_sol = ichemo
	t = tstart
	tend = t + dt
	SOL = IRKC_SET(t,y,tend,RE=rtol,AE=atol,NPDES=2,ONE_STEP=.false.)
	call IRKC(SOL,F_E,F_I)
	do ivar = 1,ODEdiff%nvars
		if (ODEdiff%vartype(ivar) == EXTRA) then
			isite = ODEdiff%extra_isite(ivar)
			allstate(ivar,ichemo) = SOL%y(2*isite-1)
		else
			isite = ODEdiff%extra_isite(ivar-1)
			allstate(ivar,ichemo) = SOL%y(2*isite)
		endif
	enddo
enddo

end subroutine

!----------------------------------------------------------------------------------
! Diffusion + decay for extracellular variables
!----------------------------------------------------------------------------------
subroutine F_E(neqn,t,y,dy)
integer :: neqn
real(REAL_KIND) ::  t, y(neqn), dy(neqn)
integer :: iex, isite, ivar, ichemo, k, kv(7)
real(REAL_KIND) :: DX2, decay_rate, dc1, dc6, cbnd, dCdiff, dCsum, val, v(7)

! need ichemo
ichemo = ichemo_sol
!write(*,*) 'F_E: ',ichemo,neqn
!write(*,'(10f7.3)') y
DX2 = DELTA_X*DELTA_X
decay_rate = chemo(ichemo)%decay_rate
dc1 = chemo(ichemo)%diff_coef/DX2
dc6 = 6*dc1 + decay_rate
cbnd = BdryConc(ichemo,t_simulation)	
dy = 0
!$omp parallel do private(isite, ivar, k, kv, v, dCdiff, dCsum, val ) default(shared) schedule(static)
do iex = 1,neqn-1,2
	isite = (iex + 1)/2
	ivar = ODEdiff%isite_extra(isite)
    dCsum = 0
    do k = 1,7
		kv(k) = ODEdiff%icoef(ivar,k)
	    if (k == 1) then
		    dCdiff = -dc6
	    else
		    dCdiff = dc1
	    endif
	    if (kv(k) == 0) then
		    val = cbnd
	    else
		    val = y(kv(k))
	    endif
	    v(k) = val
	    dCsum = dCsum + dCdiff*val
    enddo
	dy(iex) = dCsum - decay_rate*y(iex)
!	write(*,'(8i6)') isite,kv
!	write(*,'(i4,7f10.2)') iex,v
enddo
end subroutine

!----------------------------------------------------------------------------------
! Reactions + membrane transfer + decay for intracellular
! Membrane transfer for extracellular (decay handled in F_E)
! This version is for a single constituent (ichemo_sol), all others held constant
! => NPDES = 2 (1=extra, 2 = intra)
! Note: index i is the starting index of the grid-point in the vector of length neqn
!----------------------------------------------------------------------------------
subroutine F_I(i,npdes,t,y,dy,want_jac,jac)
integer :: i, npdes
real(REAL_KIND) :: t, y(npdes), dy(npdes),jac(npdes,npdes)
logical :: want_jac
integer :: isite, iex, kcell, ichemo, ict, Ng
real(REAL_KIND) :: Kmemb, membrane_flux, Vin, Vex, C, metab, dmetab, decay_rate, dCreact, Cin(MAX_CHEMO)

isite = (i+npdes-1)/npdes
iex = ODEdiff%isite_extra(isite)	! extra variable index
kcell = ODEdiff%cell_index(iex)	! cell index
if (kcell == 0) then	! vacant site
	dy = 0
	if (want_jac) then
		jac = 0
	endif
	return
endif
ichemo = ichemo_sol
!write(*,*) 'F_I: ',ichemo,y
!write(*,*) i,isite,ivar,kcell
!if (y(1) == 0 .or. y(2) == 0) stop
Cin = allstate(iex+1,:)		! could use iin = %isite_intra
decay_rate = chemo(ichemo)%decay_rate
Vex = Vextra	! use the constant value for now
Vin = Vsite - Vex
Kmemb = chemo(ichemo)%membrane_diff*Vsite
membrane_flux = Kmemb*(y(1) - y(2))
dy(1) = -membrane_flux/Vex
!if (dy(1) > 0) then
!	write(*,'(a,i4,3e12.3)') 'dy(1) > 0: ',isite,dy(1),y
!endif
ict = cell_list(kcell)%celltype
C = y(2)
!SN30K_metabolized = (SN30K%Kmet0(ict) > 0)	
if (ichemo == OXYGEN) then
	metab = O2_metab(C)
	dCreact = (-metab*chemo(OXYGEN)%max_cell_rate*1.0e6 + membrane_flux)/Vin	! convert mass rate (mol/s) to concentration rate (mM/s)
elseif (ichemo == GLUCOSE) then
	Ng = chemo(GLUCOSE)%Hill_N
	metab = C**Ng/(chemo(GLUCOSE)%MM_C0**Ng + C**Ng)
	dCreact = (-metab*chemo(ichemo)%max_cell_rate*1.0e6 + membrane_flux)/Vin	! convert mass rate (mol/s) to concentration rate (mM/s)
elseif (ichemo == TRACER) then
	dCreact = membrane_flux/Vin
elseif (ichemo == SN30000) then
    if ((SN30K%Kmet0(ict) > 0) .and. C > 0) then
		dCreact = -(SN30K%C1(ict) + SN30K%C2(ict)*SN30K%KO2(ict)/(SN30K%KO2(ict) + Cin(OXYGEN)))*SN30K%Kmet0(ict)*C
	else
		dCreact = 0
	endif
	dCreact = dCreact + membrane_flux/Vin
elseif (ichemo == SN30000_METAB) then
	if ((SN30K%Kmet0(ict) > 0) .and. Cin(SN30000) > 0) then
		dCreact = (SN30K%C1(ict) + SN30K%C2(ict)*SN30K%KO2(ict)/(SN30K%KO2(ict) + Cin(OXYGEN)))*SN30K%Kmet0(ict)*Cin(SN30000)
	else
		dCreact = 0
	endif
	dCreact = dCreact + membrane_flux/Vin
endif
dy(2) = dCreact - C*decay_rate
!if (dy(2) > 0) then
!	write(*,'(a,i4,3e12.3)') 'dy(2) > 0: ',isite,dy(2),y
!endif
!write(*,'(2i6,2e12.3)') ichemo,isite,dy
if (want_jac) then
	jac(1,1) = -Kmemb/Vex
	jac(1,2) = Kmemb/Vex
	jac(2,1) = Kmemb/Vin
	if (ichemo == OXYGEN) then
		dmetab = (O2_metab(1.01*C) - metab)/(0.01*C)
		jac(2,2) = -decay_rate - Kmemb/Vin - dmetab*chemo(OXYGEN)%max_cell_rate*1.0e6/Vin
	elseif (ichemo == GLUCOSE) then
		dmetab = ((1.01*C)**Ng/(chemo(GLUCOSE)%MM_C0**Ng + (1.01*C)**Ng) - metab)/(0.01*C)
		jac(2,2) = -decay_rate - Kmemb/Vin - dmetab*chemo(ichemo)%max_cell_rate*1.0e6/Vin
	elseif (ichemo == TRACER) then
		jac(2,2) = -decay_rate - Kmemb/Vin
	elseif (ichemo == SN30000) then
		jac(2,2) = -decay_rate - Kmemb/Vin	
		if ((SN30K%Kmet0(ict) > 0) .and. C > 0) then
			jac(2,2) = jac(2,2) - (SN30K%C1(ict) + SN30K%C2(ict)*SN30K%KO2(ict)/(SN30K%KO2(ict) + Cin(OXYGEN)))*SN30K%Kmet0(ict)
		endif
	elseif (ichemo == SN30000_METAB) then
		jac(2,2) = -decay_rate - Kmemb/Vin	
		if ((SN30K%Kmet0(ict) > 0) .and. C > 0) then
			jac(2,2) = jac(2,2) - (SN30K%C1(ict) + SN30K%C2(ict)*SN30K%KO2(ict)/(SN30K%KO2(ict) + Cin(OXYGEN)))*SN30K%Kmet0(ict)
		endif
	endif
!	write(*,'(a,i4,4e12.3)') 'jac: ',ichemo,jac
endif		
end subroutine

!----------------------------------------------------------------------------------
! Intracellular reactions, membrane diffusion, and decay.
! The number of variables, neqn = 2*nchemo, and the constituents are given by
! chemomap(:).  The first nchemo constituents are extracellular, the last intracellular.
! icase = i, the allstate variable index, for access to cell parameters
! NOT USED
!----------------------------------------------------------------------------------
subroutine f_rkc_intra(neqn,t,y,dydt,icase)
integer :: neqn, icase
real(REAL_KIND) :: t, y(neqn), dydt(neqn)
integer :: i, k, ic, ichemo, kcell, ict
real(REAL_KIND) :: decay_rate, membrane_diff, dCreact, Vc
real(REAL_KIND) :: Cin(MAX_CHEMO), Cex(MAX_CHEMO), C, Cox, Csn
logical :: bnd, metabolised
real(REAL_KIND) :: metab, dMdt

Vc = Vsite_cm3 - Vextra_cm3
do ic = 1,nchemo
	ichemo = chemomap(ic)
	Cex(ic) = y(ic)
	Cin(ic) = y(ic+nchemo)
	if (ichemo == OXYGEN) then
		Cox = Cin(ic)
	elseif (ichemo == SN30000) then
		Csn = Cin(ic)
	endif
enddo
do ic = 1,nchemo
	ichemo = chemomap(ic)
	decay_rate = chemo(ichemo)%decay_rate
!	decay_rate = 0
	membrane_diff = chemo(ichemo)%membrane_diff
!	membrane_diff = 0
	! extracellular
	dydt(ic) = -membrane_diff*(Cex(ic) - Cin(ic)) - decay_rate*Cex(ic)
	! intracellular
    kcell = ODEdiff%cell_index(icase)	! for access to cell-specific parameters
    ict = cell_list(kcell)%celltype
    metabolised = (SN30K%Kmet0(ict) > 0)
	C = Cin(ic)
	dCreact = 0
	if (ichemo == OXYGEN) then
		if (C > ODEdiff%C1_soft) then
			metab = (C-ODEdiff%deltaC_soft)/(chemo(OXYGEN)%MM_C0 + C - ODEdiff%deltaC_soft)
		elseif (C > 0) then
			metab = ODEdiff%k_soft*C*C
		else
			metab = 0
			write(*,*) 'metab = 0'
			stop
		endif
		metab = 0.001*metab
		dCreact = -metab*chemo(ichemo)%max_cell_rate*1.0e6/Vc	! convert mass rate (mol/s) to concentration rate (mM/s)
	elseif (ichemo == GLUCOSE) then
		metab = C/(chemo(GLUCOSE)%MM_C0 + C)
		dCreact = -metab*chemo(ichemo)%max_cell_rate*1.0e6/Vc	! convert mass rate (mol/s) to concentration rate (mM/s)
	elseif (ichemo == TRACER) then
		dCreact = 0
	elseif (ichemo == SN30000) then
	    if (metabolised .and. C > 0) then
			dCreact = -(1 - SN30K%C2(ict) + SN30K%C2(ict)*SN30K%KO2(ict)/(SN30K%KO2(ict) + Cox))*SN30K%Kmet0(ict)*C
		endif
	elseif (ichemo == SN30000_METAB) then
		if (metabolised .and. Csn > 0) then
			dCreact = (1 - SN30K%C2(ict) + SN30K%C2(ict)*SN30K%KO2(ict)/(SN30K%KO2(ict) + Cox))*SN30K%Kmet0(ict)*Csn
		endif
	endif
	dCreact = dCreact + membrane_diff*(Cex(ic) - Cin(ic))	
    dydt(ic+nchemo) = dCreact - Cin(ic)*decay_rate
enddo
if (icase == -1) then
	write(*,'(4e12.3)') dydt(1:4)
endif
end subroutine

!----------------------------------------------------------------------------------
! Reactions here + cross-membrane diffusion
! Should metab depend on cell volume, stage in cell cycle?
! NOT USED now
!----------------------------------------------------------------------------------
subroutine IntraReact(ichemo,Cin,Cex,vol,dCreact)
integer :: ichemo
real(REAL_KIND) :: Cin(:), Cex, vol, dCreact
real(REAL_KIND) :: metab, C, Kex, dCex
integer :: ict=1

C = Cin(ichemo)
dCreact = 0
if (ichemo == OXYGEN) then
    if (C > ODEdiff%C1_soft) then
	    metab = (C-ODEdiff%deltaC_soft)/(chemo(OXYGEN)%MM_C0 + C - ODEdiff%deltaC_soft)
    elseif (C > 0) then
	    metab = ODEdiff%k_soft*C*C
    else
	    metab = 0
	    write(*,*) 'metab = 0'
    endif
    dCreact = -metab*chemo(ichemo)%max_cell_rate*1.0e6/vol	! convert mass rate (mol/s) to concentration rate (mM/s)
elseif (ichemo == GLUCOSE) then
	metab = C/(chemo(GLUCOSE)%MM_C0 + C)
    dCreact = -metab*chemo(ichemo)%max_cell_rate*1.0e6/vol	! convert mass rate (mol/s) to concentration rate (mM/s)
elseif (ichemo == TRACER) then
	dCreact = 0
elseif (ichemo == SN30000) then
    dCreact = -(1 - SN30K%C2(ict) + SN30K%C2(ict)*SN30K%KO2(ict)/(SN30K%KO2(ict) + Cin(OXYGEN)))*SN30K%Kmet0(ict)*C
elseif (ichemo == SN30000_METAB) then
    dCreact = (1 - SN30K%C2(ict) + SN30K%C2(ict)*SN30K%KO2(ict)/(SN30K%KO2(ict) + Cin(OXYGEN)))*SN30K%Kmet0(ict)*Cin(SN30000)
endif
dCreact = dCreact + chemo(ichemo)%membrane_diff*(Cex - C)*Vsite_cm3
end subroutine

!----------------------------------------------------------------------------------------
! A new site is to be added to the blob.  The site is chosen so as to maintain the
! ellipsoidal shape.  For each current bdry site, P, imagine a line drawn from the blob centre,
! O, through the site, and determine the point Q where this line intersects the surface of the
! ellipsoid with a = aRadius, b = bRadius.  The distance from the site to the intersection
! point is calculated, and the bdry site that maximizes this distance (the site most inside the
! ellipsoid is chosen as the starting point of the search for a site to add.
! If r = site - Centre, we can parametrise the line by t.r, where t is a scalar:
! x = t.r(1), y = t.r(2), z = t.r(3)
! and t = 1 gives the point P.
! On the ellipsoid with radii a and b:
! x^2/a^2 + (y^2 + z^2)/b^2 = 1
! which implies that t^2[r(1)^2/a^2 + (r(2)^2 + r(3)^2)/b^2] = 1
! from which we obtain t^2.
! Next the 'outside' neighbours of the chosen bdry site are examined, and the one that
! is closest to the centre is selected as the new 'inside' blob site.
! The square of the length of OP is dp^2 = r(1)^2 + r(2)^2 + r(3)^2, and the square of the length
! of OQ is dq^2 = t^2.dp^2, dq = t.dp, and the distance PQ = (t-1).dp.  We want to choose the
! bdry site for which PQ is a maximum.
! Finally the bdry sites in the vicinity are tested to see if any need to be removed from
! the bdrylist. (This is why it makes sense to use a linked list).
! Removal of bdry sites is done by FixBdrylist.
! NOT USED
!----------------------------------------------------------------------------------------
subroutine AddSite(ok)
logical :: ok
type (boundary_type), pointer :: bdry
integer :: site(3), psite(3), indx(2), k, kmin, x, y, z
real(REAL_KIND) :: r(3), x2, y2, z2, dp, t, tmax, dpq, dmax, dmin
real(REAL_KIND) :: cmin, cmax

!write(*,*) 'AddSite'
dmax = -1.0e10
bdry => bdrylist
do while ( associated ( bdry )) 
    site = bdry%site
    r = site - Centre
    x2 = r(1)**2
    y2 = r(2)**2
    z2 = r(3)**2
    dp = sqrt(x2 + y2 + z2)
    t = sqrt(1/(x2/Radius**2 + y2/Radius**2 + z2/Radius**2))
    dpq = (t-1)*dp
    if (dpq > dmax) then
        psite = site
        dmax = dpq
        tmax = t
    endif
    bdry => bdry%next
enddo
if (.not.isbdry(psite))then
    write(logmsg,*) 'psite is not a bdry site'
	call logger(logmsg)
    stop
endif

! Now we need to look at all 26 neighbours of psite, to determine which 'outside' neighbour
! is closest to O.
dmin = 1.0e10
kmin = 0
do k = 1,27
	if (k == 14) cycle
	site = psite + jumpvec(:,k)
	indx = occupancy(site(1),site(2),site(3))%indx
    if (indx(1) /= OUTSIDE_TAG) cycle
    r = site - Centre
    x2 = r(1)**2
    y2 = r(2)**2
    z2 = r(3)**2
    dp = sqrt(x2 + y2 + z2)
    if (dp < dmin) then
        dmin = dp
        kmin = k
    endif
enddo    
if (kmin == 0) then
    write(logmsg,*) 'Error: no outside neighbours of bdry site'
	call logger(logmsg)
    ok = .false.
    return
endif

! This is the site to convert to 'inside'
site = psite + jumpvec(:,kmin)
occupancy(site(1),site(2),site(3))%indx = 0
Nsites = Nsites + 1
call SetRadius(Nsites)
if (isbdry(site)) then   ! add it to the bdrylist
    allocate(bdry)
    bdry%site = site
!    bdry%chemo_influx = .false.
    nullify(bdry%next)
    call bdrylist_insert(bdry,bdrylist)
    occupancy(site(1),site(2),site(3))%bdry => bdry
!    call SetBdryConcs(site)
else
    write(logmsg,*) 'Added site is not bdry: ',site,occupancy(site(1),site(2),site(3))%indx
	call logger(logmsg)
    stop
endif
ok = .true.
return
end subroutine

!----------------------------------------------------------------------------------------
! The logic is similar to AddSite.  We look for the bdry site P which is most outside
! the ellipsoid, i.e. to minimise OQ-OP, where Q is the point at
! which the line through OP intersects the ellipsoid surface.  We then choose the 'inside'
! neighbour of P that is furthest from O.
! NOT USED
!----------------------------------------------------------------------------------------
subroutine RemoveSite(ok)
logical :: ok
type (boundary_type), pointer :: bdry
integer :: site(3), psite(3), k, kmax
real(REAL_KIND) :: r(3), x2, y2, z2, dp, t, dpq, dmax, dmin

dmin = 1.0e10
psite = 0
bdry => bdrylist
do while ( associated ( bdry )) 
    site = bdry%site
    r = site - Centre
    x2 = r(1)**2
    y2 = r(2)**2
    z2 = r(3)**2
    dp = sqrt(x2 + y2 + z2)
    t = sqrt(1/(x2/Radius**2 + y2/Radius**2 + z2/Radius**2))
    dpq = (t-1)*dp
    if (dpq < dmin) then
        psite = site
        dmin = dpq
    endif
    bdry => bdry%next
enddo
if (psite(1) == 0) then
    write(logmsg,*) 'No more bdry sites: ',nbdry
	call logger(logmsg)
    ok = .false.
    return
endif
if (.not.isbdry(psite)) then
    write(logmsg,*) 'Not a bdry site: ',psite
	call logger(logmsg)
    ok = .false.
    return
endif

! Now we need to look at psite and its 26 neighbours, to determine which 'inside' neighbour
! is most distant from O.
dmax = -1.0e10
kmax = 0
do k = 1,27
	site = psite + jumpvec(:,k)
    if (occupancy(site(1),site(2),site(3))%indx(1) <= OUTSIDE_TAG) cycle   ! outside or unreachable
    r = site - Centre
    x2 = r(1)**2
    y2 = r(2)**2
    z2 = r(3)**2
    dp = sqrt(x2 + y2 + z2)
    if (dp > dmax) then
        dmax = dp
        kmax = k
    endif
enddo    
if (kmax == 0) then
    write(logmsg,*) 'Error: no inside neighbours of bdry site'
	call logger(logmsg)
    ok = .false.
    return
endif
! This is the site to convert to 'outside'
site = psite + jumpvec(:,kmax)
call ClearSite(site)
occupancy(site(1),site(2),site(3))%indx = OUTSIDE_TAG
Nsites = Nsites - 1
call SetRadius(Nsites)
if (associated(occupancy(site(1),site(2),site(3))%bdry)) then   ! remove it from the bdrylist
    call bdrylist_delete(site, bdrylist)
    nullify(occupancy(site(1),site(2),site(3))%bdry)
! Need to check for a new bdry site to replace the one removed
! Look at all the neighbours
    psite = site
    do k = 1,27
	    site = psite + jumpvec(:,k)
        if (occupancy(site(1),site(2),site(3))%indx(1) <= OUTSIDE_TAG) cycle   ! outside or unreachable
        if (associated(occupancy(site(1),site(2),site(3))%bdry)) cycle   ! bdry
        if (isbdry(site)) then
            allocate(bdry)
            bdry%site = site
!            bdry%chemo_influx = .false.
            nullify(bdry%next)
            call bdrylist_insert(bdry,bdrylist)
            occupancy(site(1),site(2),site(3))%bdry => bdry
!            call SetBdryConcs(site)
        endif
    enddo
endif
ok = .true.
end subroutine

!--------------------------------------------------------------------------------
! NOT USED
!--------------------------------------------------------------------------------
subroutine AddSites(n,ok)
integer :: n
logical :: ok
integer :: k

do k = 1,n
    call AddSite(ok)
    if (.not.ok) return
    call FixBdrylist
enddo
ok = .true.
end subroutine

!--------------------------------------------------------------------------------
! NOT USED
!--------------------------------------------------------------------------------
subroutine RemoveSites(n,ok)
integer :: n
logical :: ok
integer :: k

do k = 1,n
    call RemoveSite(ok)
    if (.not.ok) return
    call FixBdrylist
enddo
ok = .true.
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine TestAddRemoveSites
integer :: i, n, k, nb
logical :: ok
type (boundary_type), pointer :: bdry

write(*,*) 'TestAddRemoveSites'
n = 10000
do i = 1,10
    write(*,*) 'Adding sites'
    call AddSites(n,ok)
    write(*,*) 'nbdry = ',nbdry
    write(*,*) 'Removing sites'
    call RemoveSites(n,ok)
    write(*,*) 'nbdry = ',nbdry
enddo
call CreateCheckList
end subroutine

!----------------------------------------------------------------------------------------- 
! Initialise medium concentrations, etc.
! NOT USED
!-----------------------------------------------------------------------------------------
subroutine SetupMedium1
integer :: ichemo
real(REAL_KIND) :: V, V0, R1

call SetRadius(Nsites)
R1 = Radius*DELTA_X			! cm
V0 = medium_volume0			! cm3
V = V0 - (4./3.)*PI*R1**3	! cm3
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	chemo(ichemo)%medium_Cext = chemo(ichemo)%bdry_conc
	chemo(ichemo)%medium_Cbnd = chemo(ichemo)%bdry_conc
	chemo(ichemo)%medium_M = V*chemo(ichemo)%bdry_conc
	chemo(ichemo)%medium_U = 0
	write(*,'(a,i4,2e12.4)') 'SetupMedium: ',ichemo,chemo(ichemo)%medium_Cext,chemo(ichemo)%medium_M
enddo
end subroutine

!----------------------------------------------------------------------------------
! The medium concentrations are updated explicitly, assuming a sphere with boundary
! concentrations equal to the mean extracellular concentrations of boundary sites.
! Note that concentrations of O2 and glucose are not varied.
! NOT USED
!----------------------------------------------------------------------------------
subroutine UpdateMedium_old(ntvars,state,dt)
integer :: ntvars
real(REAL_KIND) :: dt, state(:,:)
integer :: nb, i, k
real(REAL_KIND) :: Csurface(MAX_CHEMO), F(MAX_CHEMO), area, C_A
logical :: bnd

if (.not.chemo(DRUG_A)%used .and. .not.chemo(DRUG_B)%used) return
! First need the spheroid radius
call SetRadius(Nsites)
! Now compute the mean boundary site concentrations Cbnd(:)
nb = 0
Csurface = 0
do i = 1,ntvars
	if (ODEdiff%vartype(i) /= EXTRA) cycle
	bnd = .false.
	do k = 1,7
		if (ODEdiff%icoef(i,k) < 0) then
			bnd = .true.
			exit
		endif
	enddo
	if (bnd) then
		nb = nb + 1
		Csurface = Csurface + state(i,:)
	endif
enddo
Csurface = Csurface/nb
area = 4*PI*Radius*Radius*DELTA_X*DELTA_X
F(:) = area*chemo(:)%diff_coef*(Csurface(:) - chemo(:)%bdry_conc)/DELTA_X
C_A = chemo(DRUG_A)%bdry_conc
chemo(DRUG_A:MAX_CHEMO)%bdry_conc = (chemo(DRUG_A:MAX_CHEMO)%bdry_conc*medium_volume + F(DRUG_A:MAX_CHEMO)*dt)/medium_volume
chemo(DRUG_A:MAX_CHEMO)%bdry_conc = chemo(DRUG_A:MAX_CHEMO)%bdry_conc*(1 - dt*chemo(DRUG_A:MAX_CHEMO)%decay_rate)
end subroutine

!-----------------------------------------------------------------------------------------
! Reduce medium volume and concentrations to account for the growth of the spheroid
! by one site.
! THIS WILL BE SUPERCEDED
!-----------------------------------------------------------------------------------------
subroutine AdjustMedium1
real(REAL_KIND) :: total(MAX_CHEMO)

total(DRUG_A:MAX_CHEMO) = chemo(DRUG_A:MAX_CHEMO)%bdry_conc*medium_volume
medium_volume = medium_volume - Vsite_cm3
chemo(DRUG_A:MAX_CHEMO)%bdry_conc = total(DRUG_A:MAX_CHEMO)/medium_volume
end subroutine


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
			TPZ%membrane_diff_in(im) = TPZ%membrane_diff_in(im)*Vsite_cm3/60		! /min -> /sec
			read(nfcell,*) TPZ%membrane_diff_out(im)
			TPZ%membrane_diff_out(im) = TPZ%membrane_diff_out(im)*Vsite_cm3/60	! /min -> /sec
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
			DNB%membrane_diff_in(im) = DNB%membrane_diff_in(im)*Vsite_cm3/60		! /min -> /sec
			read(nfcell,*) DNB%membrane_diff_out(im)
			DNB%membrane_diff_out(im) = DNB%membrane_diff_out(im)*Vsite_cm3/60	! /min -> /sec
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

#if 0
!-----------------------------------------------------------------------------------------
! Cells can be tagged to die, or finally die of anoxia, or they can be tagged for death 
! at division time if the drug is effective.
!-----------------------------------------------------------------------------------------
subroutine CellDeath1(dt,ok)
real(REAL_KIND) :: dt
logical :: ok
integer :: kcell, ict, nlist0, site(3), i, im, ityp, kpar=0 
real(REAL_KIND) :: C_O2, kmet, Kd, dMdt, kill_prob, tnow
logical :: use_TPZ_DRUG, use_DNB_DRUG

!call logger('CellDeath')
ok = .true.
use_TPZ_DRUG = chemo(TPZ_DRUG)%used
use_DNB_DRUG = chemo(DNB_DRUG)%used
tnow = istep*DELTA_T	! seconds
nlist0 = nlist
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	ityp = cell_list(kcell)%celltype
	call getO2conc(kcell,C_O2)
	if (cell_list(kcell)%anoxia_tag) then
!		write(logmsg,*) 'anoxia_tag: ',kcell,cell_list(kcell)%state,tnow,cell_list(kcell)%t_anoxia_die
!		call logger(logmsg)
		if (tnow >= cell_list(kcell)%t_anoxia_die) then
!			call logger('cell dies')
			call CellDies(kcell)
			Nanoxia_dead(ityp) = Nanoxia_dead(ityp) + 1
!			if (cell_list(kcell)%drugA_tag) then
!				NdrugA_tag(ityp) = NdrugA_tag(ityp) - 1
!			endif
!			if (cell_list(kcell)%drugB_tag) then
!				NdrugB_tag(ityp) = NdrugB_tag(ityp) - 1
!			endif
!			if (cell_list(kcell)%radiation_tag) then
!				Nradiation_tag(ityp) = Nradiation_tag(ityp) - 1
!			endif
			cycle
		endif
	else
		if (C_O2 < ANOXIA_THRESHOLD) then
			cell_list(kcell)%t_hypoxic = cell_list(kcell)%t_hypoxic + dt
			if (cell_list(kcell)%t_hypoxic > t_anoxic_limit) then
				cell_list(kcell)%anoxia_tag = .true.						! tagged to die later
				cell_list(kcell)%t_anoxia_die = tnow + anoxia_death_delay	! time that the cell will die
				Nanoxia_tag(ityp) = Nanoxia_tag(ityp) + 1
			endif
		else
			cell_list(kcell)%t_hypoxic = 0
		endif
	endif
	if (use_TPZ_DRUG .and. .not.cell_list(kcell)%drugA_tag) then	
		ict = cell_list(kcell)%celltype
		Kd = TPZ%Kd(ict)
	    kmet = (1 - TPZ%C2(ict,0) + TPZ%C2(ict,0)*TPZ%KO2(ict,0)/(TPZ%KO2(ict,0) + C_O2))*TPZ%Kmet0(ict,0)
	    dMdt = kmet*cell_list(kcell)%conc(TPZ_DRUG)
	    if (TPZ%kill_model(ict) == 1) then
		    kill_prob = Kd*dMdt*dt
	    elseif (TPZ%kill_model(ict) == 2) then
		    kill_prob = Kd*dMdt*cell_list(kcell)%conc(TPZ_DRUG)*dt
	    elseif (TPZ%kill_model(ict) == 3) then
		    kill_prob = Kd*dMdt**2*dt
	    elseif (TPZ%kill_model(ict) == 4) then
		    kill_prob = Kd*cell_list(kcell)%conc(TPZ_DRUG)*dt
	    elseif (TPZ%kill_model(ict) == 5) then
		    kill_prob = Kd*cell_list(kcell)%conc(TPZ_DRUG)**2*dt
		endif
	    if (par_uni(kpar) < kill_prob) then
            cell_list(kcell)%drugA_tag = .true.
            NdrugA_tag(ityp) = NdrugA_tag(ityp) + 1
!            write(nflog,'(a,2i6)') 'TPZ tagged: ',kcell,ict
		endif
	endif
	if (use_DNB_DRUG .and. .not.cell_list(kcell)%drugB_tag) then
		ict = cell_list(kcell)%celltype
!	    kmet = (1 - TPZ%C2(ict,0) + TPZ%C2(ict,0)*TPZ%KO2(ict,0)/(TPZ%KO2(ict,0) + C_O2))*TPZ%Kmet0(ict,0)
!	    dMdt = kmet*cell_list(kcell)%conc(TPZ_DRUG)
	    if (DNB%kill_model(ict,1) < 4 .or. DNB%kill_model(ict,2) < 4) then
			write(logmsg,*) 'Error: CellDeath: DNB kill model must be 4 or 5, not: ',DNB%kill_model(ict,1),DNB%kill_model(ict,2)
			call logger(logmsg)
			ok = .false.
			return
		endif
		kill_prob = 0
		do im = 1,2
			if (DNB%kill_model(ict,im) == 4) then
				kill_prob = kill_prob + DNB%Kd(ict,im)*cell_list(kcell)%conc(DNB_DRUG + im)*dt
			elseif (DNB%kill_model(ict,im) == 5) then
				kill_prob = kill_prob + DNB%Kd(ict,im)*(cell_list(kcell)%conc(DNB_DRUG + im)**2)*dt
			endif
		enddo
	    if (par_uni(kpar) < kill_prob) then
            cell_list(kcell)%drugB_tag = .true.
            NdrugB_tag(ityp) = NdrugB_tag(ityp) + 1
            write(nflog,'(a,2i6)') 'DNB tagged: ',kcell,ict
		endif
	endif
enddo
end subroutine
#endif


!-----------------------------------------------------------------------------------------
! Live cells = Ncells
! Drug deaths = Ndrugdead
! Hypoxia deaths = Ndead - Ndrugdead
! Total tagged for drug death on division = NdrugA_tag, NdrugB_tag
! Current tagged = Ntodie - Ntagdead 
!-----------------------------------------------------------------------------------------
subroutine get_summary1(summaryData,i_hypoxia_cutoff,i_growth_cutoff) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*), i_hypoxia_cutoff,i_growth_cutoff
integer :: Nviable(MAX_CELLTYPES), Nlive(MAX_CELLTYPES), plate_eff_10(MAX_CELLTYPES)
integer :: diam_um, vol_mm3_1000, nhypoxic(3), ngrowth(3), hypoxic_percent_10, growth_percent_10, necrotic_percent_10,  npmm3, &
    medium_oxygen_1000, medium_glucose_1000, medium_drug_1000(2), &
    bdry_oxygen_1000, bdry_glucose_1000, bdry_drug_1000(2)
integer :: TNanoxia_dead, TNradiation_dead, TNdrug_dead(2),  &
           Ntagged_anoxia(MAX_CELLTYPES), Ntagged_radiation(MAX_CELLTYPES), Ntagged_drug(2,MAX_CELLTYPES), &
           TNtagged_anoxia, TNtagged_radiation, TNtagged_drug(2)
integer :: Tplate_eff_10   
integer :: ityp
real(REAL_KIND) :: diam_cm, vol_cm3, vol_mm3, hour, plate_eff(MAX_CELLTYPES), necrotic_fraction
real(REAL_KIND) :: cmedium(MAX_CHEMO), cbdry(MAX_CHEMO)

if (i_hypoxia_cutoff == 0) stop

hour = istep*DELTA_T/3600.
call getDiamVol(diam_cm,vol_cm3)
vol_mm3 = vol_cm3*1000				! volume in mm^3
vol_mm3_1000 = vol_mm3*1000			! 1000 * volume in mm^3
diam_um = diam_cm*10000
npmm3 = Ncells/vol_mm3

Ntagged_anoxia(:) = Nanoxia_tag(:)			! number currently tagged by anoxia
Ntagged_radiation(:) = Nradiation_tag(:)	! number currently tagged by radiation
Ntagged_drug(1,:) = Ndrug_tag(1,:)			! number currently tagged by drugA
Ntagged_drug(2,:) = Ndrug_tag(2,:)			! number currently tagged by drugA

TNtagged_anoxia = sum(Ntagged_anoxia(1:Ncelltypes))
TNtagged_radiation = sum(Ntagged_radiation(1:Ncelltypes))
TNtagged_drug(1) = sum(Ntagged_drug(1,1:Ncelltypes))
TNtagged_drug(2) = sum(Ntagged_drug(2,1:Ncelltypes))

TNanoxia_dead = sum(Nanoxia_dead(1:Ncelltypes))
TNradiation_dead = sum(Nradiation_dead(1:Ncelltypes))
TNdrug_dead(1) = sum(Ndrug_dead(1,1:Ncelltypes))
TNdrug_dead(2) = sum(Ndrug_dead(2,1:Ncelltypes))

call getHypoxicCount(nhypoxic)
hypoxic_percent_10 = (1000*nhypoxic(i_hypoxia_cutoff))/Ncells	! 10* %hypoxic
call getGrowthCount(ngrowth)
growth_percent_10 = (1000*ngrowth(i_growth_cutoff))/Ncells
!if (TNanoxia_dead > 0) then
	call getNecroticFraction(necrotic_fraction,vol_cm3)
!else
!	necrotic_fraction = 0
!endif
necrotic_percent_10 = 1000*necrotic_fraction
call getNviable(Nviable, Nlive)
do ityp = 1,Ncelltypes
	if (Nlive(ityp) > 0) then
		plate_eff(ityp) = real(Nviable(ityp))/Nlive(ityp)
	else
		plate_eff(ityp) = 0
	endif
enddo
plate_eff_10 = 1000*plate_eff
Tplate_eff_10 = 0
do ityp = 1,Ncelltypes
	Tplate_eff_10 = Tplate_eff_10 + plate_eff_10(ityp)*celltype_fraction(ityp)
enddo
call getMediumConc(cmedium, cbdry)
medium_oxygen_1000 = cmedium(OXYGEN)*1000
medium_glucose_1000 = cmedium(GLUCOSE)*1000
medium_drug_1000(1) = cmedium(DRUG_A)*1000
medium_drug_1000(2) = cmedium(DRUG_B)*1000
bdry_oxygen_1000 = cbdry(OXYGEN)*1000
bdry_glucose_1000 = cbdry(GLUCOSE)*1000
bdry_drug_1000(1) = cbdry(DRUG_A)*1000
bdry_drug_1000(2) = cbdry(DRUG_B)*1000

summaryData(1:25) = [ istep, Ncells, TNanoxia_dead, TNdrug_dead(1), TNdrug_dead(2), TNradiation_dead, &
    TNtagged_anoxia, TNtagged_drug(1), TNtagged_drug(2), TNtagged_radiation, &
	diam_um, vol_mm3_1000, hypoxic_percent_10, growth_percent_10, necrotic_percent_10, Tplate_eff_10, npmm3, &
	medium_oxygen_1000, medium_glucose_1000, medium_drug_1000(1), medium_drug_1000(2), &
	bdry_oxygen_1000, bdry_glucose_1000, bdry_drug_1000(1), bdry_drug_1000(2) ]
write(nfres,'(a,a,2a12,i8,2e12.4,19i7,17e12.4)') trim(header),' ',gui_run_version, dll_run_version, &
	istep, hour, vol_mm3, diam_um, Ncells_type(1:2), &
    Nanoxia_dead(1:2), Ndrug_dead(1,1:2), &
    Ndrug_dead(2,1:2), Nradiation_dead(1:2), &
    Ntagged_anoxia(1:2), Ntagged_drug(1,1:2), &
    Ntagged_drug(2,1:2), Ntagged_radiation(1:2), &
	nhypoxic(:)/real(Ncells), ngrowth(:)/real(Ncells), &
	necrotic_fraction, plate_eff(1:2), &
	cmedium(OXYGEN), cmedium(GLUCOSE), cmedium(DRUG_A), cmedium(DRUG_B), &
	cbdry(OXYGEN), cbdry(GLUCOSE), cbdry(DRUG_A), cbdry(DRUG_B)
		
call sum_dMdt(GLUCOSE)

if (diam_count_limit > LIMIT_THRESHOLD) then
	if (Ncells > diam_count_limit) limit_stop = .true.
elseif (diam_count_limit > 0) then
	if (diam_um > diam_count_limit) limit_stop = .true.
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Cell growth, death and division are handled here.  Division occurs when cell volume 
! exceeds the divide volume. 
! As the cell grows we need to adjust both Cin and Cex to maintain mass conservation.
!-----------------------------------------------------------------------------------------
subroutine CellGrowth1(dt,ok)
real(REAL_KIND) :: dt
logical :: ok
integer :: kcell, nlist0, site(3), ityp, idrug, kpar=0
integer :: divide_list(10000), ndivide, i
real(REAL_KIND) :: tnow, C_O2, C_glucose, metab, metab_O2, metab_glucose, dVdt, vol0, r_mean(2), c_rate(2), R
real(REAL_KIND) :: Vin_0, Vex_0, dV
real(REAL_KIND) :: Cin_0(MAX_CHEMO), Cex_0(MAX_CHEMO)
character*(20) :: msg
logical :: drugkilled, glucose_growth
integer :: C_option = 1
type(cell_type), pointer :: cp

!call logger('CellGrowth')
ok = .true.
nlist0 = nlist
tnow = istep*DELTA_T
c_rate(1:2) = log(2.0)/divide_time_mean(1:2)		! Note: to randomise divide time need to use random number, not mean!
r_mean(1:2) = Vdivide0/(2*divide_time_mean(1:2))
glucose_growth = chemo(GLUCOSE)%controls_growth
ndivide = 0
do kcell = 1,nlist0
cp => cell_list(kcell)
	if (cell_list(kcell)%state == DEAD) cycle
	ityp = cell_list(kcell)%celltype
	C_O2 = cell_list(kcell)%conc(OXYGEN)
	C_glucose = cell_list(kcell)%conc(GLUCOSE)
	metab_O2 = O2_metab(C_O2)
	metab_glucose = glucose_metab(C_glucose)
	if (glucose_growth) then
		metab = metab_O2*metab_glucose
	else
		metab = metab_O2
	endif
!	if (use_constant_divide_volume) then
!		dVdt = metab*Vdivide0/(2*cell_list(kcell)%divide_time)
!	else
!		if (use_V_dependence) then
!			dVdt = c_rate(ityp)*cell_list(kcell)%volume*metab
!		else
!			dVdt = r_mean(ityp)*metab
!		endif
!	endif
	dVdt = get_dVdt(cp,metab)
	if (suppress_growth) then	! for checking solvers
		dVdt = 0
	endif
	site = cell_list(kcell)%site
	Cin_0 = cell_list(kcell)%conc
	Cex_0 = occupancy(site(1),site(2),site(3))%C
	cell_list(kcell)%dVdt = dVdt
	Vin_0 = cell_list(kcell)%volume*Vcell_cm3	! cm^3
	Vex_0 = Vsite_cm3 - Vin_0					! cm^3
	dV = dVdt*dt*Vcell_cm3						! cm^3
	cell_list(kcell)%volume = (Vin_0 + dV)/Vcell_cm3
	if (C_option == 1) then
		! Calculation based on transfer of an extracellular volume dV with constituents, i.e. holding extracellular concentrations constant
		cell_list(kcell)%conc = (Vin_0*Cin_0 + dV*Cex_0)/(Vin_0 + dV)
		occupancy(site(1),site(2),site(3))%C = (Vex_0*Cex_0 - dV*Cex_0)/(Vex_0 - dV)	! = Cex_0
	elseif (C_option == 2) then
		! Calculation based on change in volumes without mass transfer of constituents
		cell_list(kcell)%conc = Vin_0*Cin_0/(Vin_0 + dV)
		occupancy(site(1),site(2),site(3))%C = Vex_0*Cex_0/(Vex_0 - dV)
	endif
	if (cell_list(kcell)%volume > cell_list(kcell)%divide_volume) then	! time to divide
		if (cell_list(kcell)%radiation_tag) then
			R = par_uni(kpar)
!			if (cell_list(kcell)%generation > 1) then
!				write(*,'(a,2i6,2f8.4)') 'Radiation-tagged cell: generation, R, p_death: ',kcell,cell_list(kcell)%generation,R,cell_list(kcell)%p_rad_death
!			endif
			if (R < cell_list(kcell)%p_rad_death) then
				call CellDies(kcell)
				Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
!				do idrug = 1,ndrugs_used
!					if (cell_list(kcell)%drug_tag(idrug)) then
!						Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) - 1
!					endif
!				enddo
!				if (cell_list(kcell)%anoxia_tag) then
!					Nanoxia_tag(ityp) = Nanoxia_tag(ityp) - 1
!				endif
				cycle
			endif
		endif
		drugkilled = .false.
		do idrug = 1,ndrugs_used
			if (cell_list(kcell)%drug_tag(idrug)) then
				R = par_uni(kpar)
				if (R < cell_list(kcell)%p_drug_death(idrug)) then
					call CellDies(kcell)
					Ndrug_dead(idrug,ityp) = Ndrug_dead(idrug,ityp) + 1
					drugkilled = .true.
					exit
				endif
			endif
		enddo
		if (drugkilled) cycle
	    ndivide = ndivide + 1
	    divide_list(ndivide) = kcell
	endif
enddo
do i = 1,ndivide
    kcell = divide_list(i)
	kcell_dividing = kcell
	call CellDivider(kcell, ok)
	if (.not.ok) return
enddo
end subroutine