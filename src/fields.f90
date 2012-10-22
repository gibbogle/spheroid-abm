! Chemokine concentration fields

module fields
use global
use chemokine
use bdry_linked_list
use boundary
use ode_diffuse
implicit none
save

type Carray_type
	type(chemokine_type), pointer :: Cptr
end type

contains

!----------------------------------------------------------------------------------------
! Convert halflife in hours to a decay rate /sec
!----------------------------------------------------------------------------------------
real function DecayRate(halflife)
real(REAL_KIND) :: halflife

DecayRate = log(2.0)/(halflife*60*60)
end function

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine SetupChemo
integer :: ichemo

chemo(OXYGEN)%name = 'Oxygen'
chemo(OXYGEN)%used = .true.
chemo(OXYGEN)%use_secretion = .false.
chemo(OXYGEN)%bdry_rate = 0
chemo(OXYGEN)%bdry_conc = 100
chemo(OXYGEN)%diff_coef = 2.0e-5	! units = cm^2s^-1
chemo(OXYGEN)%halflife = 0.01	! hours

chemo(GLUCOSE)%name = 'Glucose'
chemo(GLUCOSE)%used = .true.
chemo(GLUCOSE)%use_secretion = .false.
chemo(GLUCOSE)%bdry_rate = 0
chemo(GLUCOSE)%bdry_conc = 50
chemo(GLUCOSE)%diff_coef = 2.0e-6	! units = cm^2s^-1
chemo(GLUCOSE)%halflife = 0.1	! hours

chemo(TRACER)%name = 'Tracer'
chemo(TRACER)%used = .true.
chemo(TRACER)%use_secretion = .false.
chemo(TRACER)%bdry_rate = 0
chemo(TRACER)%bdry_conc = 200
chemo(TRACER)%diff_coef = 2.0e-6	! units?
chemo(TRACER)%halflife = 1.0	! hours
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	chemo(ichemo)%decay_rate = DecayRate(chemo(ichemo)%halflife)
	write(*,*) 'decay_rate: ',ichemo,chemo(ichemo)%decay_rate
enddo
call AllocateConcArrays
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine AllocateConcArrays
integer :: ic

write(logmsg,*) 'AllocateConcArrays'
call logger(logmsg)
do ic = 1,MAX_CHEMO
	if (chemo(ic)%used) then
		allocate(chemo(ic)%conc(NX,NY,NZ))
		allocate(chemo(ic)%grad(3,NX,NY,NZ))
		chemo(ic)%conc = chemo(ic)%bdry_conc	
	endif
enddo
end subroutine

!----------------------------------------------------------------------------------------
! Set the concentrations at a site when it is on the boundary
!----------------------------------------------------------------------------------------
subroutine SetBdryConcs(site)
integer :: site(3)
integer :: ichemo

do ichemo = 1,MAX_CHEMO
	if (chemo(ichemo)%used) then
		chemo(ichemo)%conc(site(1),site(2),site(3)) = chemo(ichemo)%bdry_conc
	endif
enddo
end subroutine

!----------------------------------------------------------------------------------------
! Reset the concentrations at a site that is no longer on the boundary
!----------------------------------------------------------------------------------------
subroutine ResetConcs(site0)
integer :: site0(3)
integer :: site(3), ichemo, k, n
real :: csum(MAX_CHEMO)

csum = 0
n = 0
do k = 1,27
	if (k == 14) cycle
	site = site0 + jumpvec(:,k)
	if (occupancy(site(1),site(2),site(3))%indx(1) == OUTSIDE_TAG) cycle
!	if (associated(occupancy(site(1),site(2),site(3))%bdry)) cycle
	n = n+1
	do ichemo = 1,MAX_CHEMO
		if (chemo(ichemo)%used) then
			csum(ichemo) = csum(ichemo) + chemo(ichemo)%conc(site(1),site(2),site(3))
		endif
	enddo
enddo

do ichemo = 1,MAX_CHEMO
	if (chemo(ichemo)%used) then
		chemo(ichemo)%conc(site0(1),site0(2),site0(3)) = csum(ichemo)/n
	endif
enddo
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine ShowConcs
integer :: x, y, z, k, x1, site(3), ichemo
real :: C(MAX_CHEMO,100)
logical :: start

y = Centre(2) 
z = Centre(3)
start = .false.
k = 0
do x = blobrange(1,1),blobrange(1,2)
	if (occupancy(x,y,z)%indx(1) == OUTSIDE_TAG) then
		if (start) then
			exit
		else
			cycle
		endif
	endif
	if (.not.start) x1 = x
	start = .true.
	k = k+1
	do ichemo = 1,MAX_CHEMO
		if (.not.chemo(ichemo)%used) cycle
	    C(ichemo,k) = chemo(ichemo)%conc(x,y,z)
	enddo
enddo
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
    write(*,'(i4,2x,10f7.2)') istep/60,C(ichemo,1:10)
enddo
!if (C(1) < 100) then
!	write(logmsg,*) 'First conc < 100: istep: ',istep, C(1),occupancy(x1,y,z)%indx(1)
!	call logger(logmsg)
!	site = (/x1,y,z/)
!	write(logmsg,'(a,L)') 'isbdry: ',isbdry(site)
!	call logger(logmsg)
!	stop
!endif
end subroutine

!----------------------------------------------------------------------------------------
! Advance the dynamic solution for chemokine concentration fields through time interval dtstep.
! For now only treat the case of NOT use_ODE_diffusion, because use_ODE_diffusion is a bit
! more complicated (should store derivatives).
! This version solves for all the consitituents at once - needed for interactions.
!----------------------------------------------------------------------------------------
subroutine UpdateFields(dtstep)
real :: dtstep
type(chemokine_type), pointer :: Cptr
!type(Carray_type) :: Carray(MAX_CHEMO)
integer :: ichemo, it, nt = 4
real :: Kdiffusion(MAX_CHEMO), Kdecay(MAX_CHEMO), dt, dCmax_par(0:20), dCmax
integer :: z1, z2, n, kpar
integer :: x, y, z
real, allocatable :: C_par(:,:,:,:)

x = Centre(1)
y = Centre(2)
z = Centre(3)
!write(*,*) 'UpdateFields: conc: ',chemo(OXYGEN)%conc(x,y,z)
dCmax_par = 0
dt = dtstep/nt
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
!	write(*,*) 'UpdateFields: ichemo: ',ichemo
!	Carray(ichemo)%cptr => chemo(ichemo)
	Kdiffusion(ichemo) = chemo(ichemo)%diff_coef
	Kdecay(ichemo) = chemo(ichemo)%decay_rate
enddo

do it = 1,nt
	!$omp parallel do private(ichemo,z1,z2,n,C_par,dt,dCmax)
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
		allocate(C_par(MAX_CHEMO,NX,NY,n))
		call par_evolve(Kdiffusion,Kdecay,C_par,z1,z2,dt,dCmax,kpar)
		do ichemo = 1,MAX_CHEMO
			if (.not.chemo(ichemo)%used) cycle
			chemo(ichemo)%conc(:,:,z1:z2) = C_par(ichemo,:,:,1:n)
		enddo
		deallocate(C_par)
		dCmax_par(kpar) = dCmax
	enddo
enddo
if (mod(istep,60) == 0) then
	write(*,'(a,8f8.4)') 'dCmax: ',dCmax_par(0:Mnodes-1)
endif
!do ichemo = 1,MAX_CHEMO
!	if (.not.chemo(ichemo)%used) cycle
!	Cptr => chemo(ichemo)
!	call gradient(Cptr%conc,Cptr%grad)
!enddo
end subroutine

!----------------------------------------------------------------------------------------
! Updates concentrations through a timestep dt, solving the diffusion-decay eqtn by a 
! simple explicit method.  This assumes concentration boundary conditions.
!----------------------------------------------------------------------------------------
subroutine par_evolve(Kdiffusion,Kdecay,Ctemp,z1,z2,dt,dCmax,kpar)
integer :: z1, z2, kpar
!real :: C(:,:,:)
real :: Ctemp(:,:,:,:)
real :: Kdiffusion(:), Kdecay(:), dt, dCmax, Cnew
real :: sum, dV, C0(MAX_CHEMO), dMdt(MAX_CHEMO), dCdt(MAX_CHEMO)
integer :: x, y, z, zpar, xx, yy, zz, nb, k, indx(2), ichemo
logical :: source_site

!write(*,*) 'par_evolve: ',ichemo,kpar,z1,z2,dt
dCmax = 0
dV = DELTA_X**3
do zpar = 1,z2-z1+1
	z = zpar + z1-1
    do y = blobrange(2,1),blobrange(2,2)
        do x = blobrange(1,1),blobrange(1,2)
		    indx = occupancy(x,y,z)%indx
			if (.not.ChemoRegion(indx)) cycle	! site not accessible by chemokines
            source_site = .false.
            if (associated(occupancy(x,y,z)%bdry)) then
                ! Check for chemo bdry site - no change to the concentration at such a site
                do ichemo = 1,MAX_CHEMO
					if (.not.chemo(ichemo)%used) cycle
		            chemo(ichemo)%conc(x,y,z) = chemo(ichemo)%bdry_conc
			        Ctemp(ichemo,x,y,zpar) = chemo(ichemo)%bdry_conc
				enddo
		        source_site = .true.
            endif
            if (.not.source_site) then
                dMdt = 0
                C0 = 0
				do ichemo = 1,MAX_CHEMO
					if (.not.chemo(ichemo)%used) cycle
!					C0 = C(x,y,z)
					C0(ichemo) = chemo(ichemo)%conc(x,y,z)
					sum = 0
					nb = 0
					do k = 1,6
						xx = x + neumann(1,k)
						yy = y + neumann(2,k)
						zz = z + neumann(3,k)
						if (outside_xyz(xx,yy,zz)) cycle
						indx = occupancy(xx,yy,zz)%indx
						if (.not.ChemoRegion(indx)) cycle	! site not accessible by chemokines
						nb = nb + 1
!						sum = sum + C(xx,yy,zz)
						sum = sum + chemo(ichemo)%conc(xx,yy,zz)
					enddo
					dMdt(ichemo) = Kdiffusion(ichemo)*DELTA_X*(sum - nb*C0(ichemo)) - Kdecay(ichemo)*C0(ichemo)*dV ! + influx(x,y,z)
				enddo
				! reactions between chems go here, changing dMdt
				dCdt = 0
				call reactions(C0,dCdt)
				dMdt = dMdt + dCdt*dV
                do ichemo = 1,MAX_CHEMO
					if (.not.chemo(ichemo)%used) cycle
					Cnew = (C0(ichemo)*dV + dMdt(ichemo)*dt)/dV
					Ctemp(ichemo,x,y,zpar) = Cnew
					if (ichemo == 1) then
						dCmax = max(dCmax,abs(Cnew-C0(ichemo)))
					endif
				enddo
			endif
		enddo
	enddo
enddo
end subroutine

!----------------------------------------------------------------------------------------
! The processes that need to be accounted for are:
! oxygen uptake by the cell (if there is one
! 
!----------------------------------------------------------------------------------------
subroutine reactions(C0,dCdt)
real :: C0(:), dCdt(:)

dCdt = 0
end subroutine

!----------------------------------------------------------------------------------------
! Advance the dynamic solution for chemokine concentration fields through time interval dtstep.
! For now only treat the case of NOT use_ODE_diffusion, because use_ODE_diffusion is a bit
! more complicated (should store derivatives).
! This version solves for the consitituents one at a time.
!----------------------------------------------------------------------------------------
subroutine UpdateFields1(dtstep)
real :: dtstep
type(chemokine_type), pointer :: Cptr
integer :: ichemo, it, nt = 4
real(REAL_KIND) :: Kdiffusion, Kdecay, dt
integer :: z1, z2, n, kpar
integer :: x, y, z
real(REAL_KIND), allocatable :: C_par(:,:,:)

x = Centre(1)
y = Centre(2)
z = 10
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
			call par_evolve1(ichemo,Cptr%conc,Kdiffusion,Kdecay,C_par,z1,z2,dt,kpar)
			Cptr%conc(:,:,z1:z2) = C_par(:,:,1:n)
!			write(*,*) C_par(x,y,z)
			deallocate(C_par)
		enddo
	enddo
!	write(*,*) 'call gradient'
	call gradient(Cptr%conc,Cptr%grad)
enddo
!stop
end subroutine

!----------------------------------------------------------------------------------------
! Updates concentrations through a timestep dt, solving the diffusion-decay eqtn by a 
! simple explicit method.  This assumes concentration boundary conditions.
!----------------------------------------------------------------------------------------
subroutine par_evolve1(ichemo,C,Kdiffusion,Kdecay,Ctemp,z1,z2,dt,kpar)
integer :: ichemo, z1, z2, kpar
real(REAL_KIND) :: C(:,:,:), Ctemp(:,:,:)
real (REAL_KIND) :: Kdiffusion, Kdecay, dt
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
!	                if (ichemo == i .and. occupancy(x,y,z)%bdry%chemo_influx(i)) then
! We now assume that the whole boundary is a source, at one concentration
	                if (ichemo == i) then
		                C(x,y,z) = chemo(i)%bdry_conc
			            Ctemp(x,y,zpar) = C(x,y,z)
				        source_site = .true.
					endif
				enddo
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
!----------------------------------------------------------------------------------------
subroutine gradient(C,grad)
real(REAL_KIND) :: C(:,:,:), grad(:,:,:,:)
integer :: x, y, z, xx, yy, zz, x1, x2, y1, y2, z1, z2, i, k, indx(2)
real(REAL_KIND) :: g(3)
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
		        do j = 1,3
		            k = k+1
    			    if (chemo(ic)%used) then
                        gradient_array(k) = chemo(ic)%grad(j,x,y,z)
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
rad(1) = Radius
rad(2) = Radius
rad(3) = Radius
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
rad(1) = Radius
rad(2) = Radius
rad(3) = Radius
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
	            do j = 1,3
	                k = k+1
			        if (chemo(ic)%used) then
                        gradient_array(k) = chemo(ic)%grad(j,x,y,z)
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

!if (indx(1) == OUTSIDE_TAG .or. (.not.USE_CELL_SITES .and. indx(1) < 0)) then
if (indx(1) == OUTSIDE_TAG) then
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

end module