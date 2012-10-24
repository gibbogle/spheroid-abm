! Solving 3D diffusion-decay eqtn as an ODE system.
! 
! Each chemokine has a list of sites with a specified rate of secretion into each site.
! The list needs to be updated when the boundary changes, or when FDCs move and/or the number
! changes.

module ode_diffuse

use chemokine
use rkf45s

implicit none

real(REAL_KIND), parameter :: BASE_CHEMOKINE_SECRETION = 0.012

integer :: ivdbug

contains

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
!----------------------------------------------------------------------------------
subroutine SetupODEDiffusion
integer :: x, y, z, i, site(3), ifdc, ichemo, k, nc
real(REAL_KIND) :: DX, DX2, c(MAX_CHEMO,7)
real(REAL_KIND) :: secretion
logical :: left, right

call logger('Set diffusion parameters')
if (.not.allocated(ODEdiff%ivar)) then
	allocate(ODEdiff%ivar(NX,NY,NZ))
endif
if (.not.allocated(ODEdiff%varsite)) then
	allocate(ODEdiff%varsite(NX*NY*NZ,3))
endif
DX = 1.0
DX2 = DX*DX
ODEdiff%ivar = 0
i = 0
do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
			if (occupancy(x,y,z)%indx(1) >= 0) then
				i = i+1
				ODEdiff%ivar(x,y,z) = i
				ODEdiff%varsite(i,:) = (/x,y,z/)
			endif
		enddo
	enddo
enddo
ODEdiff%nvars = i

if (.not.use_ODE_diffusion) return

if (.not.allocated(ODEdiff%icoef)) then
	allocate(ODEdiff%icoef(ODEdiff%nvars,7))
endif
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	if (allocated(chemo(ichemo)%coef)) deallocate(chemo(ichemo)%coef)
	allocate(chemo(ichemo)%coef(ODEdiff%nvars,7))
enddo

do i = 1,ODEdiff%nvars
	c = 0
	nc = 0
	site = ODEdiff%varsite(i,:)
	x = site(1)
	y = site(2)
	z = site(3)
	do ichemo = 1,MAX_CHEMO
		c(ichemo,1) = -chemo(ichemo)%decay_rate
	enddo
	ODEdiff%icoef(i,1) = i
	
	left = .true.
	right = .true.
	if (x==1) then
		left = .false.
	elseif (ODEdiff%ivar(x-1,y,z) == 0) then
		left = .false.
	endif
	if (x==NX) then 
		right = .false.
	elseif (ODEdiff%ivar(x+1,y,z) == 0) then
		right = .false.
	endif
	if (left) then
		nc = nc + 1
		ODEdiff%icoef(i,2) = ODEdiff%ivar(x-1,y,z)
		do ichemo = 1,MAX_CHEMO
			c(ichemo,2) = chemo(ichemo)%diff_coef/DX2
		enddo
	endif
	if (right) then
		nc = nc + 1
		ODEdiff%icoef(i,3) = ODEdiff%ivar(x+1,y,z)
		do ichemo = 1,MAX_CHEMO
			c(ichemo,3) = chemo(ichemo)%diff_coef/DX2
		enddo
	endif
	left = .true.
	right = .true.
	if (y==1) then
		left = .false.
	elseif (ODEdiff%ivar(x,y-1,z) == 0) then
		left = .false.
	endif
	if (y==NY) then 
		right = .false.
	elseif (ODEdiff%ivar(x,y+1,z) == 0) then
		right = .false.
	endif
	if (left) then
		nc = nc + 1
		ODEdiff%icoef(i,4) = ODEdiff%ivar(x,y-1,z)
		do ichemo = 1,MAX_CHEMO
			c(ichemo,4) = chemo(ichemo)%diff_coef/DX2
		enddo
	endif
	if (right) then
		nc = nc + 1
		ODEdiff%icoef(i,5) = ODEdiff%ivar(x,y+1,z)
		do ichemo = 1,MAX_CHEMO
			c(ichemo,5) = chemo(ichemo)%diff_coef/DX2
		enddo
	endif
	left = .true.
	right = .true.
	if (z==1) then
		left = .false.
	elseif (ODEdiff%ivar(x,y,z-1) == 0) then
		left = .false.
	endif
	if (z==NZ) then 
		right = .false.
	elseif (ODEdiff%ivar(x,y,z+1) == 0) then
		right = .false.
	endif
	if (left) then
		nc = nc + 1
		ODEdiff%icoef(i,6) = ODEdiff%ivar(x,y,z-1)
		do ichemo = 1,MAX_CHEMO
			c(ichemo,6) = chemo(ichemo)%diff_coef/DX2
		enddo
	endif
	if (right) then
		nc = nc + 1
		ODEdiff%icoef(i,7) = ODEdiff%ivar(x,y,z+1)
		do ichemo = 1,MAX_CHEMO
			c(ichemo,7) = chemo(ichemo)%diff_coef/DX2
		enddo
	endif
	
	do ichemo = 1,MAX_CHEMO
		c(ichemo,1) = c(ichemo,1) - nc*chemo(ichemo)%diff_coef/DX2
		if (chemo(ichemo)%used) then
			chemo(ichemo)%coef(i,:) = c(ichemo,:) 
		endif
	enddo
!	ODEdiff%ncoef(i) = nc
	
!	if (i == ivdbug) then
!		write(*,*) 'ivdbug: ',nc, ODEdiff%varsite(ivdbug,:)
!		write(*,*) ODEdiff%ivar(x-1,y,z),ODEdiff%ivar(x+1,y,z)
!		write(*,*) ODEdiff%ivar(x,y-1,z),ODEdiff%ivar(x,y+1,z)
!		write(*,*) ODEdiff%ivar(x,y,z-1),ODEdiff%ivar(x,y,z+1)
!		write(*,'(a,7i6)') 'icoef: ',ODEdiff%icoef(i,:)
!		write(*,'(a,7f8.4)') 'coef: ',chemo(1)%coef(i,:)
!	endif
enddo
end subroutine

!----------------------------------------------------------------------------------
! Simple diffusion-decay for a single constituent (ODEdiff%ichemo)
!----------------------------------------------------------------------------------
subroutine deriv(t,v,dv,icase)
real(REAL_KIND) :: t, v(*), dv(*)
integer :: icase
integer :: i, k, kv, n, ichemo
real(REAL_KIND) :: csum, dc, DX2, val
logical :: bnd, dbug

!ichemo = ODEdiff%ichemo
ichemo = icase
DX2 = DELTA_X*DELTA_X
n = ODEdiff%nvars
!if (t < 1.0) write(*,*) icase,t
do i = 1,n
	csum = 0
	do k = 1,7
		if (ODEdiff%icoef(i,k) /= 0) then	! interior
			bnd = .false.
		else								! boundary
			bnd = .true.
		endif
		if (k == 1) then
			dc = -chemo(ichemo)%decay_rate - 6*chemo(ichemo)%diff_coef/DX2	! ODEdiff%ncoef(i)
		else
			dc = chemo(ichemo)%diff_coef/DX2
		endif
		if (bnd) then
			val = chemo(ichemo)%bdry_conc
		else
			kv = ODEdiff%icoef(i,k)
			val = v(kv)
		endif
		csum = csum + dc*val
	enddo
	dv(i) = csum
enddo
end subroutine

!----------------------------------------------------------------------------------
! Now all chemo constituents share the same FD template
!----------------------------------------------------------------------------------
subroutine deriv_all(t,v,dv)
real(REAL_KIND) :: t, v(*), dv(*)
integer :: i, k, kv, n, x, y, z, ichemo, idbug, ifdc, site(3), dx, dy, dz, ic, nf_FDC, nf_MRC
real(REAL_KIND) :: csum(MAX_CHEMO), s, vtemp, ctemp, dc, DX2, val
logical :: bnd, dbug

DX2 = DELTA_X*DELTA_X
n = ODEdiff%nvars
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
				val = chemo(ichemo)%bdry_conc
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
!	site = ODEdiff%varsite(i,:)
!	if (ODEdiff%ichemo /= CXCL13) then
!		if (associated(occupancy(site(1),site(2),site(3))%bdry)) then
!			! Check for chemo bdry site 
!			do ic = 1,MAX_CHEMO
!				if (ODEdiff%ichemo == ic .and. occupancy(site(1),site(2),site(3))%bdry%chemo_influx(ic)) then
!					if (chemo(ic)%use_secretion) then
!						dv(i) = dv(i) + chemo(ic)%bdry_rate
!					else
!						dv(i) = 0
!					endif
!				endif
!			enddo
!		endif
!	else
!		nf_FDC = occupancy(site(1),site(2),site(3))%FDC_nbdry 
!		nf_MRC = occupancy(site(1),site(2),site(3))%MRC_nbdry 
!		if (nf_FDC + nf_MRC > 0) then
!			if (chemo(CXCL13)%use_secretion) then
!				dv(i) = dv(i) + (nf_FDC+nf_MRC)*chemo(CXCL13)%bdry_rate
!			else
!				dv(i) = 0
!			endif
!		endif
!    endif
!enddo

! Secretion by FDCs
!do ifdc = 1,NFDC
!	site = FDC_list(ifdc)%site
!	do dx = -2,2
!		x = site(1)+dx
!		do dy = -2,2
!			y = site(2)+dy
!			do dz = -2,2
!				z = site(3)+dz
!				if (dx*dx+dy*dy+dz*dz > 4.1) cycle
!				if (occupancy(x,y,z)%indx(1) >= 0) then
!					i = ODEdiff%ivar(x,y,z)
!					dv(i) = dv(i) + FDC_list(ifdc)%secretion
!				endif
!			enddo
!		enddo
!	enddo
!enddo

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

write(logmsg,*) 'SolveSteadystate_B: nvars: ',ODEdiff%nvars
call logger(logmsg)
allocate(state(ODEdiff%nvars))
allocate(prev_state(ODEdiff%nvars))
allocate(statep(ODEdiff%nvars))
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
!			call r4_rkf45 ( deriv, ODEdiff%nvars, state, statep, tstart, tend, relerr, abserr, flag )
		else
!			call r8_rkf45 ( deriv, ODEdiff%nvars, state, statep, tstart, tend, relerr, abserr, flag )
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
				if (occupancy(x,y,z)%indx(1) < 0) cycle
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
do i = 1,ODEdiff%nvars
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
! In this version the diffusion/decay of each constituent is solved by a separate
! OMP thread.  Obviously this requires at least as many CPUs as there are constituents.
! Note that this required modifications to the way r4_rkf45 handles SAVEd variables,
! to avoid collisions between different threads.
!----------------------------------------------------------------------------------
subroutine TestODEDiffusion
integer flag, i, j, k, ichemo
real(REAL_KIND) :: tstart, tend, relerr, abserr
real(REAL_KIND), allocatable :: allstate(:,:), allstatep(:,:)
real(REAL_KIND), allocatable :: state(:), statep(:)
real(REAL_KIND), parameter :: dt = 20.0
real(REAL_KIND) :: res(10)
integer :: x, y, z, dx, dy, dz, xmid, ymid, zmid, ibnd, k0, site(3), flag_par(MAX_CHEMO), nchemo, icase
real(REAL_KIND) :: grad(3)
real(REAL_KIND) :: amp, r, ctemp(10)
real(REAL_KIND) :: t1, t2
integer :: nt = 180
logical :: ok

write(logmsg,*) 'TestODEDiffusion: nvars: ',ODEdiff%nvars
call logger(logmsg)
nchemo = MAX_CHEMO
!allocate(prev_state(ODEdiff%nvars))
allocate(allstate(ODEdiff%nvars,nchemo))
allocate(allstatep(ODEdiff%nvars,nchemo))
xmid = NX/2
ymid = NY/2
zmid = NZ/2
k0 = ODEdiff%ivar(xmid,ymid,zmid)
tstart = 0
do ichemo = 1,nchemo
	if (.not.chemo(ichemo)%used) cycle
	call InitState(ichemo,allstate(:,ichemo))
	allstatep(:,ichemo) = 0
	call deriv(tstart,allstate(:,ichemo),allstatep(:,ichemo),ichemo)
!	write(*,'(10f7.3)') allstate(k0:k0+9,ichemo)
enddo

!	ODEdiff%ichemo = ichemo
!	if (.not.chemo(ichemo)%used) cycle

!	call InitState(ichemo,state)
!	statep = 0
!	call deriv(tstart,state,statep)
!	write(*,*) 'did initial deriv'
	
	t1 = wtime()
!	prev_state = 0
	flag_par = 1
	do k = 1,nt

		!$omp parallel do private(tstart, tend, state, statep, flag, relerr, abserr, site, icase)
		do ichemo = 1,nchemo
			if (.not.chemo(ichemo)%used) cycle
!			ODEdiff%ichemo = ichemo
			allocate(state(ODEdiff%nvars))
			allocate(statep(ODEdiff%nvars))
			state = allstate(:,ichemo)
			statep = allstatep(:,ichemo)
			if (ichemo == 1) then
				write(*,'(i6,10f7.3)') k,state(k0:k0+9)
			endif
			if (k == 1) then
				abserr = sqrt ( epsilon ( abserr ) )/2
				relerr = sqrt ( epsilon ( relerr ) )/2
			endif
			tstart = (k-1)*dt
			tend = tstart + dt
			flag = flag_par(ichemo)
!			write(*,*) ichemo, flag, ODEdiff%nvars, tstart, tend
			if (REAL_KIND == 4) then	
				icase = ichemo
				call r4_rkf45 ( deriv, ODEdiff%nvars, state, statep, tstart, tend, relerr, abserr, flag, icase )
			else
	!			call r8_rkf45 ( deriv, ODEdiff%nvars, state, statep, tstart, tend, relerr, abserr, flag )
			endif
			if (flag /= 2) then
				write(logmsg,*) 'Bad flag: ',flag
				call logger(logmsg)
				call r4_rkf45 ( deriv, ODEdiff%nvars, state, statep, tstart, tend, relerr, abserr, flag, icase )
			endif
			flag = 2
			flag_par(ichemo) = flag
	!		do i = 1,ODEdiff%nvars
	!			site = ODEdiff%varsite(i,:)
	!			chemo(ichemo)%conc(site(1),site(2),site(3)) = state(i)
	!		enddo
			allstate(:,ichemo) = state
			allstatep(:,ichemo) = statep
			deallocate(state)
			deallocate(statep)
		enddo
		
!		call CheckConvergence(state,prev_state,ok)
!		if (ok) exit
!		prev_state = state
	enddo

	t2 = wtime()
	write(*,'(a,f10.1)') 'Time: ',t2-t1
	write(*,*) 'site:'
	if (nchemo >= 2) then
		do k = k0-20,k0+20
			site = ODEdiff%varsite(k,:)
			write(*,'(4i6,2f8.3)') k,site,allstate(k,1),allstate(k,2)
		enddo
	elseif (nchemo == 1) then
		do k = k0-20,k0+20
			site = ODEdiff%varsite(k,:)
			write(*,'(4i6,2f8.3)') k,site,allstate(k,1)
		enddo
	endif

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
! Initialise the state vector to the current concentrations
!----------------------------------------------------------------------------------
subroutine InitState(ichemo,state)
integer :: ichemo
real(REAL_KIND) :: state(*)
integer :: x, y, z, i, site(3), nz
real(REAL_KIND) :: smin, smax

write(logmsg,*) 'InitState: ',chemo(ichemo)%name
call logger(logmsg)
smin = 1.0e10
smax = -smin
do i = 1,ODEdiff%nvars
	site = ODEdiff%varsite(i,:)
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
do i = 1,ODEdiff%nvars
	site = ODEdiff%varsite(i,:)
	allstate(i,ichemo) = chemo(ichemo)%conc(site(1),site(2),site(3))
enddo
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

end module

