! Cancer cell state development

module cellstate
use global
use boundary
use fields
use ode_diffuse
implicit none

!real(REAL_KIND), parameter :: Vdivide0 = 1.6
!real(REAL_KIND), parameter :: dVdivide = 0.05
integer :: kcell_dividing = 0

contains

!-----------------------------------------------------------------------------------------
! Need to initialize site and cell concentrations when a cell divides and when there is
! cell death.
!-----------------------------------------------------------------------------------------
subroutine GrowCells(dose,dt,ok)
real(REAL_KIND) :: dose, dt
logical :: ok

!call logger('GrowCells')
ok = .true.
if (use_radiation .and. dose > 0) then
	call Irradiation(dose, ok)
	if (.not.ok) return
endif
!if (use_division) then
	call CellGrowth(dt,ok)
	if (.not.ok) return
!endif
if (use_death) then
	call CellDeath(dt,ok)
	if (.not.ok) return
endif
if (use_migration) then
	call CellMigration(ok)
	if (.not.ok) return
endif
end subroutine

!-----------------------------------------------------------------------------------------
! The O2 concentration to use with cell kcell is either the intracellular concentration,
! or is use_extracellular_O2, the corresponding extracellular concentration
!-----------------------------------------------------------------------------------------
subroutine getO2conc(kcell, C_O2)
integer :: kcell
real(REAL_KIND) :: C_O2
integer :: iv, site(3)
real(REAL_KIND) :: tnow

if (use_extracellular_O2) then
	iv = cell_list(kcell)%iv
	if (iv < 1) then
!		write(logmsg,*) 'getO2conc: ',kcell,site,iv
!		call logger(logmsg)
		tnow = istep*DELTA_T
		C_O2 = BdryConc(OXYGEN,tnow)	! assume that this is a site at the boundary
	else
		C_O2 = allstate(iv-1,OXYGEN)
	endif
else
	C_O2 = cell_list(kcell)%conc(OXYGEN)
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Irradiate cells with dose.
!-----------------------------------------------------------------------------------------
subroutine Irradiation(dose,ok)
real(REAL_KIND) :: dose
logical :: ok
integer :: kcell, site(3), iv, kpar=0
real(REAL_KIND) :: C_O2, OER_alpha_d, OER_beta_d, expon, kill_prob, R

ok = .true.
LQ%OER_am = 2.5
LQ%OER_bm = 3.0
LQ%alpha_H = 0.0473
LQ%beta_H = 0.0017
LQ%K_ms = 4.3e-3	! mM
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	if (cell_list(kcell)%radiation_tag) cycle	! we do not tag twice (yet)
	call getO2conc(kcell,C_O2)
	OER_alpha_d = dose*(LQ%OER_am*C_O2 + LQ%K_ms)/(C_O2 + LQ%K_ms)
	OER_beta_d = dose*(LQ%OER_bm*C_O2 + LQ%K_ms)/(C_O2 + LQ%K_ms)
	expon = LQ%alpha_H*OER_alpha_d + LQ%beta_H*OER_alpha_d**2
	kill_prob = 1 - exp(-expon)
!	write(*,'(i6,3e12.4)') kcell,C_O2,expon,kill_prob
	R = par_uni(kpar)
	if (R < kill_prob) then
		cell_list(kcell)%radiation_tag = .true.
		Nradiation_tag = Nradiation_tag + 1
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Cells move to preferable nearby sites.
! For now this is turned off - need to formulate a sensible interpretation of "preferable"
!-----------------------------------------------------------------------------------------
subroutine CellMigration(ok)
logical :: ok
integer :: kcell, j, indx, site0(3), site(3), jmax
real(REAL_KIND) :: C0(MAX_CHEMO), C(MAX_CHEMO), v0, v, vmax, d0, d

call logger('CellMigration is not yet implemented')
ok = .false.
return

ok = .true.
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	site0 = cell_list(kcell)%site
	C0 = occupancy(site0(1),site0(2),site0(3))%C(:)
	v0 = SiteValue(C0)
	d0 = cdistance(site0)
	jmax = 0
	vmax = -1.0e10
	do j = 1,27
		if (j == 14) cycle
		site = site0 + jumpvec(:,j)
		indx = occupancy(site(1),site(2),site(3))%indx(1)
!		if (indx < -100) then	! necrotic site
		if (indx == 0) then	!	vacant site
			C = occupancy(site(1),site(2),site(3))%C(:)
			v = SiteValue(C)
			d = cdistance(site)
			if (d > d0 .and. v > v0) then
				vmax = v
				jmax = j
			endif
		endif
	enddo
	if (jmax > 0) then
		site = site0 + jumpvec(:,jmax)
		indx = occupancy(site(1),site(2),site(3))%indx(1)
		cell_list(kcell)%site = site
		occupancy(site(1),site(2),site(3))%indx(1) = kcell
		occupancy(site0(1),site0(2),site0(3))%indx(1) = indx
!		write(*,'(i2,2f8.4,3i4,f6.1,4x,3i4,f6.1)') jmax,v0,vmax,site0,d0,site,d
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! A measure of the attractiveness of a site with concentrations C(:)
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function SiteValue(C)
real(REAL_KIND) :: C(:)

SiteValue = C(OXYGEN)
end function

!-----------------------------------------------------------------------------------------
! Cells can be tagged to die, or finally die of anoxia, or they can be tagged for death 
! at division time if the drug is effective.
!-----------------------------------------------------------------------------------------
subroutine CellDeath(dt,ok)
real(REAL_KIND) :: dt
logical :: ok
integer :: kcell, ict, nlist0, site(3), i, kpar=0 
real(REAL_KIND) :: C_O2, kmet, Kd, dMdt, pdeath, tnow
logical :: use_SN30000

!call logger('CellDeath')
ok = .true.
if (chemo(DRUG_A)%used .and. DRUG_A == SN30000) then
    use_SN30000 = .true.
else
    use_SN30000 = .false.
endif
tnow = istep*DELTA_T	! seconds
nlist0 = nlist
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
!	C_O2 = cell_list(kcell)%conc(OXYGEN)
	call getO2conc(kcell,C_O2)
	if (cell_list(kcell)%anoxia_tag) then
!		write(logmsg,*) 'anoxia_tag: ',kcell,cell_list(kcell)%state,tnow,cell_list(kcell)%t_anoxia_die
!		call logger(logmsg)
		if (tnow >= cell_list(kcell)%t_anoxia_die) then
!			call logger('cell dies')
			call CellDies(kcell)
			Nanoxia_dead = Nanoxia_dead + 1
			cycle
		endif
	else
		if (C_O2 < ANOXIA_THRESHOLD) then
			cell_list(kcell)%t_hypoxic = cell_list(kcell)%t_hypoxic + dt
			if (cell_list(kcell)%t_hypoxic > t_anoxic_limit) then
				cell_list(kcell)%anoxia_tag = .true.						! tagged to die later
				cell_list(kcell)%t_anoxia_die = tnow + anoxia_death_delay	! time that the cell will die
				Nanoxia_tag = Nanoxia_tag + 1
			endif
		endif
	endif
	if (use_SN30000) then
		ict = cell_list(kcell)%celltype
		Kd = SN30K%Kd(ict)
	    kmet = (SN30K%C1(ict) + SN30K%C2(ict)*SN30K%KO2(ict)/(SN30K%KO2(ict) + C_O2))*SN30K%Kmet0(ict)
	    dMdt = kmet*cell_list(kcell)%conc(DRUG_A)
	    if (SN30K%kill_model(ict) == 1) then
		    pdeath = Kd*dMdt*dt
	    elseif (SN30K%kill_model(ict) == 2) then
		    pdeath = Kd*dMdt*cell_list(kcell)%conc(DRUG_A)*dt
	    elseif (SN30K%kill_model(ict) == 3) then
		    pdeath = Kd*dMdt**2*dt
		endif
!	    write(*,'(4f10.5)') kmet,cell_list(kcell)%conc(DRUG_A),dMdt,pdeath
	    if (par_uni(kpar) < pdeath) then
            cell_list(kcell)%drug_tag = .true.
            Ndrug_tag = Ndrug_tag + 1
		endif
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_CellDies
integer :: kcell, i, kpar=0

do i = 1,10
	kcell = random_int(1,nlist,kpar)
	call CellDies(kcell)
enddo
stop
end subroutine

!-----------------------------------------------------------------------------------------
! If the dying cell site is less than a specified fraction f_migrate of the blob radius,
! the site migrates towards the blob centre.
! %indx -> 0
! If the site is on the boundary, it is removed from the boundary list, and %indx -> OUTSIDE_TAG
! The cell contents should be released into the site.
!-----------------------------------------------------------------------------------------
subroutine CellDies(kcell)
integer :: kcell
integer :: site(3)
real(REAL_KIND) :: V

cell_list(kcell)%state = DEAD
cell_list(kcell)%exists = .false.
Ncells = Ncells - 1
site = cell_list(kcell)%site
occupancy(site(1),site(2),site(3))%indx(1) = 0
if (associated(occupancy(site(1),site(2),site(3))%bdry)) then
	call bdrylist_delete(site,bdrylist)
    nullify(occupancy(site(1),site(2),site(3))%bdry)
	occupancy(site(1),site(2),site(3))%indx(1) = OUTSIDE_TAG
	Nsites = Nsites - 1
	bdry_changed = .true.
	call OutsideNeighbours(site)
	call AddToMedium(kcell,site)
else
	V = cell_list(kcell)%volume*Vcell
	occupancy(site(1),site(2),site(3))%C = ((Vsite - V)*occupancy(site(1),site(2),site(3))%C + V*cell_list(kcell)%conc)/Vsite
endif
!call NecroticMigration(site)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine AddToMedium(kcell,site)
integer :: kcell, site(3)
integer :: ic
real(REAL_KIND) :: V, Cex(MAX_CHEMO), Cin(MAX_CHEMO)

Cex = occupancy(site(1),site(2),site(3))%C
Cin = cell_list(kcell)%conc
V = cell_list(kcell)%volume*Vcell
do ic = 1,MAX_CHEMO
	if (.not.chemo(ic)%used) cycle
	chemo(ic)%medium_M = chemo(ic)%medium_M + V*Cin(ic) + (Vsite - V)*Cex(ic)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine RemoveFromMedium
integer :: ic

do ic = 1,MAX_CHEMO
	if (.not.chemo(ic)%used) cycle
	chemo(ic)%medium_M = chemo(ic)%medium_M - Vsite*chemo(ic)%medium_Cbnd
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! When a bdry cell dies, need to check the neighbours to see if a site is now made outside.
! To be made outside a site must have %indx(1) = 0 and at least one Neumann neighbour outside.
!-----------------------------------------------------------------------------------------
subroutine OutsideNeighbours(site0)
integer :: site0(3)
integer :: k, j, x, y, z, xx, yy, zz
logical :: isout, done

done = .false.
do while(.not.done)
	done = .true.
	do k = 1,27
		if (k == 14) cycle
		x = site0(1) + jumpvec(1,k)
		y = site0(2) + jumpvec(2,k)
		z = site0(3) + jumpvec(3,k)
		if (outside_xyz(x,y,z)) cycle
		if (occupancy(x,y,z)%indx(1) == 0) then
			isout = .false.
			do j = 1,6
				xx = x + neumann(1,j)
				yy = y + neumann(2,j)
				zz = z + neumann(3,j)
				if (outside_xyz(xx,yy,zz)) cycle
				if (occupancy(xx,yy,zz)%indx(1) == OUTSIDE_TAG) then
					isout = .true.
					exit
				endif
			enddo
			if (isout) then
				done = .false.
				occupancy(x,y,z)%indx(1) = OUTSIDE_TAG
			endif
		endif
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! A necrotic site migrates towards the blob centre, stopping when another necrotic 
! site is reached
!-----------------------------------------------------------------------------------------
subroutine NecroticMigration(site0)
integer :: site0(3)
integer :: site1(3), site2(3), site(3), j, jmin, kcell, tmp_indx
real(REAL_KIND) :: d1, d2, dmin

!write(logmsg,*) 'NecroticMigration: site0: ',site0
!call logger(logmsg)
site1 = site0
do
	d1 = cdistance(site1)
	dmin = 1.0e10
	jmin = 0
	do j = 1,27
		if (j == 14) cycle
		site = site1 + jumpvec(:,j)
!		if (occupancy(site(1),site(2),site(3))%indx(1) < 0) cycle	! do not swap with another necrotic site
		if (occupancy(site(1),site(2),site(3))%indx(1) == 0) cycle	! do not swap with another necrotic site
		d2 = cdistance(site)
		if (d2 < dmin) then
			dmin = d2
			jmin = j
		endif
	enddo
!	write(*,*) 'd1, dmin, jmin: ',d1,dmin, jmin
	if (dmin >= d1) then
		site2 = site1
		exit
	endif
	if (jmin <= 0) then
		write(*,*) 'Error: jmin: ',jmin
		stop
	endif
!	write(*,*) site1,jmin,jumpvec(:,jmin)
	site2 = site1 + jumpvec(:,jmin)
!	write(*,*) site2
	! Now swap site1 and site2
	kcell = occupancy(site2(1),site2(2),site2(3))%indx(1)
	tmp_indx = occupancy(site1(1),site1(2),site1(3))%indx(1)
	occupancy(site1(1),site1(2),site1(3))%indx(1) = kcell
	occupancy(site2(1),site2(2),site2(3))%indx(1) = tmp_indx
	cell_list(kcell)%site = site1
	site1 = site2
enddo
!write(*,*) 'NecroticMigration: site2: ',site2
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_CellDivision(ok)
logical :: ok
integer :: kcell, kpar=0
kcell = random_int(1,nlist,kpar)
call CellDivider(kcell,ok)
end subroutine

!-----------------------------------------------------------------------------------------
! Cell growth, death and division are handled here.  Division occurs when cell volume 
! exceeds the divide volume. 
! As the cell grows we need to adjust both Cin and Cex to manintain mass conservation.
!-----------------------------------------------------------------------------------------
subroutine CellGrowth(dt,ok)
real(REAL_KIND) :: dt
logical :: ok
integer :: kcell, nlist0, site(3)
integer :: divide_list(10000), ndivide, i
real(REAL_KIND) :: tnow, C_O2, metab, dVdt, vol0, r_mean, c_rate
real(REAL_KIND) :: Vin_0, Vex_0, dV
real(REAL_KIND) :: Cin_0(MAX_CHEMO), Cex_0(MAX_CHEMO)
character*(20) :: msg
logical :: C_option = 1

!call logger('CellGrowth')
ok = .true.
nlist0 = nlist
tnow = istep*DELTA_T
c_rate = log(2.0)/divide_time_mean		! Note: to randomise divide time need to use random number, not mean!
r_mean = Vdivide0/(2*divide_time_mean)
ndivide = 0
do kcell = 1,nlist0
	if (cell_list(kcell)%state == DEAD) cycle
	C_O2 = cell_list(kcell)%conc(OXYGEN)
!    if (C_O2 > ODEdiff%C1_soft) then
!	    metab = (C_O2-ODEdiff%deltaC_soft)/(chemo(OXYGEN)%MM_C0 + C_O2 - ODEdiff%deltaC_soft)
!    elseif (C_O2 > 0) then
!	    metab = ODEdiff%k_soft*C_O2*C_O2
!	else
!		metab = 0
!	endif
	metab = O2_metab(C_O2)
	if (use_V_dependence) then
		dVdt = c_rate*cell_list(kcell)%volume*metab
	else
		dVdt = r_mean*metab
	endif
!	if (istep > 1 .and. dVdt == 0) then
!		write(nflog,'(a,2i6,5e12.3)') 'dVdt: ',istep,kcell,r_mean,c_rate,C_O2,metab,dVdt
!	endif
	site = cell_list(kcell)%site
	Cin_0 = cell_list(kcell)%conc
	Cex_0 = occupancy(site(1),site(2),site(3))%C
	cell_list(kcell)%dVdt = dVdt
	Vin_0 = cell_list(kcell)%volume*Vcell	! cm^3
	Vex_0 = Vsite - Vin_0					! cm^3
	dV = dVdt*dt*Vcell						! cm^3
	cell_list(kcell)%volume = (Vin_0 + dV)/Vcell
	if (C_option == 1) then
		! Calculation based on transfer of an extracellular volume dV with constituents, i.e. holding extracellular concentrations constant
		cell_list(kcell)%conc = (Vin_0*Cin_0 + dV*Cex_0)/(Vin_0 + dV)
		occupancy(site(1),site(2),site(3))%C = (Vex_0*Cex_0 - dV*Cex_0)/(Vex_0 - dV)
	elseif (C_option == 2) then
		! Calculation based on change in volumes without mass transfer of constituents
		cell_list(kcell)%conc = Vin_0*Cin_0/(Vin_0 + dV)
		occupancy(site(1),site(2),site(3))%C = Vex_0*Cex_0/(Vex_0 - dV)
	endif
	if (cell_list(kcell)%volume > cell_list(kcell)%divide_volume) then
		if (cell_list(kcell)%radiation_tag) then
			call CellDies(kcell)
			Nradiation_dead = Nradiation_dead + 1
!			write(*,*) 'Cell died: ',kcell,
			cycle
		endif
		if (cell_list(kcell)%drug_tag) then
			call CellDies(kcell)
			Ndrug_dead = Ndrug_dead + 1
			cycle
		endif
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

!-----------------------------------------------------------------------------------------
! Note: updating of concentrations is now done in extendODEdiff()
! The dividing cell, kcell0, is at site0.
! A neighbour site site01 is chosen randomly.  This becomes the site for the daughter cell.
!-----------------------------------------------------------------------------------------
subroutine CellDivider(kcell0, ok)
integer :: kcell0
logical :: ok
integer :: kpar=0
integer :: j, k, kcell1, site0(3), site1(3), site2(3), site01(3), site(3), ichemo, nfree, bestsite(3)
integer :: npath, path(3,200)
real(REAL_KIND) :: tnow, R, v, vmax, V0, Cex(MAX_CHEMO), M0(MAX_CHEMO), M1(MAX_CHEMO), alpha(MAX_CHEMO)
logical :: freesite(27,3)
type (boundary_type), pointer :: bdry

!write(logmsg,*) 'CellDivider: ',kcell0
!call logger(logmsg)
ok = .true.
tnow = istep*DELTA_T
cell_list(kcell0)%t_divide_last = tnow
V0 = cell_list(kcell0)%volume
cell_list(kcell0)%volume = V0/2
R = par_uni(kpar)
cell_list(kcell0)%divide_volume = Vdivide0 + dVdivide*(2*R-1)
cell_list(kcell0)%M = cell_list(kcell0)%M/2
!write(logmsg,'(a,f6.1)') 'Divide time: ',tnow/3600
!call logger(logmsg)

site0 = cell_list(kcell0)%site
if (divide_option == DIVIDE_USE_CLEAR_SITE .or. &			! look for the best nearby clear site, if it exists use it
	divide_option == DIVIDE_USE_CLEAR_SITE_RANDOM) then		! make random choice from nearby clear sites
	vmax = -1.0e10
	nfree = 0
	do j = 1,27
		if (j == 14) cycle
		site01 = site0 + jumpvec(:,j)
!		if (occupancy(site01(1),site01(2),site01(3))%indx(1) < -100) then
		if (occupancy(site01(1),site01(2),site01(3))%indx(1) == 0) then
			nfree = nfree + 1
			freesite(nfree,:) = site01
			v = SiteValue(occupancy(site01(1),site01(2),site01(3))%C(:))
			if (v > vmax) then
				vmax = v
				bestsite = site01
			endif
		endif
	enddo
	if (nfree > 0) then
		if (DIVIDE_USE_CLEAR_SITE) then	! use this site for the progeny cell
			site01 = bestsite
		else	! random choice
			j = random_int(1,nfree,0)
			site01 = freesite(j,:)
		endif
		call AddCell(kcell0,kcell1,site01,ok)
		if (.not.ok) then
			call logger('Error: AddCell: vacant site')
		endif
		Nreuse = Nreuse + 1
		return
	endif
endif

!ok = .false.
!k = 0
!do while (.not.ok)
!    k = k + 1
!!	j = random_int(1,26,kpar)       ! This generates disconnected cells at the boundary
!!	if (j >= 14) j = j+1
!!	site01 = site0 + jumpvec(:,j)
!	j = random_int(1,6,kpar)
!	site01 = site0 + neumann(:,j)
!	if (k > 50 .or. occupancy(site01(1),site01(2),site01(3))%indx(1) >= -100) ok = .true.	! site01 is not necrotic
!enddo
j = random_int(1,6,kpar)
site01 = site0 + neumann(:,j)
!if (dbug) write(*,*) 'CellDivider: ',kcell0,site0,occupancy(site0(1),site0(2),site0(3))%indx
if (occupancy(site01(1),site01(2),site01(3))%indx(1) == OUTSIDE_TAG) then	! site01 is outside, use it directly
	npath = 0
	site1 = site0
elseif (bdrylist_present(site01,bdrylist)) then	! site01 is on the boundary
	npath = 1
	site1 = site01
	path(:,1) = site01
else
	call ChooseBdrysite(site01,site1)
	if (occupancy(site1(1),site1(2),site1(3))%indx(1) == 0) then
		write(*,*) 'after choose_bdrysite: site1 is VACANT: ',site1
		stop
	endif
	if (dbug) write(*,'(a,i6,9i4)') 'b4 ',kcell0,site0,site01,site1
	call SelectPath(site0,site01,site1,path,npath)
	! path(:,:) goes from site01 to site1, which is a bdry site or adjacent to a vacant site
!	if (.not.isbdry(site1)) then
!	    call logger('should be bdry or vacant, is not')
!	    stop
!	endif
endif
!if (dbug) write(*,*) 'path: ',npath
!do k = 1,npath
!	if (dbug) write(*,'(i3,2x,3i4)') path(:,k)
!enddo

if (npath > 0) then
	! Need to choose an outside or vacant site near site1
	call getOutsideSite(site1,site2)
	npath = npath+1
	path(:,npath) = site2
	! path(:,:) now goes from site01 to site2, which is an outside site next to site1
	if (occupancy(site2(1),site2(2),site2(3))%indx(1) == OUTSIDE_TAG) then
		Nsites = Nsites + 1
		call RemoveFromMedium
	endif
!	write(*,'(a,3i4,i6)') 'outside site: ',site2,occupancy(site2(1),site2(2),site2(3))%indx(1)
else
	! site01 is the destination site for the new cell
	site2 = site01
	Nsites = Nsites + 1
endif

call GetPathMass(site0,site01,path,npath,M0)
if (npath > 0) then
	call PushPath(path,npath)
	if (dbug) write(*,*) 'did push_path'
endif
call AddCell(kcell0,kcell1,site01,ok)
if (.not.ok) then
	call logger('Error: AddCell: pushed site')
	return
endif
call GetPathMass(site0,site01,path,npath,M1)
alpha = M0/M1	! scaling for concentrations on the path
call ScalePathConcentrations(site0,site01,path,npath,alpha)

!if (npath > 0) then
!	call FixPathConcentrations1(path,npath)
!endif
!! Now adjust extracellular concentrations for the parent cell site to ensure mass conservation
!Cex = occupancy(site0(1),site0(2),site0(3))%C(:)
!occupancy(site0(1),site0(2),site0(3))%C(:) = (Cex*(Vsite - V0) + occupancy(site01(1),site01(2),site01(3))%C(:)*V0/2)/(Vsite - V0/2)

call SetRadius(Nsites)
!call extendODEdiff(site2)
! Now need to fix the bdrylist.  
! site1 was on the boundary, but may no longer be.
! site2 may be now on the boundary
! First add site2
if (isbdry(site2)) then   ! add it to the bdrylist
	if (dbug) write(*,*) 'add site2 to bdrylist: ',site2
    allocate(bdry)
    bdry%site = site2
!    bdry%chemo_influx = .false.
    nullify(bdry%next)
    call bdrylist_insert(bdry,bdrylist)
    occupancy(site2(1),site2(2),site2(3))%bdry => bdry
    call SetBdryConcs(site2)
else
!    write(logmsg,'(a,3i4,i6)') 'Added site is not bdry: ',site2,occupancy(site2(1),site2(2),site2(3))%indx(1)
!	call logger(logmsg)
    call SetBdryConcs(site2)
!    stop
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
			allocate(bdry)
			bdry%site = site
		!    bdry%chemo_influx = .false.
			nullify(bdry%next)
			call bdrylist_insert(bdry,bdrylist)
			occupancy(site(1),site(2),site(3))%bdry => bdry
!		    call SetBdryConcs(site)
		endif
	else
		if (dbug) write(*,*) 'not isbdry'
		if (bdrylist_present(site,bdrylist)) then	! remove it
			if (dbug) write(*,*) 'present, remove it'
			call bdrylist_delete(site,bdrylist)
			nullify(occupancy(site(1),site(2),site(3))%bdry)
!			call ResetConcs(site)
		endif
	endif
enddo
!call AdjustMedium
if (dbug) write(*,*) 'done!'
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine GetPathMass(site0,site01,path,npath,mass)
integer :: site0(3),site01(3),path(3,200),npath
real(REAL_KIND) :: mass(:)
integer :: k, site(3), kcell, ic
real(REAL_KIND) :: V, C

do ic = 1,MAX_CHEMO
	mass(ic) = 0
	if (.not.chemo(ic)%used) cycle
	C = occupancy(site0(1),site0(2),site0(3))%C(ic)
	kcell = occupancy(site0(1),site0(2),site0(3))%indx(1)
	if (kcell > 0) then
		V = cell_list(kcell)%volume*Vcell
	else
		V = 0
	endif
	mass(ic) = mass(ic) + C*(Vsite - V)
	if (npath == 0) then
		C = occupancy(site01(1),site01(2),site01(3))%C(ic)
		kcell = occupancy(site01(1),site01(2),site01(3))%indx(1)
		if (kcell > 0) then
			V = cell_list(kcell)%volume*Vcell
		else
			V = 0
		endif
		mass(ic) = mass(ic) + C*(Vsite - V)
	else
		do k = 1, npath
			site = path(:,k)
			C = occupancy(site01(1),site01(2),site01(3))%C(ic)
			kcell = occupancy(site(1),site(2),site(3))%indx(1)
			if (kcell > 0) then
				V = cell_list(kcell)%volume*Vcell
			else
				V = 0
			endif
			mass(ic) = mass(ic) + C*(Vsite - V)
		enddo
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ScalePathConcentrations(site0,site01,path,npath,alpha)
integer :: site0(3),site01(3),path(3,200),npath
real(REAL_KIND) :: alpha(:)
integer :: k, site(3), kcell, ic

do ic = 1,MAX_CHEMO
	if (.not.chemo(ic)%used) cycle
	occupancy(site0(1),site0(2),site0(3))%C(ic) = alpha(ic)*occupancy(site0(1),site0(2),site0(3))%C(ic)
	if (npath == 0) then
		occupancy(site01(1),site01(2),site01(3))%C(ic) = alpha(ic)*occupancy(site01(1),site01(2),site01(3))%C(ic)
	else
		do k = 1, npath
			site = path(:,k)
			occupancy(site01(1),site01(2),site01(3))%C(ic) = alpha(ic)*occupancy(site01(1),site01(2),site01(3))%C(ic)
		enddo
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! The extracellular concentrations along the path are adjusted sequentially, starting from
! the last-but-one site and working backwards to the first, site01
!-----------------------------------------------------------------------------------------
subroutine FixPathConcentrations1(path,npath)
integer :: path(3,200),npath
integer :: k, site1(3), site2(3), kcell, ic
real(REAL_KIND) :: V, Cex(MAX_CHEMO)

do k = npath-1,1,-1
	site1 = path(:,k)
	kcell = occupancy(site1(1),site1(2),site1(3))%indx(1)
	V = cell_list(kcell)%volume*Vcell
	Cex = occupancy(site1(1),site1(2),site1(3))%C
	site2 = path(:,k+1)
	occupancy(site1(1),site1(2),site1(3))%C = ((Vsite - V)*Cex + V*occupancy(site2(1),site2(2),site2(3))%C)/Vsite
	do ic = 1,MAX_CHEMO
		if (.not.chemo(ic)%used) cycle
		if (occupancy(site1(1),site1(2),site1(3))%C(ic) < 0) then
			occupancy(site1(1),site1(2),site1(3))%C(ic) = max(occupancy(site2(1),site2(2),site2(3))%C(ic), Cex(ic))
		endif
	enddo
!	write(*,'(i4,4e12.4)') k,Cex(1:2),Vsite,V
!	write(*,'(4x,2e12.4)') occupancy(site2(1),site2(2),site2(3))%C(1:2)
!	write(*,'(4x,2e12.4)') occupancy(site1(1),site1(2),site1(3))%C(1:2)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine AdjustMedium

end subroutine

!-----------------------------------------------------------------------------------------
! Reduce medium volume and concentrations to account for the growth of the spheroid
! by one site.
! THIS WILL BE SUPERCEDED
!-----------------------------------------------------------------------------------------
subroutine AdjustMedium1
real(REAL_KIND) :: total(MAX_CHEMO)

total(DRUG_A:MAX_CHEMO) = chemo(DRUG_A:MAX_CHEMO)%bdry_conc*medium_volume
medium_volume = medium_volume - Vsite
chemo(DRUG_A:MAX_CHEMO)%bdry_conc = total(DRUG_A:MAX_CHEMO)/medium_volume
end subroutine

!-----------------------------------------------------------------------------------------
! Need to choose a site on the boundary in some sense near site0, and also to preserve
! the spherical shape
!-----------------------------------------------------------------------------------------
subroutine ChooseBdrysite(site0,site1)
integer :: site0(3), site1(3)
integer :: site(3), sitemin(3)
real(REAL_KIND) :: vc(3), v(3), r, rfrac, d, alpha_max, cosa_min, dmin, cosa
logical :: hit
type (boundary_type), pointer :: bdry

if (dbug) write(*,*) 'ChooseBdrysite: ',site0
vc = site0 - Centre
r = norm(vc)
vc = vc/r
rfrac = r/Radius
alpha_max = getAlphaMax(rfrac)
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
real(REAL_KIND) function getAlphaMax(r)
real(REAL_KIND) :: r
real(REAL_KIND) :: alf
real(REAL_KIND) :: r1 = 0.0, r2 = 0.8, alpha1 = PI/2, alpha2 = PI/6

if (r < r1) then
	getAlphaMax = alpha1
elseif (r > r2) then
	getAlphaMax = alpha2
else
	alf = (r-r1)/(r2-r1)
	getAlphaMax = (1-alf)*alpha1 + alf*alpha2
endif
end function


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_get_path
integer :: site0(3), site1(3), site2(3), path(3,200), npath, kpar=0
integer :: x, y, z, it, k

site0 = 0
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
	call SelectPath(site0,site1,site2,path,npath)
	write(*,*) 'path: ',npath
	do k = 1,npath
		write(*,'(i3,2x,3i4)') path(:,k)
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Find a path of sites leading from site1 to site2
! The first site in the list is site1, the last is site2
! Previously site2 was always a bdry site.  Now we stop the path when a vacant or outside site 
! is encountered.
! How do we avoid choosing a path through the site of the dividing cell? site0
!-----------------------------------------------------------------------------------------
subroutine SelectPath(site0,site1,site2,path,npath)
integer :: site0(3),site1(3), site2(3), path(3,200), npath
integer :: v(3), jump(3), site(3), k, j, jmin
real(REAL_KIND) :: r, d2, d2min
logical :: hit

if (occupancy(site2(1),site2(2),site2(3))%indx(1) == OUTSIDE_TAG) then
    call logger('site2 is OUTSIDE')
    stop
endif
if (occupancy(site2(1),site2(2),site2(3))%indx(1) == 0) then
    call logger('site2 is VACANT')
    stop
endif

k = 1
site = site1
path(:,k) = site
hit = .false.
do 
	d2min = 1.0e10
	jmin = 0
	do j = 1,27
		if (j == 14) cycle
		jump = jumpvec(:,j)
		v = site + jump
		if (v(1)==site0(1) .and. v(2)==site0(2) .and. v(3)==site0(3)) cycle
!		if (occupancy(v(1),v(2),v(3))%indx(1) == OUTSIDE_TAG) then
		if (occupancy(v(1),v(2),v(3))%indx(1) <= 0) then	! outside or vacant - we'll use this!
			site2 = site
			hit = .true.
			exit
		endif
!		if (occupancy(v(1),v(2),v(3))%indx(1) < -100) cycle		! necrotic
		v = site2 - v
		d2 = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
		if (dbug) write(*,'(8i6,2f6.1)') k,j,site+jump,v,d2,d2min
		if (d2 < d2min) then
			d2min = d2
			jmin = j
		endif
	enddo
	if (hit) exit
	if (jmin == 0) then
	    call logger('get path: stuck')
	    stop
	endif
	site = site + jumpvec(:,jmin)
	k = k+1
	if (dbug) write(*,'(11i4,f8.1)') k,site1,site2,site,jmin,d2min
	if (k==20) stop
	path(:,k) = site
	if (site(1) == site2(1) .and. site(2) == site2(2) .and. site(3) == site2(3)) exit
enddo
npath = k

end subroutine

!-----------------------------------------------------------------------------------------
! If there are any necrotic sites in the path, they are first shifted closer to the centre,
! then the path is adjusted.
!-----------------------------------------------------------------------------------------
subroutine ClearPath(path,npath)
integer :: npath, path(3,*)
logical :: clear
integer :: k, kcell, kvacant, site(3)

do
	clear = .true.
	do k = 1,npath
		site = path(:,k)
		kcell = occupancy(site(1),site(2),site(3))%indx(1)
		if (kcell < 0) then
			clear = .false.
			kvacant = k
			exit
		endif
	enddo
	if (clear) return
	
enddo			
end subroutine

!-----------------------------------------------------------------------------------------
! Need to modify the ODEdiff variables, at least %ivar(?)
!-----------------------------------------------------------------------------------------
subroutine PushPath(path,npath)
integer :: path(3,200),npath
integer :: k, site1(3), site2(3), kcell

do k = npath-1,1,-1
	site1 = path(:,k)
	kcell = occupancy(site1(1),site1(2),site1(3))%indx(1)
	if (dbug) write(*,*) k,' site1: ',site1,kcell
	site2 = path(:,k+1)
	if (dbug) write(*,*) 'site2: ',site2
	if (kcell > 0) then
    	cell_list(kcell)%site = site2
    endif
	occupancy(site2(1),site2(2),site2(3))%indx(1) = kcell
enddo
occupancy(site1(1),site1(2),site1(3))%indx = 0
end subroutine

!-----------------------------------------------------------------------------------------
! The daughter cell kcell1 is given the same characteristics as kcell0 and placed at site1.
! Random variation is introduced into %divide_volume.
!-----------------------------------------------------------------------------------------
subroutine AddCell(kcell0,kcell1,site1,ok)
integer :: kcell0, kcell1, site1(3)
logical :: ok
integer :: kpar = 0
real(REAL_KIND) :: tnow, R

ok = .true.
tnow = istep*DELTA_T
!lastID = lastID + 1
nlist = nlist + 1
if (nlist > max_nlist) then
	call logger('Dimension of cell_list() has been exceeded: increase max_nlist and rebuild')
	ok = .false.
	return
endif
Ncells = Ncells + 1
kcell1 = nlist
cell_list(kcell1)%celltype = cell_list(kcell0)%celltype
cell_list(kcell1)%state = cell_list(kcell0)%state
cell_list(kcell1)%site = site1
!cell_list(kcell1)%ID = lastID
cell_list(kcell1)%ID = cell_list(kcell0)%ID
cell_list(kcell1)%radiation_tag = .false.
cell_list(kcell1)%drug_tag = .false.
cell_list(kcell1)%anoxia_tag = .false.
cell_list(kcell1)%exists = .true.
cell_list(kcell1)%active = .true.
cell_list(kcell1)%t_divide_last = tnow
!cell_list(kcell1)%t_divide_next = tnow + DivideTime()
cell_list(kcell1)%dVdt = cell_list(kcell0)%dVdt
cell_list(kcell1)%volume = cell_list(kcell0)%volume
R = par_uni(kpar)
cell_list(kcell1)%divide_volume = Vdivide0 + dVdivide*(2*R-1)
cell_list(kcell1)%t_hypoxic = 0
cell_list(kcell1)%conc = cell_list(kcell0)%conc
cell_list(kcell1)%M = cell_list(kcell0)%M
!cell_list(kcell1)%oxygen = cell_list(kcell0)%oxygen
!cell_list(kcell1)%drug_A = cell_list(kcell0)%drug_A
!cell_list(kcell1)%drug_B = cell_list(kcell0)%drug_B
occupancy(site1(1),site1(2),site1(3))%indx(1) = kcell1
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
! The location site0 is just outside the blob.
! The aim is to find a cell near site0 that is further from the centre, and move it here.
!-----------------------------------------------------------------------------------------
subroutine Adjust(site0)
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

end module
