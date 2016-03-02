! Cancer cell state development

module cellstate
use global
use boundary
use chemokine
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
!if (use_radiation .and. dose > 0) then
if (dose > 0) then
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

if (use_extracellular_O2 .and. istep > 1) then		! fix 30/04/2015
	iv = cell_list(kcell)%ivin
	if (iv < 1) then
		C_O2 = cell_list(kcell)%conc(OXYGEN)		! use this until %ivin is set in SetupODEDiff
!		write(logmsg,*) 'getO2conc: ',istep,kcell,site,iv
!		call logger(logmsg)
!		stop
!		tnow = istep*DELTA_T
!		C_O2 = BdryConc(OXYGEN,tnow)	! assume that this is a site at the boundary
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
integer :: kcell, site(3), iv, ityp, idrug, im, ichemo, kpar=0
real(REAL_KIND) :: C_O2, SER, p_death, p_recovery, R, kill_prob, tnow
real(REAL_KIND) :: Cs							! concentration of radiosensitising drug
real(REAL_KIND) :: SER_max0, SER_Km, SER_KO2	! SER parameters of the drug
real(REAL_KIND) :: SERmax						! max sensitisation at the drug concentration

ok = .true.
call logger('Irradiation')
tnow = istep*DELTA_T	! seconds
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	if (cell_list(kcell)%radiation_tag) cycle	! we do not tag twice (yet)
	ityp = cell_list(kcell)%celltype
	call getO2conc(kcell,C_O2)
	! Compute sensitisation SER
	SER = 1.0
	do idrug = 1,Ndrugs_used
		ichemo = 4 + 3*(idrug-1)
		if (.not.chemo(ichemo)%present) cycle
		do im = 0,2
			ichemo = 4 + 3*(idrug-1) + im
			if (drug(idrug)%sensitises(ityp,im)) then
				Cs = cell_list(kcell)%conc(ichemo)	! concentration of drug/metabolite in the cell
				SER_max0 = drug(idrug)%SER_max(ityp,im)
				SER_Km = drug(idrug)%SER_Km(ityp,im)
				SER_KO2 = drug(idrug)%SER_KO2(ityp,im)
				SERmax = (Cs*SER_max0 + SER_Km)/(Cs + SER_Km)
				SER = SER*(C_O2 + SER_KO2*SERmax)/(C_O2 + SER_KO2)
			endif
		enddo
	enddo
	call get_kill_probs(ityp,dose,C_O2,SER,p_recovery,p_death)
	kill_prob = 1 - p_recovery
	R = par_uni(kpar)
	if (R < kill_prob) then
		cell_list(kcell)%radiation_tag = .true.
		Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
		cell_list(kcell)%p_rad_death = p_death
		if (LQ(ityp)%growth_delay_N > 0) then
			cell_list(kcell)%growth_delay = .true.
			cell_list(kcell)%dt_delay = LQ(ityp)%growth_delay_factor*dose
			cell_list(kcell)%N_delayed_cycles_left = LQ(ityp)%growth_delay_N
		else
			cell_list(kcell)%growth_delay = .false.
		endif
	elseif (use_radiation_growth_delay_all .and. LQ(ityp)%growth_delay_N > 0) then
!		radiation_dosed = .true.
		cell_list(kcell)%growth_delay = .true.
		cell_list(kcell)%dt_delay = LQ(ityp)%growth_delay_factor*dose
		cell_list(kcell)%N_delayed_cycles_left = LQ(ityp)%growth_delay_N
	else
		cell_list(kcell)%growth_delay = .false.
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! A cell that receives a dose of radiation either recovers completely before reaching 
! mitosis or retains damage that has a probability of causing cell death during mitosis.
! A damaged cell that does not die at this point passes the damage on to the progeny cells.
! The probability of complete recovery = p_recovery = p_r
! The probability of death for a damaged cell at mitosis = p_death = p_d
! To ensure that the short-term death probability is consistent with the previous
! LQ formulation, we require p_d(1-p_r) = kill_prob as previously calculated.
! If p_d is determined (currently it is fixed), then 1-p_r = kill_prob/p_d,
! therefore p_r = 1 - kill_prob/p_d
!-----------------------------------------------------------------------------------------
subroutine get_kill_probs(ityp,dose,C_O2,SER,p_recovery,p_death)
integer :: ityp
real(REAL_KIND) :: dose, C_O2, SER, p_recovery, p_death
real(REAL_KIND) :: OER_alpha_d, OER_beta_d, expon, kill_prob_orig

OER_alpha_d = dose*(LQ(ityp)%OER_am*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)
OER_beta_d = dose*(LQ(ityp)%OER_bm*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)
!expon_err = LQ(ityp)%alpha_H*OER_alpha_d + LQ(ityp)%beta_H*OER_alpha_d**2	! 07/08/2015

OER_alpha_d = OER_alpha_d*SER
OER_beta_d = OER_beta_d*SER

expon = LQ(ityp)%alpha_H*OER_alpha_d + LQ(ityp)%beta_H*OER_beta_d**2
p_recovery = exp(-expon)	! = SF
p_death = LQ(ityp)%death_prob

!kill_prob_orig = 1 - exp(-expon)
!call get_pdeath(ityp,dose,C_O2,p_death)
!p_recovery = 1 - kill_prob_orig/p_death
!p_recovery = 1 - kill_prob_orig
!write(nflog,'(a,4e12.3)') 'get_kill_probs: OER_alpha_d,OER_beta_d,expon,expon_err: ',OER_alpha_d,OER_beta_d,expon,expon_err
!write(nflog,'(a,2e12.3)') 'kill_prob, kill_prob_err: ',kill_prob_orig,1 - exp(-expon_err)
!write(nflog,'(a,e12.3)') 'p_recovery: ',p_recovery
end subroutine

!-----------------------------------------------------------------------------------------
! This is the probability of death at time of division of cell that received a dose of 
! radiation and did not recover.
! In general one would expect this to depend of damage, i.e. on dose and C_O2, but
! for now it is a constant for a cell type
!-----------------------------------------------------------------------------------------
subroutine get_pdeath(ityp,dose,C_O2,p_death)
integer :: ityp
real(REAL_KIND) :: dose, C_O2, p_death

p_death = LQ(ityp)%death_prob
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
integer :: kcell, nlist0, site(3), i, ichemo, idrug, im, ityp, killmodel, kpar=0 
real(REAL_KIND) :: C_O2, Cdrug, n_O2, kmet, Kd, dMdt, kill_prob, dkill_prob, death_prob, tnow
type(drug_type), pointer :: dp

!call logger('CellDeath')
ok = .true.
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
			cycle
		endif
	else
		if (C_O2 < ANOXIA_THRESHOLD) then
			cell_list(kcell)%t_hypoxic = cell_list(kcell)%t_hypoxic + dt
			if (cell_list(kcell)%t_hypoxic > t_anoxic_limit) then
				cell_list(kcell)%t_anoxia_die = tnow + anoxia_death_delay	! time that the cell will die
				if (.not.cell_list(kcell)%anoxia_tag) then	! don't retag
					Nanoxia_tag(ityp) = Nanoxia_tag(ityp) + 1
				endif
				cell_list(kcell)%anoxia_tag = .true.						! tagged to die later
			endif
		else
			cell_list(kcell)%t_hypoxic = 0
		endif
	endif
	
	do idrug = 1,ndrugs_used	
		dp => drug(idrug)
		ichemo = TRACER + 1 + 3*(idrug-1)	
		kill_prob = 0
		death_prob = 0
		do im = 0,2
			if (.not.dp%kills(ityp,im)) cycle
			killmodel = dp%kill_model(ityp,im)		! could use %drugclass to separate kill modes
			Cdrug = cell_list(kcell)%conc(ichemo + im)
			Kd = dp%Kd(ityp,im)
			n_O2 = dp%n_O2(ityp,im)
			kmet = (1 - dp%C2(ityp,im) + dp%C2(ityp,im)*dp%KO2(ityp,im)**n_O2/(dp%KO2(ityp,im)**n_O2 + C_O2**n_O2))*dp%Kmet0(ityp,im)
			dMdt = kmet*Cdrug
			call getDrugKillProb(killmodel,Kd,dMdt,Cdrug,dt,dkill_prob)
			kill_prob = kill_prob + dkill_prob
			death_prob = max(death_prob,dp%death_prob(ityp,im))
		enddo
	    if (.not.cell_list(kcell)%drug_tag(idrug) .and. par_uni(kpar) < kill_prob) then		! don't tag more than once
			cell_list(kcell)%p_drug_death(idrug) = death_prob
			cell_list(kcell)%drug_tag(idrug) = .true.
            Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) + 1
		endif
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getDrugKillProb(kill_model,Kd,dMdt,Cdrug,dt,dkill_prob)
integer :: kill_model
real(REAL_KIND) :: Kd, dMdt, Cdrug, dt, dkill_prob
real(REAL_KIND) :: SF, SF1, dtstep, kill_prob, c
integer :: Nsteps, istep

!Nsteps = dt + 0.5
!dtstep = 1.0
!SF1 = 1.0
!do istep = 1,Nsteps
!	if (kill_model == 1) then
!		kill_prob = Kd*dMdt*dtstep
!	elseif (kill_model == 2) then
!		kill_prob = Kd*dMdt*Cdrug*dtstep
!	elseif (kill_model == 3) then
!		kill_prob = Kd*dMdt**2*dtstep
!	elseif (kill_model == 4) then
!		kill_prob = Kd*Cdrug*dtstep
!	elseif (kill_model == 5) then
!		kill_prob = Kd*(Cdrug**2)*dtstep
!	endif
!    kill_prob = min(kill_prob,1.0)
!    SF1 = SF1*(1 - kill_prob)
!enddo

if (kill_model == 1) then
	c = Kd*dMdt
elseif (kill_model == 2) then
	c = Kd*dMdt*Cdrug
elseif (kill_model == 3) then
	c = Kd*dMdt**2
elseif (kill_model == 4) then
	c = Kd*Cdrug
elseif (kill_model == 5) then
	c = Kd*Cdrug**2
endif
SF = exp(-c*dt)
!write(nflog,'(a,2e12.3)') 'SF1, SF2: ',SF1, SF2
dkill_prob = 1 - SF
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
integer :: site(3), ityp, idrug
real(REAL_KIND) :: V

cell_list(kcell)%state = DEAD
cell_list(kcell)%exists = .false.
ityp = cell_list(kcell)%celltype
Ncells = Ncells - 1
Ncells_type(ityp) = Ncells_type(ityp) - 1
if (cell_list(kcell)%anoxia_tag) then
	Nanoxia_tag(ityp) = Nanoxia_tag(ityp) - 1
endif
do idrug = 1,ndrugs_used
	if (cell_list(kcell)%drug_tag(idrug)) then
		Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) - 1
	endif
enddo
if (cell_list(kcell)%radiation_tag) then
	Nradiation_tag(ityp) = Nradiation_tag(ityp) - 1
endif

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
	V = cell_list(kcell)%volume*Vcell_cm3
	occupancy(site(1),site(2),site(3))%C = ((Vsite_cm3 - V)*occupancy(site(1),site(2),site(3))%C + V*cell_list(kcell)%conc)/Vsite_cm3
endif
!call NecroticMigration(site)
ngaps = ngaps + 1
if (ngaps > max_ngaps) then
    write(logmsg,'(a,i6,i6)') 'CellDies: ngaps > max_ngaps: ',ngaps,max_ngaps
    call logger(logmsg)
    stop
endif
gaplist(ngaps) = kcell
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine AddToMedium(kcell,site)
integer :: kcell, site(3)
integer :: ichemo
real(REAL_KIND) :: V, Cex(MAX_CHEMO), Cin(MAX_CHEMO)

Cex = occupancy(site(1),site(2),site(3))%C
Cin = cell_list(kcell)%conc
V = cell_list(kcell)%volume*Vcell_cm3
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	chemo(ichemo)%medium_M = chemo(ichemo)%medium_M + V*Cin(ichemo) + (Vsite_cm3 - V)*Cex(ichemo)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine RemoveFromMedium
integer :: ichemo

do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	if (chemo(ichemo)%constant) cycle
	chemo(ichemo)%medium_M = chemo(ichemo)%medium_M - Vsite_cm3*chemo(ichemo)%medium_Cbnd
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
	if (dmin >= d1) then
		site2 = site1
		exit
	endif
	if (jmin <= 0) then
		write(nflog,*) 'Error: jmin: ',jmin
		stop
	endif
	site2 = site1 + jumpvec(:,jmin)
	! Now swap site1 and site2
	kcell = occupancy(site2(1),site2(2),site2(3))%indx(1)
	tmp_indx = occupancy(site1(1),site1(2),site1(3))%indx(1)
	occupancy(site1(1),site1(2),site1(3))%indx(1) = kcell
	occupancy(site2(1),site2(2),site2(3))%indx(1) = tmp_indx
	cell_list(kcell)%site = site1
	site1 = site2
enddo
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
! As the cell grows we need to adjust both Cin and Cex to maintain mass conservation.
! GROWTH DELAY
! When a cell has received a dose of radiation (or possibly drug - not yet considered)
! the cycle time is increased by an amount that depends on the dose.  The delay may be
! transmitted to progeny cells.
! 
!-----------------------------------------------------------------------------------------
subroutine CellGrowth(dt,ok)
real(REAL_KIND) :: dt
logical :: ok
integer :: kcell, nlist0, site(3), ityp, idrug, kpar=0
integer :: divide_list(10000), ndivide, i
real(REAL_KIND) :: tnow, C_O2, C_glucose, metab, metab_O2, metab_glucose, dVdt, vol0, R		!r_mean(2), c_rate(2)
real(REAL_KIND) :: r_mean(MAX_CELLTYPES), c_rate(MAX_CELLTYPES)
real(REAL_KIND) :: Vin_0, Vex_0, dV
real(REAL_KIND) :: Cin_0(MAX_CHEMO), Cex_0(MAX_CHEMO)
character*(20) :: msg
logical :: drugkilled, glucose_growth, first_cycle
integer :: C_option = 1

!call logger('CellGrowth')
ok = .true.
nlist0 = nlist
tnow = istep*DELTA_T
c_rate(1:2) = log(2.0)/divide_time_mean(1:2)		! Note: to randomise divide time need to use random number, not mean!
r_mean(1:2) = Vdivide0/(2*divide_time_mean(1:2))
glucose_growth = chemo(GLUCOSE)%controls_growth
ndivide = 0
do kcell = 1,nlist0
	if (cell_list(kcell)%state == DEAD) cycle
	ityp = cell_list(kcell)%celltype
!	c_rate = log(2.0)/divide_time_mean(ityp)
!	r_mean = Vdivide0/(2*divide_time_mean(ityp))
!	if (cell_list(kcell)%growth_delay) then
!		first_cycle = (cell_list(kcell)%N_delayed_cycles_left == LQ(ityp)%growth_delay_N)	! cell has not divided since hit
!		if (use_growth_suppression) then
!			if (tnow < cell_list(kcell)%t_growth_delay_end) cycle	! growth suppression
!		else
!			tdelay = cell_list(kcell)%dt_delay
!			c_rate = log(2.0)/(divide_time_mean(ityp) + tdelay)		
!			r_mean = Vdivide0/(2*(divide_time_mean(ityp) + tdelay))	! growth rate reduction
!		endif
!	endif
	if (cell_list(kcell)%volume < cell_list(kcell)%divide_volume) then	! still growing - not delayed
		C_O2 = cell_list(kcell)%conc(OXYGEN)
		C_glucose = cell_list(kcell)%conc(GLUCOSE)
		metab_O2 = O2_metab(C_O2)
		metab_glucose = glucose_metab(C_glucose)
		if (glucose_growth) then
			metab = metab_O2*metab_glucose
		else
			metab = metab_O2
		endif
		if (use_V_dependence) then
			dVdt = c_rate(ityp)*cell_list(kcell)%volume*metab
		else
			dVdt = r_mean(ityp)*metab
		endif
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
	endif
	
	if (cell_list(kcell)%volume > cell_list(kcell)%divide_volume) then	! time to divide
		! try letting tagged cell die only after end of G2/M delay
!		if (cell_list(kcell)%radiation_tag) then
!			R = par_uni(kpar)
!			if (R < cell_list(kcell)%p_rad_death) then
!				call CellDies(kcell)
!				Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
!				cycle
!			endif
		!
!			if (cell_list(kcell)%growth_delay) then
!				if (.not.cell_list(kcell)%G2_M) then	! this is the first time to divide
!					cell_list(kcell)%t_growth_delay_end = tnow + cell_list(kcell)%dt_delay
!					cell_list(kcell)%G2_M = .true.
!				endif
!				if (tnow > cell_list(kcell)%t_growth_delay_end) then
!					cell_list(kcell)%G2_M = .false.
!				else
!					cycle
!				endif
!			endif
!		endif
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
!		if (cell_list(kcell)%radiation_tag .and..not.cell_list(kcell)%G2_M) then
!			if (cell_list(kcell)%growth_delay) then
!				if (.not.cell_list(kcell)%G2_M) then	! this is the first time to divide
!					cell_list(kcell)%t_growth_delay_end = tnow + cell_list(kcell)%dt_delay
!					cell_list(kcell)%G2_M = .true.
!					write(*,'(a,i6,2f8.0)') 't_growth_delay_end: ',kcell,cell_list(kcell)%t_growth_delay_end,cell_list(kcell)%dt_delay
!				endif
!				if (tnow > cell_list(kcell)%t_growth_delay_end) then
!					write(*,*) 'tnow > t_growth_delay_end: ', kcell
!					cell_list(kcell)%G2_M = .false.
!				else
!					cycle
!				endif
!			endif
!		endif

		if (cell_list(kcell)%growth_delay) then
			if (cell_list(kcell)%G2_M) then
				if (tnow > cell_list(kcell)%t_growth_delay_end) then
					cell_list(kcell)%G2_M = .false.
				else
					cycle
				endif
			else
				cell_list(kcell)%t_growth_delay_end = tnow + cell_list(kcell)%dt_delay
				cell_list(kcell)%G2_M = .true.
				cycle
			endif
		endif
		! try moving death prob test to here
		if (cell_list(kcell)%radiation_tag) then
			R = par_uni(kpar)
			if (R < cell_list(kcell)%p_rad_death) then
				call CellDies(kcell)
				Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
				cycle
			endif
		endif
		!
!		elseif (radiation_dosed) then	! non-tagged cells also experience growth delay
!			if (cell_list(kcell)%G2_M) then
!				if (tnow > cell_list(kcell)%t_growth_delay_end) then
!					cell_list(kcell)%G2_M = .false.
!				else
!					cycle
!				endif
!			else
!				cell_list(kcell)%t_growth_delay_end = tnow + cell_list(kcell)%dt_delay
!				cell_list(kcell)%G2_M = .true.
!				cycle
!			endif		
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

!call logger('CellGrowth')
ok = .true.
nlist0 = nlist
tnow = istep*DELTA_T
c_rate(1:2) = log(2.0)/divide_time_mean(1:2)		! Note: to randomise divide time need to use random number, not mean!
r_mean(1:2) = Vdivide0/(2*divide_time_mean(1:2))
glucose_growth = chemo(GLUCOSE)%controls_growth
ndivide = 0
do kcell = 1,nlist0
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
	if (use_V_dependence) then
		dVdt = c_rate(ityp)*cell_list(kcell)%volume*metab
	else
		dVdt = r_mean(ityp)*metab
	endif
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

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine SetInitialGrowthRate
integer :: kcell, ityp
real(REAL_KIND) :: c_rate(2), r_mean(2), C_O2, C_glucose, metab, metab_O2, metab_glucose, dVdt
logical :: glucose_growth

c_rate(1:2) = log(2.0)/divide_time_mean(1:2)
r_mean(1:2) = Vdivide0/(2*divide_time_mean(1:2))
glucose_growth = chemo(GLUCOSE)%controls_growth
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	ityp = cell_list(kcell)%celltype
	C_O2 = chemo(OXYGEN)%bdry_conc
	C_glucose = cell_list(kcell)%conc(GLUCOSE)
	metab_O2 = O2_metab(C_O2)
	metab_glucose = glucose_metab(C_glucose)
	if (glucose_growth) then
		metab = metab_O2*metab_glucose
	else
		metab = metab_O2
	endif
	if (use_V_dependence) then
		dVdt = c_rate(ityp)*cell_list(kcell)%volume*metab
	else
		dVdt = r_mean(ityp)*metab
	endif
	if (suppress_growth) then	! for checking solvers
		dVdt = 0
	endif
	cell_list(kcell)%dVdt = dVdt
enddo
end subroutine

!-----------------------------------------------------------------------------------------
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
real(REAL_KIND) :: cfse0, cfse1
integer :: freesite(27,3)
type (boundary_type), pointer :: bdry

if (dbug) then
	write(logmsg,*) 'CellDivider: ',kcell0
	call logger(logmsg)
endif
ok = .true.
tnow = istep*DELTA_T
cell_list(kcell0)%t_divide_last = tnow
V0 = cell_list(kcell0)%volume
cell_list(kcell0)%volume = V0/2
cfse0 = cell_list(kcell0)%CFSE
cell_list(kcell0)%CFSE = generate_CFSE(cfse0/2)
cfse1 = cfse0 - cell_list(kcell0)%CFSE
cell_list(kcell0)%t_hypoxic = 0

R = par_uni(kpar)
cell_list(kcell0)%divide_volume = Vdivide0 + dVdivide*(2*R-1)
cell_list(kcell0)%M = cell_list(kcell0)%M/2
!write(nflog,'(a,i6,2f8.2)') 'divide: ',kcell0,cell_list(kcell0)%volume,cell_list(kcell0)%divide_volume
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
		if (divide_option == DIVIDE_USE_CLEAR_SITE) then	! use this site for the progeny cell
			site01 = bestsite
		else	! random choice
			j = random_int(1,nfree,0)
			site01 = freesite(j,:)
		endif
		call CloneCell(kcell0,kcell1,site01,ok)
		if (.not.ok) then
			call logger('Error: CloneCell: vacant site')
		endif
		cell_list(kcell1)%CFSE = cfse1
		Nreuse = Nreuse + 1
		return
	endif
endif

do
	j = random_int(1,6,kpar)
	site01 = site0 + neumann(:,j)
	if (site01(3) < zmin) then
		cycle
	endif
	if (occupancy(site01(1),site01(2),site01(3))%indx(1) == OUTSIDE_TAG) then	! site01 is outside, use it directly
		npath = 0
		site1 = site0
		if (site01(3) < zmin) then
			write(logmsg,*) 'CellDivider: OUTSIDE: z < zmin'
			call logger(logmsg)
			ok = .false.
			return
		endif
	elseif (bdrylist_present(site01,bdrylist)) then	! site01 is on the boundary
		npath = 1
		site1 = site01
		path(:,1) = site01
		if (site01(3) < zmin) then
			write(logmsg,*) 'CellDivider: boundary: z < zmin'
			call logger(logmsg)
			ok = .false.
			return
		endif
	else
		if (dbug) call logger('ChooseBdrySite')
		call ChooseBdrysite(site01,site1,ok)
		if (dbug) call logger('did ChooseBdrySite')
		if (.not.ok) return
		if (site01(3) < zmin) then
			write(logmsg,*) 'CellDivider: chosen: z < zmin'
			call logger(logmsg)
			ok = .false.
			return
		endif
		if (occupancy(site1(1),site1(2),site1(3))%indx(1) == 0) then
			write(logmsg,*) 'CellDivider: after choose_bdrysite: site1 is VACANT: ',site1
			call logger(logmsg)
			ok = .false.
			return
		endif
		call SelectPath(site0,site01,site1,path,npath,ok)
		if (.not.ok) then
			write(logmsg,*) 'CellDivider: SelectPath fails with site0, site01: ',site0,site01
			call logger(logmsg)
			ok = .false.
			return
		endif
		! path(:,:) goes from site01 to site1, which is a bdry site or adjacent to a vacant site
	!	if (.not.isbdry(site1)) then
	!	    call logger('should be bdry or vacant, is not')
	!	    stop
	!	endif
	endif
	exit
enddo

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
else
	! site01 is the destination site for the new cell
	site2 = site01
	Nsites = Nsites + 1
endif

call GetPathMass(site0,site01,path,npath,M0)
if (npath > 0) then
	call PushPath(path,npath)
!	if (bdry_debug) call CheckBdryList('after PushPath')
endif
call CloneCell(kcell0,kcell1,site01,ok)
if (.not.ok) then
	call logger('Error: CloneCell: pushed site')
	return
endif
cell_list(kcell1)%CFSE = cfse1

! Cell division completed
! Adjust concentrations, boundary list

call GetPathMass(site0,site01,path,npath,M1)
do ichemo = 1,MAX_CHEMO
	if (M0(ichemo) >= 0 .and. M1(ichemo) > 0) then
		alpha(ichemo) = M0(ichemo)/M1(ichemo)	! scaling for concentrations on the path
	else
		alpha(ichemo) = 1
	endif
enddo
call ScalePathConcentrations(site0,site01,path,npath,alpha)

!if (npath > 0) then
!	call FixPathConcentrations1(path,npath)
!endif
!! Now adjust extracellular concentrations for the parent cell site to ensure mass conservation
!Cex = occupancy(site0(1),site0(2),site0(3))%C(:)
!occupancy(site0(1),site0(2),site0(3))%C(:) = (Cex*(Vsite_cm3 - V0) + occupancy(site01(1),site01(2),site01(3))%C(:)*V0/2)/(Vsite_cm3 - V0/2)

call SetRadius(Nsites)
!call extendODEdiff(site2)
! Now need to fix the bdrylist.  
! site1 was on the boundary, but may no longer be.
! site2 may be now on the boundary
! First add site2
if (isbdry(site2)) then   ! add it to the bdrylist
	if (dbug) then
		write(logmsg,*) 'add site2 to bdrylist: ',site2
		call logger(logmsg)
	endif
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
endif
! Now check sites near site2 that may have changed their boundary status (including site1)
!do j = 1,6
!	site = site2 + neumann(:,j)
do j = 1,27
	if (j == 14) cycle
	site = site2 + jumpvec(:,j)
	if (isbdry(site)) then
		if (.not.bdrylist_present(site,bdrylist)) then	! add it
			allocate(bdry)
			bdry%site = site
			nullify(bdry%next)
			call bdrylist_insert(bdry,bdrylist)
			occupancy(site(1),site(2),site(3))%bdry => bdry
!			if (bdry_debug) call CheckBdryList('after bdrylist_insert')
		endif
	else
		if (bdrylist_present(site,bdrylist)) then	! remove it
			call bdrylist_delete(site,bdrylist)
			nullify(occupancy(site(1),site(2),site(3))%bdry)
!			if (bdry_debug) call CheckBdryList('after bdrylist_delete')
		endif
	endif
enddo
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
		V = cell_list(kcell)%volume*Vcell_cm3
	else
		V = 0
	endif
	mass(ic) = mass(ic) + C*(Vsite_cm3 - V)
	if (npath == 0) then
		C = occupancy(site01(1),site01(2),site01(3))%C(ic)
		kcell = occupancy(site01(1),site01(2),site01(3))%indx(1)
		if (kcell > 0) then
			V = cell_list(kcell)%volume*Vcell_cm3
		else
			V = 0
		endif
		mass(ic) = mass(ic) + C*(Vsite_cm3 - V)
	else
		do k = 1, npath
			site = path(:,k)
			C = occupancy(site01(1),site01(2),site01(3))%C(ic)
			kcell = occupancy(site(1),site(2),site(3))%indx(1)
			if (kcell > 0) then
				V = cell_list(kcell)%volume*Vcell_cm3
			else
				V = 0
			endif
			mass(ic) = mass(ic) + C*(Vsite_cm3 - V)
		enddo
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! NOTE: This is not satisfactory, because it can create a concentration increase -
! but it is not clear what should be done!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------------------------
subroutine ScalePathConcentrations(site0,site01,path,npath,alpha)
integer :: site0(3),site01(3),path(3,200),npath
real(REAL_KIND) :: alpha(:)
integer :: k, site(3), kcell, ichemo

do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	if (chemo(ichemo)%constant) cycle
	occupancy(site0(1),site0(2),site0(3))%C(ichemo) = alpha(ichemo)*occupancy(site0(1),site0(2),site0(3))%C(ichemo)
	if (npath == 0) then
		occupancy(site01(1),site01(2),site01(3))%C(ichemo) = alpha(ichemo)*occupancy(site01(1),site01(2),site01(3))%C(ichemo)
	else
		do k = 1, npath
			site01 = path(:,k)
			occupancy(site01(1),site01(2),site01(3))%C(ichemo) = alpha(ichemo)*occupancy(site01(1),site01(2),site01(3))%C(ichemo)
			if (isnan(occupancy(site01(1),site01(2),site01(3))%C(ichemo))) then
				write(nflog,*) 'ScalePathConcentrations: isnan: ichemo,site: ',ichemo,site01,alpha(ichemo)
			endif
			if (chemo(ichemo)%constant .and. occupancy(site01(1),site01(2),site01(3))%C(ichemo) /= chemo(ichemo)%bdry_conc) then
				write(nflog,*) 'ScalePathConcentrations: constant but C != bdry_conc: ', &
					ichemo,occupancy(site01(1),site01(2),site01(3))%C(ichemo),chemo(ichemo)%bdry_conc
				stop
			endif
			
		enddo
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! The extracellular concentrations along the path are adjusted sequentially, starting from
! the last-but-one site and working backwards to the first, site01
! NOT USED
!-----------------------------------------------------------------------------------------
subroutine FixPathConcentrations1(path,npath)
integer :: path(3,200),npath
integer :: k, site1(3), site2(3), kcell, ic
real(REAL_KIND) :: V, Cex(MAX_CHEMO)

do k = npath-1,1,-1
	site1 = path(:,k)
	kcell = occupancy(site1(1),site1(2),site1(3))%indx(1)
	V = cell_list(kcell)%volume*Vcell_cm3
	Cex = occupancy(site1(1),site1(2),site1(3))%C
	site2 = path(:,k+1)
	occupancy(site1(1),site1(2),site1(3))%C = ((Vsite_cm3 - V)*Cex + V*occupancy(site2(1),site2(2),site2(3))%C)/Vsite_cm3
	do ic = 1,MAX_CHEMO
		if (.not.chemo(ic)%used) cycle
		if (occupancy(site1(1),site1(2),site1(3))%C(ic) < 0) then
			occupancy(site1(1),site1(2),site1(3))%C(ic) = max(occupancy(site2(1),site2(2),site2(3))%C(ic), Cex(ic))
		endif
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Need to choose a site on the boundary in some sense near site0, and also to preserve
! the required shape.
! alpha_max is calculated as a function of the fractional distance of site0 from Centre,
! i.e. of r/Radius.
! This determines the subset of boundary sites that are examined.
! The criterion for a site to be a candidate is:
! the angle between v = site - Centre and vc = site0 - Centre must be < alpha_max
! In fact the decision is made based on cos(angle) > cos(alpha_max).
! From this set of sites we choose the one that departs most from the desired boundary
! in the negative sense.  In the case of a sphere this means the site with the least
! distance from Centre, but for a squashed sphere we need a different way to choose.
! 
!-----------------------------------------------------------------------------------------
subroutine ChooseBdrysite(site0,site1,ok)
integer :: site0(3), site1(3)
logical :: ok
integer :: site(3), sitemin(3)
real(REAL_KIND) :: vc(3), v(3), r, rfrac, d, alpha_max, cosa_min, dmin, cosa
real(REAL_KIND) :: z, r2, sin2, cos2, d2, dd, dsq, dsqmax, tempCentre(3)
logical :: hit
type (boundary_type), pointer :: bdry

tempCentre = Centre
if (is_dropped .and. site0(3) < Centre(3)) then
	tempCentre(3) = (site0(3) + 2*Centre(3))/3
endif
vc = site0 - tempCentre
r = norm(vc)
vc = vc/r
rfrac = r/Radius
alpha_max = getAlphaMax(rfrac)
cosa_min = cos(alpha_max)
dmin = 1.0e10
dsqmax = -1.0e10
hit = .false.
bdry => bdrylist
do while ( associated ( bdry )) 
    site = bdry%site
    v = site - tempCentre
    d = norm(v)	
	cosa = dot_product(v,vc)/d
	if (cosa > cosa_min) then	! Note: in squashed case we need to worry about cells near the surface
		hit = .true.
		if (is_dropped) then
			z = site0(3) + cdrop - zmin
			if (z < 2*bdrop*Radius) then
				cos2 = (1 - (z)/(bdrop*Radius))**2
				sin2 = 1 - cos2
			else
				r2 = (site0(1)-Centre(1))**2 + (site0(2)-Centre(2))**2
				sin2 = r2/(adrop*Radius)**2
				cos2 = 1 - sin2
			endif
			d2 = adrop*adrop*sin2 + bdrop*bdrop*cos2
			if (d2 < 0) then
				write(nflog,*) 'd2 < 0: cos2, sin2: ',cos2,sin2
				write(nflog,*) 'site0(3),cdrop,zmin: ',site0(3),cdrop,zmin,(site0(3) + cdrop - zmin)/(bdrop*Radius)
				stop
			endif
			dd = Radius*sqrt(d2)	! this is the desired distance for this z
			dsq = dd - d	! could use dd/d or dd-d
			if (dsq > dsqmax) then
				dsqmax = dsq
				sitemin = site
			endif
		else	! this is OK for a sphere, not for the squashed sphere
			if (d < dmin) then
				dmin = d
				sitemin = site
			endif
		endif
	endif
    bdry => bdry%next
enddo
if (.not.hit) then
	write(logmsg,*) 'Error: choose_bdrysite: no candidate bdry site'
	call logger(logmsg)
	ok = .false.
	return
endif
site1 = sitemin
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
! r is the fractional distance from the sphere centre = (distance from centre)/Radius
! The parameter alphamax varies: alpha1 -> alpha2 as r: r1 -> r2
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function getAlphaMax(r)
real(REAL_KIND) :: r
real(REAL_KIND) :: alf, alphamax
real(REAL_KIND) :: r1 = 0.0, r2 = 0.8, alpha1 = PI/2, alpha2 = PI/6

if (r < r1) then
	alphamax = alpha1
elseif (r > r2) then
	alphamax = alpha2
else
	alf = (r-r1)/(r2-r1)
	alphamax = (1-alf)*alpha1 + alf*alpha2
endif
getAlphaMax = alphamax
end function


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_get_path
integer :: site0(3), site1(3), site2(3), path(3,200), npath, kpar=0
integer :: x, y, z, it, k
logical :: ok

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
	call SelectPath(site0,site1,site2,path,npath,ok)
	write(nflog,*) 'path: ',npath
	do k = 1,npath
		write(nflog,'(i3,2x,3i4)') path(:,k)
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
subroutine SelectPath(site0,site1,site2,path,npath,ok)
integer :: site0(3),site1(3), site2(3), path(3,200), npath
logical :: ok
integer :: v(3), jump(3), site(3), k, j, jmin, indx
integer, parameter :: maxits = 100
real(REAL_KIND) :: r, d2, d2min
logical :: hit

if (occupancy(site2(1),site2(2),site2(3))%indx(1) <= OUTSIDE_TAG) then
    call logger('SelectPath: site2 is OUTSIDE or UNREACHABLE')
    ok = .false.
    return
endif
if (occupancy(site2(1),site2(2),site2(3))%indx(1) == 0) then
    call logger('SelectPath: site2 is VACANT')
    ok = .false.
    return
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
		indx = occupancy(v(1),v(2),v(3))%indx(1)
		if (indx == 0 .or. indx == OUTSIDE_TAG) then	! outside or vacant - we'll use this!
			site2 = site
			hit = .true.
			exit
		endif
		v = site2 - v
		d2 = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
		if (d2 < d2min) then
			d2min = d2
			jmin = j
		endif
	enddo
	if (hit) exit
	if (jmin == 0) then
	    call logger('get path: stuck')
		write(logmsg,*) 'SelectPath: stuck for site0,site1,site2: ',site0,site1,site2
		call logger(logmsg)
	    ok = .false.
		return
	endif
	site = site + jumpvec(:,jmin)
	k = k+1
	if (k > maxits) then
		write(logmsg,*) 'SelectPath: iteration limit exceeded: ',maxits,'  for site0,site1,site2: ',site0,site1,site2
		call logger(logmsg)
	    ok = .false.
		return
	endif
	path(:,k) = site
	if (site(1) == site2(1) .and. site(2) == site2(2) .and. site(3) == site2(3)) exit
enddo
npath = k
ok = .true.
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
	site2 = path(:,k+1)
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
! The concentrations of constituents must be halved NOT TRUE!!!!!!!!!!!! 
!-----------------------------------------------------------------------------------------
subroutine CloneCell(kcell0,kcell1,site1,ok)
integer :: kcell0, kcell1, site1(3), ityp, idrug
logical :: ok
integer :: kpar = 0
real(REAL_KIND) :: tnow, R

ok = .true.
tnow = istep*DELTA_T

if (use_gaplist .and. ngaps > 0) then
    kcell1 = gaplist(ngaps)
    ngaps = ngaps - 1
else
    nlist = nlist + 1
    if (nlist > max_nlist) then
		write(logmsg,*) 'Error: Dimension of cell_list() has been exceeded: increase max_nlist and rebuild'
		call logger(logmsg)
		ok = .false.
		return
	endif
    kcell1 = nlist
endif

!nlist = nlist + 1
!if (nlist > max_nlist) then
!	call logger('Dimension of cell_list() has been exceeded: increase max_nlist and rebuild')
!	ok = .false.
!	return
!endif
!kcell1 = nlist

ityp = cell_list(kcell0)%celltype
Ncells = Ncells + 1
Ncells_type(ityp) = Ncells_type(ityp) + 1
!if (cell_list(kcell0)%generation > 1 .and. cell_list(kcell0)%radiation_tag) then 
!	write(*,*) 'radiation-tagged cell gen > 1 divides: ',kcell0,cell_list(kcell0)%generation 
!	stop
!endif
cell_list(kcell0)%generation = cell_list(kcell0)%generation + 1
if (cell_list(kcell0)%growth_delay) then
	cell_list(kcell0)%N_delayed_cycles_left = cell_list(kcell0)%N_delayed_cycles_left - 1
	cell_list(kcell0)%growth_delay = (cell_list(kcell0)%N_delayed_cycles_left > 0)
endif
cell_list(kcell0)%G2_M = .false.
cell_list(kcell1)%celltype = cell_list(kcell0)%celltype
cell_list(kcell1)%state = cell_list(kcell0)%state
cell_list(kcell1)%generation = cell_list(kcell0)%generation
cell_list(kcell1)%site = site1
!cell_list(kcell1)%ID = lastID
cell_list(kcell1)%ID = cell_list(kcell0)%ID
!if (cell_list(kcell0)%ID == 582) then
!	write(*,*) 'New cell with ID=582: ',kcell1
!endif
cell_list(kcell1)%p_rad_death = cell_list(kcell0)%p_rad_death
cell_list(kcell1)%p_drug_death = cell_list(kcell0)%p_drug_death
cell_list(kcell1)%radiation_tag = cell_list(kcell0)%radiation_tag
if (cell_list(kcell1)%radiation_tag) then
	Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
endif
cell_list(kcell1)%drug_tag = cell_list(kcell0)%drug_tag
do idrug = 1,ndrugs_used
	if (cell_list(kcell1)%drug_tag(idrug)) then
		Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) + 1
	endif
enddo
cell_list(kcell1)%anoxia_tag = .false.
cell_list(kcell1)%exists = .true.
cell_list(kcell1)%active = .true.
cell_list(kcell1)%growth_delay = cell_list(kcell0)%growth_delay
if (cell_list(kcell1)%growth_delay) then
	cell_list(kcell1)%dt_delay = cell_list(kcell0)%dt_delay
	cell_list(kcell1)%N_delayed_cycles_left = cell_list(kcell0)%N_delayed_cycles_left
endif
cell_list(kcell1)%G2_M = .false.
cell_list(kcell1)%t_divide_last = tnow
cell_list(kcell1)%dVdt = cell_list(kcell0)%dVdt
cell_list(kcell1)%volume = cell_list(kcell0)%volume
R = par_uni(kpar)
cell_list(kcell1)%divide_volume = Vdivide0 + dVdivide*(2*R-1)
cell_list(kcell1)%t_hypoxic = 0
cell_list(kcell1)%conc = cell_list(kcell0)%conc
cell_list(kcell1)%Cex = cell_list(kcell0)%Cex
cell_list(kcell1)%dCdt = cell_list(kcell0)%dCdt
cell_list(kcell1)%dMdt = cell_list(kcell0)%dMdt
cell_list(kcell1)%M = cell_list(kcell0)%M
occupancy(site1(1),site1(2),site1(3))%indx(1) = kcell1
!write(nflog,'(a,i6,2f8.2)') 'cloned: ',kcell1,cell_list(kcell1)%volume,cell_list(kcell1)%divide_volume
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
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine CheckUnreachable
integer :: kcell, site(3), indx(2)

do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	site = cell_list(kcell)%site
	indx = occupancy(site(1),site(2),site(3))%indx
	if (indx(1) == UNREACHABLE_TAG .or. indx(1) == OUTSIDE_TAG) then
		write(nflog,'(a,6i6)') 'CheckUnreachable: bad indx: ',kcell,Ncells,site,indx(1)
		stop
	endif
	if (site(3) < zmin) then	
		write(nflog,'(a,6i6)') 'CheckUnreachable: bad z: ',kcell,Ncells,site,indx(1)
		stop
	endif
enddo
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
	write(nflog,*) 'Error: adjust: site is not OUTSIDE: ',site0
	stop
endif
r0 = cdistance(site0)
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
	cell_list(kcell)%site = site0
	occupancy(site0(1),site0(2),site0(3))%indx(1) = kcell
	occupancy(sitemax(1),sitemax(2),sitemax(3))%indx(1) = OUTSIDE_TAG
endif
end subroutine

end module
