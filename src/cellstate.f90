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
subroutine grow_cells(dose,dt,ok)
real(REAL_KIND) :: dose, dt
logical :: ok

!call logger('grow_cells')
ok = .true.
if (use_radiation .and. dose > 0) then
	call irradiation(dose)
endif
if (use_division) then
	call cell_division(dt,ok)
	if (.not.ok) return
endif
if (use_death) then
	call cell_death(dt)
endif
if (use_migration) then
	call cell_migration
endif
end subroutine

!-----------------------------------------------------------------------------------------
! The O2 concentration to use with cell kcell is either the intracellular concentration,
! or is use_extracellular_O2, the corresponding extracellular concentration
!-----------------------------------------------------------------------------------------
subroutine get_O2conc(kcell, C_O2)
integer :: kcell
real(REAL_KIND) :: C_O2
integer :: iv, site(3)
real(REAL_KIND) :: tnow

if (use_extracellular_O2) then
	iv = cell_list(kcell)%iv
	if (iv < 1) then
!		write(logmsg,*) 'get_O2conc: ',kcell,site,iv
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
subroutine irradiation(dose)
real(REAL_KIND) :: dose
integer :: kcell, site(3), iv, kpar=0
real(REAL_KIND) :: C_O2, OER_alpha_d, OER_beta_d, expon, kill_prob, R

LQ%OER_am = 2.5
LQ%OER_bm = 3.0
LQ%alpha_H = 0.0473
LQ%beta_H = 0.0017
LQ%K_ms = 4.3e-3	! mM
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	if (cell_list(kcell)%radiation_tag) cycle	! we do not tag twice (yet)
	call get_O2conc(kcell,C_O2)
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
!-----------------------------------------------------------------------------------------
subroutine cell_migration
integer :: kcell, j, indx, site0(3), site(3), jmax
real(REAL_KIND) :: C0(MAX_CHEMO), C(MAX_CHEMO), v0, v, vmax, d0, d

do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	site0 = cell_list(kcell)%site
	C0 = occupancy(site0(1),site0(2),site0(3))%C(:)
	v0 = site_value(C0)
	d0 = cdistance(site0)
	jmax = 0
	vmax = -1.0e10
	do j = 1,27
		if (j == 14) cycle
		site = site0 + jumpvec(:,j)
		indx = occupancy(site(1),site(2),site(3))%indx(1)
		if (indx < -100) then	! necrotic site
			C = occupancy(site(1),site(2),site(3))%C(:)
			v = site_value(C)
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
real(REAL_KIND) function site_value(C)
real(REAL_KIND) :: C(:)

site_value = C(OXYGEN)
end function

!-----------------------------------------------------------------------------------------
! Cells can be tagged to die, or finally die of anoxia, or they can be tagged for death 
! at division time if the drug is effective.
!-----------------------------------------------------------------------------------------
subroutine cell_death(dt)
real(REAL_KIND) :: dt
integer :: kcell, nlist0, site(3), i, kpar=0 
real(REAL_KIND) :: C_O2, kmet, Kd, dMdt, pdeath, tnow
logical :: use_SN30000

!call logger('cell_death')
if (chemo(DRUG_A)%used .and. DRUG_A == SN30000) then
    use_SN30000 = .true.
    Kd = SN30K%Kd
else
    use_SN30000 = .false.
endif
tnow = istep*DELTA_T	! seconds
nlist0 = nlist
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	if (cell_list(kcell)%anoxia_tag) then
!		write(logmsg,*) 'anoxia_tag: ',kcell,cell_list(kcell)%state,tnow,cell_list(kcell)%t_anoxia_die
!		call logger(logmsg)
		if (tnow >= cell_list(kcell)%t_anoxia_die) then
!			call logger('cell dies')
			call cell_dies(kcell)
			Nanoxia_dead = Nanoxia_dead + 1
			cycle
		endif
	else
!		C_O2 = cell_list(kcell)%conc(OXYGEN)
		call get_O2conc(kcell,C_O2)
		if (C_O2 < ANOXIA_FACTOR*MM_THRESHOLD) then
			cell_list(kcell)%t_hypoxic = cell_list(kcell)%t_hypoxic + dt
			if (cell_list(kcell)%t_hypoxic > t_anoxic_limit) then
				cell_list(kcell)%anoxia_tag = .true.						! tagged to die later
				cell_list(kcell)%t_anoxia_die = tnow + anoxia_death_delay	! time that the cell will die
				Nanoxia_tag = Nanoxia_tag + 1
			endif
		endif
	endif
	if (use_SN30000) then
	    kmet = (SN30K%C1 + SN30K%C2*SN30K%KO2/(SN30K%KO2 + C_O2))*SN30K%Kmet0
	    dMdt = kmet*cell_list(kcell)%conc(DRUG_A)
	    pdeath = Kd*dMdt*dt
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
subroutine test_cell_dies
integer :: kcell, i, kpar=0

do i = 1,10
	kcell = random_int(1,nlist,kpar)
	call cell_dies(kcell)
enddo
stop
end subroutine

!-----------------------------------------------------------------------------------------
! When a cell dies the site takes occupancy()%indx = -kcell.  This is a necrotic volume.
! Such necrotic sites migrate towards the blob centre.
! The cell contents should be released into the site.
!-----------------------------------------------------------------------------------------
subroutine cell_dies(kcell)
integer :: kcell
integer :: site(3)

cell_list(kcell)%state = DEAD
Ncells = Ncells - 1
!write(*,*) 'cell_dies: ',kcell,Ncells
site = cell_list(kcell)%site
occupancy(site(1),site(2),site(3))%indx(1) = -(100 + kcell)
! Adjust site concentrations?
call necrotic_migration(site)
end subroutine

!-----------------------------------------------------------------------------------------
! A necrotic site migrates towards the blob centre, stopping when another necrotic 
! site is reached
!-----------------------------------------------------------------------------------------
subroutine necrotic_migration(site0)
integer :: site0(3)
integer :: site1(3), site2(3), site(3), j, jmin, kcell, tmp_indx
real(REAL_KIND) :: d1, d2, dmin

!write(logmsg,*) 'necrotic_migration: site0: ',site0
!call logger(logmsg)
site1 = site0
do
	d1 = cdistance(site1)
	dmin = 1.0e10
	jmin = 0
	do j = 1,27
		if (j == 14) cycle
		site = site1 + jumpvec(:,j)
		if (occupancy(site(1),site(2),site(3))%indx(1) < 0) cycle	! do not swap with another necrotic site
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
!write(*,*) 'necrotic_migration: site2: ',site2
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_cell_division(ok)
logical :: ok
integer :: kcell, kpar=0
kcell = random_int(1,nlist,kpar)
call cell_divider(kcell,ok)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine cell_division(dt,ok)
real(REAL_KIND) :: dt
logical :: ok
integer :: kcell, nlist0, site(3), iv
integer :: divide_list(10000), ndivide, i
real(REAL_KIND) :: tnow, C_O2, metab, dVdt, vol0, r_mean, c_rate
character*(20) :: msg

!call logger('cell_division')
ok = .true.
nlist0 = nlist
tnow = istep*DELTA_T
c_rate = log(2.0)/divide_time_mean		! Note: to randomise divide time need to use random number, not mean!
r_mean = Vdivide0/(2*divide_time_mean)
ndivide = 0
do kcell = 1,nlist0
	if (cell_list(kcell)%state == DEAD) cycle
!	call get_O2conc(kcell,C_O2)
	C_O2 = cell_list(kcell)%conc(OXYGEN)
!	metab = max(0.0,C_O2)/(chemo(OXYGEN)%MM_C0 + C_O2)
    if (C_O2 > ODEdiff%C1_soft) then
	    metab = (C_O2-ODEdiff%deltaC_soft)/(chemo(OXYGEN)%MM_C0 + C_O2 - ODEdiff%deltaC_soft)
    elseif (C_O2 > 0) then
	    metab = ODEdiff%k_soft*C_O2*C_O2
	else
		metab = 0
	endif
	if (use_V_dependence) then
		dVdt = c_rate*cell_list(kcell)%volume*metab
	else
		dVdt = r_mean*metab
	endif
	if (istep > 1 .and. dVdt == 0) then
		write(nflog,'(a,2i6,5e12.3)') 'dVdt: ',istep,kcell,r_mean,c_rate,C_O2,metab,dVdt
	endif
	cell_list(kcell)%dVdt = dVdt
	vol0 = cell_list(kcell)%volume
	cell_list(kcell)%volume = vol0 + dVdt*dt
	cell_list(kcell)%conc = vol0*cell_list(kcell)%conc/cell_list(kcell)%volume
	if (cell_list(kcell)%volume > cell_list(kcell)%divide_volume) then
		if (cell_list(kcell)%radiation_tag) then
			call cell_dies(kcell)
			Nradiation_dead = Nradiation_dead + 1
!			write(*,*) 'Cell died: ',kcell,
			cycle
		endif
		if (cell_list(kcell)%drug_tag) then
			call cell_dies(kcell)
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
	call cell_divider(kcell, ok)
	if (.not.ok) return
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Note: updating of concentrations is now done in extendODEdiff()
! The dividing cell, kcell0, is at site0.
! A neighbour site site01 is chosen randomly.  This becomes the site for the daughter cell.
!-----------------------------------------------------------------------------------------
subroutine cell_divider(kcell0, ok)
integer :: kcell0
logical :: ok
integer :: kpar=0
integer :: j, k, kcell1, site0(3), site1(3), site2(3), site01(3), site(3), ichemo, nfree, bestsite(3)
integer :: npath, path(3,200)
real(REAL_KIND) :: tnow, R, v, vmax
logical :: freesite(27,3)
type (boundary_type), pointer :: bdry

!write(logmsg,*) 'cell_divider: ',kcell0
!call logger(logmsg)
ok = .true.
tnow = istep*DELTA_T
cell_list(kcell0)%t_divide_last = tnow
cell_list(kcell0)%volume = cell_list(kcell0)%volume/2
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
		if (occupancy(site01(1),site01(2),site01(3))%indx(1) < -100) then
			nfree = nfree + 1
			freesite(nfree,:) = site01
			v = site_value(occupancy(site01(1),site01(2),site01(3))%C(:))
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
		call add_cell(kcell0,kcell1,site01,ok)
		if (.not.ok) then
			call logger('Error: add_cell: vacant site')
		endif
		Nreuse = Nreuse + 1
		return
	endif
endif

ok = .false.
k = 0
do while (.not.ok)
    k = k + 1
!	j = random_int(1,26,kpar)       ! This generates disconnected cells at the boundary
!	if (j >= 14) j = j+1
!	site01 = site0 + jumpvec(:,j)
	j = random_int(1,6,kpar)
	site01 = site0 + neumann(:,j)
	if (k > 50 .or. occupancy(site01(1),site01(2),site01(3))%indx(1) >= -100) ok = .true.	! site01 is not necrotic
enddo
!if (dbug) write(*,*) 'cell_divider: ',kcell0,site0,occupancy(site0(1),site0(2),site0(3))%indx
if (occupancy(site01(1),site01(2),site01(3))%indx(1) == OUTSIDE_TAG) then	! site01 is outside, use it directly
	npath = 0
	site1 = site0
elseif (bdrylist_present(site01,bdrylist)) then	! site01 is on the boundary
	npath = 1
	site1 = site01
	path(:,1) = site01
else
	call choose_bdrysite(site01,site1)
	call get_path(site01,site1,path,npath)
	! path(:,:) goes from site01 to site1, which is a bdry site
	if (.not.isbdry(site1)) then
	    call logger('should be bdry, is not')
	    stop
	endif
!	call clear_path(path,npath)
endif
!if (dbug) write(*,*) 'path: ',npath
!do k = 1,npath
!	if (dbug) write(*,'(i3,2x,3i4)') path(:,k)
!enddo

if (npath > 0) then
	! Need to choose an outside site near site1
	call get_outsidesite(site1,site2)
	! path(:,:) now goes from site01 to site2, which is an outside site next to site1
	npath = npath+1
	path(:,npath) = site2
!	write(*,'(a,3i4,i6)') 'outside site: ',site2,occupancy(site2(1),site2(2),site2(3))%indx(1)
	call push_path(path,npath)
	if (dbug) write(*,*) 'did push_path'
else
	site2 = site01
endif
Nsites = Nsites + 1
call SetRadius(Nsites)
!call extendODEdiff(site2)
!call InterpolateConc(site2)
!do ichemo = 1,MAX_CHEMO
!    occupancy(site2(1),site2(2),site2(3))%C(ichemo) = chemo(ichemo)%bdry_conc
!enddo
call add_cell(kcell0,kcell1,site01,ok)
if (.not.ok) then
	call logger('Error: add_cell: pushed site')
	return
endif

! Now need to fix the bdrylist.  
! site1 was on the boundary, but may no longer be.
! site2 should be now on the boundary
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
call adjust_medium
if (dbug) write(*,*) 'done!'
end subroutine

!-----------------------------------------------------------------------------------------
! Reduce medium volume and concentrations to account for the growth of the spheroid
! by one site.
!-----------------------------------------------------------------------------------------
subroutine adjust_medium
real(REAL_KIND) :: total(MAX_CHEMO)

total(DRUG_A:MAX_CHEMO) = chemo(DRUG_A:MAX_CHEMO)%bdry_conc*medium_volume
medium_volume = medium_volume - Vsite
chemo(DRUG_A:MAX_CHEMO)%bdry_conc = total(DRUG_A:MAX_CHEMO)/medium_volume
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
!    if (occupancy(site(1),site(2),site(3))%indx(1) < -100) then    ! necrotic
!        write(*,*) 'bdry site is necrotic'
!        stop
!    endif
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
real(REAL_KIND) :: r1 = 0.0, r2 = 0.8, alpha1 = PI/2, alpha2 = PI/6

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
! The first site in the list is site1, the last is site2
!-----------------------------------------------------------------------------------------
subroutine get_path(site1,site2,path,npath)
integer :: site1(3), site2(3), path(3,200), npath
integer :: v(3), jump(3), site(3), k, j, jmin
real(REAL_KIND) :: r, d2, d2min

if (occupancy(site2(1),site2(2),site2(3))%indx(1) == OUTSIDE_TAG) then
    write(*,*) 'site2 is OUTSIDE'
    stop
endif
if (occupancy(site2(1),site2(2),site2(3))%indx(1) < -100) then
    write(*,*) 'site2 is NECROTIC'
!    stop
endif
if (dbug) write(*,'(a,3i4,2x,3i4)') 'site1, site2: ',site1, site2
k = 1
site = site1
path(:,k) = site
do 
	d2min = 1.0e10
	jmin = 0
!	do j = 1,6
!		jump = neumann(:,j)
	do j = 1,27
		if (j == 14) cycle
		jump = jumpvec(:,j)
		v = site + jump
		if (occupancy(v(1),v(2),v(3))%indx(1) == OUTSIDE_TAG) cycle
!		if (occupancy(v(1),v(2),v(3))%indx(1) < -100) cycle		! necrotic
		v = site2 - v
		d2 = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
		if (d2 < d2min) then
			d2min = d2
			jmin = j
		endif
	enddo
	if (jmin == 0) then
	    call logger('get path: stuck')
	    stop
	endif
!	site = site + neumann(:,jmin)
	site = site + jumpvec(:,jmin)
	k = k+1
!	if (k>100) write(*,'(11i4,f8.1)') k,site1,site2,site,jmin,d2min
	path(:,k) = site
!	if (occupancy(site(1),site(2),site(3))%indx(1) < -100) then		! necrotic
!		write(logmsg,*) 'Error: get_path: necrotic site in path: ',site
!		call logger(logmsg)
!		stop
!	endif
	if (site(1) == site2(1) .and. site(2) == site2(2) .and. site(3) == site2(3)) exit
enddo
npath = k

end subroutine

!-----------------------------------------------------------------------------------------
! If there are any necrotic sites in the path, they are first shifted closer to the centre,
! then the path is adjusted.
!-----------------------------------------------------------------------------------------
subroutine clear_path(path,npath)
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
! Need to deal with DEAD cells, i.e. sites that do not hold cells (necrotic region).
! Need to modify the ODEdiff variables, at least %ivar
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
!	if (kcell < 0) then
!		write(*,'(a,i6,5i4)') 'push_path: kcell<0: k,npath,site1: ',kcell,k,npath,site1
!		write(*,'(a,3i4)') 'moved to: ',site2
!		write(*,'(a,4i6)') '                  dividing cell: ',kcell_dividing,cell_list(kcell_dividing)%site
!		stop
!	endif
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
subroutine add_cell(kcell0,kcell1,site1,ok)
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

end module
