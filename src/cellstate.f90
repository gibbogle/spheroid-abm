! Cancer cell state development

module cellstate
use global
use boundary
use fields
use ode_diffuse
implicit none

real(REAL_KIND), parameter :: Vdivide = 1.6
real(REAL_KIND), parameter :: dVdivide = 0.05
real(REAL_KIND), parameter :: CO2_DEATH_THRESHOLD = 0.01
integer :: kcell_dividing = 0

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine grow_cells(dt)
real(REAL_KIND) :: dt

if (use_division) then
	call cell_division(dt)
endif
if (use_death) then
	call cell_death
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine cell_death
integer :: kcell, nlist0, site(3), i

nlist0 = nlist
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	site = cell_list(kcell)%site
	i = ODEdiff%ivar(site(1),site(2),site(3))
	if (allstate(i,OXYGEN) < CO2_DEATH_THRESHOLD) then
		call cell_dies(kcell)
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
!-----------------------------------------------------------------------------------------
subroutine cell_dies(kcell)
integer :: kcell
integer :: site(3)

!write(*,*) 'cell_dies: ',kcell
cell_list(kcell)%state = DEAD
Ncells = Ncells - 1
site = cell_list(kcell)%site
occupancy(site(1),site(2),site(3))%indx(1) = -(100 + kcell)
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

!write(*,*) 'necrotic_migration: site0: ',site0
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
subroutine test_cell_division
integer :: kcell, kpar=0
kcell = random_int(1,nlist,kpar)
call cell_divider(kcell)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine cell_division(dt)
real(REAL_KIND) :: dt
integer :: kcell, nlist0, site(3), iv
real(REAL_KIND) :: tnow, CO2, metab, dVdt, rmax

nlist0 = nlist
tnow = istep*DELTA_T
rmax = Vdivide/(2*divide_time_mean)
do kcell = 1,nlist0
	if (cell_list(kcell)%state == DEAD) cycle
	site = cell_list(kcell)%site
	iv = ODEdiff%ivar(site(1),site(2),site(3))
	CO2 = allstate(iv,OXYGEN)
	metab = max(0.0,CO2)/(chemo(OXYGEN)%MM_C0 + CO2)
	dVdt = rmax*metab
	cell_list(kcell)%volume = cell_list(kcell)%volume + dVdt*dt
	if (cell_list(kcell)%volume > cell_list(kcell)%divide_volume) then
		kcell_dividing = kcell
		call cell_divider(kcell)
	endif
!	if (cell_list(kcell)%t_divide_next <= tnow) then
!		kcell_dividing = kcell
!		call cell_divider(kcell)
!	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Note: updating of concentrations is now done in extendODEdiff()
!-----------------------------------------------------------------------------------------
subroutine cell_divider(kcell0)
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
call extendODEdiff(site2)

if (dbug) write(*,*) 'did push_path'
cell_list(kcell0)%t_divide_last = tnow
!cell_list(kcell0)%t_divide_next = tnow + DivideTime()
cell_list(kcell0)%volume = cell_list(kcell0)%volume/2
R = par_uni(kpar)
cell_list(kcell0)%divide_volume = Vdivide + dVdivide*(2*R-1)
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
!if (kcell_dividing == 3057) write(*,'(a,3i4,3f6.3)') 'site0: ',site0,v_aim
k = 1
site1 = site0
path(:,k) = site1
do 
	dmax = 0
	do j = 1,26
		if (j == 14) cycle
!		jump = neumann(:,j)
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
!if (kcell_dividing == 3057) write(*,'(i6,2x,3i3,2x,3f6.3,2x,f8.4)') j,jump,v_try,d
		if (d-del > dmax) then
			dmax = d-del
			jmax = j
		endif
	enddo
!	site1 = site1 + neumann(:,jmax)
	site1 = site1 + jumpvec(:,jmax)
	k = k+1
	path(:,k) = site1
!if (kcell_dividing == 3057) write(*,'(i3,2x,3i4,2x,i2,2x,3i3)') k,site1,jmax,neumann(:,jmax)
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
		if (occupancy(v(1),v(2),v(3))%indx(1) == OUTSIDE_TAG) cycle
		if (occupancy(v(1),v(2),v(3))%indx(1) < -100) cycle		! necrotic
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
! Need to deal with DEAD cells, i.e. sites that do not hold cells (necrotic region).
!-----------------------------------------------------------------------------------------
subroutine push_path(path,npath)
integer :: path(3,200),npath
integer :: k, site1(3), site2(3), kcell

do k = npath-1,1,-1
	site1 = path(:,k)
	kcell = occupancy(site1(1),site1(2),site1(3))%indx(1)
	if (kcell < 0) then
		write(*,'(a,4i6)') 'Error: push_path: kcell<0: site: ',kcell,site1
		write(*,'(a,4i6)') '                  dividing cell: ',kcell_dividing,cell_list(kcell_dividing)%site
		stop
	endif
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
integer :: kpar = 0
real(REAL_KIND) :: tnow, R

tnow = istep*DELTA_T
lastID = lastID + 1
nlist = nlist + 1
Ncells = Ncells + 1
Nsites = Nsites + 1
kcell1 = nlist
cell_list(kcell1)%state = cell_list(kcell0)%state
cell_list(kcell1)%site = site
cell_list(kcell1)%ID = lastID
cell_list(kcell1)%exists = .true.
cell_list(kcell1)%t_divide_last = tnow
!cell_list(kcell1)%t_divide_next = tnow + DivideTime()
cell_list(kcell1)%volume = cell_list(kcell0)%volume
R = par_uni(kpar)
cell_list(kcell1)%divide_volume = Vdivide + dVdivide*(2*R-1)
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
