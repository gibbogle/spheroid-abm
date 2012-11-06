! Cancer cell state development

module cellstate
use global
use ode_diffuse
implicit none

contains

subroutine growcells

end subroutine

!-----------------------------------------------------------------------------------------
subroutine test_cell_death
integer :: kcell, i, kpar=0

do i = 1,10
	kcell = random_int(1,nlist,kpar)
	call cell_death(kcell)
enddo
stop
end subroutine

!-----------------------------------------------------------------------------------------
! When a cell dies the site takes occupancy()%indx = -kcell.  This is a necrotic volume.
! Such necrotic sites migrate towards the blob centre.
!-----------------------------------------------------------------------------------------
subroutine cell_death(kcell)
integer :: kcell
integer :: site(3)

write(*,*) 'cell_death: ',kcell
cell_list(kcell)%state = DEAD
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
	if (dmin >= d1) exit
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
write(*,*) 'necrotic_migration: site2: ',site2
end subroutine

end module
