module boundary

use global
use bdry_linked_list
implicit none
save

integer, parameter :: MAX_BDRY = 20000

contains

!----------------------------------------------------------------------------------------
! Check to see if (x,y,z) is outside the grid
!----------------------------------------------------------------------------------------
logical function outside_xyz(x,y,z)
integer :: x, y, z
outside_xyz = .true.
if (x < 1 .or. x > NX) return
if (y < 1 .or. y > NY) return
if (z < 1 .or. z > NZ) return
outside_xyz = .false.
end function

!----------------------------------------------------------------------------------------
! A boundary site is inside and has at least one Neumann neighbour outside the blob.
!----------------------------------------------------------------------------------------
logical function isbdry(site)
integer :: site(3)
integer :: x, y, z
integer :: k, xx, yy, zz

x = site(1)
y = site(2)
z = site(3)
!if (occupancy(x,y,z)%indx(1) == OUTSIDE_TAG) then	! outside
if (occupancy(x,y,z)%indx(1) <= 0) then	! outside or unreachable or vacant
    isbdry = .false.
    return
endif
do k = 1,6
	xx = x + neumann(1,k)
	yy = y + neumann(2,k)
	zz = z + neumann(3,k)
	if (outside_xyz(xx,yy,zz)) cycle
    if (occupancy(xx,yy,zz)%indx(1) == OUTSIDE_TAG) then
        isbdry = .true.
        return
    endif
enddo
isbdry = .false.
end function

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine CreateBdrylist
integer :: x, y, z
integer :: k, site(3)
type (boundary_type), pointer :: bdry

write(nflog,*) 'CreateBdryList'
call DestroyBdrylist
!nullify(bdrylist)
nbdry = 0
do x = 1,NX
    do y = 1,NY
        do z = 1,NZ
            site = (/x,y,z/)
            if (isbdry(site)) then
                allocate(bdry)
                bdry%site = site
!                bdry%chemo_influx = .false.
                nullify(bdry%next)
                call bdrylist_insert(bdry,bdrylist)
                occupancy(x,y,z)%bdry => bdry
            endif
        enddo
    enddo
enddo
bdry_changed = .false.
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine DestroyBdrylist
integer :: x, y, z
integer :: k, site(3)
type (boundary_type), pointer :: bdry, next

bdry => bdrylist
do while ( associated ( bdry )) 
    next => bdry%next
    site = bdry%site
    nullify(occupancy(site(1),site(2),site(3))%bdry)
	nullify(bdry)
	bdry => next
enddo
nullify(bdrylist)
end subroutine

!----------------------------------------------------------------------------------------
! As well as marking bdry sites, the code checks for vacated sites that have become
! OUTSIDE.  A brute-force iteration of the procedure is used because whenever a
! vacated site is converted to OUTSIDE there is a possibility of a change to the bdry.
! Such occurrences are expected to be rare, so the computational cost is not high.
! It is important to prevent a vacated site from persisting outside the blob.
!----------------------------------------------------------------------------------------
subroutine UpdateBdrylist
integer :: x, y, z
integer :: k, site(3), Nsites_old
type (boundary_type), pointer :: bdry
logical :: done

Nsites_old = Nsites
done = .false.
do while (.not.done)
	done = .true.
	Nsites = 0
	do x = 1,NX
		do y = 1,NY
			do z = 1,NZ
				site = (/x,y,z/)
				if (occupancy(site(1),site(2),site(3))%indx(1) /= OUTSIDE_TAG) Nsites = Nsites + 1
				if (associated(occupancy(site(1),site(2),site(3))%bdry)) cycle
				if (isbdry(site)) then
					if (occupancy(site(1),site(2),site(3))%indx(1) == 0) then	! vacant site is actually outside
						write(nflog,*) 'UpdateBdryList: Vacant site was actually outside - corrected: ',site
						occupancy(site(1),site(2),site(3))%indx(1) = OUTSIDE_TAG
						done = .false.
						cycle
					endif
					allocate(bdry)
					bdry%site = site
					nullify(bdry%next)
					call bdrylist_insert(bdry,bdrylist)
					occupancy(x,y,z)%bdry => bdry
				endif
			enddo
		enddo
	enddo
enddo
bdry_changed = .false.
!if (Nsites /= Nsites_old) then
!	write(nflog,*) 'UpdateBdryList: Nsites changed: ',Nsites_old,Nsites
!endif
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine CreateChecklist
integer :: x, y, z, ncheck
integer :: k, site(3)
type (boundary_type), pointer :: bdry

write(logmsg,*) 'CreateCheckList'
call logger(logmsg)
nullify(checklist)
ncheck = 0
do x = 1,NX
    do y = 1,NY
        do z = 1,NZ
            site = (/x,y,z/)
            if (isbdry(site)) then
                ncheck = ncheck + 1
                allocate(bdry)
                bdry%site = site
!                bdry%chemo_influx = .false.
                nullify(bdry%next)
                call bdrylist_insert(bdry,checklist)
            endif
        enddo
    enddo
enddo
write(logmsg,*) 'ncheck = ',ncheck
call logger(logmsg)
end subroutine


!----------------------------------------------------------------------------------------
! Check all sites in bdrylist to ensure that they are still on the bdry, and remove any 
! that are not from bdrylist.
! When a site is removed from the list its chemokine concs and gradients are replaced by
! neighbourhood averages.
!----------------------------------------------------------------------------------------
subroutine FixBdrylist
integer :: site(3), sitelist(1000,3), k, n
type (boundary_type), pointer :: bdry

n = 0
bdry => bdrylist
do while ( associated ( bdry )) 
    site = bdry%site
    if (.not.isbdry(site)) then
        n = n+1
        sitelist(n,:) = site
    endif
    bdry => bdry%next
enddo
do k = 1,n
    site = sitelist(k,:)
    call bdrylist_delete(site,bdrylist)
    nullify(occupancy(site(1),site(2),site(3))%bdry)
    if (bdrylist_present(site, bdrylist)) then
        write(logmsg,*) 'Error: FixBdrylist: still in bdrylist: ',site
		call logger(logmsg)
        stop
    endif
!    call SetConcs(site)
enddo
end subroutine


!-----------------------------------------------------------------------------------------
! Choose an outside site or vacant site adjacent to site1
!-----------------------------------------------------------------------------------------
subroutine getoutsidesite(site1,site2)
integer :: site1(3), site2(3)
integer :: j, jmin, site(3), kcell
real(REAL_KIND) :: r, rmin

rmin = 1.0e10
do j = 1,27
	if (j == 14) cycle
!	site = site1 + neumann(:,j)
	site = site1 + jumpvec(:,j)
!	if (occupancy(site(1),site(2),site(3))%indx(1) == OUTSIDE_TAG) then
!	if (occupancy(site(1),site(2),site(3))%indx(1) <= 0) then
	kcell = occupancy(site(1),site(2),site(3))%indx(1)
	if (kcell == 0 .or. kcell == OUTSIDE_TAG) then
		r = cdistance(site)
		if (r < rmin) then
			rmin = r
			jmin = j
		endif
	endif
enddo
!site2 = site1 + neumann(:,jmin)
site2 = site1 + jumpvec(:,jmin)
end subroutine

!-----------------------------------------------------------------------------------------
! Move any cells from site to allow the site to be released (made 'outside')
!-----------------------------------------------------------------------------------------
subroutine ClearSite(csite)
integer :: csite(3)
integer :: i, k, cindx(2), kcell, site(3), indx(2), r
logical :: done

cindx = occupancy(csite(1),csite(2),csite(3))%indx
do i = 2,1,-1
    kcell = cindx(i)
    if (kcell == 0) cycle
    r = 0
    done = .false.
    do while (.not.done)
        r = r + 1
        do k = 1,27
            if (k == 14) cycle
	        site = csite + r*jumpvec(:,k)
	        if (outside_xyz(site(1),site(2),site(3))) cycle
	        indx = occupancy(site(1),site(2),site(3))%indx
            if (indx(1) <= OUTSIDE_TAG) cycle  ! outside or unreachable
            if (indx(1) == 0) then  ! use this slot
                occupancy(site(1),site(2),site(3))%indx(1) = kcell
                cell_list(kcell)%site = site
                done = .true.
                exit
            endif
        enddo
    enddo
enddo           
end subroutine

!-----------------------------------------------------------------------------------------
! Check consistency between bdrylist and occupancy%bdry
!-----------------------------------------------------------------------------------------
subroutine CheckBdryList(infomsg,ok)
character*(*) :: infomsg
logical :: ok
type (boundary_type), pointer :: bdry => null(), ocbdry => null()
integer :: x, y, z, site(3), dy, nb1, nb2, nbx, indx(2)

write(logmsg,'(a,a,i6)') 'CheckBdryList: ',infomsg,istep
call logger(logmsg)
ok = .true.
nb1 = 0
nbx = 0
bdry => bdrylist
do while ( associated ( bdry )) 
	nb1 = nb1 + 1
!	if (bdry%exit_ok) nbx = nbx + 1
    site = bdry%site
    ocbdry => occupancy(site(1),site(2),site(3))%bdry
    if (.not.associated(ocbdry,bdry)) then
		write(logmsg,'(a,a,a,3i5)') 'Error: CheckBdryList: called from: ',infomsg, ' inconsistent bdry pointers at site: ',site
		call logger(logmsg)
		ok = .false.
		return
!		stop
	endif
	indx = occupancy(site(1),site(2),site(3))%indx
	if (indx(1) <= OUTSIDE_TAG) then
		write(logmsg,'(a,a,a,3i5)') 'Error: CheckBdryList: called from: ',infomsg, ' site is outside: ',site
		call logger(logmsg)
		ok = .false.
		return
!		stop
	endif
	if (.not.isbdry(site)) then
		write(logmsg,'(a,a,a,3i5)') 'Error: CheckBdryList: called from: ',infomsg, ' site in bdrylist is NOT bdry: ',site
		call logger(logmsg)
		ok = .false.
		return
!		stop
	endif
		
    bdry => bdry%next
enddo
nb2 = 0
do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
			bdry => occupancy(x,y,z)%bdry
			if (associated(bdry)) then
				nb2 = nb2 + 1
				dy = y - blob_centre(2)
			else
				site = (/x,y,z/)
				if (isbdry(site)) then
					write(logmsg,*) 'Adding boundary site to bdrylist: ',x,y,z
					call logger(logmsg)
					! add it
					allocate(bdry)
					bdry%site = site
					nullify(bdry%next)
					call bdrylist_insert(bdry,bdrylist)
					occupancy(site(1),site(2),site(3))%bdry => bdry
				endif
			endif
		enddo
	enddo				
enddo
if (nb1 /= nb2) then
	write(logmsg,'(a,a,2i5)') 'Error: inconsistent boundary site counts: called from: ',infomsg,nb1,nb2
	call logger(logmsg)
	ok = .false.
	return
!	stop
endif
write(logmsg,*) 'bdrylist and occupancy%bdry are consistent: ',nb1,nbx
call logger(logmsg)
end subroutine

end module