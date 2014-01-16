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
if (occupancy(x,y,z)%indx(1) <= 0) then	! outside or vacant
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

nullify(bdrylist)
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
! As well as marking bdry sites, the code checks for vacated sites that have become
! OUTSIDE.  A brute-force iteration of the procedure is used because whenever a
! vacated site is converted to OUTSIDE there is a possibility of a change to the bdry.
! Such occurrences are expected to be rare, so the computational cost is not high.
! It is important to prevent a vacated site from persisting outside the blob.
!----------------------------------------------------------------------------------------
subroutine UpdateBdrylist
integer :: x, y, z
integer :: k, site(3)
type (boundary_type), pointer :: bdry
logical :: done

done = .false.
do while (.not.done)
	done = .true.
	do x = 1,NX
		do y = 1,NY
			do z = 1,NZ
				site = (/x,y,z/)
				if (associated(occupancy(site(1),site(2),site(3))%bdry)) cycle
				if (isbdry(site)) then
					if (occupancy(site(1),site(2),site(3))%indx(1) == 0) then
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
    if (occupancy(site(1),site(2),site(3))%indx(1) == OUTSIDE_TAG) cycle   ! outside
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
        if (occupancy(site(1),site(2),site(3))%indx(1) == OUTSIDE_TAG) cycle   ! outside
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

!-----------------------------------------------------------------------------------------
! Choose an outside site or vacant site adjacent to site1
!-----------------------------------------------------------------------------------------
subroutine getoutsidesite(site1,site2)
integer :: site1(3), site2(3)
integer :: j, jmin, site(3)
real(REAL_KIND) :: r, rmin

rmin = 1.0e10
do j = 1,27
	if (j == 14) cycle
!	site = site1 + neumann(:,j)
	site = site1 + jumpvec(:,j)
!	if (occupancy(site(1),site(2),site(3))%indx(1) == OUTSIDE_TAG) then
	if (occupancy(site(1),site(2),site(3))%indx(1) <= 0) then
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
            if (indx(1) == OUTSIDE_TAG) cycle  ! outside
            if (indx(1) == 0) then  ! use this slot
                occupancy(site(1),site(2),site(3))%indx(1) = kcell
                cell_list(kcell)%site = site
                done = .true.
                exit
            elseif (indx(2) == 0) then  ! use this slot
                occupancy(site(1),site(2),site(3))%indx(2) = kcell
                cell_list(kcell)%site = site
                occupancy(csite(1),csite(2),csite(3))%indx(i) = 0
                done = .true.
                exit
            endif
        enddo
    enddo
enddo           
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
! Check consistency between bdrylist and occupancy%bdry
!-----------------------------------------------------------------------------------------
subroutine CheckBdryList
type (boundary_type), pointer :: bdry => null(), ocbdry => null()
integer :: x, y, z, site(3), dy, nb1, nb2, nbx, indx(2)

nb1 = 0
nbx = 0
bdry => bdrylist
do while ( associated ( bdry )) 
	nb1 = nb1 + 1
!	if (bdry%exit_ok) nbx = nbx + 1
    site = bdry%site
    ocbdry => occupancy(site(1),site(2),site(3))%bdry
    if (.not.associated(ocbdry,bdry)) then
		write(logmsg,*) 'Error: CheckBdryList: inconsistent bdry pointers at site: ',site
		call logger(logmsg)
		stop
	endif
	indx = occupancy(site(1),site(2),site(3))%indx
	if (indx(1) == OUTSIDE_TAG) then
		write(logmsg,*) 'Error: CheckBdryList: site is outside: ',site
		call logger(logmsg)
		stop
	endif
	if (.not.isbdry(site)) then
		write(logmsg,*) 'Error: CheckBdryList: bdry site is NOT bdry: ',site
		call logger(logmsg)
		stop
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
				dy = y - Centre(2)
			else
				site = (/x,y,z/)
				if (isbdry(site)) then
					write(logmsg,*) 'Error: boundary site not in list: ',x,y,z
					call logger(logmsg)
					stop
				endif
			endif
		enddo
	enddo				
enddo
if (nb1 /= nb2) then
	write(logmsg,*) 'Error: inconsistent boundary site counts: ',nb1,nb2
	call logger(logmsg)
	stop
endif
write(logmsg,*) 'bdrylist and occupancy%bdry are consistent: ',nb1,nbx
call logger(logmsg)
end subroutine

end module