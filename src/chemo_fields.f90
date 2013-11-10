! Chemokine concentration fields

module fields
use global
use chemokine
use bdry_linked_list
use ode_diffuse
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
logical function isbdry(x,y,z)
integer :: x, y, z
integer :: k, xx, yy, zz

if (occupancy(x,y,z)%indx(1) < 0) then	! outside or DC
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
            if (isbdry(x,y,z)) then
                nbdry = nbdry + 1
                allocate(bdry)
                site = (/x,y,z/)
                bdry%site = site
                bdry%chemo_influx = .false.
                nullify(bdry%next)
                call bdrylist_insert(bdry,bdrylist)
                call AssignBdryRole(site,bdry)
                occupancy(x,y,z)%bdry => bdry
            endif
        enddo
    enddo
enddo
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
            if (isbdry(x,y,z)) then
                ncheck = ncheck + 1
                allocate(bdry)
                site = (/x,y,z/)
                bdry%site = site
                bdry%chemo_influx = .false.
                nullify(bdry%next)
                call bdrylist_insert(bdry,checklist)
                call AssignBdryRole(site,bdry)
            endif
        enddo
    enddo
enddo
write(logmsg,*) 'ncheck = ',ncheck, ' NBcells: ',NBcells
call logger(logmsg)
end subroutine

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
!----------------------------------------------------------------------------------------
subroutine AddSite(ok)
logical :: ok
type (boundary_type), pointer :: bdry
integer :: site(3), psite(3), indx(2), k, kmin, x, y, z
real :: r(3), x2, y2, z2, dp, t, tmax, dpq, dmax, dmin
real :: cmin, cmax

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
    t = sqrt(1/(x2/Radius%x**2 + y2/Radius%y**2 + z2/Radius%z**2))
    dpq = (t-1)*dp
    if (dpq > dmax) then
        psite = site
        dmax = dpq
        tmax = t
    endif
    bdry => bdry%next
enddo
if (.not.isbdry(psite(1),psite(2),psite(3)))then
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
if (isbdry(site(1),site(2),site(3))) then   ! add it to the bdrylist
    nbdry = nbdry + 1
    allocate(bdry)
    bdry%site = site
    bdry%chemo_influx = .false.
    nullify(bdry%next)
    call bdrylist_insert(bdry,bdrylist)
    call AssignBdryRole(site,bdry)
    occupancy(site(1),site(2),site(3))%bdry => bdry
    call SetBdryConcs(site)
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
    if (.not.isbdry(site(1),site(2),site(3))) then
        n = n+1
        sitelist(n,:) = site
    endif
    bdry => bdry%next
enddo
do k = 1,n
    site = sitelist(k,:)
    call bdrylist_delete(site,bdrylist)
    nullify(occupancy(site(1),site(2),site(3))%bdry)
    nbdry = nbdry - 1
    if (bdrylist_present(site, bdrylist)) then
        write(logmsg,*) 'Error: FixBdrylist: still in bdrylist: ',site
		call logger(logmsg)
        stop
    endif
    call SetConcs(site)
enddo
end subroutine

!----------------------------------------------------------------------------------------
! The logic is similar to AddSite.  We look for the bdry site P which is most outside
! the ellipsoid, i.e. to minimise OQ-OP, where Q is the point at
! which the line through OP intersects the ellipsoid surface.  We then choose the 'inside'
! neighbour of P that is furthest from O.
!----------------------------------------------------------------------------------------
subroutine RemoveSite(ok)
logical :: ok
type (boundary_type), pointer :: bdry
integer :: site(3), psite(3), k, kmax
real :: r(3), x2, y2, z2, dp, t, dpq, dmax, dmin

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
    t = sqrt(1/(x2/Radius%x**2 + y2/Radius%y**2 + z2/Radius%z**2))
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
if (.not.isbdry(psite(1),psite(2),psite(3))) then
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
    if (occupancy(site(1),site(2),site(3))%indx(1) < 0) cycle   ! outside or DC
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
    nbdry = nbdry - 1
! Need to check for a new bdry site to replace the one removed
! Look at all the neighbours
    psite = site
    do k = 1,27
	    site = psite + jumpvec(:,k)
        if (occupancy(site(1),site(2),site(3))%indx(1) < 0) cycle   ! outside or DC
        if (associated(occupancy(site(1),site(2),site(3))%bdry)) cycle   ! bdry
        if (isbdry(site(1),site(2),site(3))) then
            nbdry = nbdry + 1
            allocate(bdry)
            bdry%site = site
            bdry%chemo_influx = .false.
            nullify(bdry%next)
            call bdrylist_insert(bdry,bdrylist)
            call AssignBdryRole(site,bdry)
            occupancy(site(1),site(2),site(3))%bdry => bdry
            call SetBdryConcs(site)
        endif
    enddo
endif
ok = .true.
end subroutine

!--------------------------------------------------------------------------------
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
    write(*,*) 'nbdry = ',nbdry, ' NBcells: ',NBcells
    write(*,*) 'Removing sites'
    call RemoveSites(n,ok)
    write(*,*) 'nbdry = ',nbdry, ' NBcells: ',NBcells
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
	if (bdry%exit_ok) nbx = nbx + 1
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
	if (.not.isbdry(site(1),site(2),site(3))) then
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
				if (isbdry(x,y,z)) then
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
            if (indx(1) < 0) cycle  ! outside or DC
            if (indx(1) == 0) then  ! use this slot
                occupancy(site(1),site(2),site(3))%indx(1) = kcell
                Bcell_list(kcell)%site = site
                done = .true.
                exit
            elseif (indx(2) == 0) then  ! use this slot
                occupancy(site(1),site(2),site(3))%indx(2) = kcell
                Bcell_list(kcell)%site = site
                occupancy(csite(1),csite(2),csite(3))%indx(i) = 0
                done = .true.
                exit
            endif
        enddo
    enddo
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

end module