! To test solution of reaction-diffusion on a rectangular grid, using fgmres with ITSOL ILUK preconditioning
! This version uses two rectangular grids, coarse and fine.
! The fine grid is large enough to contain the spheroid blob completely, and has grid resolution
! such that a grid-cell can hold a few tumour cells, e.g. dxf = 20 -> 30 um
! The fine grid is embedded in a coarse grid, which has a resolution dxb that is a multiple of dxf,
! e.g. dxb = 4*dxf. Coarse grid points coincide with points in the fine grid, and in particular
! the points on the boundary of the fine grid coincides with a line of the course grid.
! The size of the coarse grid is set to match the medium volume, if not the shape.
! 
! The idea behind the solution method is that solving on the coarse grid provides the boundary
! values for solution on the fine grid.
!
! Quasi-steady-state method:
! -------------------------
! Time-dependence is handled in the coarse solution, using the IMEX 2-SBDF scheme.  The most recent 
! rates of uptake-secretion of constituents by cells are used to compute grid flux F().
! After the Cave values on the coarse grid have been determined, the boundary concentrations on
! the fine grid are estimated by interpolation, and the steady-state solution on the fine grid
! is computed.  It is this solution that provides the flux values for the next coarse grid time step.
! 
! Units:
!     time				s = seconds
!     distance			cm
!     volume			cm^3
!     mass				micromole = 10^-6 mol = mumol
!     flux				mumol/s
!     concentration		mumol/cm^3 = mM
!
! Parallel issues:
! https://software.intel.com/en-us/articles/threading-fortran-applications-for-parallel-performance-on-multi-core-systems/
!
module react_diff

use real_kind_mod
use global
use omp_lib
use sparse_map
use par_zig_mod
use chemokine
!use continuum

use, intrinsic :: iso_c_binding

implicit none

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine setup_react_diff
integer :: ix, iy, iz, ic, maxnz, ichemo
real(REAL_KIND) :: C0, total_flux
!character*(10) :: emapfile
!real(REAL_KIND), pointer :: Cprev(:,:,:), Fprev(:,:,:), Fcurr(:,:,:)
character*(10) :: bmapfile
real(REAL_KIND), pointer :: Cave_b(:,:,:), Cprev_b(:,:,:), Fprev_b(:,:,:), Fcurr_b(:,:,:)
logical :: zero
logical :: ok

write(nflog,*) 'setup_react_diff: NXB: ',NXB
!dxf = DELTA_X
!dx3 = dxf*dxf*dxf
!nrow = (NX-2)*(NY-2)*(NZ-1)		! embedded map, Dirichlet conditions on ix=1,NX, iy=1,NY, iz=NZ
!maxnz = MAX_CHEMO*nrow
!if (allocated(amap)) deallocate(amap)
!if (allocated(ja)) deallocate(ja)
!if (allocated(ia)) deallocate(ia)
!allocate(amap(maxnz,0:3))
!allocate(ja(maxnz))
!allocate(ia(nrow+1))
DXB = 120
DXB = 1.0e-4*DXB	! um -> cm
ixb0 = (1 + NXB)/2
iyb0 = (1 + NYB)/2
izb0 = 5
xb0 = (ixb0-1)*DXB
yb0 = (iyb0-1)*DXB 
zb0 = (izb0-1)*DXB
centre_b = [xb0, yb0, zb0]
dxb3 = dxb*dxb*dxb
nrow_b = NXB*NYB*NZB
maxnz = MAX_CHEMO*nrow_b
if (allocated(amap_b)) deallocate(amap_b)
if (allocated(ja_b)) deallocate(ja_b)
if (allocated(ia_b)) deallocate(ia_b)
allocate(amap_b(maxnz,0:3))
allocate(ja_b(maxnz))
allocate(ia_b(nrow_b+1))

!emapfile = ''
!write(emapfile,'(a,i2.0,a)') 'emap',NX,'.dat'
!write(nflog,*) 'setup_react_diff: ',emapfile
!call make_sparse_emap(emapfile,.true.)
!write(nflog,*) 'made emapfile: ',emapfile

bmapfile = ''
write(bmapfile,'(a,i2.0,a)') 'bmap',NXB,'.dat'
write(nflog,*) 'setup_react_diff: ',bmapfile
call make_sparse_map(bmapfile,.false.)
write(nflog,*) 'made bmapfile: ',bmapfile

call make_grid_flux_weights
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	if (ichemo > TRACER) then
		chemo(ichemo)%bdry_conc = 0
	endif
	C0 = chemo(ichemo)%bdry_conc
	Cprev_b => chemo(ichemo)%Cprev_b
	Fprev_b => chemo(ichemo)%Fprev_b
	Fcurr_b => chemo(ichemo)%Fcurr_b
	Cave_b => chemo(ichemo)%Cave_b
	Cave_b = C0
	Cprev_b = C0
	Fprev_b = 0
	call update_Cex_Cin_const(ichemo)		! initialise %dMdt with SS values
	call getF_const(ichemo,total_flux,zero)	! computes cell fluxes then estimates total flux
	Fprev_b(ixb0,iyb0,izb0) = total_flux
! Sum fine grid fluxes to initialise Fcurr_b, Fprev_b
!	Cprev => chemo(ichemo)%Cprev
!	Fprev => chemo(ichemo)%Fprev
	
!	Cprev = C0
!	Fprev = Cflux(:,:,:,ichemo)
!	call makeF_b(Fprev_b,DELTA_T,zero)
	Fcurr_b = Fprev_b
enddo
end subroutine

!-------------------------------------------------------------------------------------------
! The flux values on the coarse grid are derived from the cell flux values by
! appropriate summation.  It's important that the total flux is consistent.
!-------------------------------------------------------------------------------------------
subroutine makeF_b(F_b,dt,zero)
real(REAL_KIND) :: F_b(:,:,:), dt
logical :: zero
real(REAL_KIND) :: Fsum_b

zero = (Fsum_b == 0)	
end subroutine

!-------------------------------------------------------------------------------------------
! This is the version for the coarse grid.
! Use variable numbering (ix,iy) -> k = (ix-1)*NY + iy
! Since the equation is now dM/dt = ... need to derive different expression for Kr
! Need to add treatment of top boundary for O2
! Need to add decay!  For now add it explicitly in solver
!-------------------------------------------------------------------------------------------
subroutine make_csr_b(a_b, ichemo, dt, Cave_b, Cprev_b, Fcurr_b, Fprev_b, rhs, zero)
integer :: ichemo
real(REAL_KIND) :: dt, a_b(:), Cave_b(:,:,:), Cprev_b(:,:,:), Fcurr_b(:,:,:), Fprev_b(:,:,:), rhs(:)
logical :: zero
integer :: ixb, iyb, izb, k, i, krow, kcol, nc
real(REAL_KIND) :: Kdiff, Kr, Cbdry, Fsum
integer, parameter :: m = 3

zero = .true.
Kdiff = chemo(ichemo)%medium_diff_coef
Fsum = 0
krow = 0
do k = 1,nnz_b
	if (k == ia_b(krow+1)) krow = krow+1
	kcol = ja_b(k)
	if (amap_b(k,0) == 2*m) then
		Kr = dxb*dxb/Kdiff
		a_b(k) = 3*Kr/(2*dt) + 2*m	! ... note that Kdiff should depend on (ixb,iyb,izb), ultimately on # of cells
	else
		a_b(k) = amap_b(k,0)
	endif
	if (ichemo == OXYGEN) then
		if (amap_b(k,3) == NZB .and. kcol == krow-1) then		! check this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if (a_b(k) /= -2) then
				write(*,*) 'Error in OXYGEN bdry adjustment'
				stop
			endif
			a_b(k) = -1
		endif
	endif
enddo
!!$omp parallel do private(iyb, ixb, krow, Kr)
do izb = 1,NZB
	do iyb = 1,NYB
		do ixb = 1,NXB
			krow = (ixb-1)*NYB*NZB + (iyb-1)*NZB + izb
			Kr = dxb*dxb/Kdiff
			rhs(krow) = Kr*((-2*Fcurr_b(ixb,iyb,izb) + Fprev_b(ixb,iyb,izb))/dxb3 + (1./(2*dt))*(4*Cave_b(ixb,iyb,izb) - Cprev_b(ixb,iyb,izb)))
			if (rhs(krow) /= 0) zero = .false.
		enddo
	enddo
enddo
!!$omp end parallel do
if (ichemo == OXYGEN) then
	Cbdry = chemo(ichemo)%bdry_conc
	izb = NZB
	do ixb = 1,NXB
		do iyb = 1,NYB
			krow = (ixb-1)*NYB*NZB + (iyb-1)*NZB + izb
			rhs(krow) = rhs(krow) + Cbdry		
		enddo
	enddo
endif
end subroutine

!-------------------------------------------------------------------------------------------
! Compute SS cp%Cex(ichemo) and cp%conc(ichemo) from Cave_b(:,:,:)
!-------------------------------------------------------------------------------------------
subroutine update_Cex_Cin_const(ichemo)
integer :: ichemo
integer :: kcell
real(REAL_KIND) :: Kin, Kout, V
type(cell_type), pointer :: cp
real(REAL_KIND), pointer :: Cextra(:,:,:)
!write(*,*) 'update_Cex_Cin_const'

!Cextra => Caverage(:,:,:,ichemo)		! currently using the average concentration!
Cextra => chemo(ichemo)%Cave_b(:,:,:)
! This needs to use the coarse grid

Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	call extra_concs_const(kcell, Cextra, cp%Cex(ichemo))
	V = cp%volume*Vcell_cm3
	cp%conc(ichemo) = getCin_SS(kcell,ichemo,V,cp%Cex(ichemo))
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Estimates the extracellular concentrations at the location of the centre of a cell,
! from the field of Cextra at grid points
!-----------------------------------------------------------------------------------------
subroutine extra_concs_const(kcell,Cextra,conc)
integer :: kcell
real(REAL_KIND) :: Cextra(:,:,:), conc
real(REAL_KIND) :: cb(3), alfa(3)
integer :: i, ixb, iyb, izb, grid(3)
type(cell_type), pointer :: cp

cp => cell_list(kcell)
if (cp%state == DEAD) then
	write(*,*) 'Error: extra_concs_const: dead cell: ',kcell
	stop
endif
cb = Centre_b + (cp%site - Centre)*DELTA_X
ixb = cb(1)/DXB + 1
iyb = cb(2)/DXB + 1
izb = cb(3)/DXB + 1
grid = [ixb, iyb, izb]
do i = 1,3
	alfa(i) = (cb(i) - (grid(i)-1)*DXB)/DXB
enddo
conc = (1-alfa(1))*(1-alfa(2))*(1-alfa(3))*Cextra(ixb,iyb,izb)  &
        + (1-alfa(1))*alfa(2)*(1-alfa(3))*Cextra(ixb,iyb+1,izb)  &
        + (1-alfa(1))*alfa(2)*alfa(3)*Cextra(ixb,iyb+1,izb+1)  &
        + (1-alfa(1))*(1-alfa(2))*alfa(3)*Cextra(ixb,iyb,izb+1)  &
        + alfa(1)*(1-alfa(2))*(1-alfa(3))*Cextra(ixb+1,iyb,izb)  &
        + alfa(1)*alfa(2)*(1-alfa(3))*Cextra(ixb+1,iyb+1,izb)  &
        + alfa(1)*alfa(2)*alfa(3)*Cextra(ixb+1,iyb+1,izb+1)  &
        + alfa(1)*(1-alfa(2))*alfa(3)*Cextra(ixb+1,iyb,izb+1)
end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine update_Cex
integer :: kcell
real(REAL_KIND) :: cb(3)
type(cell_type), pointer :: cp

do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	cb = Centre_b + (cp%site - Centre)*DELTA_X
	call getConc(cb,cp%Cex)
enddo
end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine cell_location(cp,c)
type(cell_type), pointer :: cp
real(REAL_KIND) :: c(3)

c = 0

end subroutine

!-------------------------------------------------------------------------------------------
! Compute cp%Cex(ichemo) and cp%conc(ichemo) from Cave_b(:,:,:), 
!-------------------------------------------------------------------------------------------
subroutine update_Cex_Cin_dCdt_const(ichemo, dt)
integer :: ichemo
real(REAL_KIND) :: dt
integer :: kcell, ix, iy, iz
real(REAL_KIND) :: alfa(3), Clast, Kin, Kout, dC, dCexdt, dMdt, V
type(cell_type), pointer :: cp
real(REAL_KIND), pointer :: Cextra(:,:,:)
real(REAL_KIND), pointer :: Cprev(:,:,:)

!write(*,*) 'update_Cex_Cin_dCdt_const: ',ichemo

Cextra => chemo(ichemo)%Cave_b(:,:,:)		! currently using the average concentration!
Cprev => chemo(ichemo)%Cprev_b
Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out
!!$omp parallel do private(cp, ix, iy, iz, alfa, Clast, dCexdt)
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	call grid_interp(kcell, alfa)
	ix = cp%site(1)		! no longer valid, %site has a different meaning
	iy = cp%site(2)
	iz = cp%site(3)
	cp%Cex(ichemo) = (1-alfa(1))*(1-alfa(2))*(1-alfa(3))*Cextra(ix,iy,iz)  &
        + (1-alfa(1))*alfa(2)*(1-alfa(3))*Cextra(ix,iy+1,iz)  &
        + (1-alfa(1))*alfa(2)*alfa(3)*Cextra(ix,iy+1,iz+1)  &
        + (1-alfa(1))*(1-alfa(2))*alfa(3)*Cextra(ix,iy,iz+1)  &
        + alfa(1)*(1-alfa(2))*(1-alfa(3))*Cextra(ix+1,iy,iz)  &
        + alfa(1)*alfa(2)*(1-alfa(3))*Cextra(ix+1,iy+1,iz)  &
        + alfa(1)*alfa(2)*alfa(3)*Cextra(ix+1,iy+1,iz+1)  &
        + alfa(1)*(1-alfa(2))*alfa(3)*Cextra(ix+1,iy,iz+1)
	Clast = (1-alfa(1))*(1-alfa(2))*(1-alfa(3))*Cprev(ix,iy,iz)  &
        + (1-alfa(1))*alfa(2)*(1-alfa(3))*Cprev(ix,iy+1,iz)  &
        + (1-alfa(1))*alfa(2)*alfa(3)*Cprev(ix,iy+1,iz+1)  &
        + (1-alfa(1))*(1-alfa(2))*alfa(3)*Cprev(ix,iy,iz+1)  &
        + alfa(1)*(1-alfa(2))*(1-alfa(3))*Cprev(ix+1,iy,iz)  &
        + alfa(1)*alfa(2)*(1-alfa(3))*Cprev(ix+1,iy+1,iz)  &
        + alfa(1)*alfa(2)*alfa(3)*Cprev(ix+1,iy+1,iz+1)  &
        + alfa(1)*(1-alfa(2))*alfa(3)*Cprev(ix+1,iy,iz+1)
    dCexdt = (cp%Cex(ichemo) - Clast)/dt
	V = cp%volume*Vcell_cm3
	cp%conc(ichemo) = getCin(kcell,ichemo,V,cp%Cex(ichemo),dCexdt)
enddo
!!$omp end parallel do
end subroutine

!-----------------------------------------------------------------------------------------
! Now using the coarse grid
!-----------------------------------------------------------------------------------------
subroutine grid_interp(kcell,alfa)
integer :: kcell
real(REAL_KIND) :: alfa(3)
integer :: i, ix, iy, iz
real(REAL_KIND) :: centre(3)
type(cell_type), pointer :: cp

cp => cell_list(kcell)
if (cp%state == DEAD) then
	write(*,*) 'Error: grid_interp: dead cell: ',kcell
	stop
endif
!if (cp%nspheres == 1) then
!	centre = cp%centre(:,1)
!else
!	centre = 0.5*(cp%centre(:,1) + cp%centre(:,2))
!endif
call cell_location(cp,centre)
ix = cp%site(1)
iy = cp%site(2)
iz = cp%site(3)
do i = 1,3
	alfa(i) = (centre(i) - (cp%site(i)-1)*DXB)/DXB		! no longer valid, %site has a different meaning
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Flux contributions from a cell are accumulated at coarse grid points given by cnr(:,:)
! This assumes that the cell flux rates have already been computed.
! Flux units: mumol/s
!-----------------------------------------------------------------------------------------
subroutine grid_flux_weights(kcell,cnr,wt)
integer :: kcell, cnr(3,8)
real(REAL_KIND) :: wt(8)
integer :: k
real(REAL_KIND) :: centre(3), gridpt(3), r(3), d(8), sum
type(cell_type), pointer :: cp

cp => cell_list(kcell)
if (cp%state == DEAD) then
	write(*,*) 'Error: grid_flux_weights: dead cell: ',kcell
	stop
endif
!if (cp%nspheres == 1) then
!	centre = cp%centre(:,1)
!else
!	centre = 0.5*(cp%centre(:,1) + cp%centre(:,2))
!endif
call cell_location(cp,centre)
sum = 0
do k = 1,8
	gridpt(:) = (cnr(:,k)-1)*DXB
	r = centre - gridpt
	d(k) = max(sqrt(dot_product(r,r)), small_d)
	sum = sum + 1/d(k)
enddo
! The grid flux weights are (1/d(k))/sum.  Note that dMdt > 0 for +ve flux into the cell, 
do k = 1,8
!	Cflux(cnr(1,k),cnr(2,k),cnr(3,k),:) = Cflux(cnr(1,k),cnr(2,k),cnr(3,k),:) + cp%dMdt(:)*(1/d(k))/sum
	wt(k) = (1/d(k))/sum
enddo
end subroutine

!-------------------------------------------------------------------------------------------
! Estimate total flux values associated with each coarse grid pt, from the cell fluxes.
! This is the contribution to the field, i.e. cell uptake is negative flux ????
!-------------------------------------------------------------------------------------------
subroutine getF_const(ichemo, total_flux, zero)
integer :: ichemo
real(REAL_KIND) :: total_flux
logical :: zero
real(REAL_KIND) :: Kin, Kout
integer :: kcell
type(cell_type), pointer :: cp

!write(*,*) 'getF_const: ',ichemo,nlist

! Update Cex for each cell
Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out

! Compute cell fluxes cp%dMdt
total_flux = 0
!!$omp parallel do private(cp)
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	cp%dMdt(ichemo) = Kin*cp%Cex(ichemo) - Kout*cp%conc(ichemo)
	total_flux = total_flux + cp%dMdt(ichemo)
	if (kcell == 1) then
		write(nflog,'(a,i4,5e12.5)') 'getF_const: Cin, Cex, dMdt: ',ichemo, cp%conc(ichemo), cp%Cex(ichemo),cp%dMdt(ichemo),Kin,Kout
	endif
!	if (cp%dMdt(ichemo) < 0) then
!		write(*,'(a,2i6,3e12.3)') 'getF_const: dMdt< 0: ',ichemo,kcell,cp%conc(ichemo), cp%Cex(ichemo),cp%dMdt(ichemo)
!		stop
!	endif
enddo
!!$omp end parallel do
write(nflog,'(a,2i4,e12.3)') 'total_flux: ',istep,ichemo,total_flux
zero = (total_flux == 0)
!if (ichemo == OXYGEN) write(*,'(a,2e12.3)') 'Cex(O2), O2 flux: ',cell_list(1)%Cex(ichemo),cell_list(1)%dMdt(ichemo)
! Estimate grid pt flux values F
!call make_grid_flux(ichemo, Cflux_const, zero)
end subroutine

!-------------------------------------------------------------------------------------------
! alfa is the amount of the previous flux
!-------------------------------------------------------------------------------------------
subroutine make_grid_flux(ichemo,Cflux_const,zero)
integer :: ichemo
real(REAL_KIND) :: Cflux_const(:,:,:)
logical :: zero
real(REAL_KIND) :: Fsum_b
integer :: kcell, k, cnr(3,8)
type(cell_type), pointer :: cp
real(REAL_KIND) :: alfa = 0.7

Cflux_const = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	cnr = cp%cnr
	do k = 1,8
		Cflux_const(cnr(1,k),cnr(2,k),cnr(3,k)) = Cflux_const(cnr(1,k),cnr(2,k),cnr(3,k)) + cp%dMdt(ichemo)*cp%wt(k)
	enddo
enddo
Cflux_const(:,:,:) = (1-alfa)*Cflux_const(:,:,:) + alfa*chemo(ichemo)%Fprev_b(:,:,:)
Fsum_b = sum(Cflux_const)
zero = (Fsum_b == 0)
end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine make_grid_flux_weights
integer :: ic, ix, iy, iz, kcell
type(cell_type), pointer :: cp

do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	ix = cp%site(1)		! no longer valid, %site has a different meaning
	iy = cp%site(2)
	iz = cp%site(3)
	cp%cnr(:,1) = [ix, iy, iz]
	cp%cnr(:,2) = [ix, iy+1, iz]
	cp%cnr(:,3) = [ix, iy, iz+1]
	cp%cnr(:,4) = [ix, iy+1, iz+1]
	cp%cnr(:,5) = [ix+1, iy, iz]
	cp%cnr(:,6) = [ix+1, iy+1, iz]
	cp%cnr(:,7) = [ix+1, iy, iz+1]
	cp%cnr(:,8) = [ix+1, iy+1, iz+1]
	call grid_flux_weights(kcell, cp%cnr, cp%wt)
enddo
end subroutine

!-------------------------------------------------------------------------------------------
! This valid for steady-state with any constituent, with decay, when Cex is known
! In fact it is currently valid for oxygen and glucose.
! This is effectively the same as getCin() with dCexdt = 0
!-------------------------------------------------------------------------------------------
function getCin_SS(kcell,ichemo, Vin, Cex) result(Cin)
integer :: kcell, ichemo
real(REAL_KIND) :: Vin, Cex, Cin, Cex_t, Cin_t
real(REAL_KIND) :: Kin, Kout, Kd, Kmax, VKdecay, C0, a, b, c, D, r(3)
integer :: i, n, ictyp, idrug, im
real(REAL_KIND) :: CO2, C_parent, C_metab1
real(REAL_KIND) :: C2_0, C2_1, C2_2, KO2_0, KO2_1, KO2_2, Kmet0_0, Kmet0_1, Kmet0_2
real(REAL_KIND) :: K1, K2, Km, Vmax
type(drug_type), pointer :: dp

!write(*,*) 'getCin_SS: istep,ichemo,Vin,Cex: ',istep,ichemo,Vin,Cex
!if (Cex < 0) then
!	write(logmsg,*) 'getCin_SS: istep,ichemo,Vin,Cex: ',istep,ichemo,Vin,Cex
!	call logger(logmsg)
!	stop
!endif
if (Cex <= 0) then
	Cex = 0
	Cin = 0
	return
endif
ictyp = cell_list(kcell)%celltype
CO2 = cell_list(kcell)%conc(OXYGEN)
Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out
Kd = chemo(ichemo)%decay_rate
if (ichemo <= GLUCOSE) then
	Kmax = chemo(ichemo)%max_cell_rate
	VKdecay = Vin*Kd
	C0 = chemo(ichemo)%MM_C0
	if (chemo(ichemo)%Hill_N == 2) then
		b = C0*C0
		a = (Kmax + b*VKdecay - Kin*Cex)/(Kout+VKdecay)
		c = -b*Cex*Kin/(Kout + VKdecay)
		call cubic_roots(a,b,c,r,n)
		if (n == 1) then
			Cin = r(1)
		else
			n = 0
			do i = 1,3
				if (r(i) > 0) then
					n = n+1
					Cin = r(i)
				endif
			enddo
			if (n > 1) then
				write(*,*) 'getCin_SS: two roots > 0: ',r
				stop
			endif
		endif
	elseif (chemo(ichemo)%Hill_N == 1) then
		b = (Kmax + Kout*C0 + 2*VKdecay*C0 - Kin*Cex)/(Kout + VKdecay)
		c = -C0*Cex*Kin/(Kout + VKdecay)
		D = sqrt(b*b - 4*c)
		Cin = (D - b)/2
	endif
elseif (ichemo == TRACER) then
	
else	! parent drug or drug metabolite
	idrug = (ichemo - TRACER - 1)/3 + 1
	dp => drug(idrug)
	im = ichemo - TRACER - 1 - 3*(idrug-1)
	if (im == 0) then		! parent
		Kmet0_0 = dp%Kmet0(ictyp,0)
		C2_0 = dp%C2(ictyp,0)
		KO2_0 = dp%KO2(ictyp,0)
		Km = dp%Km(ictyp,0)
		Vmax = dp%Vmax(ictyp,0)
		K1 = (1 - C2_0 + C2_0*KO2_0/(KO2_0 + CO2))*Kmet0_0
		K2 = K1*Vmax/Kmet0_0
		if (K2 /= 0) then	!quadratic: a.x^2 + b.x + c = 0
			a = K1 + Kd + Kout/Vin
			b = a*Km + K2 - Kin*Cex/Vin
			c = -Kin*Cex*Km/Vin
			b = b/a
			c = c/a
			D = sqrt(b*b - 4*c)
			Cin = (D - b)/2
		else				! linear: a.x + b = 0
			a = K1 + Kd + Kout/Vin
			b = -Kin*Cex/Vin
			Cin = -b/a
		endif
		if (Cin > Cex) then
			write(*,'(a,2e12.3)') 'Cex,Cin: ',Cex,Cin
			stop
		endif
	elseif (im == 1) then	! metab1
		CO2 = cell_list(kcell)%conc(OXYGEN)
		C_parent = cell_list(kcell)%conc(ichemo-1)
		C2_0 = dp%C2(ictyp,0)
		KO2_0 = dp%KO2(ictyp,0)
		Kmet0_0 = dp%Kmet0(ictyp,0)
		C2_1 = dp%C2(ictyp,1)
		KO2_1 = dp%KO2(ictyp,1)
		Kmet0_1 = drug(idrug)%Kmet0(ictyp,1)
		Cin = ((1 - C2_0 + C2_0*KO2_0/(KO2_0 + CO2))*Kmet0_0*C_parent + Kin*Cex/Vin) &
		     /((1 - C2_1 + C2_1*KO2_1/(KO2_1 + CO2))*Kmet0_1 + Kd + Kout/Vin)
	elseif (im == 2) then	! metab2
		CO2 = cell_list(kcell)%conc(OXYGEN)
		C_metab1 = cell_list(kcell)%conc(ichemo-1)
		C2_1 = dp%C2(ictyp,1)
		KO2_1 = dp%KO2(ictyp,1)
		Kmet0_1 = dp%Kmet0(ictyp,1)
		C2_2 = dp%C2(ictyp,2)
		KO2_2 = dp%KO2(ictyp,2)
		Kmet0_2 = dp%Kmet0(ictyp,2)
		Cin = ((1 - C2_1 + C2_1*KO2_1/(KO2_1 + CO2))*Kmet0_1*C_metab1 + Kin*Cex/Vin) &
		     /((1 - C2_2 + C2_2*KO2_2/(KO2_2 + CO2))*Kmet0_2 + Kd + Kout/Vin)
	endif
endif
end function

!-------------------------------------------------------------------------------------------
! Use an approximation to dCin/dt derived from dCex/dt to arrive at a better estimate of Cin and flux
!-------------------------------------------------------------------------------------------
function getCin(kcell, ichemo, Vin, Cex, dCexdt) result(Cin)
integer :: kcell, ichemo
real(REAL_KIND) :: Vin, Cex, dCexdt, Cin
real(REAL_KIND) :: Kin, Kout, Kd, Kmax, VKdecay, dCdt, delta, C0, a, b, c, D, r(3)
integer :: i, n, ictyp, idrug, im
real(REAL_KIND) :: CO2, C_parent, C_metab1
real(REAL_KIND) :: C2_0, C2_1, C2_2, KO2_0, KO2_1, KO2_2, Kmet0_0, Kmet0_1, Kmet0_2
real(REAL_KIND) :: Km, Vmax, K1_0, K2_0, K1_1, K2_1, K1_2, K2_2
type(drug_type), pointer :: dp
type(cell_type), pointer :: cp


if (Cex < 0) then
	Cex = 0
	Cin = 0
	return
endif
cp => cell_list(kcell)
ictyp = cp%celltype
CO2 = cp%conc(OXYGEN)
Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out
Kd = chemo(ichemo)%decay_rate
if (ichemo <= GLUCOSE) then
	Kmax = chemo(ichemo)%max_cell_rate
	VKdecay = Vin*Kd
	C0 = chemo(ichemo)%MM_C0
	if (chemo(ichemo)%Hill_N == 2) then
		dCdt = dCexdt*(Kin/Vin)/(Kout/Vin + Kd + 2*Kmax*cp%conc(ichemo)*C0**2/(C0**2+cp%conc(ichemo)**2))
		delta = Vin*dCdt
		b = C0*C0
		a = (Kmax + b*VKdecay - (Kin*Cex-delta))/(Kout+VKdecay)
		c = -b*(Kin*Cex-delta)/(Kout + VKdecay)
		call cubic_roots(a,b,c,r,n)
		if (n == 1) then
			Cin = r(1)
		else
			n = 0
			do i = 1,3
				if (r(i) > 0) then
					n = n+1
					Cin = r(i)
				endif
			enddo
			if (n > 1) then
				write(*,*) 'getCin_SS: two roots > 0: ',r
				stop
			endif
		endif
	elseif (chemo(ichemo)%Hill_N == 1) then
		dCdt = dCexdt*(Kin/Vin)/(Kout/Vin + Kd + Kmax*C0/(C0 + cp%conc(ichemo)))
		delta = Vin*dCdt
		b = (Kmax + Kout*C0 + 2*VKdecay*C0 - (Kin*Cex-delta))/(Kout + VKdecay)
		c = -C0*(Kin*Cex-delta)/(Kout + VKdecay)
		D = sqrt(b*b - 4*c)
		Cin = (D - b)/2
	endif
elseif (ichemo == TRACER) then
	stop
else	! parent drug or drug metabolite
	idrug = (ichemo - TRACER - 1)/3 + 1
	dp => drug(idrug)
	im = ichemo - TRACER - 1 - 3*(idrug-1)
	if (im == 0) then		! parent
		Kmet0_0 = dp%Kmet0(ictyp,0)
		C2_0 = dp%C2(ictyp,0)
		KO2_0 = dp%KO2(ictyp,0)
		Km = dp%Km(ictyp,0)
		Vmax = dp%Vmax(ictyp,0)
		K1_0 = (1 - C2_0 + C2_0*KO2_0/(KO2_0 + CO2))*Kmet0_0
		K2_0 = K1_0*Vmax/Kmet0_0
		dCdt = (Kin*dCexdt/Vin)/(Kout/Vin + Kd + K1_0 + K2_0*Km/(Km + CO2)**2)
		if (K2_0 /= 0) then	!quadratic: a.x^2 + b.x + c = 0
			a = K1_0 + Kd + Kout/Vin
			b = a*Km + K2_0 - (Kin*Cex/Vin - dCdt)
			c = -Km*(Kin*Cex/Vin - dCdt)
			b = b/a
			c = c/a
			D = sqrt(b*b - 4*c)
			Cin = (D - b)/2
		else				! linear: a.x + b = 0
			a = K1_0 + Kd + Kout/Vin
			b = -(Kin*Cex/Vin - dCdt)
			Cin = -b/a
		endif
		Cin = max(Cin,0.0)
	elseif (im == 1) then	! metab1
		C_parent = cp%conc(ichemo-1)
		C2_0 = dp%C2(ictyp,0)
		KO2_0 = dp%KO2(ictyp,0)
		Kmet0_0 = dp%Kmet0(ictyp,0)
		C2_1 = dp%C2(ictyp,1)
		KO2_1 = dp%KO2(ictyp,1)
		Kmet0_1 = drug(idrug)%Kmet0(ictyp,1)
		K1_0 = (1 - C2_0 + C2_0*KO2_0/(KO2_0 + CO2))*Kmet0_0
		K1_1 = (1 - C2_1 + C2_1*KO2_1/(KO2_1 + CO2))*Kmet0_1
		dCdt = (Kin*dCexdt/Vin + K1_0*cp%dCdt(ichemo-1))/(Kout/Vin + Kd + K1_1)
		Cin = (K1_0*C_parent + (Kin*Cex/Vin - dCdt))/(K1_1 + Kd + Kout/Vin)
		Cin = max(Cin,0.0)
	elseif (im == 2) then	! metab2
		C_metab1 = cp%conc(ichemo-1)
		C2_1 = dp%C2(ictyp,1)
		KO2_1 = dp%KO2(ictyp,1)
		Kmet0_1 = dp%Kmet0(ictyp,1)
		C2_2 = dp%C2(ictyp,2)
		KO2_2 = dp%KO2(ictyp,2)
		Kmet0_2 = dp%Kmet0(ictyp,2)
		K1_1 = (1 - C2_1 + C2_1*KO2_1/(KO2_1 + CO2))*Kmet0_1
		K1_2 = (1 - C2_2 + C2_2*KO2_2/(KO2_2 + CO2))*Kmet0_2
		dCdt = (Kin*dCexdt/Vin + K1_1*cp%dCdt(ichemo-1))/(Kout/Vin + Kd + K1_2)
		Cin = (K1_1*C_metab1 + (Kin*Cex/Vin - dCdt))/(K1_2 + Kd + Kout/Vin)
		Cin = max(Cin,0.0)
	endif
endif
end function

!--------------------------------------------------------------------------------------
! Determine real roots r(:) of the cubic equation:
! x^3 + a.x^2 + b.x + c = 0
! If there is one real root, n=1 and the root is r(1)
! If there are three distinct real roots, n=3 and the roots are r(1), r(2), r(3)
! If there is a repeated root, n=2 and the single root is r(1), the repeated root is r(2)
!--------------------------------------------------------------------------------------
!subroutine cubic_roots(a, b, c, r, n)
!real(REAL_KIND) :: a, b, c, r(3)
!integer :: n
!real(REAL_KIND) :: QQ, RR, theta, R2, Q3, AA, BB
!
!QQ = (a*a - 3*b)/9
!RR = (2*a*a*a - 9*a*b + 27*c)/54
!Q3 = QQ*QQ*QQ
!R2 = RR*RR
!if (R2 < Q3) then
!	n = 3
!	theta = acos(RR/sqrt(Q3))
!	r(1) = -2*sqrt(QQ)*cos(theta/3) - a/3
!	r(2) = -2*sqrt(QQ)*cos((theta+2*PI)/3) - a/3
!	r(3) = -2*sqrt(QQ)*cos((theta-2*PI)/3) - a/3
!else
!	n = 1
!	AA = -sign(1.d0,RR)*(abs(RR) + sqrt(R2 - Q3))**(1.d0/3.d0)
!	if (AA == 0) then
!		BB = 0
!	else
!		BB = QQ/AA
!	endif
!	r(1) = AA + BB - a/3
!endif
!end subroutine

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!subroutine unpack_csr(a)
!real(REAL_KIND) :: a(:)
!integer :: k, irow, icol
!real(REAL_KIND) :: asum
!
!write(nflog,*) 'nrow: ',nrow
!allocate(afull(nrow,nrow))
!irow = 0
!do k = 1,nnz
!	if (k == ia(irow+1)) irow = irow+1
!	icol = ja(k)
!	afull(irow,icol) = a(k)
!enddo
!!do irow = 1,nrow
!!	asum = 0
!!	do icol = 1,nrow
!!		asum = asum + afull(irow,icol)*9
!!	enddo
!!!	write(*,'(i6,2f8.3)') irow,asum,rhs(irow)
!!enddo
!!write(*,*) 'rhs:'
!!write(*,'(25f7.2)') rhs(:)
!end subroutine

!--------------------------------------------------------------------------------------
! Use the solution for Cave_b to estimate the average blob boundary concentration: %medium_Cbnd
! by averaging over a sphere with radius = blob radius (cm) centred at the blob centre
! centre_b (unless this changes significantly)
! To generate a random point on the sphere, it is necessary only to generate two random numbers, z between -R and R, phi between 0 and 2 pi, each with a uniform distribution
! To find the latitude (theta) of this point, note that z=Rsin(theta), so theta=sin-1(z/R); its longitude is (surprise!) phi.
! In rectilinear coordinates, x=Rcos(theta)cos(phi), y=Rcos(theta)sin(phi), z=Rsin(theta)= (surprise!) z.
! Needless to say, (x,y,z) are not independent, but are rather constrained by x2+y2+z2=R2.
!--------------------------------------------------------------------------------------
subroutine set_bdry_conc
integer :: kpar = 0
real(REAL_KIND) :: rad, x, y, z, p(3), phi, theta, c(MAX_CHEMO), csum(MAX_CHEMO)
integer :: i, ic, ichemo, n = 100

call SetRadius(Nsites)
rad = Radius*DELTA_X
csum = 0
do i = 1,n
	z = -rad + 2*rad*par_uni(kpar)
	phi = 2*PI*par_uni(kpar)
	theta = asin(z/rad)
	x = xb0 + rad*cos(theta)*cos(phi)
	y = yb0 + rad*cos(theta)*sin(phi)
	z = zb0 + z
	p = [x, y, z]
	call getConc(p,c)
	csum = csum + c
enddo
do ic = 1,nchemo
	ichemo = chemomap(ic)
	chemo(ichemo)%medium_Cbnd = csum(ichemo)/n
!	write(*,'(a,i2,f8.4)') 'medium_Cbnd: ',ichemo,chemo(ichemo)%medium_Cbnd
enddo
end subroutine

!--------------------------------------------------------------------------------------
! Interpolate to obtain concentrations at p(:) = (x,y,z) from chemo(:)%Cave_b(:,:,:)
!--------------------------------------------------------------------------------------
subroutine getConc(cb,c)
real(REAL_KIND) :: cb(3), c(:)
integer :: ixb, iyb, izb, grid(3), i, ic, ichemo
real(REAL_KIND) :: alfa(3)
real(REAL_KIND), pointer :: Cextra(:,:,:)

ixb = cb(1)/DXB + 1
iyb = cb(2)/DXB + 1
izb = cb(3)/DXB + 1
grid = [ixb, iyb, izb]
do i = 1,3
	alfa(i) = (cb(i) - (grid(i)-1)*DXB)/DXB
enddo
c = 0
do ic = 1,nchemo
	ichemo = chemomap(ic)
	Cextra => chemo(ichemo)%Cave_b
	c(ichemo) = (1-alfa(1))*(1-alfa(2))*(1-alfa(3))*Cextra(ixb,iyb,izb)  &
			+ (1-alfa(1))*alfa(2)*(1-alfa(3))*Cextra(ixb,iyb+1,izb)  &
			+ (1-alfa(1))*alfa(2)*alfa(3)*Cextra(ixb,iyb+1,izb+1)  &
			+ (1-alfa(1))*(1-alfa(2))*alfa(3)*Cextra(ixb,iyb,izb+1)  &
			+ alfa(1)*(1-alfa(2))*(1-alfa(3))*Cextra(ixb+1,iyb,izb)  &
			+ alfa(1)*alfa(2)*(1-alfa(3))*Cextra(ixb+1,iyb+1,izb)  &
			+ alfa(1)*alfa(2)*alfa(3)*Cextra(ixb+1,iyb+1,izb+1)  &
			+ alfa(1)*(1-alfa(2))*alfa(3)*Cextra(ixb+1,iyb,izb+1)
enddo
end subroutine

!--------------------------------------------------------------------------------------
! See D:\ITSOL\tests\input
!
! Time for 20 iterations with NX=1000
! ILUtype=1 33 sec
! ILUtype=2 is hoplessly slow with the parameters specified
! ILUtype=3 56 sec
! ILUtype=4 197 sec
!
! Note: Solving for each constituent in parallel with OpenMP like this is possible
!       only if the constituent reactions are independent
! For drugs, the three related constituents are solved for sequentially in the fine grid.
! First drug, then metab1 (using drug results), then metab2 (using metab1 results).
!-------------------------------------------------------------------------------------- 
subroutine diff_solver(dt, framp)
real(REAL_KIND) :: dt, framp
integer :: i, k, k1, ix, iy, iz, irow, icol, kc, ic, icc, it
integer :: ixb, iyb, izb
integer :: ichemo, ierr, nfill, iters, maxits, im_krylov
real(REAL_KIND) :: R, tol, tol_b, asum, t, Vex_curr, Vex_next, Vin_curr, Vin_next, fdecay, Csum, dCsum, msum, total_flux
real(REAL_KIND), allocatable :: x(:), rhs(:)
!real(REAL_KIND), pointer :: Cave(:,:,:), Cprev(:,:,:), Fprev(:,:,:), Fcurr(:,:,:)
real(REAL_KIND), pointer :: Cave_b(:,:,:), Cprev_b(:,:,:), Fprev_b(:,:,:), Fcurr_b(:,:,:)
real(REAL_KIND), allocatable :: a(:), a_b(:)
real(REAL_KIND) :: alfa = 0.7
real(REAL_KIND) :: dCtol = 1.0e-4
integer :: ILUtype = 1
integer :: im, im1, im2, ichemof
logical :: zeroC(MAX_CHEMO)
logical :: done
logical :: use_const = .true.

nfill = 1	! Level of fill for ILUK preconditioner
tol = 1.0d-6
tol_b = 1.0d-6
im_krylov = 60	! dimension of Krylov subspace in (outer) FGMRES
maxits = 100

! Compute all steady-state grid point fluxes in advance from Cextra(:,:,:,:): Cflux(:,:,:,:)

!$omp parallel do private(Cave_b, Cprev_b, Fprev_b, Fcurr_b, a_b, x, rhs, ix, iy, iz, ixb, iyb, izb, it, done, ichemo, icc, k, dCsum, msum, total_flux, iters, ierr)
do ic = 1,nchemo
	ichemo = chemomap(ic)
	if (chemo(ichemo)%constant) cycle
!	write(*,'(a,i2)') 'coarse grid: ichemo: ',ichemo
!	write(nflog,'(a,i2)') 'coarse grid: ichemo: ',ichemo
	ichemo_curr = ichemo
	icc = ichemo - 1
	allocate(rhs(nrow_b))
	allocate(x(nrow_b))
	allocate(a_b(MAX_CHEMO*nrow_b))
	Cave_b => chemo(ichemo)%Cave_b
	Cprev_b => chemo(ichemo)%Cprev_b
	Fprev_b => chemo(ichemo)%Fprev_b
	Fcurr_b => chemo(ichemo)%Fcurr_b
!	write(*,*) 'Cave_b:'
!	write(*,'(5e15.6)') Cave_b(NXB/2,NYB/2,:)
		
	Fprev_b = Fcurr_b
	call getF_const(ichemo,total_flux,zeroC(ichemo))
	Fcurr_b(ixb0,iyb0,izb0) = total_flux*framp
	call make_csr_b(a_b, ichemo, dt, Cave_b, Cprev_b, Fcurr_b, Fprev_b, rhs, zeroC(ichemo))		! coarse grid

	! Solve Cave_b(t+dt) on coarse grid
	!----------------------------------
	call itsol_create_matrix(icc,nrow_b,nnz_b,a_b,ja_b,ia_b,ierr)
	!write(nflog,*) 'itsol_create_matrix: ierr: ',ierr
		
	if (ILUtype == 1) then
		call itsol_create_precond_ILUK(icc,nfill,ierr)
	!	write(nflog,*) 'itsol_create_precond_ILUK: ierr: ',ierr 
	elseif (ILUtype == 2) then
		call itsol_create_precond_VBILUK(icc,nfill,ierr)
	!	write(*,*) 'itsol_create_precond_VBILUK: ierr: ',ierr 
	elseif (ILUtype == 3) then
		call itsol_create_precond_ILUT(icc,nfill,tol_b,ierr)
	!	write(*,*) 'itsol_create_precond_ILUT: ierr: ',ierr 
	elseif (ILUtype == 4) then
		call itsol_create_precond_ARMS(icc,nfill,tol_b,ierr)
	!	write(*,*) 'itsol_create_precond_ARMS: ierr: ',ierr 
	endif

	do izb = 1,NZB
		do iyb = 1,NYB
			do ixb = 1,NXB
				k = (ixb-1)*NYB*NZB + (iyb-1)*NZB + izb
				x(k) = Cave_b(ixb,iyb,izb)		! initial guess
			enddo
		enddo
	enddo
	if (.not.zeroC(ichemo)) then
	!	write(nflog,*) 'call itsol_solve_fgmr_ILU'
		call itsol_solve_fgmr_ILU(icc, rhs, x, im_krylov, maxits, tol_b, iters, ierr)
	!	write(nflog,*) 'itsol_solve_fgmr_ILU: Cave_b: ierr, iters: ',ierr,iters
	else
		write(nflog,*) 'no solve, zeroC: ',ichemo
	endif
	call itsol_free_precond_ILU(icc, ierr)
!	write(nflog,*) 'did itsol_free_precond_ILU'
	call itsol_free_matrix(icc, ierr)
!	write(nflog,*) 'did itsol_free_matrix'

	Cprev_b = Cave_b
	fdecay = 1 - chemo(ichemo)%decay_rate*dt
	msum = 0
	do izb = 1,NZB
		do iyb = 1,NYB
			do ixb = 1,NXB
				k = (ixb-1)*NYB*NZB + (iyb-1)*NZB + izb
				Cave_b(ixb,iyb,izb) = fdecay*x(k)
				msum = msum + x(k)*dxb3		! this sums the mass of constituent in mumols
				if (x(k) < 0) then
					write(nflog,*) 'Cave_b < 0: ',ixb,iyb,izb,x(k)
				endif 
			enddo
		enddo
	enddo
	! interpolate Cave_b on blob boundary
	deallocate(a_b, x, rhs)
enddo
!$omp end parallel do

end subroutine

end module

