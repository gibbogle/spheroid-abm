! Transferring results to the GUI

module transfer

use global
use chemokine
use envelope
use ode_diffuse
use, intrinsic :: iso_c_binding

#include "../src/version.h"

implicit none

integer, parameter :: X_AXIS = 1
integer, parameter :: Y_AXIS = 2
integer, parameter :: Z_AXIS = 3

type, bind(C) :: celldata_type
	integer(c_int) :: tag
	real(c_double) :: radius
	real(c_double) :: centre(3)
	integer(c_int) :: celltype
	integer(c_int) :: status
end type

type, bind(C) :: fielddata_type
    integer(c_int) :: NX, NY, NZ, NCONST
    real(c_double) :: DX
    type(c_ptr) :: Conc_ptr   ! Cslice(NX,NY,NZ,NCONST)
    integer(c_int) :: ncells
    type(c_ptr) :: cell_ptr
end type

!type(celldata_type) :: cdata(4000)

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_DLL_build_version(version_array,array_len) BIND(C) 
!DEC$ ATTRIBUTES DLLEXPORT :: get_dll_build_version
use, intrinsic :: iso_c_binding
character(c_char) :: version_array(12)
integer(c_int) :: array_len
integer :: k

dll_version = DLL_BUILD_VERSION
gui_version = GUI_BUILD_VERSION
!write(nflog,*) 'get_DLL_build_version: ',dll_version
do k = 1,12
	version_array(k) = dll_version(k:k)
!	write(nflog,'(i2,a,a)') k,' ',version_array(k)
	if (version_array(k) == ' ') then
		version_array(k) = char(0)
		array_len = k
		exit
	endif
enddo
!write(nflog,*) 'array_len: ',array_len
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_dimensions(NX_dim, NY_dim, NZ_dim, nsteps_dim, deltat, maxchemo, nextra, cused, dfraction, deltax) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dimensions
use, intrinsic :: iso_c_binding
integer(c_int) :: NX_dim,NY_dim,NZ_dim,nsteps_dim, maxchemo, nextra
real(c_double) :: deltat, dfraction, deltax
logical(c_bool) :: cused(*)
integer :: ichemo

write(nflog,*) 'get_dimensions: MAX_CHEMO: ',MAX_CHEMO
NX_dim = NX
NY_dim = NY
NZ_dim = NZ
nsteps_dim = nsteps
deltat = DELTA_T
deltax = DELTA_X
maxchemo = MAX_CHEMO
nextra = N_EXTRA
do ichemo = 1,MAX_CHEMO
	cused(ichemo+1) = chemo(ichemo)%used
enddo
cused(1) = .true.			! CFSE
cused(MAX_CHEMO+2) = .true.	! Growth rate
dfraction = 2*cell_radius/DELTA_X
end subroutine

!--------------------------------------------------------------------------------
! TC = tumour cell
!--------------------------------------------------------------------------------
subroutine get_scene(nTC_list,TC_list) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scene
use, intrinsic :: iso_c_binding
integer(c_int) :: nTC_list, TC_list(*)
integer :: k, kc, kcell, site(3), j, jb, colour, nlive
integer :: col(3)
integer :: x, y, z
real(REAL_KIND) :: vol, r
integer :: itcstate, ctype, stage, region, highlight
integer :: last_id1, last_id2
logical :: ok, highlighting
integer, parameter :: axis_centre = -2	! identifies the spheroid centre
integer, parameter :: axis_end    = -3	! identifies the spheroid extent in 5 directions
integer, parameter :: axis_bottom = -4	! identifies the spheroid extent in the -Y direction, i.e. bottom surface
integer, parameter :: ninfo = 7			! the size of the info package for a cell (number of integers)
integer, parameter :: nax = 6			! number of points used to delineate the spheroid

highlighting = (show_progeny /= 0)
nTC_list = 0

k = 0
last_id1 = 0
if (.false.) then
	! Need some markers to delineate the follicle extent.  These nax "cells" are used to convey (the follicle centre
	! and) the approximate ellipsoidal blob limits in the 3 axis directions.
	do k = 1,nax
		select case (k)
	!	case (1)
	!		x = Centre(1) + 0.5
	!		y = Centre(2) + 0.5
	!		z = Centre(3) + 0.5
	!		site = (/x, y, z/)
	!		itcstate = axis_centre
		case (1)
	!		x = Centre(1) - Radius%x - 2
			x = blobrange(1,1) - 1
			y = blob_centre(2) + 0.5
			z = blob_centre(3) + 0.5
			site = (/x, y, z/)
			itcstate = axis_end
		case (2)
	!		x = Centre(1) + Radius%x + 2
			x = blobrange(1,2) + 1
			y = blob_centre(2) + 0.5
			z = blob_centre(3) + 0.5
			site = (/x, y, z/)
			itcstate = axis_end
		case (3)
			x = blob_centre(1) + 0.5
	!		y = Centre(2) - Radius%y - 2
			y = blobrange(2,1) - 1
			z = blob_centre(3) + 0.5
			site = (/x, y, z/)
			itcstate = axis_bottom
		case (4)
			x = blob_centre(1) + 0.5
	!		y = Centre(2) + Radius%y + 2
			y = blobrange(2,2) + 1
			z = blob_centre(3) + 0.5
			site = (/x, y, z/)
			itcstate = axis_end
		case (5)
			x = blob_centre(1) + 0.5
			y = blob_centre(2) + 0.5
	!		z = Centre(3) - Radius%z - 2
			z = blobrange(3,1) - 1
			site = (/x, y, z/)
			itcstate = axis_end
		case (6)
			x = blob_centre(1) + 0.5
			y = blob_centre(2) + 0.5
	!		z = Centre(3) + Radius%z + 2
			z = blobrange(3,2) + 1
			site = (/x, y, z/)
			itcstate = axis_end
		end select

		j = ninfo*(k-1)
		TC_list(j+1) = k-1
		TC_list(j+2:j+4) = site
		TC_list(j+5) = itcstate
		TC_list(j+6) = 1
		last_id1 = k-1
	enddo
	k = last_id1 + 1
endif

! Cells
nlive = 0
do kcell = 1,nlist
!	if (idbug /= 0 .and. cell_list(kcell)%ID /= idbug) cycle
	if (cell_list(kcell)%exists) then
		k = k+1
		j = ninfo*(k-1)
		site = cell_list(kcell)%site
	    if (highlighting .and. cell_list(kcell)%ID == show_progeny) then
	        highlight = 1
	    else
	        highlight = 0
	    endif
	    if (use_celltype_colour) then
			colour = cell_list(kcell)%celltype
		else
			call CellColour(kcell,highlight,col)
			colour = rgb(col)
		endif
		TC_list(j+1) = kcell + last_id1
		TC_list(j+2:j+4) = site
		TC_list(j+5) = colour
		vol = cell_list(kcell)%volume*Vcell_cm3		! cell volume in cm3
		r = (3*vol/(4*PI))**(1./3.)					! cell radius in cm
		TC_list(j+6) = 200*r/DELTA_X				! 100*diameter as fraction of DELTA_X
!		TC_list(j+6) = (cell_list(kcell)%volume)**(1./3.)	! diameter: 0.928 - 1.17
		TC_list(j+7) = highlight
		last_id2 = kcell + last_id1
		nlive = nlive + 1
	endif
enddo
!nTC_list = last_id2
nTC_list = k
!write(logmsg,*) 'nlive: ',nlive
!call logger(logmsg)
end subroutine

!-----------------------------------------------------------------------------------------
! Rendered cell colour may depend on stage, state, receptor expression level.
! col(:) = (r,g,b)
!-----------------------------------------------------------------------------------------
subroutine CellColour(kcell,highlight,col)
integer :: kcell, highlight, col(3)
integer :: stage, status
integer, parameter :: WHITE(3) = (/255,255,255/)
integer, parameter :: RED(3) = (/255,0,0/)
integer, parameter :: GREEN(3) = (/0,255,0/)
integer, parameter :: BLUE(3) = (/0,0,255/)
integer, parameter :: DEEPRED(3) = (/200,0,0/)
integer, parameter :: DEEPBLUE(3) = (/30,20,255/)
integer, parameter :: DEEPGREEN(3) = (/0,150,0/)
integer, parameter :: LIGHTRED(3) = (/255,70,90/)
integer, parameter :: LIGHTBLUE(3) = (/0,200,255/)
integer, parameter :: LIGHTGREEN(3) = (/50,255,150/)
integer, parameter :: DEEPORANGE(3) = (/240,70,0/)
integer, parameter :: LIGHTORANGE(3) = (/255,130,0/)
integer, parameter :: YELLOW(3) = (/255,255,0/)
integer, parameter :: DEEPPURPLE(3) = (/180,180,30/)
integer, parameter :: LIGHTPURPLE(3) = (/230,230,100/)
integer, parameter :: DEEPBROWN(3) = (/130,70,0/)
integer, parameter :: LIGHTBROWN(3) = (/200,100,0/)
integer, parameter :: GRAY(3) = (/128,128,128/)

integer, parameter :: Qt_white = 3
integer, parameter :: Qt_black = 2
integer, parameter :: Qt_red = 7
integer, parameter :: Qt_darkRed = 13
integer, parameter :: Qt_green = 8
integer, parameter :: Qt_darkGreen = 14
integer, parameter :: Qt_blue = 9
integer, parameter :: Qt_darkBlue = 15
integer, parameter :: Qt_cyan = 10
integer, parameter :: Qt_darkCyan = 16
integer, parameter :: Qt_magenta = 11
integer, parameter :: Qt_darkMagenta = 17
integer, parameter :: Qt_yellow = 12
integer, parameter :: Qt_darkYellow = 18
integer, parameter :: Qt_gray = 5
integer, parameter :: Qt_darkGray = 4
integer, parameter :: Qt_lightGray = 6

if (highlight == 0) then
    col = LIGHTORANGE
else
    col = LIGHTRED
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Pack the colours (r,g,b) into an integer.
!-----------------------------------------------------------------------------------------
integer function rgb(col)
integer :: col(3)

rgb = ishft(col(1),16) + ishft(col(2),8) + col(3)
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getMediumConc(cmedium, cbdry)
real(REAL_KIND) :: cmedium(:), cbdry(:)

cmedium(:) = chemo(:)%medium_Cext
cbdry(:) = chemo(:)%medium_Cbnd
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!subroutine getNecroticFraction(necrotic_fraction,totvol_cm3)
!real(REAL_KIND) :: necrotic_fraction, totvol_cm3	! vol_cm3 not used here, needed in scell
!real(REAL_KIND) :: cellvol_cm3, dvol
!!necrotic_fraction = (Nsites-Ncells)/real(Nsites)
!cellvol_cm3 = Ncells*DELTA_X**3
!dvol = totvol_cm3-cellvol_cm3
!necrotic_fraction = dvol/totvol_cm3
!if (necrotic_fraction < 0.005) necrotic_fraction = 0
!end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getDiamVol(diam_cm,vol_cm3)
real(REAL_KIND) :: diam_cm, vol_cm3
diam_cm = 2*DELTA_X*blob_radius
!vol_cm3 = Vsite_cm3*Nsites			! total volume in cm^3
vol_cm3 = Vsite_cm3*blob_volume		! total volume in cm^3
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
function rint(i) result(r)
integer :: i
real(REAL_KIND) :: r
r = i
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_summary(summaryData,i_hypoxia_cutoff,i_growth_cutoff) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
real(c_double) :: summaryData(*)
integer(c_int) :: i_hypoxia_cutoff,i_growth_cutoff
integer :: Nviable(MAX_CELLTYPES), Nlive(MAX_CELLTYPES)
integer :: nhypoxic(3), nclonohypoxic(3), ngrowth(3), &
    medium_oxygen, medium_glucose, medium_lactate, medium_drug(2,0:2), &
    bdry_oxygen, bdry_glucose, bdry_lactate, bdry_drug(2,0:2)
integer :: TNanoxia_dead, TNaglucosia_dead, TNradiation_dead, TNdrug_dead(2),  TNviable, &
           Ntagged_anoxia(MAX_CELLTYPES), Ntagged_aglucosia(MAX_CELLTYPES), Ntagged_radiation(MAX_CELLTYPES), &
           Ntagged_drug(2,MAX_CELLTYPES), &
           TNtagged_anoxia, TNtagged_aglucosia, TNtagged_radiation, TNtagged_drug(2)
integer :: ityp, i, im, idrug
real(REAL_KIND) :: diam_um, hypoxic_percent, clonohypoxic_percent, growth_percent, necrotic_percent,  npmm3, Tplate_eff
real(REAL_KIND) :: diam_cm, vol_cm3, vol_mm3, hour, necrotic_fraction, doubling_time, plate_eff(MAX_CELLTYPES)
real(REAL_KIND) :: cmedium(MAX_CHEMO), cbdry(MAX_CHEMO), volume_cm3(5), maxarea(5), diameter_um(5)
real(REAL_KIND) :: r_G, r_P, r_A, r_I, P_utilisation

hour = istep*DELTA_T/3600.
call getDiamVol(diam_cm,vol_cm3)
vol_mm3 = vol_cm3*1000				! volume in mm^3
diam_um = diam_cm*10000
npmm3 = Ncells/vol_mm3

Ntagged_anoxia(:) = Nanoxia_tag(:)			! number currently tagged by anoxia
Ntagged_aglucosia(:) = Naglucosia_tag(:)	! number currently tagged by aglucosia
Ntagged_radiation(:) = Nradiation_tag(:)	! number currently tagged by radiation
Ntagged_drug(1,:) = Ndrug_tag(1,:)			! number currently tagged by drugA
Ntagged_drug(2,:) = Ndrug_tag(2,:)			! number currently tagged by drugA

TNtagged_anoxia = sum(Ntagged_anoxia(1:Ncelltypes))
TNtagged_aglucosia = sum(Ntagged_aglucosia(1:Ncelltypes))
TNtagged_radiation = sum(Ntagged_radiation(1:Ncelltypes))
TNtagged_drug(1) = sum(Ntagged_drug(1,1:Ncelltypes))
TNtagged_drug(2) = sum(Ntagged_drug(2,1:Ncelltypes))

TNanoxia_dead = sum(Nanoxia_dead(1:Ncelltypes))
TNaglucosia_dead = sum(Naglucosia_dead(1:Ncelltypes))
TNradiation_dead = sum(Nradiation_dead(1:Ncelltypes))
TNdrug_dead(1) = sum(Ndrug_dead(1,1:Ncelltypes))
TNdrug_dead(2) = sum(Ndrug_dead(2,1:Ncelltypes))

call getNviable(Nviable, Nlive)
TNviable = sum(Nviable(1:Ncelltypes))

call getHypoxicCount(nhypoxic)
hypoxic_percent = (100.*nhypoxic(i_hypoxia_cutoff))/Ncells
call getClonoHypoxicCount(nclonohypoxic)
clonohypoxic_percent = (100.*nclonohypoxic(i_hypoxia_cutoff))/TNviable
call getGrowthCount(ngrowth)
growth_percent = (100.*ngrowth(i_growth_cutoff))/Ncells
!call getNecroticFraction(necrotic_fraction,vol_cm3)
!necrotic_percent = 100.*necrotic_fraction
do ityp = 1,Ncelltypes
	if (Nlive(ityp) > 0) then
		plate_eff(ityp) = real(Nviable(ityp))/Nlive(ityp)
	else
		plate_eff(ityp) = 0
	endif
enddo
plate_eff = 100.*plate_eff
Tplate_eff = 0
do ityp = 1,Ncelltypes
	Tplate_eff = Tplate_eff + plate_eff(ityp)*celltype_fraction(ityp)
enddo

call getMediumConc(cmedium, cbdry)

call getVolumes(volume_cm3,maxarea)
do i = 1,5
	diameter_um(i) = 10000*(6*volume_cm3(i)/PI)**(1./3)
enddo
necrotic_fraction = volume_cm3(1)/volume_cm3(5)
necrotic_percent = 100.*necrotic_fraction

!if (ndoublings > 0) then
!    doubling_time = doubling_time_sum/(3600*ndoublings)
!else
!    doubling_time = 0
!endif

summaryData(1:36) = [ rint(istep), rint(Ncells), rint(TNanoxia_dead), rint(TNaglucosia_dead), rint(TNdrug_dead(1)), rint(TNdrug_dead(2)), rint(TNradiation_dead), &
    rint(TNtagged_anoxia), rint(TNtagged_aglucosia), rint(TNtagged_drug(1)), rint(TNtagged_drug(2)), rint(TNtagged_radiation), &
	diam_um, vol_mm3, hypoxic_percent, clonohypoxic_percent, growth_percent, necrotic_percent, Tplate_eff, npmm3, &
	cmedium(OXYGEN), cmedium(GLUCOSE), cmedium(DRUG_A:DRUG_A+2), cmedium(DRUG_B:DRUG_B+2), &
	cbdry(OXYGEN), cbdry(GLUCOSE), cbdry(DRUG_A:DRUG_A+2), cbdry(DRUG_B:DRUG_B+2) ]
!	doubling_time, r_G, r_P, r_A, r_I, rint(ndoublings), P_utilisation ]
write(nfres,'(a,a,2a12,i8,7e12.4,22i7,29e12.4)') trim(header),' ',gui_run_version, dll_run_version, &
	istep, hour, vol_mm3, diameter_um, Ncells_type(1:2), &
    Nanoxia_dead(1:2), Naglucosia_dead(1:2), Ndrug_dead(1,1:2), &
    Ndrug_dead(2,1:2), Nradiation_dead(1:2), &
    Ntagged_anoxia(1:2), Ntagged_aglucosia(1:2), Ntagged_drug(1,1:2), &
    Ntagged_drug(2,1:2), Ntagged_radiation(1:2), &
	nhypoxic(:)/real(Ncells), nclonohypoxic(:)/real(TNviable), ngrowth(:)/real(Ncells), &
	necrotic_fraction, plate_eff(1:2), &
	cmedium(OXYGEN), cmedium(GLUCOSE), cmedium(DRUG_A:DRUG_A+2), cmedium(DRUG_B:DRUG_B+2), &
	cbdry(OXYGEN), cbdry(GLUCOSE), cbdry(DRUG_A:DRUG_A+2), cbdry(DRUG_B:DRUG_B+2)
!	doubling_time, r_G, r_P, r_A, r_I, real(ndoublings), P_utilisation

	
call sum_dMdt(GLUCOSE)

if (diam_count_limit > LIMIT_THRESHOLD) then
	if (Ncells > diam_count_limit) limit_stop = .true.
elseif (diam_count_limit > 0) then
	if (diam_um > diam_count_limit) limit_stop = .true.
endif

!ndoublings = 0
!doubling_time_sum = 0

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_summary1(summaryData,i_hypoxia_cutoff,i_growth_cutoff) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*), i_hypoxia_cutoff,i_growth_cutoff
integer :: Nviable(MAX_CELLTYPES), Nlive(MAX_CELLTYPES), plate_eff_10(MAX_CELLTYPES)
integer :: diam_um, vol_mm3_1000, nhypoxic(3), nclonohypoxic(3), ngrowth(3), &
    hypoxic_percent_10, clonohypoxic_percent_10, growth_percent_10, necrotic_percent_10,  npmm3, &
    medium_oxygen_1000, medium_glucose_1000, medium_drug_1000(2), &
    bdry_oxygen_1000, bdry_glucose_1000, bdry_drug_1000(2)
integer :: TNanoxia_dead, TNaglucosia_dead, TNradiation_dead, TNdrug_dead(2),  TNviable, &
           Ntagged_anoxia(MAX_CELLTYPES), Ntagged_aglucosia(MAX_CELLTYPES), Ntagged_radiation(MAX_CELLTYPES), &
           Ntagged_drug(2,MAX_CELLTYPES), &
           TNtagged_anoxia, TNtagged_aglucosia, TNtagged_radiation, TNtagged_drug(2)
integer :: Tplate_eff_10   
integer :: ityp
real(REAL_KIND) :: diam_cm, vol_cm3, vol_mm3, hour, plate_eff(MAX_CELLTYPES), necrotic_fraction
real(REAL_KIND) :: cmedium(MAX_CHEMO), cbdry(MAX_CHEMO)

hour = istep*DELTA_T/3600.
call getDiamVol(diam_cm,vol_cm3)
vol_mm3 = vol_cm3*1000				! volume in mm^3
vol_mm3_1000 = vol_mm3*1000			! 1000 * volume in mm^3
diam_um = diam_cm*10000
npmm3 = Ncells/vol_mm3

Ntagged_anoxia(:) = Nanoxia_tag(:)			! number currently tagged by anoxia
Ntagged_aglucosia(:) = Naglucosia_tag(:)	! number currently tagged by aglucosia
Ntagged_radiation(:) = Nradiation_tag(:)	! number currently tagged by radiation
Ntagged_drug(1,:) = Ndrug_tag(1,:)			! number currently tagged by drugA
Ntagged_drug(2,:) = Ndrug_tag(2,:)			! number currently tagged by drugA

TNtagged_anoxia = sum(Ntagged_anoxia(1:Ncelltypes))
TNtagged_aglucosia = sum(Ntagged_aglucosia(1:Ncelltypes))
TNtagged_radiation = sum(Ntagged_radiation(1:Ncelltypes))
TNtagged_drug(1) = sum(Ntagged_drug(1,1:Ncelltypes))
TNtagged_drug(2) = sum(Ntagged_drug(2,1:Ncelltypes))

TNanoxia_dead = sum(Nanoxia_dead(1:Ncelltypes))
TNaglucosia_dead = sum(Naglucosia_dead(1:Ncelltypes))
TNradiation_dead = sum(Nradiation_dead(1:Ncelltypes))
TNdrug_dead(1) = sum(Ndrug_dead(1,1:Ncelltypes))
TNdrug_dead(2) = sum(Ndrug_dead(2,1:Ncelltypes))

call getNviable(Nviable, Nlive)
TNviable = sum(Nviable(1:Ncelltypes))

call getHypoxicCount(nhypoxic)
hypoxic_percent_10 = (1000.*nhypoxic(i_hypoxia_cutoff))/Ncells
call getClonoHypoxicCount(nclonohypoxic)
clonohypoxic_percent_10 = (1000.*nclonohypoxic(i_hypoxia_cutoff))/TNviable
call getGrowthCount(ngrowth)
growth_percent_10 = (1000.*ngrowth(i_growth_cutoff))/Ncells
!call getNecroticFraction(necrotic_fraction,vol_cm3)
!necrotic_percent_10 = 1000.*necrotic_fraction
do ityp = 1,Ncelltypes
	if (Nlive(ityp) > 0) then
		plate_eff(ityp) = real(Nviable(ityp))/Nlive(ityp)
	else
		plate_eff(ityp) = 0
	endif
enddo
plate_eff_10 = 1000.*plate_eff
Tplate_eff_10 = 0
do ityp = 1,Ncelltypes
	Tplate_eff_10 = Tplate_eff_10 + plate_eff_10(ityp)*celltype_fraction(ityp)
enddo
call getMediumConc(cmedium, cbdry)
medium_oxygen_1000 = cmedium(OXYGEN)*1000.
medium_glucose_1000 = cmedium(GLUCOSE)*1000.
medium_drug_1000(1) = cmedium(DRUG_A)*1000.
medium_drug_1000(2) = cmedium(DRUG_B)*1000.
bdry_oxygen_1000 = cbdry(OXYGEN)*1000.
bdry_glucose_1000 = cbdry(GLUCOSE)*1000.
bdry_drug_1000(1) = cbdry(DRUG_A)*1000.
bdry_drug_1000(2) = cbdry(DRUG_B)*1000.

summaryData(1:28) = [ istep, Ncells, TNanoxia_dead, TNaglucosia_dead, TNdrug_dead(1), TNdrug_dead(2), TNradiation_dead, &
    TNtagged_anoxia, TNtagged_aglucosia, TNtagged_drug(1), TNtagged_drug(2), TNtagged_radiation, &
	diam_um, vol_mm3_1000, hypoxic_percent_10, clonohypoxic_percent_10, growth_percent_10, necrotic_percent_10, &
	Tplate_eff_10, npmm3, &
	medium_oxygen_1000, medium_glucose_1000, medium_drug_1000(1), medium_drug_1000(2), &
	bdry_oxygen_1000, bdry_glucose_1000, bdry_drug_1000(1), bdry_drug_1000(2) ]
write(nfres,'(a,a,2a12,i8,2e12.4,23i7,20e12.4)') trim(header),' ',gui_run_version, dll_run_version, &
	istep, hour, vol_mm3, diam_um, Ncells_type(1:2), &
    Nanoxia_dead(1:2), Naglucosia_dead(1:2), Ndrug_dead(1,1:2), &
    Ndrug_dead(2,1:2), Nradiation_dead(1:2), &
    Ntagged_anoxia(1:2), Ntagged_aglucosia(1:2), Ntagged_drug(1,1:2), &
    Ntagged_drug(2,1:2), Ntagged_radiation(1:2), &
	nhypoxic(:)/real(Ncells), nclonohypoxic(:)/real(TNviable), ngrowth(:)/real(Ncells), &
	necrotic_fraction, plate_eff(1:2), &
	cmedium(OXYGEN), cmedium(GLUCOSE), cmedium(DRUG_A), cmedium(DRUG_B), &
	cbdry(OXYGEN), cbdry(GLUCOSE), cbdry(DRUG_A), cbdry(DRUG_B)
		
call sum_dMdt(GLUCOSE)

if (diam_count_limit > LIMIT_THRESHOLD) then
	if (Ncells > diam_count_limit) limit_stop = .true.
elseif (diam_count_limit > 0) then
	if (diam_um > diam_count_limit) limit_stop = .true.
endif
if (Ncelltypes == 2 .and. celltype_fraction(2) > 0) then
	call get_clumpiness(2)
endif
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getHypoxicCount(nhypoxic)
integer :: nhypoxic(3)
integer :: kcell, i

nhypoxic = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	do i = 1,3
		if (cell_list(kcell)%conc(OXYGEN) < O2cutoff(i)) nhypoxic(i) = nhypoxic(i) + 1
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getClonoHypoxicCount(nclonohypoxic)
integer :: nclonohypoxic(3)
integer :: kcell, i, idrug
logical :: tagged

nclonohypoxic = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	if (cell_list(kcell)%anoxia_tag) cycle
	if (cell_list(kcell)%aglucosia_tag) cycle
	if (cell_list(kcell)%radiation_tag) cycle
	tagged = .false.
	do idrug = 1,MAX_DRUGTYPES
	    if (cell_list(kcell)%drug_tag(idrug)) tagged = .true.
	enddo
	if (tagged) cycle
	do i = 1,3
		if (cell_list(kcell)%conc(OXYGEN) < O2cutoff(i)) then
			nclonohypoxic(i) = nclonohypoxic(i) + 1
		endif
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
! Need to compare growth rate with a fraction of average growth rate
!--------------------------------------------------------------------------------
subroutine getGrowthCount(ngrowth)
integer :: ngrowth(3)
integer :: kcell, i, ityp
real(REAL_KIND) :: r_mean(2)

r_mean(1:2) = Vdivide0/(2*divide_time_mean(1:2))
ngrowth = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	ityp = cell_list(kcell)%celltype
	do i = 1,3
		if (cell_list(kcell)%dVdt < growthcutoff(i)*r_mean(ityp)) ngrowth(i) = ngrowth(i) + 1
	enddo
enddo

end subroutine

!--------------------------------------------------------------------------------
! Compute total uptake rate for a constituent
!--------------------------------------------------------------------------------
subroutine sum_dMdt(ichemo)
integer :: ichemo
integer :: kcell, Nh, Nc
real(REAL_KIND) :: C, metab, dMdt, asum, msum, Csum

if (ichemo > GLUCOSE) then
	write(*,*) 'Error: sum_dMdt: only for oxygen and glucose'
	stop
endif
Nh = chemo(ichemo)%Hill_N
asum = 0
Csum = 0
msum = 0
Nc = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	Nc = Nc + 1
	C = cell_list(kcell)%conc(ichemo)
	Csum = Csum + C
	metab = C**Nh/(chemo(ichemo)%MM_C0**Nh + C**Nh)
	msum = msum + metab
	dMdt = metab*chemo(ichemo)%max_cell_rate 
	asum = asum + dMdt
enddo
total_dMdt = total_dMdt + asum
!write(*,'(a,2i6,2e12.3)') 'sum_dMdt: ',ichemo,Nc,asum,total_dMdt*3600 
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getNviable(Nviable, Nlive)
integer :: Nviable(:), Nlive(:)
integer :: kcell, ityp, idrug
logical :: tag

Nviable = 0
Nlive = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
    ityp = cell_list(kcell)%celltype
    Nlive(ityp) = Nlive(ityp) + 1
	if (cell_list(kcell)%anoxia_tag .or. cell_list(kcell)%aglucosia_tag .or. cell_list(kcell)%radiation_tag) cycle
    tag = .false.
    do idrug = 1,ndrugs_used
		if (cell_list(kcell)%drug_tag(idrug)) tag = .true.
	enddo
	if (tag) cycle
	Nviable(ityp) = Nviable(ityp) + 1
enddo	
end subroutine

!--------------------------------------------------------------------------------
! Note: axis = 0,1,2
! Now CFSE is first in the list (0)
!--------------------------------------------------------------------------------
subroutine get_fieldinfo(nxx, axis, fraction, ns, nc, cused, res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_fieldinfo
use, intrinsic :: iso_c_binding
integer(c_int) :: nxx, axis, ns, nc, cused(0:*), res
real(c_double) :: fraction
integer rng(3,2), ichemo, kcell, x, y, z

!call logger('get_fieldinfo')
if (.not.allocated(occupancy)) then
	res = -1
	return
endif
res = 0
nxx = NX
nc = MAX_CHEMO
do ichemo = 1,MAX_CHEMO
    if (chemo(ichemo)%used) then
        cused(ichemo) = 1
    else
        cused(ichemo) = 0
    endif
enddo
cused(CFSE) = 1
cused(GROWTH_RATE) = 1
cused(CELL_VOLUME) = 1
cused(O2_BY_VOL) = 1
rng(:,1) = blob_centre(:) - (adrop*blob_radius + 2)
rng(:,2) = blob_centre(:) + (adrop*blob_radius + 2)
rng(axis,:) = blob_centre(axis) + fraction*blob_radius
!write(logmsg,*) 'Centre, Radius, axis, fraction: ',Centre, Radius, axis, fraction
!call logger(logmsg)
!write(logmsg,*) 'rng: ',rng
!call logger(logmsg)
ns = 0
do z = rng(3,1),rng(3,2)
    do y = rng(2,1),rng(2,2)
        do x = rng(1,1),rng(1,2)
            kcell = occupancy(x,y,z)%indx(1)
            if (kcell <= OUTSIDE_TAG) cycle
            ns = ns + 1
        enddo
    enddo
enddo
!write(logmsg,*) 'did get_fieldinfo: ns: ',ns
!call logger(logmsg)
end subroutine

!-----------------------------------------------------------------------------------------
! Distances are still all cm 
! Added ixyz to pass the appropriate Caverage array index
!-----------------------------------------------------------------------------------------
subroutine new_get_fielddata(axis, fraction, fdata, ixyz, res) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: new_get_fielddata
use, intrinsic :: iso_c_binding
integer(c_int) :: axis, ixyz, res
real(c_double) :: fraction
type (fielddata_type)    :: fdata
type(celldata_type) :: cdata(4000)
type(cell_type), pointer :: cp
real(REAL_KIND) :: x, y, z, rad, csum(3), bcentre(3), dx, p(3)
integer :: ichemo, ix, iy, iz, k, i, i1, i2, kcell, nc, nlump, ii
logical :: in

integer :: ixb,iyb,izb,grid(3)
real(REAL_KIND) :: cb(3), alfa(3)

nlump = 1

!write(nflog,*) 'new_get_fielddata: ',axis,fraction,ixyz
dx = nlump*DELTA_X
fdata%NX = NX/nlump
fdata%NY = NY/nlump
fdata%NZ = NZ/nlump
fdata%NCONST = MAX_CHEMO
fdata%DX = dx
fdata%conc_ptr = c_loc(Cslice)
fdata%cell_ptr = c_loc(cdata)

! Find blob centre
csum = 0
nc = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	nc = nc+1
	csum = csum + cp%site
!	write(nflog,'(a,i6,3f8.4)') 'cell centre: ',kcell,cp%centre(:,1) 
enddo
bcentre = csum*DELTA_X/nc
!write(nflog,'(a,3f8.4)') 'actual blobcentre: ',bcentre

! Set up concentration values for the slice
! start by fixing on central plane, ignoring fraction
! Note that in Qt, scene->addRect(x,y,w,w) places the square (width w) with the top left corner at (x,y)
! On a z-slice (fixed iz) Cslice(ix,iy) fills a square at x = ix*dx, y = iy*dx
! for ix=0,NX-1, iy=0,NY-1 (0-based indexing in Qt)
! Equivalent to ix=1,NX, iy=1,NY on the Fortran side.  I.e. Fortran Cslice(1,1) -> Qt Cslice(0,0) -> (x,y) = (0,0)
! 

nc = 0
if (axis == X_AXIS) then
	x = bcentre(1)	
	ix = x/dx		! approx?
	ixyz = ix
	do iy = 1,fdata%NY
	    do iz = 1,fdata%NZ
	        if (occupancy(ix,iy,iz)%indx(1) <= OUTSIDE_TAG) then
	            in = .false.
	            p(:) = [ix,iy,iz]*DELTA_X + grid_offset
	            call getConc(p,Cslice(ix,iy,iz,:))
	        else
	            in = .true.
 	            i = ODEdiff%ivar(ix,iy,iz)
	            if (i < 1) then		! this solves (sort of) some startup problems with the 2D display
					do ii = -5,5
		 	            i = ODEdiff%ivar(ix,iy+ii,iz+ii)
		 	            if (i >= 1) exit
		 	        enddo
		 	        write(nflog,*) 'new_get_fielddata: ii,i: ',ii,i
			    endif
			    Cslice(ix,iy,iz,:) = allstate(i,1:MAX_CHEMO)
    	    endif
!	        if (iy == NY/2) write(nflog,'(a,i6,L2,2e12.3)') 'blob: ',iz,in,Cslice(ix,iy,iz,1:2)
            do i1 = 1-nlump,0
            do i2 = 1-nlump,0
	        kcell = occupancy(nlump*ix,nlump*iy+i1,nlump*iz+i2)%indx(1)
	            if (kcell > 0) then
	                cp => cell_list(kcell)
	                nc = nc + 1
	                rad = (3*cp%volume*Vcell_cm3/(4*PI))**(1./3)
			        cdata(nc)%radius = rad
			        cdata(nc)%centre(1:2) = DELTA_X*[cp%site(2)- 1,NZ - cp%site(3)]
                    cdata(nc)%status = getstatus(cp)
			    endif
			enddo
			enddo
	    enddo
	enddo
elseif (axis == Y_AXIS) then
	y = bcentre(2)	
	iy = y/dx		! approx?
	ixyz = iy
	do ix = 1,fdata%NX
	    do iz = 1,fdata%NZ
	        if (occupancy(ix,iy,iz)%indx(1) <= OUTSIDE_TAG) then
	            in = .false.
	            p(:) = [ix,iy,iz]*DELTA_X + grid_offset
	            call getConc(p,Cslice(ix,iy,iz,:))
	        else
	            in = .true.
 	            i = ODEdiff%ivar(ix,iy,iz)
	            if (i < 1) then
					do ii = -5,5
		 	            i = ODEdiff%ivar(ix+ii,iy,iz+ii)
		 	            if (i >= 1) exit
		 	        enddo
		 	        write(nflog,*) 'new_get_fielddata: ii,i: ',ii,i
			    endif
			    Cslice(ix,iy,iz,:) = allstate(i,1:MAX_CHEMO)
    	    endif
!	        if (ix == NX/2) write(nflog,'(a,i6,L2,2e12.3)') 'blob: ',iz,in,Cslice(ix,iy,iz,1:2)
            do i1 = 1-nlump,0
            do i2 = 1-nlump,0
	            kcell = occupancy(nlump*ix+i1,nlump*iy,nlump*iz+i2)%indx(1)
	            if (kcell > 0) then
	                cp => cell_list(kcell)
	                nc = nc + 1
	                rad = (3*cp%volume*Vcell_cm3/(4*PI))**(1./3)
			        cdata(nc)%radius = rad
			        cdata(nc)%centre(1:2) = DELTA_X*[cp%site(1) - 1, NZ - cp%site(3)]
                    cdata(nc)%status = getstatus(cp)
			    endif
			enddo
			enddo
	    enddo
	enddo
elseif (axis == Z_AXIS) then
	z = bcentre(3)	
	iz = z/dx		! approx?
	ixyz = iz
	!
	! checking why Cslice can exceed Cave_b
!    ix = fdata%NX/2
!    iy = 1
!	cb(:) = [ix,iy,iz]*DELTA_X + grid_offset
!    ixb = cb(1)/DXB + 1
!    iyb = cb(2)/DXB + 1
!    izb = cb(3)/DXB + 1
!    grid = [ixb, iyb, izb]
!    do i = 1,3
!	    alfa(i) = (cb(i) - (grid(i)-1)*DXB)/DXB
!    enddo
!    write(nflog,'(a,6i4,3f8.3)') 'new_get_fielddata: ix,iy,iz,grid,alfa: ',ix,iy,iz,grid,alfa
!    ix = fdata%NX/2
!    iy = fdata%NY
!	cb(:) = [ix,iy,iz]*DELTA_X + grid_offset
!    ixb = cb(1)/DXB + 1
!    iyb = cb(2)/DXB + 1
!    izb = cb(3)/DXB + 1
!    grid = [ixb, iyb, izb]
!    do i = 1,3
!	    alfa(i) = (cb(i) - (grid(i)-1)*DXB)/DXB
!    enddo
!    write(nflog,'(a,6i4,3f8.3)') 'new_get_fielddata: ix,iy,iz,grid,alfa: ',ix,iy,iz,grid,alfa
    
	do iy = 1,fdata%NY
	    do ix = 1,fdata%NX
	        ! The field outside the blob is estimated by interpolation from the coarse grid solution.
	        ! Inside the blob the lattice solution is used.
	        ! Need to be sure of mapping (ix,iy,iz) -> (x,y,z) in coarse grid axes 
	        if (occupancy(ix,iy,iz)%indx(1) <= OUTSIDE_TAG) then
	            in = .false.
	            p(:) = [ix,iy,iz]*DELTA_X + grid_offset
	            call getConc(p,Cslice(ix,iy,iz,:))
	        else
	            in = .true.
 	            i = ODEdiff%ivar(ix,iy,iz)
	            if (i < 1) then
!				    res = -1
!				    write(nflog,*) 'new_get_fielddata: Z_AXIS: ixyz,i,res: ',ixyz,i,res
!				    write(nflog,*) 'ix,iy,iz: ',ix,iy,iz
!				    return
!					Cslice(ix,iy,iz,:) = 0
					do ii = -5,5
		 	            i = ODEdiff%ivar(ix+ii,iy+ii,iz)
		 	            if (i >= 1) exit
		 	        enddo
		 	        write(nflog,*) 'new_get_fielddata: ii,i: ',ii,i
			    endif
			    Cslice(ix,iy,iz,:) = allstate(i,1:MAX_CHEMO)
    	    endif
!	        if (iy == NY/2) write(nflog,'(a,3i6,L2,2e12.3)') 'Cslice: O2: ',ix,iy,iz,in,Cslice(ix,iy,iz,1)
            do i1 = 1-nlump,0
            do i2 = 1-nlump,0
	            kcell = occupancy(nlump*ix+i1,nlump*iy+i2,nlump*iz)%indx(1)
	            if (kcell > 0) then
	                cp => cell_list(kcell)
	                nc = nc + 1
	                rad = (3*cp%volume*Vcell_cm3/(4*PI))**(1./3)
!	                write(nflog,'(a,i6,2e12.3)') 'rad, DELTA_X: ',kcell,rad,DELTA_X
			        cdata(nc)%radius = rad
!			        cdata(nc)%centre(1:2) = DELTA_X*[cp%site(1)-1.5,cp%site(2)-1.5]
			        cdata(nc)%centre(1:2) = DELTA_X*[cp%site(1)-1,cp%site(2)-1]
                    cdata(nc)%status = getstatus(cp)
			    endif
			 enddo
			 enddo
	    enddo
	enddo
!	write(nflog,*) 'Cslice: O2:' 
!	write(nflog,'(10e12.3)') Cslice(NX/2,:,iz,1)

endif
ixyz = ixyz-1   ! to use in C with 0-based indexing
write(nflog,*) 'axis: ',axis,' nc: ',nc, ' NX: ',fdata%NX,' dx: ',fdata%dx
fdata%ncells = nc
res = 0

!nc = 0
!do kcell = 1,nlist
!	cp => cell_list(kcell)
!	if (cp%state == DEAD) cycle
!	do i = 1,cp%nspheres
!		c = cp%centre(:,i)
!		r = cp%radius(i)
!		if (axis == X_AXIS) then
!			if (c(1) + r < x .or. c(1) - r > x) cycle
!			rad = sqrt(r**2-(c(1)-x)**2)
!			nc = nc + 1
!			cdata(nc)%radius = rad
!			cdata(nc)%centre(1:2) = [c(2),(NX-1.1)*DELTA_X - c(3)]		! 2D coordinates
!			cdata(nc)%status = getstatus(cp)
!		elseif (axis == Y_AXIS) then
!			if (c(2) + r < y .or. c(2) - r > y) cycle
!			rad = sqrt(r**2-(c(2)-y)**2)
!			nc = nc + 1
!			cdata(nc)%radius = rad
!			cdata(nc)%centre(1:2) = [c(1),(NX-1.1)*DELTA_X - c(3)]
!			cdata(nc)%status = getstatus(cp)
!		elseif (axis == Z_AXIS) then
!			if (c(3) + r < z .or. c(3) - r > z) cycle
!			rad = sqrt(r**2-(c(3)-z)**2)
!			nc = nc + 1
!			cdata(nc)%radius = rad
!			cdata(nc)%centre(1:2) = [c(1),c(2)]		! always use centre(1:2) to store the 2D coordinates
!			cdata(nc)%status = getstatus(cp)
!		endif
!		if (cp%anoxia_tag) then
!			cdata(nc)%status = 2	! tagged to die of anoxia
!		elseif (cp%radiation_tag) then
!			cdata(nc)%status = 10
!			write(nflog,*) 'Tagged to die from radiation: ',kcell
!		elseif (cp%drug_tag(1)) then
!			cdata(nc)%status = 11
!		elseif (cp%drug_tag(1)) then
!			cdata(nc)%status = 12
!		elseif (cp%Cin(OXYGEN) < hypoxia_threshold) then
!			cdata(nc)%status = 1	! radiobiological hypoxia
!		elseif (cp%mitosis > 0) then
!			cdata(nc)%status = 3	! in mitosis
!		else
!			cdata(nc)%status = 0
!		endif
!	enddo
!enddo
!write(nflog,*) 'axis: ',axis,' nc: ',nc
!fdata%ncells = nc
!res = 0
end subroutine

!--------------------------------------------------------------------------------
! We want to evaluate Cslice(ix,iy,iz,:) at the centre of the lattice site (ix,iy,iz)
! Note that
!--------------------------------------------------------------------------------
subroutine fillCslice(ch,ixyz,val,dx)
character :: ch
integer :: ixyz
real(REAL_KIND) :: val, dx
integer :: ix, iy, iz
real(REAL_KIND) :: x, y, z

if (ch == 'X') then
    ix = ixyz
    x = val
    
elseif (ch == 'Y') then
    iy = ixyz
    y = val

elseif (ch == 'Z') then
    iz = ixyz
    z = val

endif
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
function getstatus(cp) result(status)
type(cell_type), pointer :: cp
integer :: status

if (cp%anoxia_tag) then
	status = 2	! tagged to die of anoxia
elseif (cp%aglucosia_tag) then
	status = 4	! tagged to die of aglucosia
elseif (cp%radiation_tag) then
	status = 10
elseif (cp%drug_tag(1)) then
	status = 11
elseif (cp%drug_tag(1)) then
	status = 12
elseif (cp%conc(OXYGEN) < hypoxia_threshold) then
	status = 1	! radiobiological hypoxia
!elseif (cp%mitosis > 0) then
elseif (cp%volume > 0.9*cp%divide_volume) then  ! just a surrogate for mitosis
	status = 3	! in mitosis
else
	status = 0
endif
end function

!--------------------------------------------------------------------------------
! Need to transmit medium concentration data.  This could be a separate subroutine.
!--------------------------------------------------------------------------------
subroutine get_fielddata(axis, fraction, nfdata, nc, fdata, res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_fielddata
use, intrinsic :: iso_c_binding
real(c_double) :: fraction
integer(c_int) :: axis, nc, nfdata, res
type(FIELD_DATA) :: fdata(*)
integer rng(3,2), kcell, ityp, x, y, z, i, ns
real(REAL_KIND) :: growthrate, cellvolume, cfse, rmax(2), c_rate(2), r_mean(2), rmaxx

!write(logmsg,*) 'get_fielddata: nfdata, nc: ',nfdata, nc, MAX_CHEMO
!call logger(logmsg)
res = 0
if (nc > MAX_CHEMO) then
	write(logmsg,*) 'Error: get_fielddata: dimension of conc(MAX_CHEMO) not big enough!'
	call logger(logmsg)
	res = 1
	return
endif
c_rate(1:2) = log(2.0)/divide_time_mean(1:2)
r_mean(1:2) = Vdivide0/(2*divide_time_mean(1:2))
if (use_V_dependence) then
	rmax = c_rate*Vdivide0
else
	rmax = r_mean
endif

rng(:,1) = blob_centre(:) - (adrop*blob_radius + 2)
rng(:,2) = blob_centre(:) + (adrop*blob_radius + 2)
rng(axis,:) = blob_centre(axis) + fraction*blob_radius
ns = 0
do z = rng(3,1),rng(3,2)
    do y = rng(2,1),rng(2,2)
        do x = rng(1,1),rng(1,2)
            kcell = occupancy(x,y,z)%indx(1)
            if (kcell <= OUTSIDE_TAG) cycle
!            write(nflog,*) x,y,z,kcell
            ns = ns + 1
	        i = ODEdiff%ivar(x,y,z)
	        if (i < 1) then
				res = -1
				return
			endif
            fdata(ns)%site = (/x,y,z/)
            if (kcell > 0) then
	            fdata(ns)%state = 0		! set state based on: hypoxic, tagged to die, ready to divide? 
	            fdata(ns)%state = getCellState(kcell)
                fdata(ns)%volume = cell_list(kcell)%volume
                call getExtraVars(kcell,growthrate,cellvolume,cfse)
				ityp = cell_list(kcell)%celltype 
                rmaxx = rmax(ityp)
            else
	            fdata(ns)%state = -1
                fdata(ns)%volume = 0
                growthrate = 0
                cellvolume = 0
                cfse = 0
                rmaxx = 1
            endif
!            fdata(ns)%dVdt = growthrate
			fdata(ns)%conc(0) = cfse
            fdata(ns)%conc(1:MAX_CHEMO) = allstate(i,1:MAX_CHEMO)			! cell_list(kcell)%conc(1:MAX_CHEMO)
            fdata(ns)%conc(GROWTH_RATE) = min(1.0,growthrate/rmaxx)			! fraction of max growth rate
            fdata(ns)%conc(CELL_VOLUME) = Vcell_pL*cellvolume				! Note: = fdata(ns)%volume, therefore redundant
            fdata(ns)%conc(O2_BY_VOL) = allstate(i,OXYGEN)*Vcell_pL*cellvolume	
        enddo
    enddo
enddo
!write(logmsg,*) 'did get_fielddata: ns: ',ns
!call logger(logmsg)
if (ns /= nfdata) then
    write(logmsg,*) 'Error: inconsistent nsites: ',nfdata, ns
    call logger(logmsg)
    res = 2
    return
endif

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
integer function getCellState(kcell) 
integer :: kcell
type(cell_type), pointer :: cp

cp => cell_list(kcell)
if (cp%anoxia_tag) then
	getCellState = 2	! tagged to die of anoxia
elseif (cp%aglucosia_tag) then
	getCellState = 4	! tagged to die of aglucosia
elseif (cp%radiation_tag) then
	getCellState = 10
elseif (cp%drug_tag(1)) then
	getCellState = 11
elseif (cp%drug_tag(2)) then
	getCellState = 12
elseif (cp%conc(OXYGEN) < hypoxia_threshold) then
	getCellState = 1	! radiobiological hypoxia
elseif (cp%volume/cp%divide_volume > 0.95) then
	getCellState = 3	! in mitosis
else
	getCellState = 0
endif
end function

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getExtraVars(kcell,growthrate,cellvolume,cfse)
integer :: kcell
real(REAL_KIND) :: growthrate, cellvolume, cfse

cfse = cell_list(kcell)%CFSE
growthrate = cell_list(kcell)%dVdt
cellvolume = cell_list(kcell)%volume
end subroutine

!--------------------------------------------------------------------------------
! Returns all the extracellular concentrations along a line through the blob centre.
! Together with CFSE, growth rate (dVdt), cell volume,...
! Store the constituent profiles one after the other.
! Drop the extra end points
!--------------------------------------------------------------------------------
subroutine get_concdata(nvars, ns, dx, ex_conc) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_concdata
use, intrinsic :: iso_c_binding
integer(c_int) :: nvars, ns
real(c_double) :: dx, ex_conc(0:*)
real(REAL_KIND) :: cbnd, cmin = 1.0e-7
integer :: rng(3,2), i, ic, k, ichemo, kcell, x, y, z, x1, x2, offset

!call logger('get_concdata')
nvars = 1 + MAX_CHEMO + N_EXTRA
dx = DELTA_X
rng(:,1) = blob_centre(:) - (adrop*blob_radius + 2)
rng(:,2) = blob_centre(:) + (adrop*blob_radius + 2)
!rng(axis,:) = Centre(axis) + fraction*Radius
y = blob_centre(2) + 0.5
z = blob_centre(3) + 0.5

! First need to establish the range of x that is inside the blob: (x1,x2)
x1 = 0
x2 = 0
do x = rng(1,1),rng(1,2)
	kcell = occupancy(x,y,z)%indx(1)
	if (kcell <= OUTSIDE_TAG) then
		if (x1 == 0) then
			cycle
		else
			exit
		endif
	elseif (x1 == 0) then
		x1 = x
	endif
	x2 = x
enddo
if (x2 == 0) then
	call logger('Blob is destroyed!')
	return
endif

ns = x2 - x1 + 1 
do ichemo = 0,nvars-1
	offset = ichemo*ns
	k = offset - 1
	do x = x1, x2
		k = k + 1
		kcell = occupancy(x,y,z)%indx(1)
		if (kcell <= OUTSIDE_TAG) then
			ex_conc(k) = 0
			cycle
		endif
		i = ODEdiff%ivar(x,y,z)
        if (ichemo == 0) then	! CFSE
			if (kcell > 0) then
				ex_conc(k) = cell_list(kcell)%CFSE
			else
				ex_conc(k) = 0
			endif
       elseif (ichemo <= MAX_CHEMO) then
			if (chemo(ichemo)%used) then
				if (i > 0) then
					ex_conc(k) = allstate(i,ichemo)	
				else
					ex_conc(k) = 0
				endif
			else
				ex_conc(k) = 0
			endif
        elseif (ichemo == GROWTH_RATE) then	
			if (kcell > 0) then
				ex_conc(k) = cell_list(kcell)%dVdt
			else
				ex_conc(k) = 0
			endif
        elseif (ichemo == CELL_VOLUME) then	
			if (kcell > 0) then
				ex_conc(k) = Vcell_pL*cell_list(kcell)%volume
			else
				ex_conc(k) = 0
			endif
        elseif (ichemo == O2_BY_VOL) then	
			if (kcell > 0) then
				ex_conc(k) = allstate(i,OXYGEN)*Vcell_pL*cell_list(kcell)%volume
			else
				ex_conc(k) = 0
			endif
		endif
    enddo
enddo

ichemo = GLUCOSE
offset = ichemo*ns
!write(*,'(10f7.3)') ex_conc(offset:offset+ns-1)
end subroutine

!--------------------------------------------------------------------------------
! Returns all the intracellular concentrations along a line through the blob centre.
! Together with CFSE, growth rate (dVdt), cell volume,...
! Store the constituent profiles one after the other.
!--------------------------------------------------------------------------------
subroutine get_IC_concdata(nvars, ns, dx, ic_conc) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ic_concdata
use, intrinsic :: iso_c_binding
integer(c_int) :: nvars, ns
real(c_double) :: dx, ic_conc(0:*)
real(REAL_KIND) :: cbnd, cmin = 1.0e-7
integer :: rng(3,2), i, ic, k, ichemo, kcell, x, y, z, x1, x2, offset
type(cell_type), pointer :: cp

!call logger('get_IC_concdata')
nvars = 1 + MAX_CHEMO + N_EXTRA
dx = DELTA_X
rng(:,1) = blob_centre(:) - (adrop*blob_radius + 2)
rng(:,2) = blob_centre(:) + (adrop*blob_radius + 2)
!rng(axis,:) = Centre(axis) + fraction*Radius
y = blob_centre(2) + 0.5
z = blob_centre(3) + 0.5

! First need to establish the range of x that is inside the blob: (x1,x2)
x1 = 0
x2 = 0
do x = rng(1,1),rng(1,2)
	kcell = occupancy(x,y,z)%indx(1)
	if (kcell <= OUTSIDE_TAG) then
		if (x1 == 0) then
			cycle
		else
			exit
		endif
	elseif (x1 == 0) then
		x1 = x
	endif
	x2 = x
enddo
if (x2 == 0) then
	call logger('Blob is destroyed!')
	return
endif

ns = x2 - x1 + 1 
do ichemo = 0,nvars-1
	offset = ichemo*ns
	k = offset - 1
	do x = x1, x2
		k = k + 1
		kcell = occupancy(x,y,z)%indx(1)
!		if (kcell <= OUTSIDE_TAG) then
		if (kcell <= 0) then
			ic_conc(k) = 0
			cycle
		else
		    cp => cell_list(kcell)
		endif
		i = ODEdiff%ivar(x,y,z)
        if (ichemo == 0) then	! CFSE
			if (kcell > 0) then
				ic_conc(k) = cp%CFSE
			else
				ic_conc(k) = 0
			endif
       elseif (ichemo <= MAX_CHEMO) then
			if (chemo(ichemo)%used) then
!				if (i > 0) then
!					ex_conc(k) = allstate(i,ichemo)
                if (kcell > 0) then
                    ic_conc(k) = cp%conc(ichemo)
				else
					ic_conc(k) = 0
				endif
			else
				ic_conc(k) = 0
			endif
        elseif (ichemo == GROWTH_RATE) then	
			if (kcell > 0) then
				ic_conc(k) = cp%dVdt
			else
				ic_conc(k) = 0
			endif
        elseif (ichemo == CELL_VOLUME) then	
			if (kcell > 0) then
				ic_conc(k) = Vcell_pL*cp%volume
			else
				ic_conc(k) = 0
			endif
        elseif (ichemo == O2_BY_VOL) then	
			if (kcell > 0) then
!				ex_conc(k) = allstate(i,OXYGEN)*Vcell_pL*cell_list(kcell)%volume
				ic_conc(k) = cp%conc(OXYGEN)*Vcell_pL*cp%volume
			else
				ic_conc(k) = 0
			endif
		endif
    enddo
enddo

ichemo = GLUCOSE
offset = ichemo*ns
!write(*,'(10f7.3)') ic_conc(offset:offset+ns-1)
end subroutine

!--------------------------------------------------------------------------------
! Returns the distribution of cell volume.
! nv is passed from the GUI
! Min divide volume = Vdivide0 - dVdivide
! therefore Vmin = (Vdivide0 - dVdivide)/2
! Max divide volume = Vmax = Vdivide0 + dVdivide
! dv = (Vmax - Vmin)/nv
! v0 = Vmin + dv/2
!--------------------------------------------------------------------------------
subroutine get_volprob(nv, v0, dv, prob) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_volprob
use, intrinsic :: iso_c_binding
integer(c_int) :: nv
real(c_double) :: v0, dv, prob(*)
integer :: n, kcell, k
real(REAL_KIND) :: v, Vmin, Vmax

!call logger('get_volprob')
Vmin = (Vdivide0 - dVdivide)/2
Vmax = Vdivide0 + dVdivide
dv = (Vmax - Vmin)/nv
v0 = Vmin + dv/2
prob(1:nv) = 0
n = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	v = cell_list(kcell)%volume
	k = (v - Vmin)/dv + 1
	k = min(k,nv)
	prob(k) = prob(k) + 1
	n = n+1
enddo
prob(1:nv) = prob(1:nv)/n
end subroutine

!--------------------------------------------------------------------------------
! Returns the distribution of intracellular O2 level
! nv is passed from the GUI
!--------------------------------------------------------------------------------
subroutine get_oxyprob(nv, v0, dv, prob) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_oxyprob
use, intrinsic :: iso_c_binding
integer(c_int) :: nv
real(c_double) :: v0, dv, prob(*)
integer :: n, kcell, k
real(REAL_KIND) :: v, Vmin, Vmax, O2max

!call logger('get_oxyprob')
Vmin = 0
Vmax = chemo(OXYGEN)%bdry_conc
v0 = Vmin
dv = (Vmax - Vmin)/nv
prob(1:nv) = 0
n = 0
O2max = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	v = cell_list(kcell)%conc(OXYGEN)
	O2max = max(O2max,v)
	k = (v - Vmin)/dv + 1
	k = min(k,nv)
	k = max(k,1)
	prob(k) = prob(k) + 1
	n = n+1
enddo
prob(1:nv) = prob(1:nv)/n
end subroutine

!--------------------------------------------------------------------------------
! Returns the distribution of intracellular constituent concentration
! nv is passed from the GUI
!--------------------------------------------------------------------------------
subroutine get_concprob(ichemo, nv, v0, dv, prob) !BIND(C)
!!DEC$ ATTRIBUTES DLLEXPORT :: get_concprob
use, intrinsic :: iso_c_binding
integer(c_int) :: ichemo, nv
real(c_double) :: v0, dv, prob(*)
integer :: n, kcell, k
real(REAL_KIND) :: v, Vmin, Vmax

!call logger('get_oxyprob')
Vmin = 0
!Vmax = chemo(ichemo)%bdry_conc
Vmax = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	v = cell_list(kcell)%conc(ichemo)
	Vmax = max(Vmax,v)
enddo
VMax = 1.05*Vmax;
v0 = Vmin
dv = (Vmax - Vmin)/nv
prob(1:nv) = 0
n = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	v = cell_list(kcell)%conc(ichemo)
	k = (v - Vmin)/dv + 1
	k = min(k,nv)
	k = max(k,1)
	prob(k) = prob(k) + 1
	n = n+1
enddo
prob(1:nv) = prob(1:nv)/n
end subroutine

!--------------------------------------------------------------------------------
!subroutine get_distdata(nv, ddata) BIND(C)
!!DEC$ ATTRIBUTES DLLEXPORT :: get_distdata
!integer(c_int) :: nv
!type(dist_data) :: ddata(0:*)
!integer :: ichemo
!
!do ichemo = 0,MAX_CHEMO+N_EXTRA
!	if (.not.ddata(ichemo).used) continue
!	if (ichemo  == 0) then	! CFSE
!	elseif (ichemo <= MAX_CHEMO) then
!	elseif (ichemo == GROWTH_RATE) then
!	elseif (ichemo == CELL_VOLUME) then
!	endif
!enddo
!
!end subroutine

!-----------------------------------------------------------------------------------------
! Get number of live cells
!-----------------------------------------------------------------------------------------
subroutine get_nFACS(n) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nfacs
use, intrinsic :: iso_c_binding
integer(c_int) :: n
integer :: k, kcell

n = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	n = n+1
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_FACS(facs_data) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_facs
use, intrinsic :: iso_c_binding
real(c_double) :: val, facs_data(*)
integer :: k, kcell, iextra, ichemo, ivar, nvars, var_index(32)
real(REAL_KIND) :: cfse_min

nvars = 1	! CFSE
var_index(nvars) = 0
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	nvars = nvars + 1
	var_index(nvars) = ichemo
enddo
do iextra = 1,N_EXTRA-1
	nvars = nvars + 1
	var_index(nvars) = MAX_CHEMO + iextra
enddo
cfse_min = 1.0e20
k = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cell_list(kcell)%CFSE
			cfse_min = min(val,cfse_min)
		elseif (ichemo <= MAX_CHEMO) then
			val = cell_list(kcell)%conc(ichemo)
		elseif (ichemo == GROWTH_RATE) then
			val = cell_list(kcell)%dVdt
		elseif (ichemo == CELL_VOLUME) then
			val = cell_list(kcell)%volume
		elseif (ichemo == O2_BY_VOL) then
			val = cell_list(kcell)%volume*cell_list(kcell)%conc(OXYGEN)
		endif
		k = k+1
		facs_data(k) = val
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! nhisto is the number of histogram boxes
! vmax(ivar) is the maximum value for variable ivar
! Probably need to adjust vmax() to a roundish value
!
! Compute 3 distributions: 1 = both cell types
!                          2 = type 1
!                          3 = type 2
! Stack three cases in vmax() and histo_data()
!-----------------------------------------------------------------------------------------
subroutine get_histo(nhisto, histo_data, vmin, vmax, histo_data_log, vmin_log, vmax_log) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_histo
use, intrinsic :: iso_c_binding
integer(c_int),value :: nhisto
real(c_double) :: vmin(*), vmax(*), histo_data(*)
real(c_double) :: vmin_log(*), vmax_log(*), histo_data_log(*)
real(REAL_KIND) :: val, val_log
integer :: n(3), i, ih, k, kcell, ict, ichemo, ivar, nvars, var_index(32)
integer,allocatable :: cnt(:,:,:)
real(REAL_KIND),allocatable :: dv(:,:), valmin(:,:), valmax(:,:)
integer,allocatable :: cnt_log(:,:,:)
real(REAL_KIND),allocatable :: dv_log(:,:), valmin_log(:,:), valmax_log(:,:)
!real(REAL_KIND) :: vmin_log(100), vmax_log(100)
!real(REAL_KIND),allocatable :: histo_data_log(:)

nvars = 1	! CFSE
var_index(nvars) = 0
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	nvars = nvars + 1
	var_index(nvars) = ichemo
enddo
nvars = nvars + 1
var_index(nvars) = GROWTH_RATE
nvars = nvars + 1
var_index(nvars) = CELL_VOLUME
nvars = nvars + 1
var_index(nvars) = O2_BY_VOL

allocate(cnt(3,nvars,nhisto))
allocate(dv(3,nvars))
allocate(valmin(3,nvars))
allocate(valmax(3,nvars))
allocate(cnt_log(3,nvars,nhisto))
allocate(dv_log(3,nvars))
allocate(valmin_log(3,nvars))
allocate(valmax_log(3,nvars))
!allocate(histo_data_log(10000))
cnt = 0
valmin = 1.0e10
valmax = -1.0e10
cnt_log = 0
valmin_log = 1.0e10
valmax_log = -1.0e10
n = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	ict = cell_list(kcell)%celltype
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cell_list(kcell)%CFSE
		elseif (ichemo <= MAX_CHEMO) then
			val = cell_list(kcell)%conc(ichemo)
		elseif (ichemo == GROWTH_RATE) then
			val = cell_list(kcell)%dVdt
		elseif (ichemo == CELL_VOLUME) then
			val = Vcell_pL*cell_list(kcell)%volume
		elseif (ichemo == O2_BY_VOL) then
			val = cell_list(kcell)%conc(OXYGEN)*Vcell_pL*cell_list(kcell)%volume
		endif
		valmax(ict+1,ivar) = max(valmax(ict+1,ivar),val)	! cell type 1 or 2
		valmax(1,ivar) = max(valmax(1,ivar),val)			! both
		if (val <= 1.0e-8) then
			val_log = -8
		else
			val_log = log10(val)
		endif
		valmin_log(ict+1,ivar) = min(valmin_log(ict+1,ivar),val_log)	! cell type 1 or 2
		valmin_log(1,ivar) = min(valmin_log(1,ivar),val_log)			! both
		valmax_log(ict+1,ivar) = max(valmax_log(ict+1,ivar),val_log)	! cell type 1 or 2
		valmax_log(1,ivar) = max(valmax_log(1,ivar),val_log)			! both
	enddo
	n(ict+1) = n(ict+1) + 1
	n(1) = n(1) + 1
enddo
do ivar = 1,nvars
	ichemo = var_index(ivar)
	if (ichemo == CELL_VOLUME) then
		valmin(:,ivar) = Vcell_pL*0.8
		valmin_log(:,ivar) = log10(Vcell_pL*0.8)
	else
		valmin(:,ivar) = 0
	endif
enddo

dv = (valmax - valmin)/nhisto
!write(nflog,*) 'dv'
!write(nflog,'(e12.3)') dv
dv_log = (valmax_log - valmin_log)/nhisto
!write(nflog,*) 'dv_log'
!write(nflog,'(e12.3)') dv_log
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	ict = cell_list(kcell)%celltype
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cell_list(kcell)%CFSE
		elseif (ichemo <= MAX_CHEMO) then
			val = cell_list(kcell)%conc(ichemo)
		elseif (ichemo == GROWTH_RATE) then
			val = cell_list(kcell)%dVdt
		elseif (ichemo == CELL_VOLUME) then
			val = Vcell_pL*cell_list(kcell)%volume
		elseif (ichemo == O2_BY_VOL) then
			val = cell_list(kcell)%conc(OXYGEN)*Vcell_pL*cell_list(kcell)%volume
		endif
		k = (val-valmin(1,ivar))/dv(1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt(1,ivar,k) = cnt(1,ivar,k) + 1
		k = (val-valmin(ict+1,ivar))/dv(ict+1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt(ict+1,ivar,k) = cnt(ict+1,ivar,k) + 1
		if (val <= 1.0e-8) then
			val_log = -8
		else
			val_log = log10(val)
		endif
		k = (val_log-valmin_log(1,ivar))/dv_log(1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt_log(1,ivar,k) = cnt_log(1,ivar,k) + 1
		k = (val_log-valmin_log(ict+1,ivar))/dv_log(ict+1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt_log(ict+1,ivar,k) = cnt_log(ict+1,ivar,k) + 1
	enddo
enddo

do i = 1,3
	if (n(i) == 0) then
		vmin((i-1)*nvars+1:i*nvars) = 0
		vmax((i-1)*nvars+1:i*nvars) = 0
		histo_data((i-1)*nvars*nhisto+1:i*nhisto*nvars) = 0
		vmin_log((i-1)*nvars+1:i*nvars) = 0
		vmax_log((i-1)*nvars+1:i*nvars) = 0
		histo_data_log((i-1)*nvars*nhisto+1:i*nhisto*nvars) = 0
	else
		do ivar = 1,nvars
			vmin((i-1)*nvars+ivar) = valmin(i,ivar)
			vmax((i-1)*nvars+ivar) = valmax(i,ivar)
			do ih = 1,nhisto
				k = (i-1)*nvars*nhisto + (ivar-1)*nhisto + ih
				histo_data(k) = (100.*cnt(i,ivar,ih))/n(i)
			enddo
			vmin_log((i-1)*nvars+ivar) = valmin_log(i,ivar)
			vmax_log((i-1)*nvars+ivar) = valmax_log(i,ivar)
			do ih = 1,nhisto
				k = (i-1)*nvars*nhisto + (ivar-1)*nhisto + ih
				histo_data_log(k) = (100.*cnt_log(i,ivar,ih))/n(i)
			enddo
		enddo
	endif
enddo
deallocate(cnt)
deallocate(dv)
deallocate(valmin)
deallocate(valmax)
deallocate(cnt_log)
deallocate(dv_log)
deallocate(valmin_log)
deallocate(valmax_log)
!deallocate(histo_data_log)
end subroutine

!--------------------------------------------------------------------------------
! This version computes concentrations on a line through the blob centre.
!--------------------------------------------------------------------------------
subroutine WriteProfileData1
integer :: ns
real(REAL_KIND) :: dx
real(REAL_KIND), allocatable :: ex_conc(:,:)
real(REAL_KIND) :: cbnd, cmin = 1.0e-7
integer :: rng(3,2), i, k, ichemo, kcell, x, y, z, ic, nc, kmax
character*(16) :: title(1:MAX_CHEMO+N_EXTRA)
character*(128) :: filename
character*(6) :: mintag

dx = DELTA_X
rng(:,1) = blob_centre(:) - (adrop*blob_radius + 2)
rng(:,2) = blob_centre(:) + (adrop*blob_radius + 2)
kmax = rng(1,2)-rng(1,1)+ 3
allocate(ex_conc(1:MAX_CHEMO+N_EXTRA,kmax))
y = blob_centre(2) + 0.5
z = blob_centre(3) + 0.5
ic = 0
do ichemo = 0,MAX_CHEMO+N_EXTRA
	if (ichemo == CFSE) then
		if (.not.chemo(ichemo)%used) cycle
		ic = ic + 1
		title(ic) = 'CFSE'
	    cbnd = 0
	elseif (ichemo <= MAX_CHEMO) then
		if (.not.chemo(ichemo)%used) cycle
		ic = ic + 1
		title(ic) = chemo(ichemo)%name
	    cbnd = BdryConc(ichemo)
	elseif (ichemo == GROWTH_RATE) then
		ic = ic + 1
		title(ic) = 'Growth_rate'
		cbnd = 0
	elseif (ichemo == CELL_VOLUME) then
		ic = ic + 1
		title(ic) = 'Cell_volume'
		cbnd = 0
	elseif (ichemo == O2_BY_VOL) then
		ic = ic + 1
		title(ic) = 'Cell_O2xVol'
		cbnd = 0
	endif
	k = 1
    ex_conc(ic,k) = cbnd
	do x = rng(1,1),rng(1,2)
		kcell = occupancy(x,y,z)%indx(1)
		if (kcell <= OUTSIDE_TAG) cycle
		i = ODEdiff%ivar(x,y,z)
		k = k+1
		if (ichemo == CFSE) then
			if (kcell > 0) then
				ex_conc(ic,k) = cell_list(kcell)%cfse
			else
				ex_conc(ic,k) = 0
			endif
		elseif (ichemo <= MAX_CHEMO) then
			if (i > 0) then
				ex_conc(ic,k) = allstate(i,ichemo)
			else
				ex_conc(ic,k) = 0
			endif
		elseif (ichemo == GROWTH_RATE) then
			if (kcell > 0) then
				ex_conc(ic,k) = cell_list(kcell)%dVdt
			else
				ex_conc(ic,k) = 0
			endif
		elseif (ichemo == CELL_VOLUME) then
			if (kcell > 0) then
				ex_conc(ic,k) = Vcell_pL*cell_list(kcell)%volume
			else
				ex_conc(ic,k) = 0
			endif
		elseif (ichemo == O2_BY_VOL) then
			if (kcell > 0) then
				ex_conc(ic,k) = allstate(i,OXYGEN)*Vcell_pL*cell_list(kcell)%volume
			else
				ex_conc(ic,k) = 0
			endif
		endif
	enddo
	k = k+1
    ex_conc(ic,k) = cbnd	! Add concentration at the boundary  
enddo
ns = k
nc = ic
! Remove time tag from the filename for download from NeSI
!write(mintag,'(i6)') int(istep*DELTA_T/60)
filename = saveprofile%filebase
!filename = trim(filename)//'_'
!filename = trim(filename)//trim(adjustl(mintag))
!filename = trim(filename)//'min.dat'
open(nfprofile,file=filename,status='replace')
write(nfprofile,'(a,a)') 'GUI version: ',gui_run_version
write(nfprofile,'(a,a)') 'DLL version: ',dll_run_version
write(nfprofile,'(i6,a)') int(istep*DELTA_T/60),' minutes'
write(nfprofile,'(i6,a)') ns,' sites'
write(nfprofile,'(f6.2,a)') 1000*DELTA_X,' dx (um)'
write(nfprofile,'(32a16)') title(1:nc)
do k = 1,ns
	write(nfprofile,'(32(e12.3,4x))') ex_conc(1:nc,k)
enddo
close(nfprofile)
deallocate(ex_conc)
end subroutine

!--------------------------------------------------------------------------------
! This version computes concentrations in a spherical shell of thickness dr.
! Note that for a cell at (ix,iy,iz):
! x = (ix-0.5)*dx, y = (iy-0.5)*dx, z = (iz-0.5)*dx
!--------------------------------------------------------------------------------
subroutine WriteProfileData
integer :: ns
real(REAL_KIND), allocatable :: ex_conc(:,:)
real(REAL_KIND) :: dr = 20.e-4	! 20um -> cm
real(REAL_KIND) :: dx, xyz0(3), dxyz(3), r2, r
integer :: ntot, ir, nr, ichemo, kcell, site(3)
integer :: i, ic, nc
integer, parameter :: max_shells = 100
integer :: cnt(max_shells)
integer :: icmap(0:MAX_CHEMO+N_EXTRA)		! maps ichemo -> ic
character*(16) :: title(1:MAX_CHEMO+N_EXTRA)
character*(128) :: filename
character*(6) :: mintag
type(cell_type), pointer :: cp

!write(nflog,*) 'WriteProfileData: MAX_CHEMO,N_EXTRA: ',MAX_CHEMO,N_EXTRA
! First find the centre of the blob
dx = DELTA_X
xyz0 = 0
ntot = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	ntot = ntot+1
	xyz0 = xyz0 + (cp%site - 0.5)*dx
enddo
if (ntot == 0) return
xyz0 = xyz0/ntot	! this is the blob centre

! Set up icmap
ic = 0
do ichemo = 0,MAX_CHEMO+N_EXTRA
	if (ichemo == CFSE) then
		if (.not.chemo(ichemo)%used) cycle
		ic = ic + 1
		title(ic) = 'CFSE'
	    icmap(ichemo) = ic
	elseif (ichemo <= MAX_CHEMO) then
		if (.not.chemo(ichemo)%used) cycle
		ic = ic + 1
		title(ic) = chemo(ichemo)%name
	    icmap(ichemo) = ic
	elseif (ichemo == GROWTH_RATE) then
		ic = ic + 1
		title(ic) = 'Growth_rate'
	    icmap(ichemo) = ic
	elseif (ichemo == CELL_VOLUME) then
		ic = ic + 1
		title(ic) = 'Cell_volume'
	    icmap(ichemo) = ic
!	elseif (ichemo == O2_BY_VOL) then
!		ic = ic + 1
!		title(ic) = 'Cell_O2xVol'
!	    icmap(ichemo) = ic
	endif
enddo
nc = ic
allocate(ex_conc(nc,max_shells))
ex_conc = 0

! Now look at shells at dr spacing. 
nr = 0
cnt = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	dxyz = (cp%site - 0.5)*dx - xyz0
	r2 = 0
	do i = 1,3
		r2 = r2 + dxyz(i)**2
	enddo
	r = sqrt(r2)
	ir = r/dr + 1
	nr = max(ir,nr)
	cnt(ir) = cnt(ir) + 1
	do ichemo = 0,MAX_CHEMO+N_EXTRA
		ic = icmap(ichemo)
		if (ichemo == CFSE .and. chemo(ichemo)%used) then
			ex_conc(ic,ir) = ex_conc(ic,ir) + cp%cfse
		elseif (ichemo <= MAX_CHEMO .and. chemo(ichemo)%used) then
!			if (cp%conc(ichemo) > 10) then
!			write(*,'(a,3i6,e12.3)') 'bad conc: ',kcell,ichemo,ic,cp%conc(ichemo)
!			stop
!			endif
			ex_conc(ic,ir) = ex_conc(ic,ir) + cp%conc(ichemo)
		elseif (ichemo == GROWTH_RATE) then
			ex_conc(ic,ir) = ex_conc(ic,ir) + cp%dVdt
		elseif (ichemo == CELL_VOLUME) then
			ex_conc(ic,ir) = ex_conc(ic,ir) + Vcell_pL*cp%volume
		endif
	enddo
enddo

! Average
do ir = 1,nr
	do ic = 1,nc
		ex_conc(ic,ir) = ex_conc(ic,ir)/cnt(ir)
	enddo
enddo			
	
! Remove time tag from the filename for download from NeSI
write(mintag,'(i6)') int(istep*DELTA_T/60)
filename = saveprofile%filebase
filename = trim(filename)//'_'
filename = trim(filename)//trim(adjustl(mintag))
filename = trim(filename)//'min.dat'
open(nfprofile,file=filename,status='replace')
write(nfprofile,'(a,a)') 'GUI version: ',gui_run_version
write(nfprofile,'(a,a)') 'DLL version: ',dll_run_version
write(nfprofile,'(i6,a)') int(istep*DELTA_T/60),' minutes'
write(nfprofile,'(i6,a)') nr,' shells'
write(nfprofile,'(f6.2,a)') 10000*dr,' dr (um)'
write(nfprofile,'(32a16)') title(1:nc)
do ir = 1,nr
	write(nfprofile,'(32(e12.3,4x))') ex_conc(1:nc,ir)
enddo
close(nfprofile)
deallocate(ex_conc)
!write(nflog,*) 'did WriteProfileData: ',filename
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine WriteSliceData
character*(128) :: filename
character*(6) :: mintag
integer, parameter :: NSMAX = 200
integer :: iz, nslices
real(REAL_KIND) :: dzslice
real(REAL_KIND) :: rad(NSMAX,2)

! Remove time tag from the filename for download from NeSI
write(mintag,'(i6)') int(istep*DELTA_T/60)
filename = saveslice%filebase
filename = trim(filename)//'_'
filename = trim(filename)//trim(adjustl(mintag))
filename = trim(filename)//'min.dat'
open(nfslice,file=filename,status='replace')
write(nfslice,'(a,a)') 'GUI version: ',gui_run_version
write(nfslice,'(a,a)') 'DLL version: ',dll_run_version
write(nfslice,'(i6,a)') int(istep*DELTA_T/60),' minutes'
call getSlices(nslices,dzslice,NSMAX,rad)
write(nfslice,'(e12.3,a)') dzslice,'  slice thickness (cm)'
write(nfslice,'(i4,a)') nslices,'  number of slices'
do iz = 1,nslices
    write(nfslice,'(i4,3e12.3)') iz,rad(iz,:),rad(iz,2)-rad(iz,1)
enddo
close(nfslice)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine WriteFACSData
character*(128) :: filename
character*(6) :: mintag
!integer :: ic, nc, ichem(2)
type(cell_type), pointer :: cp
real(REAL_KIND) :: hour, val, facs_data(32)
integer :: k, kcell, iextra, ichemo, ivar, nvars, var_index(32)
character*(24) :: var_name(32)

hour = istep*DELTA_T/3600.
nvars = 1	! CFSE
var_index(nvars) = 0
var_name(nvars) = 'CFSE'
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	nvars = nvars + 1
	var_index(nvars) = ichemo
	var_name(nvars) = chemo(ichemo)%name
enddo
do iextra = 1,N_EXTRA-1
	nvars = nvars + 1
	var_index(nvars) = MAX_CHEMO + iextra
	if (iextra == 1) then
		var_name(nvars) = 'GROWTH_RATE'
	elseif (iextra == 2) then
		var_name(nvars) = 'CELL_VOLUME'
	elseif (iextra == 3) then
		var_name(nvars) = 'O2_BY_VOL'
	endif
enddo

!nc = 0
!if (chemo(DRUG_A)%present) then
!	nc = nc + 1
!	ichem(nc) = DRUG_A
!endif
!if (chemo(DRUG_B)%present) then
!	nc = nc + 1
!	ichem(nc) = DRUG_B
!endif
! Remove time tag from the filename for download from NeSI
write(mintag,'(i6)') int(istep*DELTA_T/60)
filename = saveFACS%filebase
filename = trim(filename)//'_'
filename = trim(filename)//trim(adjustl(mintag))
filename = trim(filename)//'min.dat'
open(nfFACS,file=filename,status='replace')
write(nfFACS,'(a,a,2a12,a,i8,a,f8.2)') trim(header),' ',gui_run_version, dll_run_version, ' step: ',istep, ' hour: ',hour
do ivar = 1,nvars
	write(nfFACS,'(3x,a,$)') trim(var_name(ivar))
enddo
write(nfFACS,*)
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
!	write(nfFACS,'(2e14.5)') (cp%conc(ichem(ic)),ic=1,nc)
	k = 0
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cell_list(kcell)%CFSE
		elseif (ichemo <= MAX_CHEMO) then
			val = cell_list(kcell)%conc(ichemo)
		elseif (ichemo == GROWTH_RATE) then
			val = cell_list(kcell)%dVdt
		elseif (ichemo == CELL_VOLUME) then
			val = cell_list(kcell)%volume
		elseif (ichemo == O2_BY_VOL) then
			val = cell_list(kcell)%volume*cell_list(kcell)%conc(OXYGEN)
		endif
		k = k+1
		facs_data(k) = val
	enddo
	write(nfFACS,'(32e14.5)') facs_data(1:nvars)
enddo
close(nfFACS)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_constituents(nvars,cvar_index,nvarlen,name_array,narraylen) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_constituents
use, intrinsic :: iso_c_binding
character(c_char) :: name_array(0:*)
integer(c_int) :: nvars, cvar_index(0:*), nvarlen, narraylen
integer :: ivar, k, ichemo
character*(24) :: name
character(c_char) :: c

write(nflog,*) 'get_constituents'
nvarlen = 24
ivar = 0
k = ivar*nvarlen
cvar_index(ivar) = 0	! CFSE
name = 'CFSE'
call copyname(name,name_array(k),nvarlen)
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	ivar = ivar + 1
	k = ivar*nvarlen
	cvar_index(ivar) = ichemo
	name = chemo(ichemo)%name
	write(nflog,*) 'get_constituents: ',ichemo,name
	call copyname(name,name_array(k),nvarlen)
enddo
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = GROWTH_RATE
name = 'Growth rate'
call copyname(name,name_array(k),nvarlen)
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = CELL_VOLUME
name = 'Cell volume'
call copyname(name,name_array(k),nvarlen)
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = O2_BY_VOL
name = 'Cell O2xVol'
call copyname(name,name_array(k),nvarlen)
nvars = ivar + 1
write(nflog,*) 'did get_constituents'
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_clumpiness(ityp)
integer :: ityp
integer :: ntypcells, kcell0, kcell, k0, k, n, ix, iy, iz, site0(3), site(3)
real(REAL_KIND) :: df, fraction, dist(0:8)
type(cell_type), pointer :: cp

ntypcells = 0
fraction = 0
dist = 0
do kcell0 = 1,nlist
	cp => cell_list(kcell0)
	if (cp%state == DEAD) cycle
	if (cp%celltype == ityp) then
		ntypcells = ntypcells + 1
		n = 0
		site0 = cell_list(kcell0)%site
		do ix = -1,1,2
			do iy = -1,1,2
				do iz = -1,1,2
					site = site0 + [ix,iy,iz]
					kcell = occupancy(site(1),site(2),site(3))%indx(1)
					if (kcell > 0) then
						if (cell_list(kcell)%celltype == ityp) then
							n = n+1
						endif
					endif
				enddo
			enddo
		enddo
		df = n/8.
		fraction = fraction + df
		dist(n) = dist(n) + 1
	endif
enddo

dist = dist/ntypcells
write(nflog,*)
write(nflog,'(a)') 'Clumpiness:'
write(nflog,'(a,i1,f8.3)') 'Actual fraction of activator cells: celltype: ',ityp,real(ntypcells)/Ncells
write(nflog,'(a,f8.3)') 'Average fraction of neighbour activator cells: ',fraction/ntypcells
write(nflog,'(a,9f8.4)') 'Probability of n = 0,..8: ',dist(0:8)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine copyname(name,name_array,n)
character*(*) :: name
character :: name_array(*)
integer :: n
integer :: k

do k = 1,n
	name_array(k) = name(k:k)
enddo
end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_string(bufptr) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_string
use, intrinsic :: iso_c_binding
type(c_ptr) :: bufptr
character(c_char) :: buf(1024)
integer :: buflen
character*(1024), save :: string

string = 'A test string'
buflen = len(trim(string))
!write(*,*) 'buflen: ',buflen
string(buflen+1:buflen+1) = char(0)
bufptr = c_loc(string)
end subroutine

!--------------------------------------------------------------------------------
! Returns all the extracellular concentrations along a line through the blob centre.
! Together with CFSE, growth rate (dVdt), cell volume,...
!--------------------------------------------------------------------------------
subroutine old_get_concdata(ns, dx, ex_conc) BIND(C)
!!DEC$ ATTRIBUTES DLLEXPORT :: get_concdata
use, intrinsic :: iso_c_binding
integer(c_int) :: ns
real(c_double) :: dx, ex_conc(*)
real(REAL_KIND) :: cbnd, cmin = 1.0e-7
integer :: rng(3,2), i, ic, k, ichemo, kcell, x, y, z, nvars

nvars = 1 + MAX_CHEMO + N_EXTRA
!call logger('get_concdata')
dx = DELTA_X
rng(:,1) = blob_centre(:) - (adrop*blob_radius + 2)
rng(:,2) = blob_centre(:) + (adrop*blob_radius + 2)
!rng(axis,:) = Centre(axis) + fraction*Radius
y = blob_centre(2) + 0.5
z = blob_centre(3) + 0.5
ns = 1
do x = rng(1,1),rng(1,2)
    kcell = occupancy(x,y,z)%indx(1)
    if (kcell <= OUTSIDE_TAG) cycle
	i = ODEdiff%ivar(x,y,z)
    ns = ns + 1
    do ichemo = 1,nvars
        k = (ns-1)*nvars + ic
        if (ichemo == 0) then	! CFSE
			if (kcell > 0) then
				ex_conc(k) = cell_list(kcell)%CFSE
			else
				ex_conc(k) = 0
			endif
       elseif (ichemo <= MAX_CHEMO) then
			if (chemo(ichemo)%used) then
				if (i > 0) then
					ex_conc(k) = allstate(i,ichemo)	
				else
					ex_conc(k) = 0
				endif
			else
				ex_conc(k) = 0
			endif
        elseif (ichemo == MAX_CHEMO+1) then	! growth rate
			if (kcell > 0) then
				ex_conc(k) = cell_list(kcell)%dVdt
			else
				ex_conc(k) = 0
			endif
! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,ns,kcell,k,ex_conc(k)
        elseif (ichemo == MAX_CHEMO+2) then	! cell volume
			if (kcell > 0) then
				ex_conc(k) = cell_list(kcell)%volume
			else
				ex_conc(k) = 0
			endif
! 			write(nflog,'(a,4i6,f8.3)') 'Cell volume: ',x,ns,i,k,ex_conc(k)
		endif
    enddo
enddo
! Add concentrations at the two boundaries 
! At ns=1, at at ns=ns+1
ns = ns+1
do ic = 1,nvars
	ichemo = ic - 1
	if (ichemo == 0) then	! CFSE
		k = ic	
		ex_conc(k) = 0
        k = (ns-1)*nvars + ic
        ex_conc(k) = 0	
	elseif (ichemo <= MAX_CHEMO) then
		k = ic
		if (chemo(ichemo)%used) then
			ex_conc(k) = BdryConc(ichemo)
		else
			ex_conc(k) = 0
		endif      
		k = (ns-1)*nvars + ic
		if (chemo(ichemo)%used) then
			ex_conc(k) = BdryConc(ichemo)
		else
			ex_conc(k) = 0
		endif      
    elseif (ichemo == MAX_CHEMO+1) then	! growth rate
		k = ic
		ex_conc(k) = 0
! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,1,i,k,ex_conc(k)
        k = (ns-1)*nvars + ic
        ex_conc(k) = 0
! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,ns,i,k,ex_conc(k)
    elseif (ichemo == MAX_CHEMO+2) then	! cell volume
		k = ic
		ex_conc(k) = 0
! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,1,i,k,ex_conc(k)
        k = (ns-1)*nvars + ic
        ex_conc(k) = 0
! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,ns,i,k,ex_conc(k)
    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! EC values should be average in the medium.
! DRUG_A = TPZ_DRUG e.g. SN30000
! DRUG_B = DNB_DRUG e.g. PR104A 
!-----------------------------------------------------------------------------------------
subroutine get_values(nvars,varID,ysim)
!DEC$ ATTRIBUTES DLLEXPORT :: get_values
integer :: nvars
character*(24) :: varID(nvars)
real(REAL_KIND) :: ysim(nvars)
integer :: ivar, ityp
integer :: Nviable(MAX_CELLTYPES), Nlive(MAX_CELLTYPES)
real(REAL_KIND) :: plate_eff(MAX_CELLTYPES)

call MakeMediumCave
do ivar = 1,nvars
	if (varID(ivar) == 'OXYGEN_EC') then
		ysim(ivar) = chemo(OXYGEN)%medium_Cave
	elseif (varID(ivar) == 'GLUCOSE_EC') then
		ysim(ivar) = chemo(GLUCOSE)%medium_Cave
!	elseif (varID(ivar) == 'LACTATE_EC') then
!		ysim(ivar) = chemo(LACTATE)%medium_Cave
	elseif (varID(ivar) == 'SN30000_EC') then
		ysim(ivar) = chemo(DRUG_A)%medium_Cave
	elseif (varID(ivar) == 'SN30000_METAB1_EC') then
		ysim(ivar) = chemo(DRUG_A+1)%medium_Cave
	elseif (varID(ivar) == 'SN30000_METAB2_EC') then
		ysim(ivar) = chemo(DRUG_A+2)%medium_Cave
	elseif (varID(ivar) == 'PR104A_EC') then
		ysim(ivar) = chemo(DRUG_B)%medium_Cave
	elseif (varID(ivar) == 'PR104A_METAB1_EC') then
		ysim(ivar) = chemo(DRUG_B+1)%medium_Cave
	elseif (varID(ivar) == 'PR104A_METAB2_EC') then
		ysim(ivar) = chemo(DRUG_B+2)%medium_Cave
	elseif (varID(ivar) == 'OXYGEN_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+OXYGEN)
		ysim(ivar) = 0
	elseif (varID(ivar) == 'GLUCOSE_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+GLUCOSE)
		ysim(ivar) = 0
!	elseif (varID(ivar) == 'LACTATE_IC') then
!!		ysim(ivar) = Caverage(MAX_CHEMO+LACTATE)
!		ysim(ivar) = 0
	elseif (varID(ivar) == 'SN30000_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+DRUG_A)
		ysim(ivar) = 0
	elseif (varID(ivar) == 'SN30000_METAB1_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+DRUG_A+1)
		ysim(ivar) = 0
	elseif (varID(ivar) == 'SN30000_METAB2_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+DRUG_A+2)
		ysim(ivar) = 0
	elseif (varID(ivar) == 'PR104A_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+DRUG_A)
		ysim(ivar) = 0
	elseif (varID(ivar) == 'PR104A_METAB1_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+DRUG_A+1)
		ysim(ivar) = 0
	elseif (varID(ivar) == 'PR104A_METAB2_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+DRUG_A+2)
		ysim(ivar) = 0
	elseif (varID(ivar) == 'NCELLS') then
		ysim(ivar) = Ncells			! for now, total live cells
	elseif (varID(ivar) == 'PE') then
		call getNviable(Nviable, Nlive)
		do ityp = 1,Ncelltypes
			if (Nlive(ityp) > 0) then
				plate_eff(ityp) = real(Nviable(ityp))/Nlive(ityp)
			else
				plate_eff(ityp) = 0
			endif
		enddo
		ysim(ivar) = plate_eff(1)	! for now, just type 1 cells
	elseif (varID(ivar) == 'RADIATION') then
		ysim(ivar) = -1
	else
		write(*,*) 'varID is not in the list of possible IDs: ',varID(ivar)
		stop
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine MakeMediumCave
integer :: ichemo

do ichemo = 1,MAX_CHEMO
	if (chemo(ichemo)%used) then
!		chemo(ichemo)%medium_Cave = sum(chemo(ichemo)%Cave_b)/(NXB*NYB*NZB)
		chemo(ichemo)%medium_Cave = chemo(ichemo)%medium_Cext
	endif
enddo
end subroutine
end module
