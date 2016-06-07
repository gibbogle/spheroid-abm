! Transferring results to the GUI

module transfer

use global
use chemokine
use envelope
use ode_diffuse

#include "../src/version.h"

implicit none

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
subroutine getNecroticFraction(necrotic_fraction,totvol_cm3)
real(REAL_KIND) :: necrotic_fraction, totvol_cm3	! vol_cm3 not used here, needed in scell
real(REAL_KIND) :: cellvol_cm3, dvol
!necrotic_fraction = (Nsites-Ncells)/real(Nsites)
cellvol_cm3 = Ncells*DELTA_X**3
dvol = totvol_cm3-cellvol_cm3
necrotic_fraction = dvol/totvol_cm3
if (necrotic_fraction < 0.005) necrotic_fraction = 0
end subroutine

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
subroutine get_summary(summaryData,i_hypoxia_cutoff,i_growth_cutoff) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*), i_hypoxia_cutoff,i_growth_cutoff
integer :: Nviable(MAX_CELLTYPES), Nlive(MAX_CELLTYPES), plate_eff_10(MAX_CELLTYPES)
integer :: diam_um, vol_mm3_1000, nhypoxic(3), nclonohypoxic(3), ngrowth(3), &
    hypoxic_percent_10, clonohypoxic_percent_10, growth_percent_10, necrotic_percent_10,  npmm3, &
    medium_oxygen_1000, medium_glucose_1000, medium_drug_1000(2), &
    bdry_oxygen_1000, bdry_glucose_1000, bdry_drug_1000(2)
integer :: TNanoxia_dead, TNradiation_dead, TNdrug_dead(2),  TNviable, &
           Ntagged_anoxia(MAX_CELLTYPES), Ntagged_radiation(MAX_CELLTYPES), Ntagged_drug(2,MAX_CELLTYPES), &
           TNtagged_anoxia, TNtagged_radiation, TNtagged_drug(2)
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
Ntagged_radiation(:) = Nradiation_tag(:)	! number currently tagged by radiation
Ntagged_drug(1,:) = Ndrug_tag(1,:)			! number currently tagged by drugA
Ntagged_drug(2,:) = Ndrug_tag(2,:)			! number currently tagged by drugA

TNtagged_anoxia = sum(Ntagged_anoxia(1:Ncelltypes))
TNtagged_radiation = sum(Ntagged_radiation(1:Ncelltypes))
TNtagged_drug(1) = sum(Ntagged_drug(1,1:Ncelltypes))
TNtagged_drug(2) = sum(Ntagged_drug(2,1:Ncelltypes))

TNanoxia_dead = sum(Nanoxia_dead(1:Ncelltypes))
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
!if (TNanoxia_dead > 0) then
	call getNecroticFraction(necrotic_fraction,vol_cm3)
!else
!	necrotic_fraction = 0
!endif
necrotic_percent_10 = 1000.*necrotic_fraction
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

summaryData(1:26) = [ istep, Ncells, TNanoxia_dead, TNdrug_dead(1), TNdrug_dead(2), TNradiation_dead, &
    TNtagged_anoxia, TNtagged_drug(1), TNtagged_drug(2), TNtagged_radiation, &
	diam_um, vol_mm3_1000, hypoxic_percent_10, clonohypoxic_percent_10, growth_percent_10, necrotic_percent_10, &
	Tplate_eff_10, npmm3, &
	medium_oxygen_1000, medium_glucose_1000, medium_drug_1000(1), medium_drug_1000(2), &
	bdry_oxygen_1000, bdry_glucose_1000, bdry_drug_1000(1), bdry_drug_1000(2) ]
write(nfres,'(2a12,i8,2e12.4,19i7,20e12.4)') gui_run_version, dll_run_version, &
	istep, hour, vol_mm3, diam_um, Ncells_type(1:2), &
    Nanoxia_dead(1:2), Ndrug_dead(1,1:2), &
    Ndrug_dead(2,1:2), Nradiation_dead(1:2), &
    Ntagged_anoxia(1:2), Ntagged_drug(1,1:2), &
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

end subroutine

!-----------------------------------------------------------------------------------------
! Live cells = Ncells
! Drug deaths = Ndrugdead
! Hypoxia deaths = Ndead - Ndrugdead
! Total tagged for drug death on division = NdrugA_tag, NdrugB_tag
! Current tagged = Ntodie - Ntagdead 
!-----------------------------------------------------------------------------------------
subroutine get_summary1(summaryData,i_hypoxia_cutoff,i_growth_cutoff) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*), i_hypoxia_cutoff,i_growth_cutoff
integer :: Nviable(MAX_CELLTYPES), Nlive(MAX_CELLTYPES), plate_eff_10(MAX_CELLTYPES)
integer :: diam_um, vol_mm3_1000, nhypoxic(3), ngrowth(3), hypoxic_percent_10, growth_percent_10, necrotic_percent_10,  npmm3, &
    medium_oxygen_1000, medium_glucose_1000, medium_drug_1000(2), &
    bdry_oxygen_1000, bdry_glucose_1000, bdry_drug_1000(2)
integer :: TNanoxia_dead, TNradiation_dead, TNdrug_dead(2),  &
           Ntagged_anoxia(MAX_CELLTYPES), Ntagged_radiation(MAX_CELLTYPES), Ntagged_drug(2,MAX_CELLTYPES), &
           TNtagged_anoxia, TNtagged_radiation, TNtagged_drug(2)
integer :: Tplate_eff_10   
integer :: ityp
real(REAL_KIND) :: diam_cm, vol_cm3, vol_mm3, hour, plate_eff(MAX_CELLTYPES), necrotic_fraction
real(REAL_KIND) :: cmedium(MAX_CHEMO), cbdry(MAX_CHEMO)

if (i_hypoxia_cutoff == 0) stop

hour = istep*DELTA_T/3600.
call getDiamVol(diam_cm,vol_cm3)
vol_mm3 = vol_cm3*1000				! volume in mm^3
vol_mm3_1000 = vol_mm3*1000			! 1000 * volume in mm^3
diam_um = diam_cm*10000
npmm3 = Ncells/vol_mm3

Ntagged_anoxia(:) = Nanoxia_tag(:)			! number currently tagged by anoxia
Ntagged_radiation(:) = Nradiation_tag(:)	! number currently tagged by radiation
Ntagged_drug(1,:) = Ndrug_tag(1,:)			! number currently tagged by drugA
Ntagged_drug(2,:) = Ndrug_tag(2,:)			! number currently tagged by drugA

TNtagged_anoxia = sum(Ntagged_anoxia(1:Ncelltypes))
TNtagged_radiation = sum(Ntagged_radiation(1:Ncelltypes))
TNtagged_drug(1) = sum(Ntagged_drug(1,1:Ncelltypes))
TNtagged_drug(2) = sum(Ntagged_drug(2,1:Ncelltypes))

TNanoxia_dead = sum(Nanoxia_dead(1:Ncelltypes))
TNradiation_dead = sum(Nradiation_dead(1:Ncelltypes))
TNdrug_dead(1) = sum(Ndrug_dead(1,1:Ncelltypes))
TNdrug_dead(2) = sum(Ndrug_dead(2,1:Ncelltypes))

call getHypoxicCount(nhypoxic)
hypoxic_percent_10 = (1000*nhypoxic(i_hypoxia_cutoff))/Ncells	! 10* %hypoxic
call getGrowthCount(ngrowth)
growth_percent_10 = (1000*ngrowth(i_growth_cutoff))/Ncells
!if (TNanoxia_dead > 0) then
	call getNecroticFraction(necrotic_fraction,vol_cm3)
!else
!	necrotic_fraction = 0
!endif
necrotic_percent_10 = 1000*necrotic_fraction
call getNviable(Nviable, Nlive)
do ityp = 1,Ncelltypes
	if (Nlive(ityp) > 0) then
		plate_eff(ityp) = real(Nviable(ityp))/Nlive(ityp)
	else
		plate_eff(ityp) = 0
	endif
enddo
plate_eff_10 = 1000*plate_eff
Tplate_eff_10 = 0
do ityp = 1,Ncelltypes
	Tplate_eff_10 = Tplate_eff_10 + plate_eff_10(ityp)*celltype_fraction(ityp)
enddo
call getMediumConc(cmedium, cbdry)
medium_oxygen_1000 = cmedium(OXYGEN)*1000
medium_glucose_1000 = cmedium(GLUCOSE)*1000
medium_drug_1000(1) = cmedium(DRUG_A)*1000
medium_drug_1000(2) = cmedium(DRUG_B)*1000
bdry_oxygen_1000 = cbdry(OXYGEN)*1000
bdry_glucose_1000 = cbdry(GLUCOSE)*1000
bdry_drug_1000(1) = cbdry(DRUG_A)*1000
bdry_drug_1000(2) = cbdry(DRUG_B)*1000

summaryData(1:25) = [ istep, Ncells, TNanoxia_dead, TNdrug_dead(1), TNdrug_dead(2), TNradiation_dead, &
    TNtagged_anoxia, TNtagged_drug(1), TNtagged_drug(2), TNtagged_radiation, &
	diam_um, vol_mm3_1000, hypoxic_percent_10, growth_percent_10, necrotic_percent_10, Tplate_eff_10, npmm3, &
	medium_oxygen_1000, medium_glucose_1000, medium_drug_1000(1), medium_drug_1000(2), &
	bdry_oxygen_1000, bdry_glucose_1000, bdry_drug_1000(1), bdry_drug_1000(2) ]
write(nfres,'(a,a,2a12,i8,2e12.4,19i7,17e12.4)') trim(header),' ',gui_run_version, dll_run_version, &
	istep, hour, vol_mm3, diam_um, Ncells_type(1:2), &
    Nanoxia_dead(1:2), Ndrug_dead(1,1:2), &
    Ndrug_dead(2,1:2), Nradiation_dead(1:2), &
    Ntagged_anoxia(1:2), Ntagged_drug(1,1:2), &
    Ntagged_drug(2,1:2), Ntagged_radiation(1:2), &
	nhypoxic(:)/real(Ncells), ngrowth(:)/real(Ncells), &
	necrotic_fraction, plate_eff(1:2), &
	cmedium(OXYGEN), cmedium(GLUCOSE), cmedium(DRUG_A), cmedium(DRUG_B), &
	cbdry(OXYGEN), cbdry(GLUCOSE), cbdry(DRUG_A), cbdry(DRUG_B)
		
call sum_dMdt(GLUCOSE)

if (diam_count_limit > LIMIT_THRESHOLD) then
	if (Ncells > diam_count_limit) limit_stop = .true.
elseif (diam_count_limit > 0) then
	if (diam_um > diam_count_limit) limit_stop = .true.
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
	if (cell_list(kcell)%anoxia_tag .or. cell_list(kcell)%radiation_tag) cycle
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
real(REAL_KIND) :: cbnd, cmin = 1.0e-6
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
! Returns all the extracellular concentrations along a line through the blob centre.
! Together with CFSE, growth rate (dVdt), cell volume,...
! Store the constituent profiles one after the other.
! Drop the extra end points
!--------------------------------------------------------------------------------
subroutine get_IC_concdata(nvars, ns, dx, ic_conc) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ic_concdata
use, intrinsic :: iso_c_binding
integer(c_int) :: nvars, ns
real(c_double) :: dx, ic_conc(0:*)
real(REAL_KIND) :: cbnd, cmin = 1.0e-6
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
		if (kcell <= OUTSIDE_TAG) then
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
!	k = k+1
!	facs_data(k) = cell_list(kcell)%CFSE
!	k = k+1
!	facs_data(k) = cell_list(kcell)%dVdt
!	k = k+1
!	facs_data(k) = cell_list(kcell)%conc(OXYGEN)
!	if (cell_list(kcell)%conc(OXYGEN) <= 0.00001 .or. cell_list(kcell)%dVdt < 2.0e-6) then
!		write(nflog,'(2i6,2e12.3)') istep,kcell,cell_list(kcell)%dVdt,cell_list(kcell)%conc(OXYGEN)
!	endif
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
real(REAL_KIND) :: cbnd, cmin = 1.0e-6
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
real(REAL_KIND) :: cbnd, cmin = 1.0e-6
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
 			write(nflog,'(a,4i6,f8.3)') 'Cell volume: ',x,ns,i,k,ex_conc(k)
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

end module
