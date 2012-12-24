! Chemokine data

module chemokine

use global

implicit none

!type receptor_type
!	character(8) :: name
!	logical :: used
!	integer :: chemokine
!	integer :: sign
!	real(REAL_KIND) :: strength
!	real(REAL_KIND) :: level(5)
!	real(REAL_KIND) :: saturation_threshold
!	real(REAL_KIND) :: refractory_time
!end type
!type(receptor_type), target :: receptor(MAX_RECEPTOR)

type chemokine_type
	character(8) :: name
	logical :: used
	logical :: use_secretion
	real(REAL_KIND) :: bdry_rate
	real(REAL_KIND) :: bdry_conc
	real(REAL_KIND) :: diff_coef
	real(REAL_KIND) :: cell_diff
	real(REAL_KIND) :: halflife
	real(REAL_KIND) :: decay_rate
	logical :: bdry_decay
	real(REAL_KIND) :: bdry_halflife
	real(REAL_KIND) :: bdry_decay_rate
	real(REAL_KIND) :: max_cell_rate		! Vmax
	real(REAL_KIND) :: MM_C0			! Km
	real(REAL_KIND), allocatable :: coef(:,:)
	real(REAL_KIND), allocatable :: conc(:,:,:)
	real(REAL_KIND), allocatable :: grad(:,:,:,:)
end type
type(chemokine_type), target :: chemo(MAX_CHEMO)

type ODEdiff_type
	integer :: ichemo
	integer :: nextra
	integer :: nintra
	integer :: nvars
	integer, allocatable :: ivar(:,:,:)
	integer, allocatable :: varsite(:,:)
	integer, allocatable :: icoef(:,:)
	integer, allocatable :: vartype(:)
	integer, allocatable :: cell_index(:)
!	integer, allocatable :: ncoef(:)
	real(REAL_KIND) :: deltaC
	real(REAL_KIND) :: k
	real(REAL_KIND) :: C1
end type
type(ODEdiff_type) :: ODEdiff

end module