!==============================================================================!
  subroutine Allocate_Real_Cell(Work, Grid, n)
!------------------------------------------------------------------------------!
!>  Allocates memory for real-typed working arrays associated with cells
!>  in the Work object.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work     !! parent; the singleton Work object
  type(Grid_Type)  :: Grid(:)  !! grids on which the Work will be used
  integer          :: n        !! number of real cell arrays
!-----------------------------------[Locals]-----------------------------------!
  integer :: nc, nb, i
!==============================================================================!

  if(n .eq. 0) return

  ! Get number of cells and boundary cells
  nc = maxval(Grid(1:size(Grid)) % n_cells)
  nb = maxval(Grid(1:size(Grid)) % n_bnd_cells)

  allocate(Work % r_cell(n))

  do i = 1, n
    allocate(Work % r_cell(i) % array(-nb:nc))
    Work % r_cell(i) % array(:) = 0.0
  end do

  Work % last_r_cell = 0

  end subroutine
