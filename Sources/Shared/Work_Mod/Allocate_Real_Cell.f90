!==============================================================================!
  subroutine Allocate_Real_Cell(Work, Grid, n)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work
  type(Grid_Type)  :: Grid(:)
  integer          :: n     ! number of real cell arrays
!-----------------------------------[Locals]-----------------------------------!
  integer :: nc, nb, i
!==============================================================================!

  if(n .eq. 0) return

  ! Get number of cells and boundary cells
  nc = maxval(Grid(1:size(Grid)) % n_cells)
  nb = maxval(Grid(1:size(Grid)) % n_bnd_cells)

  allocate(Work % r_cell(n))

  do i = 1, n
    allocate(Work % r_cell(i) % ptr(-nb:nc))
    Work % r_cell(i) % ptr(:) = 0.0
  end do

  Work % last_r_cell = 0

  end subroutine
