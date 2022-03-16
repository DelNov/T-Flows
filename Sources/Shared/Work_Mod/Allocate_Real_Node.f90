!==============================================================================!
  subroutine Allocate_Real_Node(Work, Grid, n)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work
  type(Grid_Type)  :: Grid(:)
  integer          :: n     ! number of real cell arrays
!-----------------------------------[Locals]-----------------------------------!
  integer :: nn, i
!==============================================================================!

  if(n .eq. 0) return

  ! Get number of cells and boundary cells
  nn = maxval(Grid(1:size(Grid)) % n_nodes)

  allocate(Work % r_node(n))

  do i = 1, n
    allocate(Work % r_node(i) % ptr(1:nn))
    Work % r_node(i) % ptr(:) = 0.0
  end do

  Work % last_r_node = 0

  end subroutine
