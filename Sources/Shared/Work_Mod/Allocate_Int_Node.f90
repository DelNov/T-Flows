!==============================================================================!
  subroutine Allocate_Int_Node(Work, Grid, n)
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

  allocate(Work % i_node(n))

  do i = 1, n
    allocate(Work % i_node(i) % ptr(1:nn))
    Work % i_node(i) % ptr(:) = 0.0
  end do

  Work % last_i_node = 0

  end subroutine
