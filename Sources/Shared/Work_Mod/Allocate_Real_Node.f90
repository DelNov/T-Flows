!==============================================================================!
  subroutine Allocate_Real_Node(Work, Grid, n)
!------------------------------------------------------------------------------!
!>  Allocates memory for integer-typed working arrays associated with nodes
!>  in the Work object.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work     !! parent; the singleton Work object
  type(Grid_Type)  :: Grid(:)  !! grids on which the Work will be used
  integer          :: n        !! number of real node arrays
!-----------------------------------[Locals]-----------------------------------!
  integer :: nn, i
!==============================================================================!

  if(n .eq. 0) return

  ! Get number of nodes
  nn = maxval(Grid(1:size(Grid)) % n_nodes)

  allocate(Work % r_node(n))

  do i = 1, n
    allocate(Work % r_node(i) % array(1:nn))
    Work % r_node(i) % array(:) = 0.0
  end do

  Work % last_r_node = 0

  end subroutine
