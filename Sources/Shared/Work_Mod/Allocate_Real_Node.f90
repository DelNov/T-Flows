!==============================================================================!
  subroutine Allocate_Real_Node(Work, r)
!------------------------------------------------------------------------------!
!>  Allocates memory for integer-typed working arrays associated with nodes
!>  in the Work object.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work     !! parent; the singleton Work object
  integer          :: r        !! number of real node arrays
!==============================================================================!

  if(r .eq. 0) return

  allocate(Work % r_node(r) % array(1 : Work % max_nn))
  Work % r_node(r) % array(:) = 0.0

  call Gpu % Vector_Real_Copy_To_Device(lbound(Work % r_node(r) % array, 1),  &
                                        ubound(Work % r_node(r) % array, 1),  &
                                               Work % r_node(r) % array)

  end subroutine
