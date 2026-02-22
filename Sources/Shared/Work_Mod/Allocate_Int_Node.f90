!==============================================================================!
  subroutine Allocate_Int_Node(Work, r)
!------------------------------------------------------------------------------!
!>  Allocates memory for integer-typed working arrays associated with nodes
!>  in the Work object.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work     !! parent; the singleton Work object
  integer          :: r        !! number of integer node arrays
!==============================================================================!

  if(r .eq. 0) return

  allocate(Work % i_node(r) % array(1 : Work % max_nn))
  Work % i_node(r) % array(:) = 0

  call Gpu % Vector_Int_Copy_To_Device(lbound(Work % i_node(r) % array, 1),  &
                                       ubound(Work % i_node(r) % array, 1),  &
                                              Work % i_node(r) % array)

  end subroutine
