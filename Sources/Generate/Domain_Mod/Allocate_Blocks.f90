!==============================================================================!
  subroutine Allocate_Blocks(Dom, n)
!------------------------------------------------------------------------------!
!>  Sets the number of blocks in the domain and allocates memory for them.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Domain_Type)  :: Dom  !! domain under consideration
  integer, intent(in) :: n    !! number of blocks
!==============================================================================!

  Dom % n_blocks = n
  allocate(Dom % blocks(n))

  end subroutine
