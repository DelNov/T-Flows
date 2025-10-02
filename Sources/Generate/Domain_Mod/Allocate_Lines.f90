!==============================================================================!
  subroutine Allocate_Lines(Dom, n)
!------------------------------------------------------------------------------!
!>  Sets the number of lines in the domain and allocates memory for them.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Domain_Type)  :: Dom  !! domain under consideration
  integer, intent(in) :: n    !! number of lines in the .dom file
!==============================================================================!

  Dom % n_lines = n
  allocate(Dom % lines(n))

  end subroutine
