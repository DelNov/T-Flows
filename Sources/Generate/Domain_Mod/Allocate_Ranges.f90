!==============================================================================!
  subroutine Allocate_Ranges(Dom, n)
!------------------------------------------------------------------------------!
!>  Sets the number of ranges (for defining boundary conditions and materials)
!>  in the domain and allocates memory for these ranges.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Domain_Type)  :: Dom  !! domain under consideration
  integer, intent(in) :: n    !! number of ranges (boundary conditions
                              !! or materials)
!==============================================================================!

  Dom % n_ranges = n
  allocate(Dom % ranges(n))

  end subroutine
