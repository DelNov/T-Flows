!==============================================================================!
  subroutine Allocate_Points(Dom, n)
!------------------------------------------------------------------------------!
!>  Sets the number of points and allocates memory for point storage.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Domain_Type)  :: Dom  !! domain under consideration
  integer, intent(in) :: n    !! number of points in the .dom file
!==============================================================================!

  Dom % n_points = n
  allocate(Dom % points(n))

  end subroutine
