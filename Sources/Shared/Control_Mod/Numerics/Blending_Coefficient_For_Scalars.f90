!==============================================================================!
  subroutine Blending_Coefficient_For_Scalars(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('BLENDING_COEFFICIENT_FOR_SCALARS',  &
                                 1.0, val, verbose)

  end subroutine
